#include "global_defs.h"
#include "element_definitions.h"
#include "petsc_citcoms.h"


/* return ||V||^2 */
double global_v_norm2_PETSc( struct All_variables *E,  Vec v )
{
    int i, m, d;
    int eqn1, eqn2, eqn3;
    double prod = 0.0, temp = 0.0;
    PetscErrorCode ierr;

    PetscScalar *V;
    ierr = VecGetArray( v, &V ); CHKERRQ( ierr );
    for (m=1; m<=E->sphere.caps_per_proc; m++)
        for (i=1; i<=E->lmesh.nno; i++) {
            eqn1 = E->id[m][i].doff[1];
            eqn2 = E->id[m][i].doff[2];
            eqn3 = E->id[m][i].doff[3];
            /* L2 norm  */
            temp += (V[eqn1] * V[eqn1] +
                     V[eqn2] * V[eqn2] +
                     V[eqn3] * V[eqn3]) * E->NMass[m][i];
        }
    ierr = VecRestoreArray( v, &V ); CHKERRQ( ierr );

    MPI_Allreduce(&temp, &prod, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);

    return (prod/E->mesh.volume);
}


/* return ||P||^2 */
double global_p_norm2_PETSc( struct All_variables *E,  Vec p )
{
    int i, m;
    double prod = 0.0, temp = 0.0;
    PetscErrorCode ierr;

    PetscScalar *P;
    ierr = VecGetArray( p, &P ); CHKERRQ( ierr );

    for (m=1; m<=E->sphere.caps_per_proc; m++)
        for (i=1; i<=E->lmesh.npno; i++) {
            /* L2 norm */
            temp += P[i] * P[i] * E->eco[m][i].area;
        }
    ierr = VecRestoreArray( p, &P ); CHKERRQ( ierr );

    MPI_Allreduce(&temp, &prod, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);

    return (prod/E->mesh.volume);
}

/* return ||A||^2, where A_i is \int{div(u) d\Omega_i} */
double global_div_norm2_PETSc( struct All_variables *E,  Vec a )
{
    int i, m;
    double prod = 0.0, temp = 0.0;
    PetscErrorCode ierr;

    PetscScalar *A;
    ierr = VecGetArray( a, &A ); CHKERRQ( ierr );


    for (m=1; m<=E->sphere.caps_per_proc; m++)
        for (i=1; i<=E->lmesh.npno; i++) {
            /* L2 norm of div(u) */
            temp += A[i] * A[i] / E->eco[m][i].area;

            /* L1 norm */
            /*temp += fabs(A[i]);*/
        }
    ierr = VecRestoreArray( a, &A ); CHKERRQ( ierr );

    MPI_Allreduce(&temp, &prod, 1, MPI_DOUBLE, MPI_SUM, E->parallel.world);

    return (prod/E->mesh.volume);
}

/* =====================================================
   Assemble grad(rho_ref*ez)*V element by element.
   Note that the storage is not zero'd before assembling.
   =====================================================  */

PetscErrorCode assemble_c_u_PETSc( struct All_variables *E,
                         Vec U, Vec result, int level )
//                  double **U, double **result, int level)
{
    int e,j1,j2,j3,p,a,b,m;

    const int nel = E->lmesh.NEL[level];
    const int ends = enodes[E->mesh.nsd];
    const int dims = E->mesh.nsd;
    const int npno = E->lmesh.NPNO[level];

    PetscErrorCode ierr;
    PetscScalar *U_temp, *result_temp;

    ierr = VecGetArray( U, &U_temp ); CHKERRQ( ierr );
    ierr = VecGetArray( result, &result_temp ); CHKERRQ( ierr );

  for( m = 1; m <= E->sphere.caps_per_proc; m++ ) {
    for(a=1;a<=ends;a++) {
      p = (a-1)*dims;
      for(e=1;e<=nel;e++) {
        b = E->IEN[level][m][e].node[a];
        j1= E->ID[level][m][b].doff[1];
        j2= E->ID[level][m][b].doff[2];
        j3= E->ID[level][m][b].doff[3];

        result_temp[e]  += E->elt_c[level][m][e].c[p  ][0] * U_temp[j1]
                         + E->elt_c[level][m][e].c[p+1][0] * U_temp[j2]
                         + E->elt_c[level][m][e].c[p+2][0] * U_temp[j3];
      }
    }
  }
  ierr = VecRestoreArray( U, &U_temp ); CHKERRQ( ierr );
  ierr = VecRestoreArray( result, &result_temp ); CHKERRQ( ierr );

  PetscFunctionReturn(0);
}

void strip_bcs_from_residual_PETSc( 
    struct All_variables *E, Vec Res, int level )
{
  int i, m, low, high;
  VecGetOwnershipRange( Res, &low, &high );
 
  for( m = 1; m <= E->sphere.caps_per_proc; m++ ) {
    if( E->num_zero_resid[level][m] ) {
      for( i = 1; i <= E->num_zero_resid[level][m]; i++ ) {
	      VecSetValue( Res, E->zero_resid[level][m][i]+low, 0.0, INSERT_VALUES );
      }
    }
  }
  VecAssemblyBegin( Res );
  VecAssemblyEnd( Res );
}

PetscErrorCode initial_vel_residual_PETSc( struct All_variables *E,
                                 Vec V, Vec P, Vec F,
                                 double acc )
{
    int neq = E->lmesh.neq;
    int lev = E->mesh.levmax;
    int npnp = E->lmesh.npno;
    int i, m, valid;
    PetscErrorCode ierr;

    Vec u1;
    ierr = VecCreateMPI( PETSC_COMM_WORLD, neq+1, PETSC_DECIDE, &u1 ); 
    CHKERRQ( ierr );

    /* F = F - grad(P) - K*V */
    // u1 = grad(P) i.e. G*P
    ierr = MatMult( E->G, P, u1 ); CHKERRQ( ierr );
    // F = F - u1
    ierr = VecAXPY( F, -1.0, u1 ); CHKERRQ( ierr ); 
    // u1 = del2(V) i.e. K*V
    ierr = MatMult( E->K, V, u1 ); CHKERRQ( ierr );
    // F = F - u1
    ierr = VecAXPY( F, -1.0, u1 ); CHKERRQ( ierr ); 

    strip_bcs_from_residual_PETSc(E, F, lev);

    /* solve K*u1 = F for u1 */
    //ierr = KSPSetTolerances( ... );
    ierr = KSPSolve( E->ksp, F, u1 ); CHKERRQ( ierr );

    strip_bcs_from_residual_PETSc(E, u1, lev);

    /* V = V + u1 */
    ierr = VecAXPY( V, 1.0, u1 ); CHKERRQ( ierr );
  PetscFunctionReturn(0);
}

PetscErrorCode PC_Apply_MultiGrid( PC pc, Vec x, Vec y )
{

  PetscErrorCode ierr;
  struct MultiGrid_PC *ctx;
  ierr = PCShellGetContext( pc, (void **)&ctx ); CHKERRQ( ierr );
  int count, valid;
  double residual;
  int m, i;


  int low, high;
  VecGetOwnershipRange( x, &low, &high );

  VecAssemblyBegin(x);
  VecAssemblyEnd(x);
  for( m=1; m<=ctx->E->sphere.caps_per_proc; ++m ) {
    for( i=0; i<ctx->nno; ++i ) {
      PetscInt ix[] = {i+low};
      ierr = VecGetValues( x, 1, ix, &ctx->RR[m][i] );
      CHKERRQ( ierr );
    }
  }

  /* initialize the space for the solution */
  for( i = 0; i < ctx->nno; i++ )
      ctx->V[1][i] = 0.0;

  count = 0;

  do {
    residual = multi_grid( ctx->E, ctx->V, ctx->RR, ctx->acc, ctx->level );
    valid = (residual < ctx->acc) ? 1 : 0;
    count++;
  } while ( (!valid) && (count < ctx->max_vel_iterations) );
  ctx->status = residual;

  VecGetOwnershipRange( y, &low, &high );
  for( i = 0; i < ctx->nno; i++ ) {
    ierr = VecSetValue( y, i+low, ctx->V[1][i], INSERT_VALUES ); 
    CHKERRQ( ierr );
  }
  VecAssemblyBegin(y);
  VecAssemblyEnd(y);

  PetscFunctionReturn(0);
}

PetscErrorCode MatShellMult_del2_u( Mat K, Vec U, Vec KU )
{
  // K is (neq+1) x (neq+1)
  // U is (neq+1)
  // KU is (neq+1)
  int i, j, neq;
  PetscErrorCode ierr;
  PetscScalar *KData, *KUData;
  struct MatMultShell *ctx;
  MatShellGetContext( K, (void **)&ctx );
  neq = ctx->iSize; // ctx->iSize SHOULD be the same as ctx->oSize
#if 0
  int low, high;

  VecAssemblyBegin(U);
  VecAssemblyEnd(U);
  VecGetOwnershipRange( U, &low, &high );
  for( i = 1; i <= ctx->E->sphere.caps_per_proc; i++ ) {
    for( j = 0; j <= neq; j++ ) {
      PetscInt ix[] = {j+low};
      ierr = VecGetValues( U, 1, ix, &ctx->iData[i][j] );
      CHKERRQ( ierr );
    }
  }
#endif
  ierr = VecGetArray(K, &KData); CHKERRQ(ierr);
  for(j = 0; j <=neq; j++)
    ctx->iData[1][j] = KData[j];
  ierr = VecRestoreArray(K, &KData); CHKERRQ(ierr);
  // actual CitcomS operation
  assemble_del2_u( ctx->E, ctx->iData, ctx->oData, ctx->level, 1 );
  ierr = VecGetArray(KU, &KUData); CHKERRQ(ierr);
  for(j = 0; j <= neq; j++)
    KUData[j] = ctx->oData[1][j];
  ierr = VecRestoreArray(KU, &KUData); CHKERRQ(ierr);
#if 0 
  VecGetOwnershipRange( KU, &low, &high );
  for( i=1; i <= ctx->E->sphere.caps_per_proc; ++i ) {
    for( j = 0; j <= neq; j++ ) {
      ierr = VecSetValue( KU, j+low, ctx->oData[i][j], INSERT_VALUES );
      CHKERRQ( ierr );
    }
  }
  VecAssemblyBegin( KU );
  VecAssemblyEnd( KU );
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode MatShellMult_grad_p( Mat G, Vec P, Vec GP )
{
  // G is (neq+1) x (nel)
  // P is (nel) x 1
  // GP is (neq+1) x 1 of which first neq entries (0:neq-1) are actual values
  int i, j, neq, nel;
  PetscErrorCode ierr;
  PetscScalar *PData, *GPData;
  struct MatMultShell *ctx;
  MatShellGetContext( G, (void **)&ctx );
  nel = ctx->iSize;
  neq = ctx->oSize;
#if 0  
  int low, high;

  VecGetOwnershipRange( P, &low, &high );

  VecAssemblyBegin( P );
  VecAssemblyEnd( P );
  for( i = 1; i <= ctx->E->sphere.caps_per_proc; i++ ) {
    for( j = 0; j < nel; j++ ) {
      PetscInt ix[] = {j+low};
      ierr = VecGetValues( P, 1, ix, &ctx->iData[i][j+1] );
      CHKERRQ( ierr );
    }
  }
#endif
  ierr = VecGetArray(P, &PData); CHKERRQ(ierr);
  for(j = 0; j < nel; j++)
    ctx->iData[1][j+1] = PData[j];
  ierr = VecRestoreArray(P, &PData); CHKERRQ(ierr);


  // actual CitcomS operation
  assemble_grad_p( ctx->E, ctx->iData, ctx->oData, ctx->level );

  ierr = VecGetArray(GP, &GPData); CHKERRQ(ierr);
  for(j = 0; j < neq; j++)
    GPData[j] = ctx->oData[1][j];
  ierr = VecRestoreArray(GP, &GPData); CHKERRQ(ierr);
#if 0
  VecGetOwnershipRange( GP, &low, &high );
  for( i = 1; i <= ctx->E->sphere.caps_per_proc; i++ ) {
    for( j = 0; j < neq; j++ ) {
      ierr = VecSetValue( GP, j+low, ctx->oData[i][j], INSERT_VALUES );
      CHKERRQ( ierr );
    }
  }
  VecAssemblyBegin( GP );
  VecAssemblyEnd( GP );
#endif
  PetscFunctionReturn(0);
}

PetscErrorCode MatShellMult_div_u( Mat D, Vec U, Vec DU )
{
  // D is nel x (neq+1)
  // U is (neq+1) x 1, of which first neq values (0:neq-1) are actual values
  // DU is nel x 1
  int i, j, neq, nel;
  PetscErrorCode ierr;
  struct MatMultShell *ctx;
  MatShellGetContext( D, (void **)&ctx );
  neq = ctx->iSize;
  nel = ctx->oSize;

  int low, high;
  VecGetOwnershipRange( U, &low, &high );

  VecAssemblyBegin( U );
  VecAssemblyEnd( U );
  for( i = 1; i <= ctx->E->sphere.caps_per_proc; i++ ) {
    for( j = 0; j < neq; j++ ) {
      PetscInt ix[] = {j+low};
      ierr = VecGetValues( U, 1, ix, &ctx->iData[i][j] );
      CHKERRQ( ierr );
    }
  }

  // actual CitcomS operation
  assemble_div_u( ctx->E, ctx->iData, ctx->oData, ctx->level );

  VecGetOwnershipRange( DU, &low, &high );
  for( i = 1; i <= ctx->E->sphere.caps_per_proc; i++ ) {
    for( j = 0; j < nel; j++ ) {
      ierr = VecSetValue( DU, j+low, ctx->oData[i][j+1], INSERT_VALUES );
      CHKERRQ( ierr );
    }
  }
  VecAssemblyBegin( DU );
  VecAssemblyEnd( DU );

  PetscFunctionReturn(0);
}

PetscErrorCode MatShellMult_div_rho_u( Mat DC, Vec U, Vec DU )
{
  // DC is nel x (neq+1)
  // U is (neq+1) x 1, of which first neq values (0:neq-1) are actual values
  // DU is nel x 1
  int i, j, neq, nel;
  PetscErrorCode ierr;
  struct MatMultShell *ctx;
  MatShellGetContext( DC, (void **)&ctx );
  neq = ctx->iSize;
  nel = ctx->oSize;

  int low, high;
  VecGetOwnershipRange( U, &low, &high );

  VecAssemblyBegin( U );
  VecAssemblyEnd( U );
  for( i = 1; i <= ctx->E->sphere.caps_per_proc; i++ ) {
    for( j = 0; j < neq; j++ ) {
      PetscInt ix[] = {j+low};
      ierr = VecGetValues( U, 1, ix, &ctx->iData[i][j] );
      CHKERRQ( ierr );
    }
  }

  // actual CitcomS operation
  assemble_div_rho_u( ctx->E, ctx->iData, ctx->oData, ctx->level );

  VecGetOwnershipRange( DU, &low, &high );
  for( i = 1; i <= ctx->E->sphere.caps_per_proc; i++ ) {
    for( j = 0; j < nel; j++ ) {
      ierr = VecSetValue( DU, j+low, ctx->oData[i][j+1], INSERT_VALUES );
      CHKERRQ( ierr );
    }
  }
  VecAssemblyBegin( DU );
  VecAssemblyEnd( DU );

  PetscFunctionReturn(0);
}
