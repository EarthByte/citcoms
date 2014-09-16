/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */


/*   Functions which solve for the velocity and pressure fields using Uzawa-type iteration loop.  */

#include <math.h>
#include <string.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"
#include <stdlib.h>

#ifdef USE_PETSC
#include "petsc_citcoms.h"
#endif

void myerror(struct All_variables *,char *);

static void solve_Ahat_p_fhat(struct All_variables *E,
                              double **V, double **P, double **F,
                              double imp, int *steps_max);
static void solve_Ahat_p_fhat_CG(struct All_variables *E,
                                 double **V, double **P, double **F,
                                 double imp, int *steps_max);
static void solve_Ahat_p_fhat_BiCG(struct All_variables *E,
                                    double **V, double **P, double **F,
                                    double imp, int *steps_max);
static void solve_Ahat_p_fhat_iterCG(struct All_variables *E,
                                      double **V, double **P, double **F,
                                      double imp, int *steps_max);

#ifdef USE_PETSC
static PetscErrorCode solve_Ahat_p_fhat_PETSc_Schur(struct All_variables *E,
    double **V, double **P, double **F, double imp, int *steps_max);

static PetscErrorCode solve_Ahat_p_fhat_CG_PETSc(struct All_variables *E, 
    double **V, double **P, double **F, double imp, int *steps_max);

static PetscErrorCode solve_Ahat_p_fhat_BiCG_PETSc(struct All_variables *E,
    double **V, double **P, double **F, double imp, int *steps_max);
#endif

static void initial_vel_residual(struct All_variables *E,
                                 double **V, double **P, double **F,
                                 double imp);


/* Master loop for pressure and (hence) velocity field */

void solve_constrained_flow_iterative(E)
     struct All_variables *E;

{
    void v_from_vector();
    void v_from_vector_pseudo_surf();
    void p_to_nodes();

    int cycles;

    cycles=E->control.p_iterations;

    /* Solve for velocity and pressure, correct for bc's */

    solve_Ahat_p_fhat(E,E->U,E->P,E->F,E->control.accuracy,&cycles);

    if(E->control.pseudo_free_surf)
        v_from_vector_pseudo_surf(E);
    else
        v_from_vector(E);

    p_to_nodes(E,E->P,E->NP,E->mesh.levmax);

    return;
}
/* ========================================================================= */

static double momentum_eqn_residual(struct All_variables *E,
                                    double **V, double **P, double **F)
{
    /* Compute the norm of (F - grad(P) - K*V)
     * This norm is ~= E->monitor.momentum_residual */
    void assemble_del2_u();
    void assemble_grad_p();
    void strip_bcs_from_residual();
    double global_v_norm2();

    int i, m;
    double *r1[NCS], *r2[NCS];
    double res;
    const int lev = E->mesh.levmax;
    const int neq = E->lmesh.neq;

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        r1[CPPR] = malloc((neq+1)*sizeof(double));
        r2[CPPR] = malloc((neq+1)*sizeof(double));
    }

    /* r2 = F - grad(P) - K*V */
    assemble_grad_p(E, P, E->u1, lev);
    assemble_del2_u(E, V, r1, lev, 1);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            r2[CPPR][i] = F[CPPR][i] - E->u1[CPPR][i] - r1[CPPR][i];

    strip_bcs_from_residual(E, r2, lev);

    res = sqrt(global_v_norm2(E, r2));

    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        free(r1[CPPR]);
        free(r2[CPPR]);
    }
    return(res);
}


static void print_convergence_progress(struct All_variables *E,
                                       int count, double time0,
                                       double v_norm, double p_norm,
                                       double dv, double dp,
                                       double div)
{
    double CPU_time0(), t;
    t = CPU_time0() - time0;

    fprintf(E->fp, "(%03d) %5.1f s v=%e p=%e "
            "div/v=%.2e dv/v=%.2e dp/p=%.2e step %d\n",
            count, t, v_norm, p_norm, div, dv, dp,
            E->monitor.solution_cycles);
    fprintf(stderr, "(%03d) %5.1f s v=%e p=%e "
            "div/v=%.2e dv/v=%.2e dp/p=%.2e step %d\n",
            count, t, v_norm, p_norm, div, dv, dp,
            E->monitor.solution_cycles);
    fflush(stderr);
    return;
}

static void print_convergence_progress_schur(struct All_variables *E, 
                                             int inner, 
                                             int outer, 
                                             double time0,
                                             double vnorm, double pnorm)
{
  double t = CPU_time0() - time0;

  fprintf(E->fp, "#inner=%03d #outer=%03d time=%5.1f s v=%e p=%e step %d\n", 
      inner, outer, t, vnorm, pnorm,E->monitor.solution_cycles);
  fprintf(stderr, "#inner=%03d #outer=%03d time=%5.1f s v=%e p=%e step %d\n", 
      inner, outer, t, vnorm, pnorm,E->monitor.solution_cycles);
  fflush(stderr);
}

static int keep_iterating(struct All_variables *E,
                          double acc, int converging)
{
    const int required_converging_loops = 2;

    if(E->control.check_continuity_convergence)
        return (E->monitor.incompressibility > acc) ||
	    (converging < required_converging_loops);
    else
        return (E->monitor.incompressibility > acc) &&
	    (converging < required_converging_loops);
}

static void solve_Ahat_p_fhat(struct All_variables *E,
                               double **V, double **P, double **F,
                               double imp, int *steps_max)
{
#ifdef USE_PETSC
  if(E->control.use_petsc)
  {
    if(E->control.petsc_schur) // use Schur complement reduction
    {
      //myerror(E, "Schur complement reduction approach is not implemented yet");
      solve_Ahat_p_fhat_PETSc_Schur(E, V, P, F, imp, steps_max);
    }
    else                       // use the Uzawa algorithm
    {
      if(E->control.inv_gruneisen == 0)
        solve_Ahat_p_fhat_CG_PETSc(E, V, P, F, imp, steps_max);
      else
      {
        if(strcmp(E->control.uzawa, "cg") == 0)
          solve_Ahat_p_fhat_iterCG(E, V, P, F, imp, steps_max);
        else if(strcmp(E->control.uzawa, "bicg") == 0)
          solve_Ahat_p_fhat_BiCG_PETSc(E, V, P, F, imp, steps_max);
        else
          myerror(E, "Error: unknown Uzawa iteration\n");
      }
    }
  }
  else                         // the original non-PETSc CitcomS code
  {
#endif 
    if(E->control.inv_gruneisen == 0)
        solve_Ahat_p_fhat_CG(E, V, P, F, imp, steps_max);
    else {
        if(strcmp(E->control.uzawa, "cg") == 0)
            solve_Ahat_p_fhat_iterCG(E, V, P, F, imp, steps_max);
        else if(strcmp(E->control.uzawa, "bicg") == 0)
            solve_Ahat_p_fhat_BiCG(E, V, P, F, imp, steps_max);
        else
            myerror(E, "Error: unknown Uzawa iteration\n");
    }
#ifdef USE_PETSC
  }
#endif
}

#ifdef USE_PETSC
static PetscErrorCode solve_Ahat_p_fhat_PETSc_Schur(struct All_variables *E,
  double **V, double **P, double **F, double imp, int *steps_max)
{
  int i, npno, neq, lev, nel, N, count;
  PetscErrorCode ierr;

  Mat S;
  Vec FF, PVec, VVec, fhat, fstar, t1;
  PetscReal *P_data, *V_data;
  KSP inner_ksp, S_ksp;
  PC inner_pc, S_pc;
  PetscInt inner1, inner2, outer;

  nel = E->lmesh.nel;
  npno = E->lmesh.npno;
  neq = E->lmesh.neq;
  lev = E->mesh.levmax;

  count = *steps_max;
  double time0 = CPU_time0();

  // Create the force Vec
  ierr = VecCreateMPI( PETSC_COMM_WORLD, neq, PETSC_DECIDE, &FF ); 
  CHKERRQ(ierr);
  double *FF_data;
  ierr = VecGetArray( FF, &FF_data ); CHKERRQ( ierr );
  for( i = 0; i < neq; i++ )
    FF_data[i] = F[1][i];
  ierr = VecRestoreArray( FF, &FF_data ); CHKERRQ( ierr );

  ierr = VecDuplicate(FF, &t1); CHKERRQ(ierr);

  /* create the pressure vector and initialize it to zero */
  ierr = VecCreateMPI( PETSC_COMM_WORLD, nel, PETSC_DECIDE, &PVec ); 
  CHKERRQ(ierr);
  ierr = VecSet( PVec, 0.0 ); CHKERRQ( ierr );

  /* create the velocity vector and copy contents of V into it */
  ierr = VecCreateMPI( PETSC_COMM_WORLD, neq, PETSC_DECIDE, &VVec ); 
  CHKERRQ(ierr);
  ierr = VecGetArray( VVec, &V_data ); CHKERRQ( ierr );
  for( i = 0; i < neq; i++ )
    V_data[i] = V[1][i];
  ierr = VecRestoreArray( VVec, &V_data ); CHKERRQ( ierr );

  /*----------------------------------*/
  /* Define a Schur complement matrix */
  /*----------------------------------*/
  ierr = MatCreateSchurComplement(E->K,E->K,E->G,E->D,PETSC_NULL, &S); 
  CHKERRQ(ierr);
  ierr = MatSchurComplementGetKSP(S, &inner_ksp); CHKERRQ(ierr);
  // the following commented out line doesn't work with shell matrices
  //ierr = KSPSetType(inner_ksp, "preonly"); CHKERRQ(ierr);
  ierr = KSPGetPC(inner_ksp, &inner_pc); CHKERRQ(ierr);
  // the following commented out line doesn't work with shell matrices
  //ierr = PCSetType(inner_pc, "lu"); CHKERRQ(ierr);

  /*-------------------------------------------------*/
  /* Build the RHS of the Schur Complement Reduction */
  /*-------------------------------------------------*/
  ierr = MatGetVecs(S, PETSC_NULL, &fhat); CHKERRQ(ierr);
  ierr = KSPSolve(inner_ksp, FF, t1); CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(inner_ksp, &inner1); CHKERRQ(ierr);
  ierr = MatMult(E->D, t1, fhat); CHKERRQ(ierr);

  /*-------------------------------------------*/
  /* Build the solver for the Schur complement */
  /*-------------------------------------------*/
  ierr = KSPCreate(PETSC_COMM_WORLD, &S_ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(S_ksp, S, S); CHKERRQ(ierr);
  //ierr = KSPGetPC(S_ksp, &S_pc); CHKERRQ(ierr);
  ierr = KSPSetType(S_ksp, "cg"); CHKERRQ(ierr);
  //ierr = PCSetType(S_pc, "none"); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(S_ksp, PETSC_TRUE); CHKERRQ(ierr);

  /*--------------------*/
  /* Solve for pressure */
  /*--------------------*/
  ierr = KSPSolve(S_ksp, fhat, PVec); CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(inner_ksp, &outer); CHKERRQ(ierr);

  /*--------------------*/
  /* Solve for velocity */
  /*--------------------*/
  ierr = MatGetVecs(E->K, PETSC_NULL, &fstar); CHKERRQ(ierr);
  ierr = MatMult(E->G, PVec, fstar); CHKERRQ(ierr);
  ierr = VecAYPX(fstar, 1.0, FF); CHKERRQ(ierr);
  ierr = KSPSetInitialGuessNonzero(inner_ksp, PETSC_TRUE); CHKERRQ(ierr);
  ierr = KSPSolve(inner_ksp, fstar, VVec); CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(inner_ksp, &inner2); CHKERRQ(ierr);

  /*-----------------------------------------------*/
  /* copy the values of VVec and PVec into V and P */
  /*-----------------------------------------------*/
  ierr = VecGetArray( VVec, &V_data ); CHKERRQ( ierr );
  for( i = 0; i < neq; i++ )
    V[1][i] = V_data[i];
  ierr = VecRestoreArray( VVec, &V_data ); CHKERRQ( ierr );
  //velocities_conform_bcs( E, VVec, lev );
  
  ierr = VecGetArray( PVec, &P_data ); CHKERRQ( ierr );
  for( i = 0; i < nel; i++ )
    P[1][i] = P_data[i]; 
  ierr = VecRestoreArray( PVec, &P_data ); CHKERRQ( ierr );

	if((E->sphere.caps == 12) && (E->control.inner_remove_rigid_rotation)){
	  /* allow for removal of net rotation at each iterative step (expensive) */
	  if(E->control.pseudo_free_surf) /* move from U to V */
	    v_from_vector_pseudo_surf(E);
	  else
	    v_from_vector(E);
	  remove_rigid_rot(E);	/* correct V */
	  assign_v_to_vector(E); /* assign V to U */
	}


  /*---------------------------------------------------------------------------------------*/
  /* output v_norm and p_norm values; dvelocity and dpressure are currently not meaningful */
  /*---------------------------------------------------------------------------------------*/
  E->monitor.vdotv = global_v_norm2(E, V);
  E->monitor.pdotp = global_p_norm2(E, P);
  double v_norm = sqrt(E->monitor.vdotv);
  double p_norm = sqrt(E->monitor.pdotp);
  double dvelocity = 1.0;
  double dpressure = 1.0;
  int converging = 0;
  if (E->control.print_convergence && E->parallel.me==0)  {
    print_convergence_progress_schur(E, (inner1+inner2), outer, time0, v_norm, p_norm);
  }

  PetscFunctionReturn(0);
}
#endif

/* Solve incompressible Stokes flow using
 * conjugate gradient (CG) iterations
 */

static void solve_Ahat_p_fhat_CG(struct All_variables *E,
                                 double **V, double **P, double **FF,
                                 double imp, int *steps_max)
{
    int m, j, count, valid, lev, npno, neq;

    double *r1[NCS], *r2[NCS], *z1[NCS], *s1[NCS], *s2[NCS], *cu[NCS];
    double *F[NCS];
    double *shuffle[NCS];
    double alpha, delta, r0dotz0, r1dotz1;
    double v_res;
    double inner_imp;
    double global_pdot();
    double global_v_norm2(), global_p_norm2(), global_div_norm2();

    double time0, CPU_time0();
    double v_norm, p_norm;
    double dvelocity, dpressure;
    int converging;
    void assemble_c_u();
    void assemble_div_u();
    void assemble_del2_u();
    void assemble_grad_p();
    void strip_bcs_from_residual();
    int  solve_del2_u();
    void parallel_process_termination();
    void v_from_vector();
    void v_from_vector_pseudo_surf();
    void assign_v_to_vector();

    inner_imp = imp * E->control.inner_accuracy_scale; /* allow for different innner loop accuracy */

    npno = E->lmesh.npno;
    neq = E->lmesh.neq;
    lev = E->mesh.levmax;

    for (m=1; m<=E->sphere.caps_per_proc; m++)   {
        F[CPPR] = (double *)malloc(neq*sizeof(double));
        r1[CPPR] = (double *)malloc(npno*sizeof(double));
        r2[CPPR] = (double *)malloc(npno*sizeof(double));
        z1[CPPR] = (double *)malloc(npno*sizeof(double));
        s1[CPPR] = (double *)malloc(npno*sizeof(double));
        s2[CPPR] = (double *)malloc(npno*sizeof(double));
        cu[CPPR] = (double *)malloc(npno*sizeof(double));
    }

    time0 = CPU_time0();
    count = 0;
    v_res = E->monitor.fdotf;

    /* copy the original force vector since we need to keep it intact
       between iterations */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(j=0;j<neq;j++)
            F[CPPR][j] = FF[CPPR][j];


    /* calculate the contribution of compressibility in the continuity eqn */
    if(E->control.inv_gruneisen != 0) {
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<npno;j++)
                cu[CPPR][j] = 0.0;

        assemble_c_u(E, V, cu, lev);
    }


    /* calculate the initial velocity residual */
    /* In the compressible case, the initial guess of P might be bad.
     * Do not correct V with it. */
    if(E->control.inv_gruneisen == 0)
        initial_vel_residual(E, V, P, F, inner_imp*v_res);


    /* initial residual r1 = div(V) */
    assemble_div_u(E, V, r1, lev);


    /* add the contribution of compressibility to the initial residual */
    if(E->control.inv_gruneisen != 0)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<npno;j++) {
                r1[CPPR][j] += cu[CPPR][j];
            }

    E->monitor.vdotv = global_v_norm2(E, V);
    E->monitor.incompressibility = sqrt(global_div_norm2(E, r1)
                                        / (1e-32 + E->monitor.vdotv));

    v_norm = sqrt(E->monitor.vdotv);
    p_norm = sqrt(E->monitor.pdotp);
    dvelocity = 1.0;
    dpressure = 1.0;
    converging = 0;

    if (E->control.print_convergence && E->parallel.me==0)  {
        print_convergence_progress(E, count, time0,
                                   v_norm, p_norm,
                                   dvelocity, dpressure,
                                   E->monitor.incompressibility);
    }

  
    r0dotz0 = 0;

    while( (count < *steps_max) && keep_iterating(E, imp, converging) ) {
        /* require two consecutive converging iterations to quit the while-loop */

        /* preconditioner BPI ~= inv(K), z1 = BPI*r1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<npno; j++)
                z1[CPPR][j] = E->BPI[lev][CPPR][j+1] * r1[CPPR][j]; /* E->BPI[lev][m][j] when it is made 0-based */


        /* r1dotz1 = <r1, z1> */
        r1dotz1 = global_pdot(E, r1, z1, lev);
        assert(r1dotz1 != 0.0  /* Division by zero in head of incompressibility iteration */);

        /* update search direction */
        if(count == 0)
            for (m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=0; j<npno; j++)
                    s2[CPPR][j] = z1[CPPR][j];
        else {
            /* s2 = z1 + s1 * <r1,z1>/<r0,z0> */
            delta = r1dotz1 / r0dotz0;
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=0; j<npno; j++)
                    s2[CPPR][j] = z1[CPPR][j] + delta * s1[CPPR][j];
        }

        /* solve K*u1 = grad(s2) for u1 */
        assemble_grad_p(E, s2, F, lev);
        valid = solve_del2_u(E, E->u1, F, inner_imp*v_res, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 1\n", stderr);
            fputs("Warning: solver not converging! 1\n", E->fp);
        }
        //strip_bcs_from_residual(E, E->u1, lev);


        /* F = div(u1) */
        assemble_div_u(E, E->u1, F, lev);


        /* alpha = <r1, z1> / <s2, F> */
        alpha = r1dotz1 / global_pdot(E, s2, F, lev);


        /* r2 = r1 - alpha * div(u1) */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<npno; j++)
                r2[CPPR][j] = r1[CPPR][j] - alpha * F[CPPR][j];


        /* P = P + alpha * s2 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<npno; j++)
                P[CPPR][j] += alpha * s2[CPPR][j];

        /* V = V - alpha * u1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                V[CPPR][j] -= alpha * E->u1[CPPR][j];


        /* compute velocity and incompressibility residual */
        E->monitor.vdotv = global_v_norm2(E, V);
        E->monitor.pdotp = global_p_norm2(E, P);
        v_norm = sqrt(E->monitor.vdotv);
        p_norm = sqrt(E->monitor.pdotp);
        dvelocity = alpha * sqrt(global_v_norm2(E, E->u1) / (1e-32 + E->monitor.vdotv));
        dpressure = alpha * sqrt(global_p_norm2(E, s2) / (1e-32 + E->monitor.pdotp));

       

        assemble_div_u(E, V, z1, lev);
        if(E->control.inv_gruneisen != 0)
            for(m=1;m<=E->sphere.caps_per_proc;m++)
                for(j=0;j<npno;j++) {
                    z1[CPPR][j] += cu[CPPR][j];
            }
        E->monitor.incompressibility = sqrt(global_div_norm2(E, z1)
                                            / (1e-32 + E->monitor.vdotv));

        count++;


        if (E->control.print_convergence && E->parallel.me==0)  {
            print_convergence_progress(E, count, time0,
                                       v_norm, p_norm,
                                       dvelocity, dpressure,
                                       E->monitor.incompressibility);
        }

	if(!valid){
            /* reset consecutive converging iterations */
            converging = 0;
	}else{
            /* how many consecutive converging iterations? */
            if(E->control.check_pressure_convergence) {
                /* check dv and dp */
                if(dvelocity < imp && dpressure < imp)
                    converging++;
                else
                    converging = 0;
            }else{
                /* check dv only */
                if(dvelocity < imp)
                    converging++;
                else
                    converging = 0;
            }
	  
	}

        /* shift array pointers */
        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            shuffle[CPPR] = s1[CPPR];
            s1[CPPR] = s2[CPPR];
            s2[CPPR] = shuffle[CPPR];

            shuffle[CPPR] = r1[CPPR];
            r1[CPPR] = r2[CPPR];
            r2[CPPR] = shuffle[CPPR];
        }

        /* shift <r0, z0> = <r1, z1> */
        r0dotz0 = r1dotz1;
	if((E->sphere.caps == 12) && (E->control.inner_remove_rigid_rotation)){
	  /* allow for removal of net rotation at each iterative step
	     (expensive) */
	  if(E->control.pseudo_free_surf) /* move from U to V */
	    v_from_vector_pseudo_surf(E);
	  else
	    v_from_vector(E);
	  remove_rigid_rot(E);	/* correct V */
	  assign_v_to_vector(E); /* assign V to U */
	}
    } /* end loop for conjugate gradient */

    assemble_div_u(E, V, z1, lev);
    if(E->control.inv_gruneisen != 0)
        for(m=1;m<=E->sphere.caps_per_proc;m++)
            for(j=0;j<npno;j++) {
                z1[CPPR][j] += cu[CPPR][j];
            }


    for(m=1; m<=E->sphere.caps_per_proc; m++) {
        free((void *) F[CPPR]);
        free((void *) r1[CPPR]);
        free((void *) r2[CPPR]);
        free((void *) z1[CPPR]);
        free((void *) s1[CPPR]);
        free((void *) s2[CPPR]);
        free((void *) cu[CPPR]);
    }

    *steps_max=count;

    return;
}

#ifdef USE_PETSC
/*
 * Implementation of the Conjugate Gradient Uzawa algorithm using PETSc
 * Vec, Mat and KSPSolve
 */
static PetscErrorCode solve_Ahat_p_fhat_CG_PETSc( struct All_variables *E,
				     double **V, double **P, double **F,
				     double imp, int *steps_max )
{
  PetscErrorCode ierr;
  Vec V_k, P_k, s_1, s_2, r_1, r_2, z_1, BPI, FF, Gsk, u_k, Duk, cu;
  PetscReal alpha, delta, r_1_norm, r1dotz1, r0dotz0, s_2_dot_F;
  int i,j,count,m;
  double time0, CPU_time0();
  double v_norm, p_norm, dvelocity, dpressure;
  double inner_imp, v_res;

  int lev = E->mesh.levmax;
  int npno = E->lmesh.npno;
  int neq = E->lmesh.neq;
  int nel = E->lmesh.nel;

  inner_imp = imp * E->control.inner_accuracy_scale; 
  v_res = E->monitor.fdotf;

  time0 = CPU_time0();
  count = 0;

  // Create the force Vec
  ierr = VecCreateMPI( PETSC_COMM_WORLD, neq, PETSC_DECIDE, &FF ); 
  CHKERRQ(ierr);
  double *F_tmp;
  ierr = VecGetArray( FF, &F_tmp ); CHKERRQ( ierr );
  for( m = 1; m <= E->sphere.caps_per_proc; m++ ) {
    for( i = 0; i < neq; i++ )
      F_tmp[i] = F[CPPR][i];
  }
  ierr = VecRestoreArray( FF, &F_tmp ); CHKERRQ( ierr );

  // create the pressure vector and initialize it to zero
  ierr = VecCreateMPI( PETSC_COMM_WORLD, nel, PETSC_DECIDE, &P_k ); 
  CHKERRQ(ierr);
  ierr = VecSet( P_k, 0.0 ); CHKERRQ( ierr );

  // create the velocity vector
  ierr = VecCreateMPI( PETSC_COMM_WORLD, neq, PETSC_DECIDE, &V_k ); 
  CHKERRQ(ierr);

  // Copy the contents of V into V_k
  PetscScalar *V_k_tmp;
  ierr = VecGetArray( V_k, &V_k_tmp ); CHKERRQ( ierr );
  for( m = 1; m <= E->sphere.caps_per_proc; m++ ) {
    for( i = 0; i < neq; i++ )
      V_k_tmp[i] = V[CPPR][i];
  }
  ierr = VecRestoreArray( V_k, &V_k_tmp ); CHKERRQ( ierr );

  // PETSc bookkeeping --- create various temporary Vec objects with
  // the appropriate sizes, including the PETSc vec version of E->BPI
  // preconditioner
  ierr = VecCreateMPI( PETSC_COMM_WORLD, npno, PETSC_DECIDE, &r_1 ); 
  CHKERRQ(ierr);
  ierr = VecDuplicate( V_k, &Gsk ); CHKERRQ( ierr );
  ierr = VecDuplicate( V_k, &u_k ); CHKERRQ( ierr );
  ierr = VecDuplicate( r_1, &s_1 ); CHKERRQ( ierr );
  ierr = VecDuplicate( r_1, &s_2 ); CHKERRQ( ierr );
  ierr = VecDuplicate( r_1, &r_2 ); CHKERRQ( ierr );
  ierr = VecDuplicate( r_1, &z_1 ); CHKERRQ( ierr );
  ierr = VecDuplicate( r_1, &cu ); CHKERRQ( ierr );
  ierr = VecDuplicate( r_1, &Duk ); CHKERRQ( ierr );
  ierr = VecDuplicate( r_1, &BPI ); CHKERRQ( ierr );
  PetscReal *bpi;
  ierr = VecGetArray( BPI, &bpi ); CHKERRQ( ierr );
  for( m = 1; m <= E->sphere.caps_per_proc; m++ )
    for( j = 0; j < npno; j++ )
      bpi[j] = E->BPI[lev][CPPR][j+1];
  ierr = VecRestoreArray( BPI, &bpi ); CHKERRQ( ierr );

  /* calculate the contribution of compressibility in the continuity eqn */
  if( E->control.inv_gruneisen != 0 ) {
    ierr = VecSet( cu, 0.0 ); CHKERRQ( ierr );
    assemble_c_u_PETSc( E, V_k, cu, lev );
  }

  /* calculate the initial velocity residual */
  /* In the compressible case, the initial guess of P might be bad.
   * Do not correct V with it. */
  if( E->control.inv_gruneisen == 0 ) {
    initial_vel_residual_PETSc(E, V_k, P_k, FF, inner_imp*v_res);
  }

  /* initial residual r1 = div(V) */
  ierr = MatMult( E->D, V_k, r_1 ); CHKERRQ( ierr );

  /* add the contribution of compressibility to the initial residual */
  if( E->control.inv_gruneisen != 0 ) {
    // r_1 += cu
    ierr = VecAXPY( r_1, 1.0, cu ); CHKERRQ( ierr );
  }

  E->monitor.vdotv = global_v_norm2_PETSc( E, V_k );
  E->monitor.incompressibility = sqrt( global_div_norm2_PETSc( E, r_1 ) / (1e-32 + E->monitor.vdotv) );

  v_norm = sqrt( E->monitor.vdotv );
  p_norm = sqrt( E->monitor.pdotp );
  dvelocity = 1.0;
  dpressure = 1.0;

  if( E->control.print_convergence && E->parallel.me == 0 ) {
    print_convergence_progress( E, count, time0, 
                                v_norm, p_norm,
                                dvelocity, dpressure,
                                E->monitor.incompressibility );
  }

  
  r0dotz0 = 0;
  ierr = VecNorm( r_1, NORM_2, &r_1_norm ); CHKERRQ( ierr );

  while( (r_1_norm > E->control.petsc_uzawa_tol) && (count < *steps_max) )
    {
      
      /* preconditioner BPI ~= inv(K), z1 = BPI*r1 */
      ierr = VecPointwiseMult( z_1, BPI, r_1 ); CHKERRQ( ierr );

      /* r1dotz1 = <r1, z1> */
      ierr = VecDot( r_1, z_1, &r1dotz1 ); CHKERRQ( ierr );
      assert( r1dotz1 != 0.0  /* Division by zero in head of incompressibility
                                 iteration */);

      /* update search direction */
      if( count == 0 )
	    {
	      // s_2 = z_1
	      ierr = VecCopy( z_1, s_2 ); CHKERRQ( ierr ); // s2 = z1
	    }
      else
	    {
	      // s2 = z1 + s1 * <r1,z1>/<r0,z0>
	      delta = r1dotz1 / r0dotz0;
	      ierr = VecWAXPY( s_2, delta, s_1, z_1 ); CHKERRQ( ierr );
	    }

      // Solve K*u_k = grad(s_2) for u_k
      ierr = MatMult( E->G, s_2, Gsk ); CHKERRQ( ierr );
      ierr = KSPSolve( E->ksp, Gsk, u_k ); CHKERRQ( ierr );
      strip_bcs_from_residual_PETSc( E, u_k, lev );

      // Duk = D*u_k ( D*u_k is the same as div(u_k) )
      ierr = MatMult( E->D, u_k, Duk ); CHKERRQ( ierr );

      // alpha = <r1,z1> / <s2,F>
      ierr = VecDot( s_2, Duk, &s_2_dot_F ); CHKERRQ( ierr );
      alpha = r1dotz1 / s_2_dot_F;

      // r2 = r1 - alpha * div(u_k)
      ierr = VecWAXPY( r_2, -1.0*alpha, Duk, r_1 ); CHKERRQ( ierr );

      // P = P + alpha * s_2
      ierr = VecAXPY( P_k, 1.0*alpha, s_2 ); CHKERRQ( ierr );
      
      // V = V - alpha * u_1
      ierr = VecAXPY( V_k, -1.0*alpha, u_k ); CHKERRQ( ierr );
      //strip_bcs_from_residual_PETSc( E, V_k, E->mesh.levmax );

      /* compute velocity and incompressibility residual */
      E->monitor.vdotv = global_v_norm2_PETSc( E, V_k );
      E->monitor.pdotp = global_p_norm2_PETSc( E, P_k );
      v_norm = sqrt( E->monitor.vdotv );
      p_norm = sqrt( E->monitor.pdotp );
      dvelocity = alpha * sqrt( global_v_norm2_PETSc( E, u_k ) / (1e-32 + E->monitor.vdotv) );
      dpressure = alpha * sqrt( global_p_norm2_PETSc( E, s_2 ) / (1e-32 + E->monitor.pdotp) );

      // compute the updated value of z_1, z1 = div(V) 
      ierr = MatMult( E->D, V_k, z_1 ); CHKERRQ( ierr );
      if( E->control.inv_gruneisen != 0 )
	    {
        // z_1 += cu
        ierr = VecAXPY( z_1, 1.0, cu ); CHKERRQ( ierr );
      }

      E->monitor.incompressibility = sqrt( global_div_norm2_PETSc( E, z_1 ) / (1e-32 + E->monitor.vdotv) );

      count++;

      if( E->control.print_convergence && E->parallel.me == 0 ) {
        print_convergence_progress( E, count, time0,
                                    v_norm, p_norm,
                                    dvelocity, dpressure,
                                    E->monitor.incompressibility );
      }

      /* shift array pointers */
      ierr = VecSwap( s_2, s_1 ); CHKERRQ( ierr );
      ierr = VecSwap( r_2, r_1 ); CHKERRQ( ierr );

      /* shift <r0, z0> = <r1, z1> */
      r0dotz0 = r1dotz1;

      // recompute the norm
      ierr = VecNorm( r_1, NORM_2, &r_1_norm ); CHKERRQ( ierr );

    } /* end loop for conjugate gradient */

  // converged. now copy the converged values of V_k and P_k into V and P
  PetscReal *P_tmp, *V_tmp;
  ierr = VecGetArray( V_k, &V_tmp ); CHKERRQ( ierr );
  for( m = 1; m <= E->sphere.caps_per_proc; ++m ) {
    for( i = 0; i < neq; i++ )
      V[CPPR][i] = V_tmp[i];
  }
  ierr = VecRestoreArray( V_k, &V_tmp ); CHKERRQ( ierr );
  
  ierr = VecGetArray( P_k, &P_tmp ); CHKERRQ( ierr );
  for( m = 1; m <= E->sphere.caps_per_proc; ++m ) {
    for( i = 0; i < nel; i++ )
      P[CPPR][i+1] = P_tmp[i]; 
  }
  ierr = VecRestoreArray( P_k, &P_tmp ); CHKERRQ( ierr );

  // PETSc cleanup of all temporary Vec objects
  ierr = VecDestroy( &V_k ); CHKERRQ( ierr );
  ierr = VecDestroy( &P_k ); CHKERRQ( ierr );
  ierr = VecDestroy( &s_1 ); CHKERRQ( ierr );
  ierr = VecDestroy( &s_2 ); CHKERRQ( ierr );
  ierr = VecDestroy( &r_1 ); CHKERRQ( ierr );
  ierr = VecDestroy( &r_2 ); CHKERRQ( ierr );
  ierr = VecDestroy( &z_1 ); CHKERRQ( ierr );
  ierr = VecDestroy( &BPI ); CHKERRQ( ierr );
  ierr = VecDestroy( &FF );  CHKERRQ( ierr );
  ierr = VecDestroy( &Gsk ); CHKERRQ( ierr );
  ierr = VecDestroy( &u_k ); CHKERRQ( ierr );
  ierr = VecDestroy( &Duk ); CHKERRQ( ierr );

  *steps_max = count;

  PetscFunctionReturn(0);
}

/*
 * BiCGstab for compressible Stokes flow using PETSc Vec, Mat and KSPSolve
 */
static PetscErrorCode solve_Ahat_p_fhat_BiCG_PETSc( struct All_variables *E,
					  double **V, double **P, double **F,
					  double imp, int *steps_max )
{
  PetscErrorCode ierr;
  Vec FF, P0, p1, p2, pt, r1, r2, rt, s0, st, t0, u0, u1, V0, BPI, v0;

  PetscReal alpha, omega, beta;
  PetscReal r1_norm, r1dotrt, r0dotrt, rtdotV0, t0dots0, t0dott0;

  int i,j,k,m, count;

  double time0, CPU_time0();
  double v_norm, p_norm, inner_imp, v_res, dvelocity, dpressure;

  int lev = E->mesh.levmax;
  int npno = E->lmesh.npno;
  int neq = E->lmesh.neq;
  int nel = E->lmesh.nel;
  
  // Create the force Vec
  ierr = VecCreateMPI( PETSC_COMM_WORLD, neq, PETSC_DECIDE, &FF ); 
  CHKERRQ(ierr);
  double *F_tmp;
  ierr = VecGetArray( FF, &F_tmp ); CHKERRQ( ierr );
  for( m=1; m<=E->sphere.caps_per_proc; ++m ) {
    for( i = 0; i < neq; i++ )
      F_tmp[i] = F[CPPR][i];
  }
  ierr = VecRestoreArray( FF, &F_tmp ); CHKERRQ( ierr );

  inner_imp = imp * E->control.inner_accuracy_scale;
  time0 = CPU_time0();
  count = 0;
  v_res = E->monitor.fdotf;


  // create the pressure vector and initialize it to zero
  ierr = VecCreateMPI( PETSC_COMM_WORLD, nel, PETSC_DECIDE, &P0 ); 
  CHKERRQ(ierr);
  ierr = VecSet( P0, 0.0 ); CHKERRQ( ierr );

  // create the velocity vector
  ierr = VecCreateMPI( PETSC_COMM_WORLD, neq, PETSC_DECIDE, &V0 ); 
  CHKERRQ(ierr);

  // Copy the contents of V into V0
  PetscScalar *V0_tmp;
  ierr = VecGetArray( V0, &V0_tmp ); CHKERRQ( ierr );
  for( m = 1; m <= E->sphere.caps_per_proc; m++ ) {
    for( i = 0; i < neq; i++ )
      V0_tmp[i] = V[CPPR][i];
  }
  ierr = VecRestoreArray( V0, &V0_tmp ); CHKERRQ( ierr );


  ierr = VecDuplicate( V0, &u0 ); CHKERRQ( ierr );
  ierr = VecDuplicate( V0, &u1 ); CHKERRQ( ierr );

  /* calculate the initial velocity residual */
  initial_vel_residual_PETSc( E, V0, P0, FF, inner_imp*v_res );

  /* initial residual r1 = div(rho_ref*V) */
  ierr = VecCreateMPI( PETSC_COMM_WORLD, npno, PETSC_DECIDE, &r1 ); 

  CHKERRQ(ierr);
  ierr = MatMult( E->DC, V0, r1 ); CHKERRQ( ierr );

  E->monitor.vdotv = global_v_norm2_PETSc( E, V0 );
  E->monitor.incompressibility = sqrt( global_div_norm2_PETSc( E, r1 ) / (1e-32 + E->monitor.vdotv) );
  v_norm = sqrt( E->monitor.vdotv );
  p_norm = sqrt( E->monitor.pdotp );
  dvelocity = 1.0;
  dpressure = 1.0;

  if( E->control.print_convergence && E->parallel.me == 0 ) {
    print_convergence_progress( E, count, time0,
                                v_norm, p_norm,
                                dvelocity, dpressure,
                                E->monitor.incompressibility );
  }

  // create all the vectors for later use
  ierr = VecDuplicate( r1, &rt ); CHKERRQ( ierr );
  ierr = VecDuplicate( r1, &p1 ); CHKERRQ( ierr );
  ierr = VecDuplicate( r1, &p2 ); CHKERRQ( ierr );
  ierr = VecDuplicate( r1, &BPI ); CHKERRQ( ierr );
  ierr = VecDuplicate( r1, &pt ); CHKERRQ( ierr );
  ierr = VecDuplicate( r1, &s0 ); CHKERRQ( ierr );
  ierr = VecDuplicate( r1, &st ); CHKERRQ( ierr );
  ierr = VecDuplicate( r1, &t0 ); CHKERRQ( ierr );
  ierr = VecDuplicate( r1, &r2 ); CHKERRQ( ierr );
  ierr = VecDuplicate( r1, &v0 ); CHKERRQ( ierr );

  /* initial conjugate residual rt = r1 */
  ierr = VecCopy( r1, rt ); CHKERRQ( ierr );

  /* BPI ~= K inverse */
  PetscReal *bpi;
  ierr = VecGetArray( BPI, &bpi ); CHKERRQ( ierr );
  for( m = 1; m <= E->sphere.caps_per_proc; m++ )
    for( j = 0; j < npno; j++ )
      bpi[j] = E->BPI[lev][CPPR][j+1];
  ierr = VecRestoreArray( BPI, &bpi ); CHKERRQ( ierr );

  r0dotrt = alpha = omega = 0.0;

  ierr = VecNorm( r1, NORM_2, &r1_norm ); CHKERRQ( ierr );

  while( (r1_norm > E->control.petsc_uzawa_tol) && (count < *steps_max) )
  {

    /* r1dotrt = <r1, rt> */
    // r1dotrt = global_pdot( E, r1, rt, lev )
    ierr = VecDot( r1, rt, &r1dotrt ); CHKERRQ( ierr );

    if( r1dotrt == 0.0 ) {
      fprintf( E->fp, "BiCGstab method failed!!\n" );
      fprintf( stderr, "BiCGstab method failed!!\n" );
      parallel_process_termination();
    }

    /* update search direction */
    if( count == 0 ) {
      ierr = VecCopy( r1, p2 ); CHKERRQ( ierr );
    }
    else {
      /* p2 = r1 + <r1,rt>/<r0,rt> * alpha/omega * (p1 - omega*v0) */
      beta = (r1dotrt/r0dotrt)*(alpha/omega);
      ierr = VecAXPY( p1, -1.0*omega, v0); CHKERRQ( ierr );
      ierr = VecWAXPY( p2, beta, p1, r1 ); CHKERRQ( ierr );
    }

    /* preconditioner BPI ~= inv(K), pt = BPI*p2 */
    ierr = VecPointwiseMult( pt, BPI, p2 ); CHKERRQ( ierr );

    /* solve K*u0 = grad(pt) for u1 */
    ierr = MatMult( E->G, pt, FF ); CHKERRQ( ierr );
    ierr = KSPSolve( E->ksp, FF, u0 ); CHKERRQ( ierr );
    strip_bcs_from_residual_PETSc( E, u0, lev );

    /* v0 = div(rho_ref*u0) */
    ierr = MatMult( E->DC, u0, v0 ); CHKERRQ( ierr );
    
    /* alpha = r1dotrt / <rt, v0> */
    ierr = VecDot( rt, v0, &rtdotV0 ); CHKERRQ( ierr );
    alpha = r1dotrt / rtdotV0;

    /* s0 = r1 - alpha * v0 */
    ierr = VecWAXPY( s0, -1.0*alpha, v0, r1 ); CHKERRQ( ierr );

    /* preconditioner BPI ~= inv(K), st = BPI*s0 */
    ierr = VecPointwiseMult( st, BPI, s0 ); CHKERRQ( ierr );

    /* solve K*u1 = grad(st) for u1*/
    ierr = MatMult( E->G, st, FF ); CHKERRQ( ierr );
    ierr = KSPSolve( E->ksp, FF, u1 ); CHKERRQ( ierr );
    strip_bcs_from_residual_PETSc( E, u1, lev );

    /* t0 = div(rho_ref*u1) */
    ierr = MatMult( E->DC, u1, t0 ); CHKERRQ( ierr );

    /* omega = <t0, s0> / <t0, t0> */
    ierr = VecDot( t0, s0, &t0dots0 ); CHKERRQ( ierr );
    ierr = VecDot( t0, t0, &t0dott0 ); CHKERRQ( ierr );
    omega = t0dots0 / t0dott0;

    /* r2 = s0 - omega * t0 */
    ierr = VecWAXPY( r2, -1.0*omega, t0, s0 ); CHKERRQ( ierr );

    /* P = P + alpha * pt + omega * st */
    ierr = VecAXPBY( st, alpha, omega, pt ); CHKERRQ( ierr );
    ierr = VecAXPY( P0, 1.0, st ); CHKERRQ( ierr );

    /* V = V - alpha * u0 - omega * u1 */
    ierr = VecAXPBY( u1, alpha, omega, u0 ); CHKERRQ( ierr );
    ierr = VecAXPY( V0, -1.0, u1 ); CHKERRQ( ierr );

    /* compute velocity and incompressibility residual */
    E->monitor.vdotv = global_v_norm2_PETSc(E, V0);
    E->monitor.pdotp = global_p_norm2_PETSc(E, P0);
    v_norm = sqrt( E->monitor.vdotv );
    p_norm = sqrt( E->monitor.pdotp );
    dvelocity = sqrt( global_v_norm2_PETSc( E, u1 ) / (1e-32 + E->monitor.vdotv) );
    dpressure = sqrt( global_p_norm2_PETSc( E, st ) / (1e-32 + E->monitor.pdotp) );


    ierr = MatMult( E->DC, V0, t0 ); CHKERRQ( ierr );
    E->monitor.incompressibility = sqrt( global_div_norm2_PETSc( E, t0 )
                                        / (1e-32 + E->monitor.vdotv) );

    count++;

    if( E->control.print_convergence && E->parallel.me == 0 ) {
      print_convergence_progress( E, count, time0,
                                  v_norm, p_norm,
                                  dvelocity, dpressure,
                                  E->monitor.incompressibility );
    }

    /* shift array pointers */
    ierr = VecSwap( p1, p2 ); CHKERRQ( ierr );
    ierr = VecSwap( r1, r2 ); CHKERRQ( ierr );

    /* shift <r0, rt> = <r1, rt> */
    r0dotrt = r1dotrt;

    // recompute the norm of the residual
    ierr = VecNorm( r1, NORM_2, &r1_norm ); CHKERRQ( ierr );

  }

  // converged. now copy the converged values of V0 and P0 into V and P
  PetscReal *P_tmp, *V_tmp;
  ierr = VecGetArray( V0, &V_tmp ); CHKERRQ( ierr );
  for( m=1; m<=E->sphere.caps_per_proc; ++m ) {
    for( i = 0; i < neq; i++ )
      V[CPPR][i] = V_tmp[i]; 
  }
  ierr = VecRestoreArray( V0, &V_tmp ); CHKERRQ( ierr );
  
  ierr = VecGetArray( P0, &P_tmp ); CHKERRQ( ierr );
  for( m = 1; m < E->sphere.caps_per_proc; ++m ) {
    for( i = 0; i < nel; i++ )
      P[CPPR][i+1] = P_tmp[i]; 
  }
  ierr = VecRestoreArray( P0, &P_tmp ); CHKERRQ( ierr );


  ierr = VecDestroy( &rt ); CHKERRQ( ierr );
  ierr = VecDestroy( &p1 ); CHKERRQ( ierr );
  ierr = VecDestroy( &p2 ); CHKERRQ( ierr );
  ierr = VecDestroy( &BPI ); CHKERRQ( ierr );
  ierr = VecDestroy( &pt ); CHKERRQ( ierr );
  ierr = VecDestroy( &u0 ); CHKERRQ( ierr );
  ierr = VecDestroy( &s0 ); CHKERRQ( ierr );
  ierr = VecDestroy( &st ); CHKERRQ( ierr );
  ierr = VecDestroy( &u1 ); CHKERRQ( ierr );
  ierr = VecDestroy( &t0 ); CHKERRQ( ierr );
  ierr = VecDestroy( &r2 ); CHKERRQ( ierr );
  ierr = VecDestroy( &v0 ); CHKERRQ( ierr );
  ierr = VecDestroy( &V0 ); CHKERRQ( ierr );
  ierr = VecDestroy( &P0 ); CHKERRQ( ierr );

  *steps_max = count;

  PetscFunctionReturn(0);
}
#endif /* USE_PETSC */

/* Solve compressible Stokes flow using
 * bi-conjugate gradient stablized (BiCG-stab) iterations
 */

static void solve_Ahat_p_fhat_BiCG(struct All_variables *E,
                                   double **V, double **P, double **FF,
                                   double imp, int *steps_max)
{
    void assemble_div_rho_u();
    void assemble_del2_u();
    void assemble_grad_p();
    void strip_bcs_from_residual();
    int  solve_del2_u();
    void parallel_process_termination();

    double global_pdot();
    double global_v_norm2(), global_p_norm2(), global_div_norm2();
    double CPU_time0();

    int npno, neq;
    int m, j, count, lev;
    int valid;

    double alpha, beta, omega,inner_imp;
    double r0dotrt, r1dotrt;
    double v_norm, p_norm;
    double dvelocity, dpressure;
    int converging;

    double *F[NCS];
    double *r1[NCS], *r2[NCS], *pt[NCS], *p1[NCS], *p2[NCS];
    double *rt[NCS], *v0[NCS], *s0[NCS], *st[NCS], *t0[NCS];
    double *u0[NCS];
    double *shuffle[NCS];

    double time0, v_res;
    
    inner_imp = imp * E->control.inner_accuracy_scale; /* allow for different innner loop accuracy */

    npno = E->lmesh.npno;
    neq = E->lmesh.neq;
    lev = E->mesh.levmax;

    for (m=1; m<=E->sphere.caps_per_proc; m++)   {
        F[CPPR] = (double *)malloc(neq*sizeof(double));
        r1[CPPR] = (double *)malloc(npno*sizeof(double));
        r2[CPPR] = (double *)malloc(npno*sizeof(double));
        pt[CPPR] = (double *)malloc(npno*sizeof(double));
        p1[CPPR] = (double *)malloc(npno*sizeof(double));
        p2[CPPR] = (double *)malloc(npno*sizeof(double));
        rt[CPPR] = (double *)malloc(npno*sizeof(double));
        v0[CPPR] = (double *)malloc(npno*sizeof(double));
        s0[CPPR] = (double *)malloc(npno*sizeof(double));
        st[CPPR] = (double *)malloc(npno*sizeof(double));
        t0[CPPR] = (double *)malloc(npno*sizeof(double));

        u0[CPPR] = (double *)malloc(neq*sizeof(double));
    }

    time0 = CPU_time0();
    count = 0;
    v_res = E->monitor.fdotf;

    /* copy the original force vector since we need to keep it intact
       between iterations */
    for(m=1;m<=E->sphere.caps_per_proc;m++)
        for(j=0;j<neq;j++)
            F[CPPR][j] = FF[CPPR][j];


    /* calculate the initial velocity residual */
    initial_vel_residual(E, V, P, F, inner_imp*v_res);


    /* initial residual r1 = div(rho_ref*V) */
    assemble_div_rho_u(E, V, r1, lev);

    E->monitor.vdotv = global_v_norm2(E, V);
    E->monitor.incompressibility = sqrt(global_div_norm2(E, r1)
                                        / (1e-32 + E->monitor.vdotv));

    v_norm = sqrt(E->monitor.vdotv);
    p_norm = sqrt(E->monitor.pdotp);
    dvelocity = 1.0;
    dpressure = 1.0;
    converging = 0;


    if (E->control.print_convergence && E->parallel.me==0)  {
        print_convergence_progress(E, count, time0,
                                   v_norm, p_norm,
                                   dvelocity, dpressure,
                                   E->monitor.incompressibility);
    }


    /* initial conjugate residual rt = r1 */
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(j=0; j<npno; j++)
            rt[CPPR][j] = r1[CPPR][j];


    valid = 1;
    r0dotrt = alpha = omega = 0;

    while( (count < *steps_max) && keep_iterating(E, imp, converging) ) {
        /* require two consecutive converging iterations to quit the while-loop */

        /* r1dotrt = <r1, rt> */
        r1dotrt = global_pdot(E, r1, rt, lev);
        if(r1dotrt == 0.0) {
            /* XXX: can we resume the computation when BiCGstab failed? */
            fprintf(E->fp, "BiCGstab method failed!!\n");
            fprintf(stderr, "BiCGstab method failed!!\n");
            parallel_process_termination();
        }


        /* update search direction */
        if(count == 0)
            for (m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=0; j<npno; j++)
                    p2[CPPR][j] = r1[CPPR][j];
        else {
            /* p2 = r1 + <r1,rt>/<r0,rt> * alpha/omega * (p1 - omega*v0) */
            beta = (r1dotrt / r0dotrt) * (alpha / omega);
            for(m=1; m<=E->sphere.caps_per_proc; m++)
                for(j=0; j<npno; j++)
                    p2[CPPR][j] = r1[CPPR][j] + beta*(p1[CPPR][j] - omega*v0[CPPR][j]);
        }


        /* preconditioner BPI ~= inv(K), pt = BPI*p2 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<npno; j++)
                /* change to E->BPI[lev][m][j] after it has been made 0-based */
                pt[CPPR][j] = E->BPI[lev][CPPR][j+1] * p2[CPPR][j];


        /* solve K*u0 = grad(pt) for u1 */
        assemble_grad_p(E, pt, F, lev);
        valid = solve_del2_u(E, u0, F, inner_imp*v_res, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 1\n", stderr);
            fputs("Warning: solver not converging! 1\n", E->fp);
        }
        strip_bcs_from_residual(E, u0, lev);


        /* v0 = div(rho_ref*u0) */
        assemble_div_rho_u(E, u0, v0, lev);


        /* alpha = r1dotrt / <rt, v0> */
        alpha = r1dotrt / global_pdot(E, rt, v0, lev);


        /* s0 = r1 - alpha * v0 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<npno; j++)
                s0[CPPR][j] = r1[CPPR][j] - alpha * v0[CPPR][j];


        /* preconditioner BPI ~= inv(K), st = BPI*s0 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<npno; j++)
                /* change to E->BPI[lev][m][j] after it has been made 0-based */
                st[CPPR][j] = E->BPI[lev][CPPR][j+1] * s0[CPPR][j];


        /* solve K*u1 = grad(st) for u1 */
        assemble_grad_p(E, st, F, lev);
        valid = solve_del2_u(E, E->u1, F, inner_imp*v_res, lev);
        if(!valid && (E->parallel.me==0)) {
            fputs("Warning: solver not converging! 2\n", stderr);
            fputs("Warning: solver not converging! 2\n", E->fp);
        }
        strip_bcs_from_residual(E, E->u1, lev);


        /* t0 = div(rho_ref * u1) */
        assemble_div_rho_u(E, E->u1, t0, lev);


        /* omega = <t0, s0> / <t0, t0> */
        omega = global_pdot(E, t0, s0, lev) / global_pdot(E, t0, t0, lev);


        /* r2 = s0 - omega * t0 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<npno; j++)
                r2[CPPR][j] = s0[CPPR][j] - omega * t0[CPPR][j];


        /* P = P + alpha * pt + omega * st */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<npno; j++)
                s0[CPPR][j] = alpha * pt[CPPR][j] + omega * st[CPPR][j];

        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<npno; j++)
                P[CPPR][j] += s0[CPPR][j];


        /* V = V - alpha * u0 - omega * u1 */
        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                F[CPPR][j] = alpha * u0[CPPR][j] + omega * E->u1[CPPR][j];

        for(m=1; m<=E->sphere.caps_per_proc; m++)
            for(j=0; j<neq; j++)
                V[CPPR][j] -= F[CPPR][j];


        /* compute velocity and incompressibility residual */
        E->monitor.vdotv = global_v_norm2(E, V);
        E->monitor.pdotp = global_p_norm2(E, P);
        v_norm = sqrt(E->monitor.vdotv);
        p_norm = sqrt(E->monitor.pdotp);
        dvelocity = sqrt(global_v_norm2(E, F) / (1e-32 + E->monitor.vdotv));
        dpressure = sqrt(global_p_norm2(E, s0) / (1e-32 + E->monitor.pdotp));
	

        assemble_div_rho_u(E, V, t0, lev);
        E->monitor.incompressibility = sqrt(global_div_norm2(E, t0)
                                            / (1e-32 + E->monitor.vdotv));

        count++;

        if(E->control.print_convergence && E->parallel.me==0) {
            print_convergence_progress(E, count, time0,
                                       v_norm, p_norm,
                                       dvelocity, dpressure,
                                       E->monitor.incompressibility);
        }

	if(!valid){
            /* reset consecutive converging iterations */
            converging = 0;
	}else{
            /* how many consecutive converging iterations? */
            if(E->control.check_pressure_convergence) {
                /* check dv and dp */
                if(dvelocity < imp && dpressure < imp)
                    converging++;
                else
                    converging = 0;
            }else{
                /* check dv only */
                if(dvelocity < imp)
                    converging++;
                else
                    converging = 0;
            }
	  
	}

	/* shift array pointers */
        for(m=1; m<=E->sphere.caps_per_proc; m++) {
            shuffle[CPPR] = p1[CPPR];
            p1[CPPR] = p2[CPPR];
            p2[CPPR] = shuffle[CPPR];

            shuffle[CPPR] = r1[CPPR];
            r1[CPPR] = r2[CPPR];
            r2[CPPR] = shuffle[CPPR];
        }

        /* shift <r0, rt> = <r1, rt> */
        r0dotrt = r1dotrt;

    } /* end loop for conjugate gradient */


    for(m=1; m<=E->sphere.caps_per_proc; m++) {
    	free((void *) F[CPPR]);
        free((void *) r1[CPPR]);
        free((void *) r2[CPPR]);
        free((void *) pt[CPPR]);
        free((void *) p1[CPPR]);
        free((void *) p2[CPPR]);
        free((void *) rt[CPPR]);
        free((void *) v0[CPPR]);
        free((void *) s0[CPPR]);
        free((void *) st[CPPR]);
        free((void *) t0[CPPR]);

        free((void *) u0[CPPR]);
    }

    *steps_max=count;
}


/* Solve compressible Stokes flow using
 * conjugate gradient (CG) iterations with an outer iteration
 */

static void solve_Ahat_p_fhat_iterCG(struct All_variables *E,
                                     double **V, double **P, double **F,
                                     double imp, int *steps_max)
{
    int m, i;
    int cycles, num_of_loop;
    double relative_err_v, relative_err_p;
    double *old_v[NCS], *old_p[NCS],*diff_v[NCS],*diff_p[NCS];
    double div_res;
    const int npno = E->lmesh.npno;
    const int neq = E->lmesh.neq;
    const int lev = E->mesh.levmax;

    double global_v_norm2(),global_p_norm2();
    double global_div_norm2();
    void assemble_div_rho_u();
    
    for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    	old_v[CPPR] = (double *)malloc(neq*sizeof(double));
    	diff_v[CPPR] = (double *)malloc(neq*sizeof(double));
    	old_p[CPPR] = (double *)malloc(npno*sizeof(double));
    	diff_p[CPPR] = (double *)malloc(npno*sizeof(double));
    }

    cycles = E->control.p_iterations;

    initial_vel_residual(E, V, P, F,
                         imp * E->control.inner_accuracy_scale * E->monitor.fdotf);

    div_res = 1.0;
    relative_err_v = 1.0;
    relative_err_p = 1.0;
    num_of_loop = 0;

    while((relative_err_v >= imp || relative_err_p >= imp) &&
          (div_res > imp) &&
          (num_of_loop <= E->control.compress_iter_maxstep)) {

        for (m=1;m<=E->sphere.caps_per_proc;m++) {
            for(i=0;i<neq;i++) old_v[CPPR][i] = V[CPPR][i];
            for(i=0;i<npno;i++) old_p[CPPR][i] = P[CPPR][i];
        }
#ifdef USE_PETSC
        if(E->control.use_petsc)
          solve_Ahat_p_fhat_CG_PETSc(E, V, P, F, imp, &cycles);
        else
#endif
          solve_Ahat_p_fhat_CG(E, V, P, F, imp, &cycles);


        /* compute norm of div(rho*V) */
        assemble_div_rho_u(E, V, E->u1, lev);
        div_res = sqrt(global_div_norm2(E, E->u1) / (1e-32 + E->monitor.vdotv));

        for (m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=0;i<neq;i++) 
              diff_v[CPPR][i] = V[CPPR][i] - old_v[CPPR][i];

        relative_err_v = sqrt( global_v_norm2(E,diff_v) /
                               (1.0e-32 + E->monitor.vdotv) );

        for (m=1;m<=E->sphere.caps_per_proc;m++)
            for(i=0;i<npno;i++) 
              diff_p[CPPR][i] = P[CPPR][i] - old_p[CPPR][i];

        relative_err_p = sqrt( global_p_norm2(E,diff_p) /
                               (1.0e-32 + E->monitor.pdotp) );

        if(E->parallel.me == 0) {
            fprintf(stderr, "itercg -- div(rho*v)/v=%.2e dv/v=%.2e and dp/p=%.2e loop %d\n\n", div_res, relative_err_v, relative_err_p, num_of_loop);
            fprintf(E->fp, "itercg -- div(rho*v)/v=%.2e dv/v=%.2e and dp/p=%.2e loop %d\n\n", div_res, relative_err_v, relative_err_p, num_of_loop);
        }

        num_of_loop++;

    } /* end of while */

    for (m=1;m<=E->sphere.caps_per_proc;m++)   {
    	free((void *) old_v[CPPR]);
    	free((void *) old_p[CPPR]);
	free((void *) diff_v[CPPR]);
	free((void *) diff_p[CPPR]);
    }

    return;
}


static void initial_vel_residual(struct All_variables *E,
                                 double **V, double **P, double **F,
                                 double acc)
{
    void assemble_del2_u();
    void assemble_grad_p();
    void strip_bcs_from_residual();
    int  solve_del2_u();

    int neq = E->lmesh.neq;
    int lev = E->mesh.levmax;
    int i, m, valid;

    /* F = F - grad(P) - K*V */
    assemble_grad_p(E, P, E->u1, lev);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            F[CPPR][i] = F[CPPR][i] - E->u1[CPPR][i];

    assemble_del2_u(E, V, E->u1, lev, 1);
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            F[CPPR][i] = F[CPPR][i] - E->u1[CPPR][i];

    strip_bcs_from_residual(E, F, lev);


    /* solve K*u1 = F for u1 */
    valid = solve_del2_u(E, E->u1, F, acc, lev);
    if(!valid && (E->parallel.me==0)) {
        fputs("Warning: solver not converging! 0\n", stderr);
        fputs("Warning: solver not converging! 0\n", E->fp);
    }
    strip_bcs_from_residual(E, E->u1, lev);


    /* V = V + u1 */
    for(m=1; m<=E->sphere.caps_per_proc; m++)
        for(i=0; i<neq; i++)
            V[CPPR][i] += E->u1[CPPR][i];

    return;
}
