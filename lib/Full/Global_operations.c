#include <mpi.h>

#include <math.h>
#include <sys/types.h>
#include "element_definitions.h"
#include "global_defs.h"


MPI_Comm get_horizontal_communicator(E, layer)
     struct All_variables *E;
     int layer;
{
  static MPI_Comm MY_HORZ_COMM = MPI_COMM_NULL;
  MPI_Group world_g, horizon_g;
  int m,i,j,k,nproch;
  int *processors;
	
  if (MY_HORZ_COMM == MPI_COMM_NULL) {
    nproch = E->parallel.nprocxy*E->parallel.surf_proc_per_cap;
    processors = (int *) malloc(nproch*sizeof(int));

    k=0;
    for (m=0;m<E->sphere.caps;m++)
      for (i=0;i<E->parallel.nprocx;i++)
	for (j=0;j<E->parallel.nprocy;j++) 
	  processors[k++] = E->parallel.loc2proc_map[m][i][j][layer];
	
    MPI_Comm_group(MPI_COMM_WORLD, &world_g);
    MPI_Group_incl(world_g, nproch, processors, &horizon_g);
    MPI_Comm_create(MPI_COMM_WORLD, horizon_g, &MY_HORZ_COMM);

    MPI_Group_free(&horizon_g);
    MPI_Group_free(&world_g);
    free ((void *) processors);
  }

  return(MY_HORZ_COMM);
}	

/* ===============================================
   strips horizontal average from nodal field X. 
   Assumes orthogonal mesh, otherwise, horizontals
   aren't & another method is required.
   =============================================== */

void remove_horiz_ave(E,X,H,store_or_not)
     struct All_variables *E;
     double **X, *H;
     int store_or_not;

{
    int m,i,j,k,l,n,ln,nox,noz,noy;
    void return_horiz_ave();

    const int dims = E->mesh.nsd;

    noy = E->lmesh.noy;
    noz = E->lmesh.noz;
    nox = E->lmesh.nox;

    return_horiz_ave(E,X,H);

    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(k=1;k<=noy;k++)     
	for(j=1;j<=nox;j++)
	  for(i=1;i<=noz;i++) {
	    n = i+(j-1)*noz+(k-1)*noz*nox;
	    X[m][n] -= H[i];  
	  }
    
    return; 
}


void return_horiz_ave(E,X,H)
     struct All_variables *E;
     double **X, *H;
{ 
  const int dims = E->mesh.nsd;
  int m,i,j,k,l,d,nint,noz,nox,noy,el,elz,elx,ely,j1,j2,i1,i2,k1,k2,nproc;
  int top,lnode[5], sizeofH, noz2,iroot;
  double *Have,*temp,aa[5];
  struct Shape_function1 M;
  struct Shape_function1_dA dGamma;
  void get_global_1d_shape_fn();
  MPI_Comm get_horizontal_communicator();
  
  sizeofH = (2*E->lmesh.noz+2)*sizeof(double);
  Have = (double *)malloc(sizeofH);
  temp = (double *)malloc(sizeofH);

  noz = E->lmesh.noz;
  noy = E->lmesh.noy;
  elz = E->lmesh.elz;
  elx = E->lmesh.elx;
  ely = E->lmesh.ely;
  noz2 = 2*noz;

  for (i=1;i<=elz;i++)  {
    temp[i] = temp[i+noz] = 0.0;
    temp[i+1] = temp[i+1+noz] = 0.0;
    top = 0;
    if (i==elz) top = 1;
    for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (k=1;k<=ely;k++)
        for (j=1;j<=elx;j++)     {
          el = i + (j-1)*elz + (k-1)*elx*elz;
          get_global_1d_shape_fn(E,el,&M,&dGamma,top,m);   

          lnode[1] = E->ien[m][el].node[1];
          lnode[2] = E->ien[m][el].node[2];
          lnode[3] = E->ien[m][el].node[3];
          lnode[4] = E->ien[m][el].node[4];

          for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
            for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
              temp[i] += X[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] 
                          * dGamma.vpt[GMVGAMMA(0,nint)];
            temp[i+noz] += dGamma.vpt[GMVGAMMA(0,nint)];
            }

          if (i==elz)  {
            lnode[1] = E->ien[m][el].node[5];
            lnode[2] = E->ien[m][el].node[6];
            lnode[3] = E->ien[m][el].node[7];
            lnode[4] = E->ien[m][el].node[8];

            for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
              for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
                temp[i+1] += X[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] 
                          * dGamma.vpt[GMVGAMMA(1,nint)];
              temp[i+1+noz] += dGamma.vpt[GMVGAMMA(1,nint)];
              }

            }   /* end of if i==elz    */
          }   /* end of j  and k, and m  */
     }        /* Done for i */ 


  if (E->parallel.nproc>1)  {
    MPI_Allreduce(temp,Have,noz2+1,MPI_DOUBLE,MPI_SUM,
    			  get_horizontal_communicator(E, E->parallel.me_loc[3]));
  }
  else
    for (i=1;i<=noz2;i++)
      Have[i] = temp[i];

  for (i=1;i<=noz;i++) {
    if(Have[i+noz] != 0.0) 
       H[i] = Have[i]/Have[i+noz];
    }
 /* if (E->parallel.me==0)
    for(i=1;i<=noz;i++)
      fprintf(stderr,"area %d %d %g\n",E->parallel.me,i,Have[i+noz]);
*/
  free ((void *) Have);
  free ((void *) temp);
  
  return; 
  }
	   
void return_horiz_ave_f(E,X,H)
     struct All_variables *E;
     float **X, *H;
{ 
  const int dims = E->mesh.nsd;
  int m,i,j,k,d,nint,noz,nox,noy,el,elz,elx,ely,j1,j2,i1,i2,k1,k2,nproc;
  int top,lnode[5], sizeofH, noz2,iroot;
  float *Have,*temp,aa[5];
  struct Shape_function1 M;
  struct Shape_function1_dA dGamma;
  void get_global_1d_shape_fn();
  MPI_Comm get_horizontal_communicator();

  sizeofH = (2*E->lmesh.noz+2)*sizeof(float);
  Have = (float *)malloc(sizeofH);
  temp = (float *)malloc(sizeofH);

  noz = E->lmesh.noz;
  noy = E->lmesh.noy;
  elz = E->lmesh.elz;
  elx = E->lmesh.elx;
  ely = E->lmesh.ely;
  noz2 = 2*noz;

  for (i=1;i<=elz;i++)  {
    temp[i] = temp[i+noz] = 0.0;
    temp[i+1] = temp[i+1+noz] = 0.0;
    top = 0;
    if (i==elz) top = 1;
    for (m=1;m<=E->sphere.caps_per_proc;m++)
      for (k=1;k<=ely;k++)
        for (j=1;j<=elx;j++)     {
          el = i + (j-1)*elz + (k-1)*elx*elz;
          get_global_1d_shape_fn(E,el,&M,&dGamma,top,m);   

          lnode[1] = E->ien[m][el].node[1];
          lnode[2] = E->ien[m][el].node[2];
          lnode[3] = E->ien[m][el].node[3];
          lnode[4] = E->ien[m][el].node[4];

          for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
            for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
              temp[i] += X[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] 
                          * dGamma.vpt[GMVGAMMA(0,nint)];
            temp[i+noz] += dGamma.vpt[GMVGAMMA(0,nint)];
            }

          if (i==elz)  {
            lnode[1] = E->ien[m][el].node[5];
            lnode[2] = E->ien[m][el].node[6];
            lnode[3] = E->ien[m][el].node[7];
            lnode[4] = E->ien[m][el].node[8];

            for(nint=1;nint<=onedvpoints[E->mesh.nsd];nint++)   {
              for(d=1;d<=onedvpoints[E->mesh.nsd];d++)
                temp[i+1] += X[m][lnode[d]] * E->M.vpt[GMVINDEX(d,nint)] 
                          * dGamma.vpt[GMVGAMMA(1,nint)];
              temp[i+1+noz] += dGamma.vpt[GMVGAMMA(1,nint)];
              }

            }   /* end of if i==elz    */
          }   /* end of j  and k, and m  */
     }        /* Done for i */ 


  if (E->parallel.nproc>1)  {
    MPI_Allreduce(temp,Have,noz2+1,MPI_FLOAT,MPI_SUM,
    			  get_horizontal_communicator(E, E->parallel.me_loc[3]));
    }
  else
    for (i=1;i<=noz2;i++)
      Have[i] = temp[i];

  for (i=1;i<=noz;i++) {
    if(Have[i+noz] != 0.0) 
       H[i] = Have[i]/Have[i+noz];
    }
 /* if (E->parallel.me==0)
    for(i=1;i<=noz;i++)
      fprintf(stderr,"area %d %d %g\n",E->parallel.me,i,Have[i+noz]);
*/
  free ((void *) Have);
  free ((void *) temp);
  
  return; 
  }

float return_bulk_value(E,Z,average) 
     struct All_variables *E;
     float **Z;
     int average;

{  
    void get_global_shape_fn();
    void float_global_operation();

    double rtf[4][9];

    int n,i,j,k,el,m;
    float volume,integral,volume1,integral1;
   
    struct Shape_function GN;
    struct Shape_function_dx GNx;
    struct Shape_function_dA dOmega;
    
    const int vpts = vpoints[E->mesh.nsd];
    const int ends = enodes[E->mesh.nsd];
    const int sphere_key=1;

    volume1=0.0;
    integral1=0.0;
     
    for (m=1;m<=E->sphere.caps_per_proc;m++)
       for (el=1;el<=E->lmesh.nel;el++)  {

	  get_global_shape_fn(E,el,&GN,&GNx,&dOmega,0,sphere_key,rtf,E->mesh.levmax,m);
	
	  for(j=1;j<=vpts;j++)
	    for(i=1;i<=ends;i++) {
		n = E->ien[m][el].node[i];
		volume1 += E->N.vpt[GNVINDEX(i,j)] * dOmega.vpt[j];
		integral1 += Z[m][n] * E->N.vpt[GNVINDEX(i,j)] * dOmega.vpt[j];
                }

          }


    MPI_Allreduce(&volume1  ,&volume  ,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&integral1,&integral,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
     
    if(average && volume != 0.0)
 	   integral /= volume;

    return((float)integral);
}

/* ================================================== */
/* void sum_across_depth_sph1(E,sphc,sphs) */
/* struct All_variables *E; */
/* float *sphc,*sphs; */
/* { */
/*  int jumpp,total,j,d; */

/*  static float *sphcs,*temp; */
/*  static int been_here=0; */
/*  static int *processors,nproc; */

/*  static MPI_Comm world, horizon_p; */
/*  static MPI_Group world_g, horizon_g; */

/* if (been_here==0)  { */
/*  processors = (int *)malloc((E->parallel.nprocz+2)*sizeof(int)); */
/*  temp = (float *) malloc((E->sphere.hindice*2+3)*sizeof(float)); */
/*  sphcs = (float *) malloc((E->sphere.hindice*2+3)*sizeof(float)); */

/*  nproc = 0; */
/*  for (j=0;j<E->parallel.nprocz;j++) { */
/*    d =E->parallel.me_sph*E->parallel.nprocz+E->parallel.nprocz-1-j;  */
/*    processors[nproc] =  d; */
/*    nproc ++; */
/*    } */

/*  if (nproc>0)  { */
/*     world = MPI_COMM_WORLD; */
/*     MPI_Comm_group(world, &world_g); */
/*     MPI_Group_incl(world_g, nproc, processors, &horizon_g); */
/*     MPI_Comm_create(world, horizon_g, &horizon_p); */
/*     } */

/*  been_here++; */
/*  } */

/*  total = E->sphere.hindice*2+3; */
/*   jumpp = E->sphere.hindice; */
/*   for (j=0;j<E->sphere.hindice;j++)   { */
/*       sphcs[j] = sphc[j]; */
/*       sphcs[j+jumpp] = sphs[j]; */
/*      } */


/*  if (nproc>0)  { */

/*     MPI_Allreduce(sphcs,temp,total,MPI_FLOAT,MPI_SUM,horizon_p); */

/*     for (j=0;j<E->sphere.hindice;j++)   { */
/*       sphc[j] = temp[j]; */
/*       sphs[j] = temp[j+jumpp]; */
/*      } */

/*     } */

/* return; */
/* } */

/* ================================================== */
float find_max_horizontal(E,Tmax)
struct All_variables *E;
float Tmax;
{
 int j,d;
 static int been_here=0;
 float ttmax;
 MPI_Comm get_horizontal_communicator();

 if (E->parallel.nproc>1)
   MPI_Allreduce(&Tmax,&ttmax,1,MPI_FLOAT,MPI_MAX,
   				 get_horizontal_communicator(E, E->parallel.me_loc[3]));
 else
   ttmax = Tmax;

 return(ttmax);
 }

/* ================================================== */
void sum_across_surface(E,data,total)
struct All_variables *E;
float *data;
int total;
{
 int j,d;
 float *temp;
 MPI_Comm get_horizontal_communicator();

 if (E->parallel.nproc>1)  {

    temp = (float *)malloc((total+1)*sizeof(float));
    MPI_Allreduce(data,temp,total,MPI_FLOAT,MPI_SUM,
    			  get_horizontal_communicator(E, E->parallel.me_loc[3]));

    for (j=0;j<total;j++)   {
      data[j] = temp[j];
      }
    free((void *)temp); 
    }
		    
return;
}

/* ================================================== */
/* ================================================== */

/* ================================================== */
void sum_across_surf_sph1(E,sphc,sphs)
struct All_variables *E;
float *sphc,*sphs;
{
 int jumpp,total,j,d;
 static float *sphcs,*temp;
 static int been_here=0;
 MPI_Comm get_horizontal_communicator();
 
 if (been_here==0)  {
 temp = (float *) malloc((E->sphere.hindice*2+3)*sizeof(float));
 sphcs = (float *) malloc((E->sphere.hindice*2+3)*sizeof(float));
 been_here++;
 }

 jumpp = E->sphere.hindice;
 total = E->sphere.hindice*2+3;
 for (j=0;j<E->sphere.hindice;j++)   {
    sphcs[j] = sphc[j];
    sphcs[j+jumpp] = sphs[j];
    }

 if (E->parallel.nproc>1)  {
    MPI_Allreduce(sphcs,temp,total,MPI_FLOAT,MPI_SUM,
				  get_horizontal_communicator(E, E->parallel.me_loc[3]));
    for (j=0;j<E->sphere.hindice;j++)   {
      sphc[j] = temp[j];
      sphs[j] = temp[j+jumpp];
      }
    }
		    
return;
}

/* ================================================== */


float global_fvdot(E,A,B,lev)
   struct All_variables *E;
   float **A,**B;
   int lev;

{
  int m,i,neq;
  float prod, temp,temp1;

  neq=E->lmesh.NEQ[lev];

  temp = 0.0;
  temp1 = 0.0;
  prod = 0.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    neq=E->lmesh.NEQ[lev];
    temp1 = 0.0;
    for (i=0;i<neq;i++)
      temp += A[m][i]*B[m][i]; 

    for (i=1;i<=E->parallel.Skip_neq[lev][m];i++)
       temp1 += A[m][E->parallel.Skip_id[lev][m][i]]*B[m][E->parallel.Skip_id[lev][m][i]];

    temp -= temp1;

    }  

  MPI_Allreduce(&temp, &prod,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

  return (prod);
}


double kineticE_radial(E,A,lev)
   struct All_variables *E;
   double **A;
   int lev;

{
  int m,i,neq;
  double prod, temp,temp1;

    temp = 0.0;
    prod = 0.0;

  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    neq=E->lmesh.NEQ[lev];
    temp1 = 0.0;
    for (i=0;i<neq;i++)  
      if ((i+1)%3==0)
        temp += A[m][i]*A[m][i];

    for (i=1;i<=E->parallel.Skip_neq[lev][m];i++)
      if ((E->parallel.Skip_id[lev][m][i]+1)%3==0)
        temp1 += A[m][E->parallel.Skip_id[lev][m][i]]*A[m][E->parallel.Skip_id[lev][m][i]];

    temp -= temp1;

    }  

  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  return (prod);
}

double global_vdot(E,A,B,lev)
   struct All_variables *E;
   double **A,**B;
   int lev;

{
  int m,i,neq;
  double prod, temp,temp1;

    temp = 0.0;
    prod = 0.0;

  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    neq=E->lmesh.NEQ[lev];
    temp1 = 0.0;
    for (i=0;i<neq;i++)
      temp += A[m][i]*B[m][i]; 

    for (i=1;i<=E->parallel.Skip_neq[lev][m];i++)
       temp1 += A[m][E->parallel.Skip_id[lev][m][i]]*B[m][E->parallel.Skip_id[lev][m][i]];

    temp -= temp1;

    }  

  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  return (prod);
}


double global_pdot(E,A,B,lev)
   struct All_variables *E;
   double **A,**B;
   int lev;

{
  int i,m,npno;
  double prod, temp;

  npno=E->lmesh.NPNO[lev];

  temp = 0.0;
  prod = 0.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    npno=E->lmesh.NPNO[lev];
    for (i=1;i<=npno;i++)
      temp += A[m][i]*B[m][i]; 
    }

  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  return (prod);
  }


double global_tdot_d(E,A,B,lev)
   struct All_variables *E;
   double **A,**B;
   int lev;

{
  int i,nno,m;
  double prod, temp;

  nno=E->lmesh.NNO[lev];

  temp = 0.0;
  prod = 0.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    nno=E->lmesh.NNO[lev];
    for (i=1;i<=nno;i++)
    if (!(E->NODE[lev][m][i] & SKIP))
      temp += A[m][i]; 
    }

  MPI_Allreduce(&temp, &prod,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  return (prod);
  }

float global_tdot(E,A,B,lev)
   struct All_variables *E;
   float **A,**B;
   int lev;

{
  int i,nno,m;
  float prod, temp;


  temp = 0.0;
  prod = 0.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)  {
    nno=E->lmesh.NNO[lev];
    for (i=1;i<=nno;i++)
      if (!(E->NODE[lev][m][i] & SKIP))
        temp += A[m][i]*B[m][i]; 
    }

  MPI_Allreduce(&temp, &prod,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

  return (prod);
  }


float global_fmin(E,a)
   struct All_variables *E;
   float a;
{
  float temp;
  MPI_Allreduce(&a, &temp,1,MPI_FLOAT,MPI_MIN,MPI_COMM_WORLD);
  return (temp);
  }

double global_dmax(E,a)
   struct All_variables *E;
   double a;
{
  double temp;
  MPI_Allreduce(&a, &temp,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  return (temp);
  }


float global_fmax(E,a)
   struct All_variables *E;
   float a;
{
  float temp;
  MPI_Allreduce(&a, &temp,1,MPI_FLOAT,MPI_MAX,MPI_COMM_WORLD);
  return (temp);
  }

double Tmaxd(E,T)
  struct All_variables *E;
  double **T;
{
  double global_dmax(),temp,temp1;
  int i,m;

  temp = -10.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=E->lmesh.nno;i++)
      temp = max(T[m][i],temp);

  temp1 = global_dmax(E,temp);
  return (temp1);
  }


float Tmax(E,T)
  struct All_variables *E;
  float **T;
{
  float global_fmax(),temp,temp1;
  int i,m;

  temp = -10.0;
  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=E->lmesh.nno;i++)
      temp = max(T[m][i],temp);

  temp1 = global_fmax(E,temp);
  return (temp1);
  }


float  vnorm_nonnewt(E,dU,U,lev)
  struct All_variables *E;
  float **dU,**U;
  int lev;
{
 float temp1,temp2,dtemp,temp;
 int a,e,i,m,node;
 const int dims = E->mesh.nsd;
 const int ends = enodes[dims];
 const int nel=E->lmesh.nel;

 dtemp=0.0;
 temp=0.0;
 for (m=1;m<=E->sphere.caps_per_proc;m++)    
  for (e=1;e<=nel;e++)
   if (E->mat[m][e]==1)
     for (i=1;i<=dims;i++)
       for (a=1;a<=ends;a++) {
	 node = E->IEN[lev][m][e].node[a];
         dtemp += dU[m][ E->ID[lev][m][node].doff[i] ]*
                  dU[m][ E->ID[lev][m][node].doff[i] ];
         temp += U[m][ E->ID[lev][m][node].doff[i] ]*
                 U[m][ E->ID[lev][m][node].doff[i] ];
         }
   

  MPI_Allreduce(&dtemp, &temp2,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&temp, &temp1,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);

  temp1 = sqrt(temp2/temp1);

  return (temp1);
  }
