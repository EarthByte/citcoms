/* Routine to process the output of the finite element cycles
   and to turn them into a coherent suite  files  */


#include <stdlib.h>

#include "element_definitions.h"
#include "global_defs.h"
#include "output.h"

void output_coord(struct All_variables *);
void output_velo(struct All_variables *, int);
void output_visc_prepare(struct All_variables *, float **);
void output_visc(struct All_variables *, int);
void output_surf_botm(struct All_variables *, int);
void output_ave_r(struct All_variables *, int);
void output_mat(struct All_variables *);

void parallel_process_termination();

/*********************************************************************/


void output(struct All_variables *E, int cycles)
{

  if (cycles == 0) {
    output_coord(E);
    output_mat(E);
  }

  output_velo(E, cycles);
  output_visc(E, cycles);
  output_surf_botm(E, cycles);

  /* disable horizontal average output   by Tan2 */
  /* output_ave_r(E, cycles); */

  return;
}


FILE* output_open(char *filename)
{
  FILE *fp1;

  // if filename is empty, output to stderr.
  if (*filename) {
    fp1 = fopen(filename,"w");
    if (!fp1) {
      fprintf(stderr,"Cannot open file '%s'\n",filename);
      parallel_process_termination();
    }
  }
  else
    fp1 = stderr;

  return fp1;
}


void output_coord(struct All_variables *E)
{
  int i, j;
  char output_file[255];
  FILE *fp1;

  sprintf(output_file,"%s.coord.%d",E->control.data_file,E->parallel.me);
  fp1 = output_open(output_file);

  for(j=1;j<=E->sphere.caps_per_proc;j++)     {
    fprintf(fp1,"%3d %7d\n",j,E->lmesh.nno);
    for(i=1;i<=E->lmesh.nno;i++)
      fprintf(fp1,"%.3e %.3e %.3e\n",E->sx[j][1][i],E->sx[j][2][i],E->sx[j][3][i]);
  }

  fclose(fp1);

  return;
}


void output_visc_prepare(struct All_variables *E, float **VE)
{
  void get_ele_visc();
  void visc_from_ele_to_gint();
  void visc_from_gint_to_nodes();

  float *EV, *VN[NCS];
  const int lev = E->mesh.levmax;
  const int nsd = E->mesh.nsd;
  const int vpts = vpoints[nsd];
  int i, m;


  // Here is a bug in the original code. EV is not allocated for each
  // E->sphere.caps_per_proc. Later, when elemental viscosity is written
  // to it (in get_ele_visc()), viscosity in high cap number will overwrite
  // that in a lower cap number.
  //
  // Since current CitcomS only support 1 cap per processor, this bug won't
  // manifest itself. So, I will leave it here.
  // by Tan2 5/22/2003
  int size2 = (E->lmesh.nel+1)*sizeof(float);
  EV = (float *) malloc (size2);

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    VN[m]=(float *)malloc((1+E->lmesh.nel*vpts)*sizeof(float));
  }

  get_ele_visc(E,EV,1);

  for(i=1;i<=E->lmesh.nel;i++)
    VE[1][i]=EV[i];

  visc_from_ele_to_gint(E, VE, VN, lev);
  visc_from_gint_to_nodes(E, VN, VE, lev);

  free((void *) EV);
  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    free((void *) VN[m]);
  }

  return;
}


void output_visc(struct All_variables *E, int cycles)
{
  int i, j, m;
  char output_file[255];
  FILE *fp1;
  float *VE[NCS];

  sprintf(output_file,"%s.visc.%d.%d",E->control.data_file,E->parallel.me,cycles);
  fp1 = output_open(output_file);

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    VE[m]=(float *)malloc((1+E->lmesh.nno)*sizeof(float));
  }

  output_visc_prepare(E, VE);

  for(j=1;j<=E->sphere.caps_per_proc;j++) {
    fprintf(fp1,"%3d %7d\n",j,E->lmesh.nno);
    for(i=1;i<=E->lmesh.nno;i++)
      fprintf(fp1,"%.3e\n",VE[1][i]);
  }

  for(m=1;m<=E->sphere.caps_per_proc;m++) {
    free((void*) VE[m]);
  }

  fclose(fp1);

  return;
}



void output_velo(struct All_variables *E, int cycles)
{
  int i, j;
  char output_file[255];
  FILE *fp1;

  sprintf(output_file,"%s.velo.%d.%d",E->control.data_file,E->parallel.me,cycles);
  fp1 = output_open(output_file);

  fprintf(fp1,"%d %d %.5e\n",cycles,E->lmesh.nno,E->monitor.elapsed_time);

  for(j=1;j<=E->sphere.caps_per_proc;j++) {
    fprintf(fp1,"%3d %7d\n",j,E->lmesh.nno);
    for(i=1;i<=E->lmesh.nno;i++)
      fprintf(fp1,"%.6e %.6e %.6e %.6e\n",E->sphere.cap[j].V[1][i],E->sphere.cap[j].V[2][i],E->sphere.cap[j].V[3][i],E->T[j][i]);
  }

  fclose(fp1);

  return;
}


void output_surf_botm(struct All_variables *E, int cycles)
{
  int i, j, s;
  char output_file[255];
  FILE* fp2;

  void heat_flux();
  void get_STD_topo();

  heat_flux(E);
  get_STD_topo(E,E->slice.tpg,E->slice.tpgb,E->slice.divg,E->slice.vort,cycles);

  if (E->parallel.me_locl[3]==E->parallel.nproczl-1) {
    sprintf(output_file,"%s.surf.%d.%d",E->control.data_file,E->parallel.me,cycles);
    fp2 = output_open(output_file);

    for(j=1;j<=E->sphere.caps_per_proc;j++)  {
      fprintf(fp2,"%3d %7d\n",j,E->lmesh.nsf);
      for(i=1;i<=E->lmesh.nsf;i++)   {
	s = i*E->lmesh.noz;
        fprintf(fp2,"%.4e %.4e %.4e %.4e\n",E->slice.tpg[j][i],E->slice.shflux[j][i],E->sphere.cap[j].V[1][s],E->sphere.cap[j].V[2][s]);
      }
    }
    fclose(fp2);
  }


  if (E->parallel.me_locl[3]==0)      {
    sprintf(output_file,"%s.botm.%d.%d",E->control.data_file,E->parallel.me,cycles);
    fp2 = output_open(output_file);

    for(j=1;j<=E->sphere.caps_per_proc;j++)  {
      fprintf(fp2,"%3d %7d\n",j,E->lmesh.nsf);
      for(i=1;i<=E->lmesh.nsf;i++)  {
	s = (i-1)*E->lmesh.noz + 1;
        fprintf(fp2,"%.4e %.4e %.4e %.4e\n",E->slice.tpgb[j][i],E->slice.bhflux[j][i],E->sphere.cap[j].V[1][s],E->sphere.cap[j].V[2][s]);
      }
    }
    fclose(fp2);
  }

  return;
}


void output_ave_r(struct All_variables *E, int cycles)
{
  int j;
  char output_file[255];
  FILE* fp2;

  // compute horizontal average here....

  // only the first nproczl processors need to output
  if (E->parallel.me<E->parallel.nproczl)  {
    sprintf(output_file,"%s.ave_r.%d.%d",E->control.data_file,E->parallel.me,cycles);
    fp2 = output_open(output_file);

    for(j=1;j<=E->lmesh.noz;j++)  {
      fprintf(fp2,"%.4e %.4e %.4e %.4e\n",E->sx[1][3][j],E->Have.T[j],E->Have.V[1][j],E->Have.V[2][j]);
    }
    fclose(fp2);
  }


  return;
}



void output_mat(struct All_variables *E)
{
  int m, el;
  char output_file[255];
  FILE* fp;

  sprintf(output_file,"%s.mat.%d",E->control.data_file,E->parallel.me);
  fp = output_open(output_file);

  for (m=1;m<=E->sphere.caps_per_proc;m++)
    for(el=1;el<=E->lmesh.nel;el++)
      fprintf(fp,"%d %d %f\n", el,E->mat[m][el],E->VIP[m][el]);

  fclose(fp);

  return;
}



