
#include <math.h>

#include "global_defs.h"

//#include "age_related.h"
#include "parallel_related.h"
#include "parsing.h"
#include "lith_age.h"

float find_age_in_MY();
void lith_age_restart_conform_tbc(struct All_variables *E);


void lith_age_input(struct All_variables *E)
{
  int m = E->parallel.me;

  E->control.lith_age = 0;
  E->control.lith_age_time = 0;
  E->control.temperature_bound_adj = 0;

  input_int("lith_age",&(E->control.lith_age),"0",m);
  input_float("mantle_temp",&(E->control.lith_age_mantle_temp),"1.0",m);

  if (E->control.lith_age) {
    input_int("lith_age_time",&(E->control.lith_age_time),"0",m);
    input_string("lith_age_file",&(E->control.lith_age_file),"",m);
    input_float("lith_age_depth",&(E->control.lith_age_depth),"0.0471",m);

    input_int("temperature_bound_adj",&(E->control.temperature_bound_adj),"0",m);
    if (E->control.temperature_bound_adj) {
      input_float("depth_bound_adj",&(E->control.depth_bound_adj),"0.1570",m);
      input_float("width_bound_adj",&(E->control.width_bound_adj),"0.08727",m);
    }
  }
  return;
}

void lith_age_restart_tic(struct All_variables *E)
{
  int i, j, k, m, node, nodeg;
  int nox, noy, noz, gnox, gnoy, gnoz;
  double r1, temp;
  float age;

  char output_file[255];
  FILE *fp1;

  if(E->control.lith_age_time==1)   {
    /* if opening lithosphere age info every timestep - naming is different*/
    age=find_age_in_MY(E);
    sprintf(output_file,"%s%0.0f",E->control.lith_age_file,age);
  }
  else {     /* just open lithosphere age info here*/
    sprintf(output_file,"%s",E->control.lith_age_file);
  }

  if(E->parallel.me==0)  {
    fprintf(E->fp,"%s %s\n","Initial Lithosphere age info:",output_file);
  }

  fp1=fopen(output_file,"r");
  if (fp1 == NULL) {
    fprintf(E->fp,"(Convection.c #2) Cannot open %s\n",output_file);
    parallel_process_termination();
  }

  noy=E->lmesh.noy;
  nox=E->lmesh.nox;
  noz=E->lmesh.noz;

  gnox=E->mesh.nox;
  gnoy=E->mesh.noy;
  gnoz=E->mesh.noz;
  for(i=1;i<=gnoy;i++)
    for(j=1;j<=gnox;j++) {
      node=j+(i-1)*gnox;
      fscanf(fp1,"%f",&(E->age_t[node]));
      E->age_t[node]=E->age_t[node]/E->data.scalet;
    }
  fclose(fp1);

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=noy;i++)
      for(j=1;j<=nox;j++)
	for(k=1;k<=noz;k++)  {
	  nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	  node=k+(j-1)*noz+(i-1)*nox*noz;
	  r1=E->sx[m][3][node];
	  if( r1 >= E->sphere.ro-E->control.lith_age_depth )
	    { /* if closer than (lith_age_depth) from top */
	      temp = (E->sphere.ro-r1) *0.5 /sqrt(E->age_t[nodeg]);
	      E->T[m][node] = E->control.lith_age_mantle_temp * erf(temp);
	    }
	}

  /* modify temperature BC to be concorded with restarted T */
  lith_age_restart_conform_tbc(E);

  return;
}


void lith_age_restart_conform_tbc(struct All_variables *E)
{
  int i, j, k, m, node;
  int nox, noy, noz;
  double r1, rout, rin;
  const float e_4=1.e-4;

  noy = E->lmesh.noy;
  nox = E->lmesh.nox;
  noz = E->lmesh.noz;
  rout = E->sphere.ro;
  rin = E->sphere.ri;

  for(m=1;m<=E->sphere.caps_per_proc;m++)
    for(i=1;i<=noy;i++)
      for(j=1;j<=nox;j++)
	for(k=1;k<=noz;k++)  {
	  node=k+(j-1)*noz+(i-1)*nox*noz;
	  r1=E->sx[m][3][node];

	  if(fabs(r1-rout)>=e_4 && fabs(r1-rin)>=e_4)  {
	    E->sphere.cap[m].TB[1][node]=E->T[m][node];
	    E->sphere.cap[m].TB[2][node]=E->T[m][node];
	    E->sphere.cap[m].TB[3][node]=E->T[m][node];
	  }
	}

  return;
}



void lith_age_construct_tic(struct All_variables *E)
{
  lith_age_restart_tic(E);
  return;
}


void lith_age_read_files(struct All_variables *E, int output)
{
  FILE *fp1, *fp2;
  float age, newage1, newage2;
  char output_file1[255],output_file2[255];
  float inputage1, inputage2;
  int nox,noz,noy,nox1,noz1,noy1,lev;
  int i,j,node;
  int intage, pos_age;

  pos_age = 1;

  nox=E->mesh.nox;
  noy=E->mesh.noy;
  noz=E->mesh.noz;
  nox1=E->lmesh.nox;
  noz1=E->lmesh.noz;
  noy1=E->lmesh.noy;
  lev=E->mesh.levmax;

  age=find_age_in_MY(E);
  intage = (int) age;
  newage1 = 1.0*intage;
  newage2 = 1.0*intage + 1.0;
  if (newage1 < 0.0) {
    /* age is negative -> use age=0 for input files */
    newage1 = 0.0;
    pos_age = 0;
  }

  /* read ages for lithosphere tempperature boundary conditions */
  sprintf(output_file1,"%s%0.0f",E->control.lith_age_file,newage1);
  sprintf(output_file2,"%s%0.0f",E->control.lith_age_file,newage2);
  fp1=fopen(output_file1,"r");
  if (fp1 == NULL) {
    fprintf(E->fp,"(Problem_related #6) Cannot open %s\n",output_file1);
    parallel_process_termination();
  }
  if (pos_age) {
    fp2=fopen(output_file2,"r");
    if (fp2 == NULL) {
      fprintf(E->fp,"(Problem_related #7) Cannot open %s\n",output_file2);
    parallel_process_termination();
    }
  }
  if((E->parallel.me==0) && output) {
    fprintf(E->fp,"Age: Starting Age = %g, Elapsed time = %g, Current Age = %g\n",E->control.start_age,E->monitor.elapsed_time,age);
    fprintf(E->fp,"Age: File1 = %s\n",output_file1);
    if (pos_age)
      fprintf(E->fp,"Age: File2 = %s\n",output_file2);
    else
      fprintf(E->fp,"Age: File2 = No file inputted (negative age)\n");
  }

  for(i=1;i<=noy;i++)
    for(j=1;j<=nox;j++) {
      node=j+(i-1)*nox;
      fscanf(fp1,"%f",&inputage1);
      if (pos_age) { /* positive ages - we must interpolate */
	fscanf(fp2,"%f",&inputage2);
	E->age_t[node] = (inputage1 + (inputage2-inputage1)/(newage2-newage1)*(age-newage1))/E->data.scalet;
      }
      else { /* negative ages - don't do the interpolation */
	E->age_t[node] = inputage1;
      }
    }
  fclose(fp1);
  if (pos_age) fclose(fp2);

  return;
}


void lith_age_temperature_bound_adj(struct All_variables *E, int lv)
{
  int j,node,nno;
  float ttt2,ttt3,fff2,fff3;

  nno=E->lmesh.nno;

  ttt2=E->control.theta_min + E->control.width_bound_adj;
  ttt3=E->control.theta_max - E->control.width_bound_adj;
  fff2=E->control.fi_min + E->control.width_bound_adj;
  fff3=E->control.fi_max - E->control.width_bound_adj;


/* NOTE: To start, the relevent bits of "node" are zero. Thus, they only
get set to TBX/TBY/TBZ if the node is in one of the bounding regions.
Also note that right now, no matter which bounding region you are in,
all three get set to true. CPC 6/20/00 */

  if (E->control.temperature_bound_adj) {
    if(lv==E->mesh.gridmax)
      for(j=1;j<=E->sphere.caps_per_proc;j++)
	for(node=1;node<=E->lmesh.nno;node++)  {
	  if( ((E->sx[j][1][node]<=ttt2) && (E->sx[j][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) || ((E->sx[j][1][node]>=ttt3) && (E->sx[j][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) )
	    /* if < (width) from x bounds AND (depth) from top */
	    {
	      E->node[j][node]=E->node[j][node] | TBX;
	      E->node[j][node]=E->node[j][node] & (~FBX);
	      E->node[j][node]=E->node[j][node] | TBY;
	      E->node[j][node]=E->node[j][node] & (~FBY);
	      E->node[j][node]=E->node[j][node] | TBZ;
	      E->node[j][node]=E->node[j][node] & (~FBZ);
	    }

	  if( ((E->sx[j][2][node]<=fff2) && (E->sx[j][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) )
	    /* if fi is < (width) from side AND z is < (depth) from top */
	    {
	      E->node[j][node]=E->node[j][node] | TBX;
	      E->node[j][node]=E->node[j][node] & (~FBX);
	      E->node[j][node]=E->node[j][node] | TBY;
	      E->node[j][node]=E->node[j][node] & (~FBY);
	      E->node[j][node]=E->node[j][node] | TBZ;
	      E->node[j][node]=E->node[j][node] & (~FBZ);
	    }

	  if( ((E->sx[j][2][node]>=fff3) && (E->sx[j][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) )
	    /* if fi is < (width) from side AND z is < (depth) from top */
	    {
	      E->node[j][node]=E->node[j][node] | TBX;
	      E->node[j][node]=E->node[j][node] & (~FBX);
	      E->node[j][node]=E->node[j][node] | TBY;
	      E->node[j][node]=E->node[j][node] & (~FBY);
	      E->node[j][node]=E->node[j][node] | TBZ;
	      E->node[j][node]=E->node[j][node] & (~FBZ);
	    }

	}
  } /* end E->control.temperature_bound_adj */

  if (E->control.lith_age_time) {
    if(lv==E->mesh.gridmax)
      for(j=1;j<=E->sphere.caps_per_proc;j++)
	for(node=1;node<=E->lmesh.nno;node++)  {
	  if(E->sx[j][3][node]>=E->sphere.ro-E->control.lith_age_depth)
	    { /* if closer than (lith_age_depth) from top */
	      E->node[j][node]=E->node[j][node] | TBX;
	      E->node[j][node]=E->node[j][node] & (~FBX);
	      E->node[j][node]=E->node[j][node] | TBY;
	      E->node[j][node]=E->node[j][node] & (~FBY);
	      E->node[j][node]=E->node[j][node] | TBZ;
	      E->node[j][node]=E->node[j][node] & (~FBZ);
	    }

	}
  } /* end E->control.lith_age_time */

  return;
}


void lith_age_conform_tbc(struct All_variables *E)
{
  int m,j,node,nox,noz,noy,gnox,gnoy,gnoz,nodeg,i,k;
  float ttt2,ttt3,fff2,fff3;
  float r1,t1,f1,t0,temp;
  float e_4;
  FILE *fp1;
  char output_file[255];
  int output;

  e_4=1.e-4;
  output = 0;

  ttt2=E->control.theta_min + E->control.width_bound_adj;
  ttt3=E->control.theta_max - E->control.width_bound_adj;
  fff2=E->control.fi_min + E->control.width_bound_adj;
  fff3=E->control.fi_max - E->control.width_bound_adj;

  gnox=E->mesh.nox;
  gnoy=E->mesh.noy;
  gnoz=E->mesh.noz;
  nox=E->lmesh.nox;
  noy=E->lmesh.noy;
  noz=E->lmesh.noz;

  if(E->control.lith_age_time==1)   {
    /* to open files every timestep */
    if(E->monitor.solution_cycles==0) {
      E->age_t=(float*) malloc((gnox*gnoy+1)*sizeof(float));
      E->control.lith_age_old_cycles = E->monitor.solution_cycles;
    }
    if (E->control.lith_age_old_cycles != E->monitor.solution_cycles) {
      /*update so that output only happens once*/
      output = 1;
      E->control.lith_age_old_cycles = E->monitor.solution_cycles;
    }
    lith_age_read_files(E,output);
  }
  else {
    /* otherwise, just open for the first timestep */
    /* NOTE: This is only used if we are adjusting the boundaries */

    if(E->monitor.solution_cycles==0) {
      E->age_t=(float*) malloc((gnox*gnoy+1)*sizeof(float));
      sprintf(output_file,"%s",E->control.lith_age_file);
      fp1=fopen(output_file,"r");
      if (fp1 == NULL) {
	fprintf(E->fp,"(Boundary_conditions #1) Can't open %s\n",output_file);
	parallel_process_termination();
      }
      for(i=1;i<=gnoy;i++)
	for(j=1;j<=gnox;j++) {
	  node=j+(i-1)*gnox;
	  fscanf(fp1,"%f",&(E->age_t[node]));
	  E->age_t[node]=E->age_t[node]*E->data.scalet;
	}
      fclose(fp1);
    } /* end E->monitor.solution_cycles == 0 */
  } /* end E->control.lith_age_time == false */


  /* NOW SET THE TEMPERATURES IN THE BOUNDARY REGIONS */
  if(E->monitor.solution_cycles>1 && E->control.temperature_bound_adj) {
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=noy;i++)
	for(j=1;j<=nox;j++)
	  for(k=1;k<=noz;k++)  {
	    nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	    node=k+(j-1)*noz+(i-1)*nox*noz;
	    t1=E->sx[m][1][node];
	    f1=E->sx[m][2][node];
	    r1=E->sx[m][3][node];

	    if(fabs(r1-E->sphere.ro)>=e_4 && fabs(r1-E->sphere.ri)>=e_4)  { /* if NOT right on the boundary */
	      if( ((E->sx[m][1][node]<=ttt2) && (E->sx[m][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) || ((E->sx[m][1][node]>=ttt3) && (E->sx[m][3][node]>=E->sphere.ro-E->control.depth_bound_adj)) ) {
		/* if < (width) from x bounds AND (depth) from top */
		temp = (E->sphere.ro-r1) *0.5 /sqrt(E->age_t[nodeg]);
		t0 = E->control.lith_age_mantle_temp * erf(temp);

		/* keep the age the same! */
		E->sphere.cap[m].TB[1][node]=t0;
		E->sphere.cap[m].TB[2][node]=t0;
		E->sphere.cap[m].TB[3][node]=t0;
	      }

	      if( ((E->sx[m][2][node]<=fff2) || (E->sx[m][2][node]>=fff3)) && (E->sx[m][3][node]>=E->sphere.ro-E->control.depth_bound_adj) ) {
		/* if < (width) from y bounds AND (depth) from top */

		/* keep the age the same! */
		temp = (E->sphere.ro-r1) *0.5 /sqrt(E->age_t[nodeg]);
		t0 = E->control.lith_age_mantle_temp * erf(temp);

		E->sphere.cap[m].TB[1][node]=t0;
		E->sphere.cap[m].TB[2][node]=t0;
		E->sphere.cap[m].TB[3][node]=t0;

	      }

	    }

	  } /* end k   */

  }   /*  end of solution cycles  && temperature_bound_adj */


  /* NOW SET THE TEMPERATURES IN THE LITHOSPHERE IF CHANGING EVERY TIME STEP */
  if(E->monitor.solution_cycles>1 && E->control.lith_age_time)   {
    for(m=1;m<=E->sphere.caps_per_proc;m++)
      for(i=1;i<=noy;i++)
	for(j=1;j<=nox;j++)
	  for(k=1;k<=noz;k++)  {
	    nodeg=E->lmesh.nxs-1+j+(E->lmesh.nys+i-2)*gnox;
	    node=k+(j-1)*noz+(i-1)*nox*noz;
	    t1=E->sx[m][1][node];
	    f1=E->sx[m][2][node];
	    r1=E->sx[m][3][node];

	    if(fabs(r1-E->sphere.ro)>=e_4 && fabs(r1-E->sphere.ri)>=e_4)  { /* if NOT right on the boundary */
	      if(  E->sx[m][3][node]>=E->sphere.ro-E->control.lith_age_depth ) {
		/* if closer than (lith_age_depth) from top */

		/* set a new age from the file */
		temp = (E->sphere.ro-r1) *0.5 /sqrt(E->age_t[nodeg]);
		t0 = E->control.lith_age_mantle_temp * erf(temp);

		E->sphere.cap[m].TB[1][node]=t0;
		E->sphere.cap[m].TB[2][node]=t0;
		E->sphere.cap[m].TB[3][node]=t0;
	      }
	    }
	  }     /* end k   */
  }   /*  end of solution cycles  && lith_age_time */

  return;
}
