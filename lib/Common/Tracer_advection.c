/*

      Tracer_advection.c

      A program which initiates the distribution of tracers
      and advects those tracers in a time evolving velocity field.

      Called and used from the CitCOM finite element code.

      Written 2/96 M. Gurnis for Citcom in cartesian geometry

      Modified by Lijie in 1998 and by Vlad and Eh in (March 23) 2005 for CitcomS
*/



#include <mpi.h>
#include <math.h>
#include <sys/types.h>
#include <malloc.h>
#include <stdlib.h> /* for "system" command */

#include "element_definitions.h"
#include "global_defs.h"
#include "parsing.h"

void tracer_input(struct All_variables *E)
{
  int m=E->parallel.me;

  input_int("tracer",&(E->control.tracer),"0",m);
  input_string("tracer_file",E->control.tracer_file,"tracer.dat",m);
}


void tracer_initial_settings(E)
   struct All_variables *E;
 {
   void tracer_setup();
   void tracer_advection();

   E->problem_tracer_setup=tracer_setup;
   E->problem_tracer_advection=tracer_advection;
}


void tracer_setup(E)
    struct All_variables *E;
{
	int i,j,k,node;
	int m,c,nn;
	float idummy,xdummy,ydummy,zdummy;
        float xmin,xmax,ymin,ymax,zmin,zmax;
	FILE *fp,*fp1;
        int nox,noy,noz,gnox,gnoy,gnoz;
        char output_file[255];
	char aaa[100];

        gnox=E->mesh.nox;
        gnoy=E->mesh.noy;
        gnoz=E->mesh.noz;
        nox=E->lmesh.nox;
        noy=E->lmesh.noy;
        noz=E->lmesh.noz;

        sprintf(output_file,"%s",E->control.tracer_file);
       fp=fopen(output_file,"r");
	if (fp == NULL) {
          fprintf(E->fp,"(Tracer_advection #1) Cannot open %s\n",output_file);
          exit(8);
	}
         fscanf(fp,"%d",&(E->Tracer.NUM_TRACERS));

        E->Tracer.tracer_x=(float*) malloc((E->Tracer.NUM_TRACERS+1)*sizeof(float));
        E->Tracer.tracer_y=(float*) malloc((E->Tracer.NUM_TRACERS+1)*sizeof(float));
        E->Tracer.tracer_z=(float*) malloc((E->Tracer.NUM_TRACERS+1)*sizeof(float));
        E->Tracer.itcolor=(float*) malloc((E->Tracer.NUM_TRACERS+1)*sizeof(float));

        E->Tracer.x_space=(float*) malloc((gnox+1)*sizeof(float));
        E->Tracer.y_space=(float*) malloc((gnoy+1)*sizeof(float));
        E->Tracer.z_space=(float*) malloc((gnoz+1)*sizeof(float));


	/***comment by Vlad 03/08/2005
	checking the coor file
	***/

       fp1=fopen("coor.dat","r");
        if (fp1 == NULL) {
          c = 0;
                         }
        else             {
          c = 1;
                         }



        /***comment by Vlad 1/26/2005
        reading the mesh coordinate
        ***/
	
	if(E->parallel.me == 0) {

        if(c == 0) {

        for(m=1;m<=E->sphere.caps_per_proc;m++)  {

        for(i=1;i<=gnox;i++)
	   {
	     xmin = E->sphere.cap[m].theta[1];
	     xmax = E->sphere.cap[m].theta[3];
             E->Tracer.x_space[i]=xmin + (i-1)*(xmax-xmin)/(gnox-1);
	   }

        for(j=1;j<=gnoy;j++)
	   {
             ymin = E->sphere.cap[m].fi[1];
             ymax = E->sphere.cap[m].fi[3];
             E->Tracer.y_space[j]=ymin + (j-1)*(ymax-ymin)/(gnoy-1);
	   }

        for(k=1;k<=gnoz;k++)
	   {
	     zmin=E->sphere.ri;
             zmax=E->sphere.ro;
             E->Tracer.z_space[k]=zmin + (k-1)*(zmax-zmin)/(gnoz-1);
	   }   /* end of reading the coord */

           }   /* end of m  */

	   }   /* end of c */


        else {

        for(m=1;m<=E->sphere.caps_per_proc;m++)  {

        fscanf(fp1, "%s %d", aaa, &nn);
        for(i=1;i<=gnox;i++)  {
        fscanf(fp1, "%d %f", &nn, &(E->Tracer.x_space[i]));
                              }
        fscanf(fp1, "%s %d", aaa, &nn);
        for(j=1;j<=gnoy;j++)  {
        fscanf(fp1, "%d %f", &nn, &(E->Tracer.y_space[j]));
                              }
        fscanf(fp1, "%s %d", aaa, &nn);
        for(k=1;k<=gnoz;k++)  {
        fscanf(fp1, "%d %f", &nn, &(E->Tracer.z_space[k]));
                              }

                                                }
             }


	}

           for(i=1;i<=E->Tracer.NUM_TRACERS;i++)
	     {
	          fscanf(fp,"%f %f %f %f",
	          &idummy, &xdummy, &ydummy, &zdummy);
	          E->Tracer.itcolor[i]=idummy; E->Tracer.tracer_x[i]=xdummy;
		  E->Tracer.tracer_y[i]=ydummy;E->Tracer.tracer_z[i]=zdummy;

              }

return;

}

void tracer_advection(E)
    struct All_variables *E;
    {
      int n,loc;
      int i,j,k,m,ii,jj,kk;
      int n_x,n_y,n_z;
      int node1,node2,node3,node4,node5,node6,node7,node8;
      int nno,nox,noy,noz,gnox,gnoy,gnoz,nodeg,nodel;
      int nxs,nys,nzs,proc,nprocxl,nprocyl,nproczl;
      int iteration;
      int root;
      float x_tmp, y_tmp, z_tmp;
      float volume, tr_dx, tr_dy, tr_dz, dx, dy, dz;
      float w1,w2,w3,w4,w5,w6,w7,w8;
      float tr_v[NCS][4];
      MPI_Comm world;
      float xmin,xmax,ymin,ymax,zmin,zmax;


      world=E->parallel.world;

      gnox=E->mesh.nox;
      gnoy=E->mesh.noy;
      gnoz=E->mesh.noz;
      nox=E->lmesh.nox;
      noy=E->lmesh.noy;
      noz=E->lmesh.noz;
      nno=nox*noy*noz;


      /*** comment by Tan2 1/25/2005
	   This part of code find the bounding box of the global mesh.
	   Not checked.
      ***/

       for(m=1;m<=E->sphere.caps_per_proc;m++)  {

       if(E->sphere.cap[m].theta[1]<=E->sphere.cap[m].theta[3]) {
         xmin=E->sphere.cap[m].theta[1];
         xmax=E->sphere.cap[m].theta[3];
        }
       else {
         xmin=E->sphere.cap[m].theta[3];
         xmax=E->sphere.cap[m].theta[1];
        }

       if(E->sphere.cap[m].fi[1]<=E->sphere.cap[m].fi[3]) {
         ymin=E->sphere.cap[m].fi[1];
         ymax=E->sphere.cap[m].fi[3];

        }
       else {
         ymin=E->sphere.cap[m].fi[3];
         ymax=E->sphere.cap[m].fi[1];
        }

         zmin=E->sphere.ri;
         zmax=E->sphere.ro;

     }   /* end of m */
      
      /*** comment by Tan2 1/25/2005
           Copy the velocity array. Is it necessary?
      ***/


      for(m=1;m<=E->sphere.caps_per_proc;m++)   {
      for(i=0;i<=nno-1;i++)
         for(j=1;j<=3;j++)   {
           E->V[m][j][i]=E->sphere.cap[m].V[j][i+1];

        }
      }



    root=0;


      /*** comment by Tan2 1/25/2005
           Gather the global velocity to root processor from all processors.
           Seems the code is designed for multi-processor run, although the
           tracer advection only runs in root processor.
      ***/


      for(m=1;m<=E->sphere.caps_per_proc;m++)   {

      MPI_Gather(E->V[m][1], nno, MPI_FLOAT, E->GV1[m][1], nno, MPI_FLOAT, root, world);
      MPI_Gather(E->V[m][2], nno, MPI_FLOAT, E->GV1[m][2], nno, MPI_FLOAT, root, world);
      MPI_Gather(E->V[m][3], nno, MPI_FLOAT, E->GV1[m][3], nno, MPI_FLOAT, root, world);

      }

    loc=0;
    nprocxl=E->parallel.nprocx;
    nprocyl=E->parallel.nprocy;
    nproczl=E->parallel.nprocz;


      /*** comment by Tan2 1/25/2005
	   Combine the gathered velocity to get the global velocity field.
      ***/

 if(E->parallel.me==root)    {
 for(m=1;m<=E->sphere.caps_per_proc;m++)
  for (k=1;k<=nprocyl;k++)
    for (j=1;j<=nprocxl;j++)
      for (i=1;i<=nproczl;i++)    {

        proc = i-1 + (j-1)*nproczl+(k-1)*nprocxl*nproczl;

        nxs = (j-1)*(nox-1) + 1;
        nys = (k-1)*(noy-1) + 1;
        nzs = (i-1)*(noz-1) + 1;

        for (kk=1;kk<=noy;kk++)
          for (jj=1;jj<=nox;jj++)
            for (ii=1;ii<=noz;ii++)              {

              nodeg = ii + nzs -1 + gnoz*
                  (nxs+jj-2)
                  + (nys+kk-2)*gnox*gnoz;

              nodel = (kk-1)*nox*noz + (jj-1)*noz + ii;

              E->GV[m][1][nodeg] = E->GV1[m][1][nodel+loc-1];
              E->GV[m][2][nodeg] = E->GV1[m][2][nodel+loc-1];
              E->GV[m][3][nodeg] = E->GV1[m][3][nodel+loc-1];
              
	      }

        loc += nno;


        }

      }


      /*** comment by Tan2 1/25/2005
	   Now, root processor has velocity field of the whole mesh and starts
	   to advect the tracers.
      ***/

      if(E->parallel.me==root)   {
      for(n=1;n<=E->Tracer.NUM_TRACERS;n++) {
      n_x=0;
      n_y=0;
      n_z=0;

      /*** comment by Tan2 1/25/2005
	   itcolor is like the starting age of the tracer. Disable it temporarly.
      ***/

      if(1)   {

     /*  mid point method uses 2 iterations */
      for(iteration=1;iteration<=2;iteration++)
      {

         if(iteration ==1) {
	    x_tmp=E->Tracer.tracer_x[n];
	    y_tmp=E->Tracer.tracer_y[n];
	    z_tmp=E->Tracer.tracer_z[n];

	 }

	 /*** comment by Tan2 1/25/2005
	      Find the element that contains the tracer.

      nodex   n_x                 n_x+1
	      |           *        |
	      <----------->
	          tr_dx

	      <-------------------->
	               dx
	  ***/

	 for(i=1;i<gnox;i++) {
             if(x_tmp >= E->Tracer.x_space[i] && x_tmp <= E->Tracer.x_space[i+1]) {
		tr_dx=x_tmp-E->Tracer.x_space[i];
		dx=E->Tracer.x_space[i+1]-E->Tracer.x_space[i];
		n_x=i;
             }
	 }
	 for(j=1;j<gnoz;j++) {
             if(z_tmp >= E->Tracer.z_space[j] && z_tmp <= E->Tracer.z_space[j+1]) {
		tr_dz=z_tmp-E->Tracer.z_space[j];
		dz=E->Tracer.z_space[j+1]-E->Tracer.z_space[j];
	        n_z=j;
	     }
	 }
	 for(k=1;k<gnoy;k++) {
             if(y_tmp >= E->Tracer.y_space[k] && y_tmp <= E->Tracer.y_space[k+1]) {
		tr_dy=y_tmp-E->Tracer.y_space[k];
		dy=E->Tracer.y_space[k+1]-E->Tracer.y_space[k];
	        n_y=k;

	     }
	 }

	 /*** comment by Tan2 1/25/2005
	      Calculate shape functions from tr_dx, tr_dy, tr_dz
	      This assumes linear element
	  ***/


          /* compute volumetic weighting functions */
	    w1=tr_dx*tr_dz*tr_dy;
	    w2=(dx-tr_dx)*tr_dz*tr_dy;
	    w3=tr_dx*(dz-tr_dz)*tr_dy;
	    w4=(dx-tr_dx)*(dz-tr_dz)*tr_dy;
	    w5=tr_dx*tr_dz*(dy-tr_dy);
	    w6=(dx-tr_dx)*tr_dz*(dy-tr_dy);
	    w7=tr_dx*(dz-tr_dz)*(dy-tr_dy);
	    w8=(dx-tr_dx)*(dz-tr_dz)*(dy-tr_dy);

	    volume=dx*dz*dy;

 
	/*** comment by Tan2 1/25/2005
		 Calculate the 8 node numbers of current element
	     ***/

	 node1 = n_z + (n_x-1)*gnoz + (n_y-1)*gnoz*gnox;
	 node2 = n_z + n_x*gnoz + (n_y-1)*gnoz*gnox;
	 node3 = n_z+1  + (n_x-1)*gnoz + (n_y-1)*gnoz*gnox;
	 node4 = n_z+1  + n_x*gnoz + (n_y-1)*gnoz*gnox;
	 node5 = n_z + (n_x-1)*gnoz + n_y*gnoz*gnox;
	 node6 = n_z + n_x*gnoz + n_y*gnoz*gnox;
	 node7 = n_z+1 + (n_x-1)*gnoz + n_y*gnoz*gnox;
	 node8 = n_z+1 + n_x*gnoz + n_y*gnoz*gnox;
	

 //       printf("%d %d %d %d %d %d %d %d %d\n",E->parallel.me, node1,node2,node3,node4,node5,node6,node7,node8);
   //     printf("%d %f %f %f %f %f %f %f %f\n", E->parallel.me, E->GV[1][2][node1], E->GV[1][2][node2], E->GV[1][2][node3], E->GV[1][2][node4], E->GV[1][2][node5], E->GV[1][2][node6], E->GV[1][2][node7], E->GV[1][2][node8]);


	 /*** comment by Tan2 1/25/2005
	      Interpolate the velocity on the tracer position
	  ***/

            for(m=1;m<=E->sphere.caps_per_proc;m++)   {
              for(j=1;j<=3;j++)   {
	     tr_v[m][j]=w8*E->GV[m][j][node1]
		    +w7*E->GV[m][j][node2]
		    +w6*E->GV[m][j][node3]
		    +w5*E->GV[m][j][node4]
		    +w4*E->GV[m][j][node5]
		    +w3*E->GV[m][j][node6]
		    +w2*E->GV[m][j][node7]
		    +w1*E->GV[m][j][node8];
              tr_v[m][j]=tr_v[m][j]/volume;

	    }
                                                                                                                                                         
	      /*** comment by Tan2 1/25/2005
		   advect tracer using mid-point method (2nd order accuracy)
	       ***/

          /* mid point method */

	  if(iteration == 1) {
	     x_tmp = x_tmp + (E->advection.timestep/2.0)*tr_v[m][1]/E->Tracer.z_space[n_z];
	     y_tmp = y_tmp + (E->advection.timestep/2.0)*tr_v[m][2]/(E->Tracer.z_space[n_z]*sin(E->Tracer.x_space[n_x]));
	     z_tmp = z_tmp + (E->advection.timestep/2.0)*tr_v[m][3];
            }
	  if( iteration == 2) {
	     E->Tracer.tracer_x[n] += E->advection.timestep*tr_v[m][1]/E->Tracer.z_space[n_z];
	     E->Tracer.tracer_y[n] += E->advection.timestep*tr_v[m][2]/(E->Tracer.z_space[n_z]*sin(E->Tracer.x_space[n_x]));
	     E->Tracer.tracer_z[n] += E->advection.timestep*tr_v[m][3];

	  }


      }   /*  end of m  */


      } /* end of iteration loop */


      } /* if E->Tracer.itcolor   */

      } /* end of tracer loop */

     }



/* make sure that the new position of the tracers within the box  */

       for(n=1;n<=E->Tracer.NUM_TRACERS;n++) {
        if(E->Tracer.tracer_x[n]>xmax) E->Tracer.tracer_x[n]=2.0*xmax-E->Tracer.tracer_x[n];
        if(E->Tracer.tracer_x[n]<xmin) E->Tracer.tracer_x[n]=2.0*xmin-E->Tracer.tracer_x[n];
        if(E->Tracer.tracer_y[n]>ymax) E->Tracer.tracer_y[n]=2.0*ymax-E->Tracer.tracer_y[n];
        if(E->Tracer.tracer_y[n]<ymin) E->Tracer.tracer_y[n]=2.0*ymin-E->Tracer.tracer_y[n];
        if(E->Tracer.tracer_z[n]>zmax) E->Tracer.tracer_z[n]=2.0*zmax-E->Tracer.tracer_z[n];
        if(E->Tracer.tracer_z[n]<zmin) E->Tracer.tracer_z[n]=2.0*zmin-E->Tracer.tracer_z[n];
      }

return;
    }
