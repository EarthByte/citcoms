/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/*

  Tracer_setup.c

      A program which initiates the distribution of tracers
      and advects those tracers in a time evolving velocity field.
      Called and used from the CitCOM finite element code.
      Written 2/96 M. Gurnis for Citcom in cartesian geometry
      Modified by Lijie in 1998 and by Vlad and Eh in 2005 for the
      regional version of CitcomS. In 2003, Allen McNamara wrote the
      tracer module for the global version of CitcomS. In 2007, Eh Tan
      merged the two versions of tracer codes together.
*/

#include <math.h>
#include "global_defs.h"
#include "parsing.h"

void tracer_post_processing(struct All_variables *E);


static void predict_tracers(struct All_variables *E);
static void correct_tracers(struct All_variables *E);
void count_tracers_of_flavors(struct All_variables *E);
static void check_sum(struct All_variables *E);
int isum_tracers(struct All_variables *E);
static void put_lost_tracers(struct All_variables *E,
                             int isend[13][13], double *send[13][13]);
int icheck_that_processor_shell(struct All_variables *E,
                                int j, int nprocessor, double rad);


void tracer_input(struct All_variables *E)
{
    void full_tracer_input();
    int m=E->parallel.me;

    input_int("tracer",&(E->control.tracer),"0",m);
    if(E->control.tracer) {

        /* Initial condition, this option is ignored if E->control.restart is 1,
         *  ie. restarted from a previous run */
        /* tracer_ic_method=0 (random generated array) */
        /* tracer_ic_method=1 (all proc read the same file) */
        /* tracer_ic_method=2 (each proc reads its restart file) */
        if(E->control.restart)
            E->trace.ic_method = 2;
        else {
            input_int("tracer_ic_method",&(E->trace.ic_method),"0,0,nomax",m);

            if (E->trace.ic_method==0)
                input_int("tracers_per_element",&(E->trace.itperel),"10,0,nomax",m);
            else if (E->trace.ic_method==1)
                input_string("tracer_file",E->trace.tracer_file,"tracer.dat",m);
            else if (E->trace.ic_method==2) {
            }
            else {
                fprintf(stderr,"Sorry, tracer_ic_method only 0, 1 and 2 available\n");
                fflush(stderr);
                parallel_process_termination();
            }
        }


        /* How many flavors of tracers */
        /* If tracer_flavors > 0, each element will report the number of
         * tracers of each flavor inside it. This information can be used
         * later for many purposes. One of it is to compute composition,
         * either using absolute method or ratio method. */
        input_int("tracer_flavors",&(E->trace.nflavors),"0,0,nomax",m);


        input_int("ic_method_for_flavors",&(E->trace.ic_method_for_flavors),"0,0,nomax",m);
        if (E->trace.ic_method_for_flavors == 0)
            input_double("z_interface",&(E->trace.z_interface),"0.7",m);


        /* Advection Scheme */

        /* itracer_advection_scheme=1 (simple predictor corrector -uses only V(to)) */
        /* itracer_advection_scheme=2 (predictor-corrector - uses V(to) and V(to+dt)) */

        E->trace.itracer_advection_scheme=2;
        input_int("tracer_advection_scheme",&(E->trace.itracer_advection_scheme),
                  "2,0,nomax",m);

        if (E->trace.itracer_advection_scheme==1)
            {}
        else if (E->trace.itracer_advection_scheme==2)
            {}
        else {
            fprintf(stderr,"Sorry, only advection scheme 1 and 2 available (%d)\n",E->trace.itracer_advection_scheme);
            fflush(stderr);
            parallel_process_termination();
        }




        if(E->parallel.nprocxy == 12)
            full_tracer_input(E);


        composition_input(E);
    }

    return;
}


void tracer_initial_settings(struct All_variables *E)
{
   void full_tracer_setup();
   void full_get_velocity();
   int full_iget_element();
   void regional_tracer_setup();
   void regional_get_velocity();
   int regional_iget_element();

   if(E->parallel.nprocxy == 1) {
       E->problem_tracer_setup = regional_tracer_setup;

       E->trace.get_velocity = regional_get_velocity;
       E->trace.iget_element = regional_iget_element;
   }
   else {
       E->problem_tracer_setup = full_tracer_setup;

       E->trace.get_velocity = full_get_velocity;
       E->trace.iget_element = full_iget_element;
   }
}



/*****************************************************************************/
/* This function is the primary tracing routine called from Citcom.c         */
/* In this code, unlike the original 3D cartesian code, force is filled      */
/* during Stokes solution. No need to call thermal_buoyancy() after tracing. */


void tracer_advection(struct All_variables *E)
{

    fprintf(E->trace.fpt,"STEP %d\n",E->monitor.solution_cycles);

    /* advect tracers */
    predict_tracers(E);
    correct_tracers(E);

    /* check that the number of tracers is conserved */
    check_sum(E);

    /* count # of tracers of each flavor */
    if (E->trace.nflavors > 0)
        count_tracers_of_flavors(E);

    /* update the composition field */
    if (E->composition.on) {
        fill_composition(E);
    }

    tracer_post_processing(E);

    return;
}



/********* TRACER POST PROCESSING ****************************************/

void tracer_post_processing(struct All_variables *E)
{
    void get_bulk_composition();
    char output_file[200];

    double convection_time,tracer_time;
    double trace_fraction,total_time;


    fprintf(E->trace.fpt,"Number of times for all element search  %d\n",E->trace.istat1);

    fprintf(E->trace.fpt,"Number of tracers sent to other processors: %d\n",E->trace.istat_isend);

    fprintf(E->trace.fpt,"Number of times element columns are checked: %d \n",E->trace.istat_elements_checked);

    if (E->composition.ichemical_buoyancy==1) {
        fprintf(E->trace.fpt,"Empty elements filled with old compositional values: %d (%f percent)\n",
                E->trace.istat_iempty,(100.0*E->trace.istat_iempty)/E->lmesh.nel);
    }


    /* reset statistical counters */

    E->trace.istat_isend=0;
    E->trace.istat_iempty=0;
    E->trace.istat_elements_checked=0;
    E->trace.istat1=0;

    /* compositional and error fraction data files */
    //TODO: move
    if (E->composition.ichemical_buoyancy==1) {
        get_bulk_composition(E);
        if (E->parallel.me==0) {
            fprintf(E->fp,"composition: %e %e\n",E->monitor.elapsed_time,E->composition.bulk_composition);
            fprintf(E->fp,"composition_error_fraction: %e %e\n",E->monitor.elapsed_time,E->composition.error_fraction);

        }

    }

    fflush(E->trace.fpt);

    return;
}


/*********** PREDICT TRACERS **********************************************/
/*                                                                        */
/* This function predicts tracers performing an euler step                */
/*                                                                        */
/*                                                                        */
/* Note positions used in tracer array                                    */
/* [positions 0-5 are always fixed with current coordinates               */
/*  regardless of advection scheme].                                      */
/*  Positions 6-8 contain original Cartesian coordinates.                 */
/*  Positions 9-11 contain original Cartesian velocities.                 */
/*                                                                        */


static void predict_tracers(struct All_variables *E)
{

    int numtracers;
    int j;
    int kk;
    int nelem;

    double dt;
    double theta0,phi0,rad0;
    double x0,y0,z0;
    double theta_pred,phi_pred,rad_pred;
    double x_pred,y_pred,z_pred;
    double velocity_vector[4];

    void cart_to_sphere();
    void keep_in_sphere();
    void find_tracers();


    dt=E->advection.timestep;


    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        numtracers=E->trace.ntracers[j];

        for (kk=1;kk<=numtracers;kk++) {

            theta0=E->trace.basicq[j][0][kk];
            phi0=E->trace.basicq[j][1][kk];
            rad0=E->trace.basicq[j][2][kk];
            x0=E->trace.basicq[j][3][kk];
            y0=E->trace.basicq[j][4][kk];
            z0=E->trace.basicq[j][5][kk];

            nelem=E->trace.ielement[j][kk];
            (E->trace.get_velocity)(E,j,nelem,theta0,phi0,rad0,velocity_vector);

            x_pred=x0+velocity_vector[1]*dt;
            y_pred=y0+velocity_vector[2]*dt;
            z_pred=z0+velocity_vector[3]*dt;


            /* keep in box */

            cart_to_sphere(E,x_pred,y_pred,z_pred,&theta_pred,&phi_pred,&rad_pred);
            keep_in_sphere(E,&x_pred,&y_pred,&z_pred,&theta_pred,&phi_pred,&rad_pred);

            /* Current Coordinates are always kept in positions 0-5. */

            E->trace.basicq[j][0][kk]=theta_pred;
            E->trace.basicq[j][1][kk]=phi_pred;
            E->trace.basicq[j][2][kk]=rad_pred;
            E->trace.basicq[j][3][kk]=x_pred;
            E->trace.basicq[j][4][kk]=y_pred;
            E->trace.basicq[j][5][kk]=z_pred;

            /* Fill in original coords (positions 6-8) */

            E->trace.basicq[j][6][kk]=x0;
            E->trace.basicq[j][7][kk]=y0;
            E->trace.basicq[j][8][kk]=z0;

            /* Fill in original velocities (positions 9-11) */

            E->trace.basicq[j][9][kk]=velocity_vector[1];  /* Vx */
            E->trace.basicq[j][10][kk]=velocity_vector[2];  /* Vy */
            E->trace.basicq[j][11][kk]=velocity_vector[3];  /* Vz */


        } /* end kk, predicting tracers */
    } /* end caps */

    /* find new tracer elements and caps */

    find_tracers(E);

    return;

}


/*********** CORRECT TRACERS **********************************************/
/*                                                                        */
/* This function corrects tracers using both initial and                  */
/* predicted velocities                                                   */
/*                                                                        */
/*                                                                        */
/* Note positions used in tracer array                                    */
/* [positions 0-5 are always fixed with current coordinates               */
/*  regardless of advection scheme].                                      */
/*  Positions 6-8 contain original Cartesian coordinates.                 */
/*  Positions 9-11 contain original Cartesian velocities.                 */
/*                                                                        */


static void correct_tracers(struct All_variables *E)
{

    int j;
    int kk;
    int nelem;


    double dt;
    double x0,y0,z0;
    double theta_pred,phi_pred,rad_pred;
    double x_pred,y_pred,z_pred;
    double theta_cor,phi_cor,rad_cor;
    double x_cor,y_cor,z_cor;
    double velocity_vector[4];
    double Vx0,Vy0,Vz0;
    double Vx_pred,Vy_pred,Vz_pred;

    void cart_to_sphere();
    void keep_in_sphere();
    void find_tracers();


    dt=E->advection.timestep;


    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        for (kk=1;kk<=E->trace.ntracers[j];kk++) {

            theta_pred=E->trace.basicq[j][0][kk];
            phi_pred=E->trace.basicq[j][1][kk];
            rad_pred=E->trace.basicq[j][2][kk];
            x_pred=E->trace.basicq[j][3][kk];
            y_pred=E->trace.basicq[j][4][kk];
            z_pred=E->trace.basicq[j][5][kk];

            x0=E->trace.basicq[j][6][kk];
            y0=E->trace.basicq[j][7][kk];
            z0=E->trace.basicq[j][8][kk];

            Vx0=E->trace.basicq[j][9][kk];
            Vy0=E->trace.basicq[j][10][kk];
            Vz0=E->trace.basicq[j][11][kk];

            nelem=E->trace.ielement[j][kk];

            (E->trace.get_velocity)(E,j,nelem,theta_pred,phi_pred,rad_pred,velocity_vector);

            Vx_pred=velocity_vector[1];
            Vy_pred=velocity_vector[2];
            Vz_pred=velocity_vector[3];

            x_cor=x0 + dt * 0.5*(Vx0+Vx_pred);
            y_cor=y0 + dt * 0.5*(Vy0+Vy_pred);
            z_cor=z0 + dt * 0.5*(Vz0+Vz_pred);

            cart_to_sphere(E,x_cor,y_cor,z_cor,&theta_cor,&phi_cor,&rad_cor);
            keep_in_sphere(E,&x_cor,&y_cor,&z_cor,&theta_cor,&phi_cor,&rad_cor);

            /* Fill in Current Positions (other positions are no longer important) */

            E->trace.basicq[j][0][kk]=theta_cor;
            E->trace.basicq[j][1][kk]=phi_cor;
            E->trace.basicq[j][2][kk]=rad_cor;
            E->trace.basicq[j][3][kk]=x_cor;
            E->trace.basicq[j][4][kk]=y_cor;
            E->trace.basicq[j][5][kk]=z_cor;

        } /* end kk, correcting tracers */
    } /* end caps */

    /* find new tracer elements and caps */

    find_tracers(E);

    return;
}


/************ FIND TRACERS *************************************/
/*                                                             */
/* This function finds tracer elements and moves tracers to    */
/* other processor domains if necessary.                       */
/* Array ielement is filled with elemental values.                */

void find_tracers(struct All_variables *E)
{

    int iel;
    int kk;
    int j;
    int it;
    int iprevious_element;
    int num_tracers;

    double theta,phi,rad;
    double x,y,z;
    double time_stat1;
    double time_stat2;

    void put_away_later();
    void eject_tracer();
    void reduce_tracer_arrays();
    void lost_souls();
    void sphere_to_cart();

    time_stat1=CPU_time0();


    for (j=1;j<=E->sphere.caps_per_proc;j++) {


        /* initialize arrays and statistical counters */

        E->trace.ilater[j]=0;

        E->trace.istat1=0;
        for (kk=0;kk<=4;kk++) {
            E->trace.istat_ichoice[j][kk]=0;
        }

        //TODO: use while-loop instead of for-loop
        /* important to index by it, not kk */

        it=0;
        num_tracers=E->trace.ntracers[j];

        for (kk=1;kk<=num_tracers;kk++) {

            it++;

            theta=E->trace.basicq[j][0][it];
            phi=E->trace.basicq[j][1][it];
            rad=E->trace.basicq[j][2][it];
            x=E->trace.basicq[j][3][it];
            y=E->trace.basicq[j][4][it];
            z=E->trace.basicq[j][5][it];

            iprevious_element=E->trace.ielement[j][it];

            iel=(E->trace.iget_element)(E,j,iprevious_element,x,y,z,theta,phi,rad);
            /* debug *
            fprintf(E->trace.fpt,"BB. kk %d %d %d %d %f %f %f %f %f %f\n",kk,j,iprevious_element,iel,x,y,z,theta,phi,rad);
            fflush(E->trace.fpt);
            /**/

            E->trace.ielement[j][it]=iel;

            if (iel<0) {
                put_away_later(E,j,it);
                eject_tracer(E,j,it);
                it--;
            }

        } /* end tracers */

    } /* end j */


    /* Now take care of tracers that exited cap */

    /* REMOVE */
    /*
      parallel_process_termination();
    */

    lost_souls(E);

    /* Free later arrays */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        if (E->trace.ilater[j]>0) {
            for (kk=0;kk<=((E->trace.number_of_tracer_quantities)-1);kk++) {
                free(E->trace.rlater[j][kk]);
            }
        }
    } /* end j */


    /* Adjust Array Sizes */

    reduce_tracer_arrays(E);

    time_stat2=CPU_time0();

    fprintf(E->trace.fpt,"AA: time for find tracers: %f\n", time_stat2-time_stat1);

    return;
}

/************** LOST SOULS ****************************************************/
/*                                                                            */
/* This function is used to transport tracers to proper processor domains.    */
/* (MPI parallel)                                                             */
/*  All of the tracers that were sent to rlater arrays are destined to another */
/*  cap and sent there. Then they are raised up or down for multiple z procs.  */
/*  isend[j][n]=number of tracers this processor cap is sending to cap n        */
/*  ireceive[j][n]=number of tracers this processor cap is receiving from cap n */


void lost_souls(struct All_variables *E)
{
    int ithiscap;
    int ithatcap=1;
    int isend[13][13];
    int ireceive[13][13];
    int isize[13];
    int kk,pp,j;
    int mm;
    int numtracers;
    int icheck=0;
    int isend_position;
    int ipos,ipos2,ipos3;
    int idb;
    int idestination_proc=0;
    int isource_proc;
    int isend_z[13][3];
    int ireceive_z[13][3];
    int isum[13];
    int irad;
    int ival;
    int ithat_processor;
    int ireceive_position;
    int ihorizontal_neighbor;
    int ivertical_neighbor;
    int ilast_receiver_position;
    int it;
    int irec[13];
    int irec_position;
    int iel;
    int num_tracers;
    int isize_send;
    int isize_receive;
    int itemp_size;
    int itracers_subject_to_vertical_transport[13];

    double x,y,z;
    double theta,phi,rad;
    double *send[13][13];
    double *receive[13][13];
    double *send_z[13][3];
    double *receive_z[13][3];
    double *REC[13];

    void expand_tracer_arrays();

    int number_of_caps=12;
    int lev=E->mesh.levmax;
    int num_ngb;

    /* Note, if for some reason, the number of neighbors exceeds */
    /* 50, which is unlikely, the MPI arrays must be increased.  */
    MPI_Status status[200];
    MPI_Request request[200];
    MPI_Status status1;
    MPI_Status status2;
    int itag=1;


    parallel_process_sync(E);
    fprintf(E->trace.fpt, "Entering lost_souls()\n");


    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        E->trace.istat_isend=E->trace.ilater[j];
    }


    /* debug *
    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        for (kk=1; kk<=E->trace.istat_isend; kk++) {
            fprintf(E->trace.fpt, "tracer#=%d xx=(%g,%g,%g)\n", kk,
                    E->trace.rlater[j][0][kk],
                    E->trace.rlater[j][1][kk],
                    E->trace.rlater[j][2][kk]);
        }
        fflush(E->trace.fpt);
    }
    /**/


    /* initialize isend and ireceive */
    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        /* # of neighbors in the horizontal plane */
        num_ngb = E->parallel.TNUM_PASS[lev][j];
        isize[j]=E->trace.ilater[j]*E->trace.number_of_tracer_quantities;
        for (kk=0;kk<=num_ngb;kk++) isend[j][kk]=0;
        for (kk=0;kk<=num_ngb;kk++) ireceive[j][kk]=0;
    }

    /* Allocate Maximum Memory to Send Arrays */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        itemp_size=max(isize[j],1);

        num_ngb = E->parallel.TNUM_PASS[lev][j];
        for (kk=0;kk<=num_ngb;kk++) {
            if ((send[j][kk]=(double *)malloc(itemp_size*sizeof(double)))==NULL) {
                fprintf(E->trace.fpt,"Error(lost souls)-no memory (u389)\n");
                fflush(E->trace.fpt);
                exit(10);
            }
        }
    }


    /* debug *
    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        ithiscap=E->sphere.capid[j];
        for (kk=1;kk<=E->parallel.TNUM_PASS[lev][j];kk++) {
            ithatcap=E->parallel.PROCESSOR[lev][j].pass[kk]+1;
            fprintf(E->trace.fpt,"cap: %d proc %d TNUM: %d ithatcap: %d\n",
                    ithiscap,E->parallel.me,kk,ithatcap);

        }
        fflush(E->trace.fpt);
    }
    /**/


    /* Pre communication */
    if (E->parallel.nprocxy == 12)
        full_put_lost_tracers(E, isend, send);
    else
        regional_put_lost_tracers(E, isend, send);

    /* Send info to other processors regarding number of send tracers */

    /* idb is the request array index variable */
    /* Each send and receive has a request variable */

    idb=0;
    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        ithiscap=0;

        /* if tracer is in same cap (nprocz>1) */

        if (E->parallel.nprocz>1) {
            ireceive[j][ithiscap]=isend[j][ithiscap];
        }

        for (kk=1;kk<=E->parallel.TNUM_PASS[lev][j];kk++) {
            ithatcap=kk;

            /* if neighbor cap is in another processor, send information via MPI */

            idestination_proc=E->parallel.PROCESSOR[lev][j].pass[kk];

            idb++;
            MPI_Isend(&isend[j][ithatcap],1,MPI_INT,idestination_proc,
                      11,E->parallel.world,
                      &request[idb-1]);

        } /* end kk, number of neighbors */

    } /* end j, caps per proc */

    /* Receive tracer count info */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        for (kk=1;kk<=E->parallel.TNUM_PASS[lev][j];kk++) {
            ithatcap=kk;

            /* if neighbor cap is in another processor, receive information via MPI */

            isource_proc=E->parallel.PROCESSOR[lev][j].pass[kk];

            if (idestination_proc!=E->parallel.me) {

                idb++;
                MPI_Irecv(&ireceive[j][ithatcap],1,MPI_INT,isource_proc,
                          11,E->parallel.world,
                          &request[idb-1]);

            } /* end if */

        } /* end kk, number of neighbors */
    } /* end j, caps per proc */

    /* Wait for non-blocking calls to complete */

    MPI_Waitall(idb,request,status);

    /* Testing, should remove */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        for (kk=1;kk<=E->parallel.TNUM_PASS[lev][j];kk++) {
            isource_proc=E->parallel.PROCESSOR[lev][j].pass[kk];
            fprintf(E->trace.fpt,"j: %d send %d to cap %d\n",j,isend[j][kk],isource_proc);
            fprintf(E->trace.fpt,"j: %d rec  %d from cap %d\n",j,ireceive[j][kk],isource_proc);
        }
    }


    /* Allocate memory in receive arrays */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        num_ngb = E->parallel.TNUM_PASS[lev][j];
        for (ithatcap=1;ithatcap<=num_ngb;ithatcap++) {
            isize[j]=ireceive[j][ithatcap]*E->trace.number_of_tracer_quantities;

            itemp_size=max(1,isize[j]);

            if ((receive[j][ithatcap]=(double *)malloc(itemp_size*sizeof(double)))==NULL) {
                fprintf(E->trace.fpt,"Error(lost souls)-no memory (c721)\n");
                fflush(E->trace.fpt);
                exit(10);
            }
        }
    } /* end j */

    /* Now, send the tracers to proper caps */

    idb=0;
    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        ithiscap=0;

        /* same cap */

        if (E->parallel.nprocz>1) {

            ithatcap=ithiscap;
            isize[j]=isend[j][ithatcap]*E->trace.number_of_tracer_quantities;
            for (mm=0;mm<=(isize[j]-1);mm++) {
                receive[j][ithatcap][mm]=send[j][ithatcap][mm];
            }

        }

        /* neighbor caps */

        for (kk=1;kk<=E->parallel.TNUM_PASS[lev][j];kk++) {
            ithatcap=kk;

            idestination_proc=E->parallel.PROCESSOR[lev][j].pass[kk];

            isize[j]=isend[j][ithatcap]*E->trace.number_of_tracer_quantities;

            idb++;

            MPI_Isend(send[j][ithatcap],isize[j],MPI_DOUBLE,idestination_proc,
                      11,E->parallel.world,
                      &request[idb-1]);

        } /* end kk, number of neighbors */

    } /* end j, caps per proc */


    /* Receive tracers */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        ithiscap=0;
        for (kk=1;kk<=E->parallel.TNUM_PASS[lev][j];kk++) {
            ithatcap=kk;

            isource_proc=E->parallel.PROCESSOR[lev][j].pass[kk];

            idb++;

            isize[j]=ireceive[j][ithatcap]*E->trace.number_of_tracer_quantities;

            MPI_Irecv(receive[j][ithatcap],isize[j],MPI_DOUBLE,isource_proc,
                      11,E->parallel.world,
                      &request[idb-1]);

        } /* end kk, number of neighbors */

    } /* end j, caps per proc */

    /* Wait for non-blocking calls to complete */

    MPI_Waitall(idb,request,status);


    /* Put all received tracers in array REC[j] */
    /* This makes things more convenient.       */

    /* Sum up size of receive arrays (all tracers sent to this processor) */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        isum[j]=0;

        ithiscap=0;

        for (kk=1;kk<=E->parallel.TNUM_PASS[lev][j];kk++) {
            ithatcap=kk;
            isum[j]=isum[j]+ireceive[j][ithatcap];
        }
        if (E->parallel.nprocz>1) isum[j]=isum[j]+ireceive[j][ithiscap];

        itracers_subject_to_vertical_transport[j]=isum[j];
    }

    /* Allocate Memory for REC array */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        isize[j]=isum[j]*E->trace.number_of_tracer_quantities;
        isize[j]=max(isize[j],1);
        if ((REC[j]=(double *)malloc(isize[j]*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"Error(lost souls)-no memory (g323)\n");
            fflush(E->trace.fpt);
            exit(10);
        }
        REC[j][0]=0.0;
    }

    /* Put Received tracers in REC */


    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        irec[j]=0;

        irec_position=0;

        ithiscap=0;

        /* horizontal neighbors */

        for (ihorizontal_neighbor=1;ihorizontal_neighbor<=E->parallel.TNUM_PASS[lev][j];ihorizontal_neighbor++) {

            ithatcap=ihorizontal_neighbor;

            for (pp=1;pp<=ireceive[j][ithatcap];pp++) {
                irec[j]++;
                ipos=(pp-1)*E->trace.number_of_tracer_quantities;

                for (mm=0;mm<=(E->trace.number_of_tracer_quantities-1);mm++) {
                    ipos2=ipos+mm;
                    REC[j][irec_position]=receive[j][ithatcap][ipos2];

                    irec_position++;

                } /* end mm (cycling tracer quantities) */
            } /* end pp (cycling tracers) */
        } /* end ihorizontal_neighbors (cycling neighbors) */

        /* for tracers in the same cap (nprocz>1) */

        if (E->parallel.nprocz>1) {
            ithatcap=ithiscap;
            for (pp=1;pp<=ireceive[j][ithatcap];pp++) {
                irec[j]++;
                ipos=(pp-1)*E->trace.number_of_tracer_quantities;

                for (mm=0;mm<=(E->trace.number_of_tracer_quantities-1);mm++) {
                    ipos2=ipos+mm;

                    REC[j][irec_position]=receive[j][ithatcap][ipos2];

                    irec_position++;

                } /* end mm (cycling tracer quantities) */

            } /* end pp (cycling tracers) */

        } /* endif nproc>1 */

    } /* end j (cycling caps) */

    /* Done filling REC */



    /* VERTICAL COMMUNICATION */

    /* Note: For generality, I assume that both multiple */
    /* caps per processor as well as multiple processors */
    /* in the radial direction. These are probably       */
    /* inconsistent parameter choices for the regular    */
    /* CitcomS code.                                     */

    if (E->parallel.nprocz>1) {

        /* Allocate memory for send_z */
        /* Make send_z the size of receive array (max size) */
        /* (No dynamic reallocation of send_z necessary)    */

        for (j=1;j<=E->sphere.caps_per_proc;j++) {
            for (kk=1;kk<=E->parallel.TNUM_PASSz[lev];kk++) {
                isize[j]=itracers_subject_to_vertical_transport[j]*E->trace.number_of_tracer_quantities;
                isize[j]=max(isize[j],1);

                if ((send_z[j][kk]=(double *)malloc(isize[j]*sizeof(double)))==NULL) {
                    fprintf(E->trace.fpt,"Error(lost souls)-no memory (c721)\n");
                    fflush(E->trace.fpt);
                    exit(10);
                                }
            }
        } /* end j */


        for (j=1;j<=E->sphere.caps_per_proc;j++) {

            for (ivertical_neighbor=1;ivertical_neighbor<=E->parallel.TNUM_PASSz[lev];ivertical_neighbor++) {

                ithat_processor=E->parallel.PROCESSORz[lev].pass[ivertical_neighbor];

                /* initialize isend_z and ireceive_z array */

                isend_z[j][ivertical_neighbor]=0;
                ireceive_z[j][ivertical_neighbor]=0;

                /* sort through receive array and check radius */

                it=0;
                num_tracers=irec[j];
                for (kk=1;kk<=num_tracers;kk++) {

                    it++;

                    ireceive_position=((it-1)*E->trace.number_of_tracer_quantities);

                    irad=ireceive_position+2;

                    rad=REC[j][irad];

                    ival=icheck_that_processor_shell(E,j,ithat_processor,rad);


                    /* if tracer is in other shell, take out of receive array and give to send_z*/

                    if (ival==1) {

                        isend_z[j][ivertical_neighbor]++;

                        isend_position=(isend_z[j][ivertical_neighbor]-1)*E->trace.number_of_tracer_quantities;

                        ilast_receiver_position=(irec[j]-1)*E->trace.number_of_tracer_quantities;

                        for (mm=0;mm<=(E->trace.number_of_tracer_quantities-1);mm++) {
                            ipos=ireceive_position+mm;
                            ipos2=isend_position+mm;

                            send_z[j][ivertical_neighbor][ipos2]=REC[j][ipos];


                            /* eject tracer info from REC array, and replace with last tracer in array */

                            ipos3=ilast_receiver_position+mm;
                            REC[j][ipos]=REC[j][ipos3];

                        }

                        it--;
                        irec[j]--;

                    } /* end if ival===1 */

                    /* Otherwise, leave tracer */

                } /* end kk (cycling through tracers) */

            } /* end ivertical_neighbor */

        } /* end j */


        /* Send arrays are now filled.                         */
        /* Now send send information to vertical processor neighbor */

        for (j=1;j<=E->sphere.caps_per_proc;j++) {
            for (ivertical_neighbor=1;ivertical_neighbor<=E->parallel.TNUM_PASSz[lev];ivertical_neighbor++) {

                MPI_Sendrecv(&isend_z[j][ivertical_neighbor],1,MPI_INT,
                             E->parallel.PROCESSORz[lev].pass[ivertical_neighbor],itag,
                             &ireceive_z[j][ivertical_neighbor],1,MPI_INT,
                             E->parallel.PROCESSORz[lev].pass[ivertical_neighbor],
                             itag,E->parallel.world,&status1);

                /* for testing - remove */
                /*
                  fprintf(E->trace.fpt,"PROC: %d IVN: %d (P: %d) SEND: %d REC: %d\n",
                  E->parallel.me,ivertical_neighbor,E->parallel.PROCESSORz[lev].pass[ivertical_neighbor],
                  isend_z[j][ivertical_neighbor],ireceive_z[j][ivertical_neighbor]);
                  fflush(E->trace.fpt);
                */

            } /* end ivertical_neighbor */

        } /* end j */


        /* Allocate memory to receive_z arrays */


        for (j=1;j<=E->sphere.caps_per_proc;j++) {
            for (kk=1;kk<=E->parallel.TNUM_PASSz[lev];kk++) {
                isize[j]=ireceive_z[j][kk]*E->trace.number_of_tracer_quantities;
                isize[j]=max(isize[j],1);

                if ((receive_z[j][kk]=(double *)malloc(isize[j]*sizeof(double)))==NULL) {
                    fprintf(E->trace.fpt,"Error(lost souls)-no memory (t590)\n");
                    fflush(E->trace.fpt);
                    exit(10);
                }
            }
        } /* end j */

        /* Send Tracers */

        for (j=1;j<=E->sphere.caps_per_proc;j++) {
            for (ivertical_neighbor=1;ivertical_neighbor<=E->parallel.TNUM_PASSz[lev];ivertical_neighbor++) {
                isize_send=isend_z[j][ivertical_neighbor]*E->trace.number_of_tracer_quantities;
                isize_receive=ireceive_z[j][ivertical_neighbor]*E->trace.number_of_tracer_quantities;

                MPI_Sendrecv(send_z[j][ivertical_neighbor],isize_send,
                             MPI_DOUBLE,
                             E->parallel.PROCESSORz[lev].pass[ivertical_neighbor],itag+1,
                             receive_z[j][ivertical_neighbor],isize_receive,
                             MPI_DOUBLE,
                             E->parallel.PROCESSORz[lev].pass[ivertical_neighbor],
                             itag+1,E->parallel.world,&status2);

            }
        }

        /* Put tracers into REC array */

        /* First, reallocate memory to REC */

        for (j=1;j<=E->sphere.caps_per_proc;j++) {
            isum[j]=0;
            for (ivertical_neighbor=1;ivertical_neighbor<=E->parallel.TNUM_PASSz[lev];ivertical_neighbor++) {
                isum[j]=isum[j]+ireceive_z[j][ivertical_neighbor];
            }

            isum[j]=isum[j]+irec[j];

            isize[j]=isum[j]*E->trace.number_of_tracer_quantities;

            if (isize[j]>0) {
                if ((REC[j]=(double *)realloc(REC[j],isize[j]*sizeof(double)))==NULL) {
                    fprintf(E->trace.fpt,"Error(lost souls)-no memory (i981)\n");
                    fprintf(E->trace.fpt,"isize: %d\n",isize[j]);
                    fflush(E->trace.fpt);
                    exit(10);
                }
            }
        }


        for (j=1;j<=E->sphere.caps_per_proc;j++) {
            for (ivertical_neighbor=1;ivertical_neighbor<=E->parallel.TNUM_PASSz[lev];ivertical_neighbor++) {

                for (kk=1;kk<=ireceive_z[j][ivertical_neighbor];kk++) {
                    irec[j]++;

                    irec_position=(irec[j]-1)*E->trace.number_of_tracer_quantities;
                    ireceive_position=(kk-1)*E->trace.number_of_tracer_quantities;

                    for (mm=0;mm<=(E->trace.number_of_tracer_quantities-1);mm++) {
                        REC[j][irec_position+mm]=receive_z[j][ivertical_neighbor][ireceive_position+mm];
                    }
                }

            }
        }

        /* Free Vertical Arrays */

        for (j=1;j<=E->sphere.caps_per_proc;j++) {
            for (ivertical_neighbor=1;ivertical_neighbor<=E->parallel.TNUM_PASSz[lev];ivertical_neighbor++) {
                free(send_z[j][ivertical_neighbor]);
                free(receive_z[j][ivertical_neighbor]);
            }
        }

    } /* endif nprocz>1 */

    /* END OF VERTICAL TRANSPORT */

    /* Put away tracers */


    for (j=1;j<=E->sphere.caps_per_proc;j++) {
        for (kk=1;kk<=irec[j];kk++) {
            E->trace.ntracers[j]++;

            if (E->trace.ntracers[j]>(E->trace.max_ntracers[j]-5)) expand_tracer_arrays(E,j);

            ireceive_position=(kk-1)*E->trace.number_of_tracer_quantities;

            for (mm=0;mm<=(E->trace.number_of_basic_quantities-1);mm++) {
                ipos=ireceive_position+mm;

                E->trace.basicq[j][mm][E->trace.ntracers[j]]=REC[j][ipos];
            }
            for (mm=0;mm<=(E->trace.number_of_extra_quantities-1);mm++) {
                ipos=ireceive_position+E->trace.number_of_basic_quantities+mm;

                E->trace.extraq[j][mm][E->trace.ntracers[j]]=REC[j][ipos];
            }

            theta=E->trace.basicq[j][0][E->trace.ntracers[j]];
            phi=E->trace.basicq[j][1][E->trace.ntracers[j]];
            rad=E->trace.basicq[j][2][E->trace.ntracers[j]];
            x=E->trace.basicq[j][3][E->trace.ntracers[j]];
            y=E->trace.basicq[j][4][E->trace.ntracers[j]];
            z=E->trace.basicq[j][5][E->trace.ntracers[j]];


            iel=(E->trace.iget_element)(E,j,-99,x,y,z,theta,phi,rad);

            if (iel<1) {
                fprintf(E->trace.fpt,"Error(lost souls) - element not here?\n");
                fprintf(E->trace.fpt,"x,y,z-theta,phi,rad: %f %f %f - %f %f %f\n",x,y,z,theta,phi,rad);
                fflush(E->trace.fpt);
                exit(10);
            }

            E->trace.ielement[j][E->trace.ntracers[j]]=iel;

        }
    }

    fprintf(E->trace.fpt,"Freeing memory in lost_souls()\n");
    fflush(E->trace.fpt);
    parallel_process_sync(E);

    /* Free Arrays */

    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        ithiscap=0;

        free(REC[j]);

        free(send[j][ithiscap]);

        for (kk=1;kk<=E->parallel.TNUM_PASS[lev][j];kk++) {
            ithatcap=kk;

            free(send[j][ithatcap]);
            free(receive[j][ithatcap]);

        }

    }
    fprintf(E->trace.fpt,"Leaving lost_souls()\n");
    fflush(E->trace.fpt);

    return;
}


/***********************************************************************/
/* This function computes the number of tracers in each element.       */
/* Each tracer can be of different "flavors", which is the 0th index   */
/* of extraq. How to interprete "flavor" is left for the application.  */

void count_tracers_of_flavors(struct All_variables *E)
{

    int j, flavor, e, kk;
    int numtracers;

    for (j=1; j<=E->sphere.caps_per_proc; j++) {

        /* first zero arrays */
        for (flavor=0; flavor<E->trace.nflavors; flavor++)
            for (e=1; e<=E->lmesh.nel; e++)
                E->trace.ntracer_flavor[j][flavor][e] = 0;

        numtracers=E->trace.ntracers[j];

        /* Fill arrays */
        for (kk=1; kk<=numtracers; kk++) {
            e = E->trace.ielement[j][kk];
            flavor = E->trace.extraq[j][0][kk];
            E->trace.ntracer_flavor[j][flavor][e]++;
        }
    }

    /* debug */
    /**
    for (j=1; j<=E->sphere.caps_per_proc; j++) {
        for (e=1; e<=E->lmesh.nel; e++) {
            fprintf(E->trace.fpt, "element=%d ntracer_flaver =", e);
            for (flavor=0; flavor<E->trace.nflavors; flavor++) {
                fprintf(E->trace.fpt, " %d",
                        E->trace.ntracer_flavor[j][flavor][e]);
            }
            fprintf(E->trace.fpt, "\n");
        }
    }
    fflush(E->trace.fpt);
    /**/

    return;
}



void initialize_tracers(struct All_variables *E)
{
    void make_tracer_array();
    void read_tracer_file();
    void restart_tracers();
    int isum_tracers();

    if (E->trace.ic_method==0)
        make_tracer_array(E);
    else if (E->trace.ic_method==1)
        read_tracer_file(E);
    else if (E->trace.ic_method==2)
        restart_tracers(E);
    else {
        fprintf(E->trace.fpt,"Not ready for other inputs yet\n");
        fflush(E->trace.fpt);
        parallel_process_termination();
    }

    /* total number of tracers  */

    E->trace.ilast_tracer_count = isum_tracers(E);
    fprintf(E->trace.fpt, "Sum of Tracers: %d\n", E->trace.ilast_tracer_count);

    return;
}


/************** MAKE TRACER ARRAY ********************************/
/* Here, each processor will generate tracers somewhere          */
/* in the sphere - check if its in this cap  - then check radial */

void make_tracer_array(struct All_variables *E)
{

    int tracers_cap;
    int j;
    double processor_fraction;

    void generate_random_tracers();
    void init_tracer_flavors();

    if (E->parallel.me==0) fprintf(stderr,"Making Tracer Array\n");
    fflush(stderr);


    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        processor_fraction=E->lmesh.volume/E->mesh.volume;
        tracers_cap=E->mesh.nel*E->trace.itperel*processor_fraction;
        /*
          fprintf(stderr,"AA: proc frac: %f (%d) %d %d %f %f\n",processor_fraction,tracers_cap,E->lmesh.nel,E->parallel.nprocz, E->sx[j][3][E->lmesh.noz],E->sx[j][3][1]);
        */

        fprintf(E->trace.fpt,"\nGenerating %d Tracers\n",tracers_cap);

        generate_random_tracers(E, tracers_cap, j);



    }/* end j */


    /* Initialize tracer flavors */
    if (E->trace.nflavors) init_tracer_flavors(E);


    fprintf(stderr,"DONE Making Tracer Array (%d)\n",E->parallel.me);
    fflush(stderr);

    return;
}



void generate_random_tracers(struct All_variables *E,
                             int tracers_cap, int j)
{
    void cart_to_sphere();
    void keep_in_sphere();
    void initialize_tracer_arrays();

    int kk;
    int ival;
    int number_of_tries=0;
    int max_tries;

    double x,y,z;
    double theta,phi,rad;
    double dmin,dmax;
    double random1,random2,random3;


    initialize_tracer_arrays(E,j,tracers_cap);


    /* Tracers are placed randomly in cap */
    /* (intentionally using rand() instead of srand() )*/

    dmin=-1.0*E->sphere.ro;
    dmax=E->sphere.ro;

    while (E->trace.ntracers[j]<tracers_cap) {

        number_of_tries++;
        max_tries=500*tracers_cap;

        if (number_of_tries>max_tries) {
            fprintf(E->trace.fpt,"Error(make_tracer_array)-too many tries?\n");
            fprintf(E->trace.fpt,"%d %d %d\n",max_tries,number_of_tries,RAND_MAX);
            fflush(E->trace.fpt);
            exit(10);
        }


        random1=(1.0*rand())/(1.0*RAND_MAX);
        random2=(1.0*rand())/(1.0*RAND_MAX);
        random3=(1.0*rand())/(1.0*RAND_MAX);

        x=dmin+random1*(dmax-dmin);
        y=dmin+random2*(dmax-dmin);
        z=dmin+random3*(dmax-dmin);

        /* first check if within shell */

        cart_to_sphere(E,x,y,z,&theta,&phi,&rad);

        if (rad>=E->sx[j][3][E->lmesh.noz]) continue;
        if (rad<E->sx[j][3][1]) continue;


        /* check if in current cap */
        if (E->parallel.nprocxy==1)
            ival=regional_icheck_cap(E,0,theta,phi,rad,rad);
        else
            ival=full_icheck_cap(E,0,x,y,z,rad);

        if (ival!=1) continue;

        /* Made it, so record tracer information */

        keep_in_sphere(E,&x,&y,&z,&theta,&phi,&rad);

        E->trace.ntracers[j]++;
        kk=E->trace.ntracers[j];

        E->trace.basicq[j][0][kk]=theta;
        E->trace.basicq[j][1][kk]=phi;
        E->trace.basicq[j][2][kk]=rad;
        E->trace.basicq[j][3][kk]=x;
        E->trace.basicq[j][4][kk]=y;
        E->trace.basicq[j][5][kk]=z;

    } /* end while */

    return;
}


/******** READ TRACER ARRAY *********************************************/
/*                                                                      */
/* This function reads tracers from input file.                         */
/* All processors read the same input file, then sort out which ones    */
/* belong.                                                              */

void read_tracer_file(E)
     struct All_variables *E;
{

    char input_s[1000];

    int number_of_tracers, ncolumns;
    int kk;
    int icheck;
    int iestimate;
    int icushion;
    int i, j;

    int icheck_processor_shell();
    int isum_tracers();
    void initialize_tracer_arrays();
    void keep_in_sphere();
    void sphere_to_cart();
    void cart_to_sphere();
    void expand_tracer_arrays();

    double x,y,z;
    double theta,phi,rad;
    double extra[100];

    FILE *fptracer;

    fptracer=fopen(E->trace.tracer_file,"r");
    fprintf(E->trace.fpt,"Opening %s\n",E->trace.tracer_file);

    fgets(input_s,200,fptracer);
    sscanf(input_s,"%d %d",&number_of_tracers,&ncolumns);
    fprintf(E->trace.fpt,"%d Tracers, %d columns in file \n",
            number_of_tracers, ncolumns);

    /* some error control */
    if (E->trace.number_of_extra_quantities+3 != ncolumns) {
        fprintf(E->trace.fpt,"ERROR(read tracer file)-wrong # of columns\n");
        fflush(E->trace.fpt);
        exit(10);
    }


    /* initially size tracer arrays to number of tracers divided by processors */

    icushion=100;

    iestimate=number_of_tracers/E->parallel.nproc + icushion;

    for (j=1;j<=E->sphere.caps_per_proc;j++) {

        initialize_tracer_arrays(E,j,iestimate);

        for (kk=1;kk<=number_of_tracers;kk++) {
            fgets(input_s,200,fptracer);
            if (E->trace.number_of_extra_quantities==0) {
                sscanf(input_s,"%lf %lf %lf\n",&theta,&phi,&rad);
            }
            else if (E->trace.number_of_extra_quantities==1) {
                sscanf(input_s,"%lf %lf %lf %lf\n",&theta,&phi,&rad,&extra[0]);
            }
            /* XXX: if E->trace.number_of_extra_quantities is greater than 1 */
            /* this part has to be changed... */
            else {
                fprintf(E->trace.fpt,"ERROR(restart tracers)-huh?\n");
                fflush(E->trace.fpt);
                exit(10);
            }

            sphere_to_cart(E,theta,phi,rad,&x,&y,&z);


            /* make sure theta, phi is in range, and radius is within bounds */

            keep_in_sphere(E,&x,&y,&z,&theta,&phi,&rad);

            /* check whether tracer is within processor domain */

            icheck=1;
            if (E->parallel.nprocz>1) icheck=icheck_processor_shell(E,j,rad);
            if (icheck!=1) continue;

            if (E->parallel.nprocxy==1)
                icheck=regional_icheck_cap(E,0,theta,phi,rad,rad);
            else
                icheck=full_icheck_cap(E,0,x,y,z,rad);

            if (icheck==0) continue;

            /* if still here, tracer is in processor domain */


            E->trace.ntracers[j]++;

            if (E->trace.ntracers[j]>=(E->trace.max_ntracers[j]-5)) expand_tracer_arrays(E,j);

            E->trace.basicq[j][0][E->trace.ntracers[j]]=theta;
            E->trace.basicq[j][1][E->trace.ntracers[j]]=phi;
            E->trace.basicq[j][2][E->trace.ntracers[j]]=rad;
            E->trace.basicq[j][3][E->trace.ntracers[j]]=x;
            E->trace.basicq[j][4][E->trace.ntracers[j]]=y;
            E->trace.basicq[j][5][E->trace.ntracers[j]]=z;

            for (i=0; i<E->trace.number_of_extra_quantities; i++)
                E->trace.extraq[j][i][E->trace.ntracers[j]]=extra[i];

        } /* end kk, number of tracers */

        fprintf(E->trace.fpt,"Number of tracers in this cap is: %d\n",
                E->trace.ntracers[j]);

    } /* end j */

    fclose(fptracer);

    icheck=isum_tracers(E);

    if (icheck!=number_of_tracers) {
        fprintf(E->trace.fpt,"ERROR(read_tracer_file) - tracers != number in file\n");
        fprintf(E->trace.fpt,"Tracers in system: %d\n", icheck);
        fprintf(E->trace.fpt,"Tracers in file: %d\n", number_of_tracers);
        fflush(E->trace.fpt);
        exit(10);
    }

    return;
}


/************** RESTART TRACERS ******************************************/
/*                                                                       */
/* This function restarts tracers written from previous calculation      */
/* and the tracers are read as seperate files for each processor domain. */

void restart_tracers(E)
     struct All_variables *E;
{

    char output_file[200];
    char input_s[1000];

    int i,j,kk;
    int idum1,ncolumns;
    int numtracers;

    double rdum1;
    double theta,phi,rad;
    double extra[100];
    double x,y,z;

    void initialize_tracer_arrays();
    void sphere_to_cart();
    void keep_in_sphere();

    FILE *fp1;

    if (E->trace.number_of_extra_quantities>99) {
        fprintf(E->trace.fpt,"ERROR(restart_tracers)-increase size of extra[]\n");
        fflush(E->trace.fpt);
        parallel_process_termination();
    }

    sprintf(output_file,"%s.tracer.%d.%d",E->control.old_P_file,E->parallel.me,E->monitor.solution_cycles_init);

    if ( (fp1=fopen(output_file,"r"))==NULL) {
        fprintf(E->trace.fpt,"ERROR(restart tracers)-file not found %s\n",output_file);
        fflush(E->trace.fpt);
        exit(10);
    }

    fprintf(stderr,"Restarting Tracers from %s\n",output_file);
    fflush(stderr);


    for(j=1;j<=E->sphere.caps_per_proc;j++) {
        fgets(input_s,200,fp1);
        sscanf(input_s,"%d %d %d %lf",
               &idum1, &numtracers, &ncolumns, &rdum1);

        /* some error control */
        if (E->trace.number_of_extra_quantities+3 != ncolumns) {
            fprintf(E->trace.fpt,"ERROR(restart tracers)-wrong # of columns\n");
            fflush(E->trace.fpt);
            exit(10);
        }

        /* allocate memory for tracer arrays */

        initialize_tracer_arrays(E,j,numtracers);
        E->trace.ntracers[j]=numtracers;

        for (kk=1;kk<=numtracers;kk++) {
            fgets(input_s,200,fp1);
            if (E->trace.number_of_extra_quantities==0) {
                sscanf(input_s,"%lf %lf %lf\n",&theta,&phi,&rad);
            }
            else if (E->trace.number_of_extra_quantities==1) {
                sscanf(input_s,"%lf %lf %lf %lf\n",&theta,&phi,&rad,&extra[0]);
            }
            /* XXX: if E->trace.number_of_extra_quantities is greater than 1 */
            /* this part has to be changed... */
            else {
                fprintf(E->trace.fpt,"ERROR(restart tracers)-huh?\n");
                fflush(E->trace.fpt);
                exit(10);
            }

            sphere_to_cart(E,theta,phi,rad,&x,&y,&z);

            /* it is possible that if on phi=0 boundary, significant digits can push phi over 2pi */

            keep_in_sphere(E,&x,&y,&z,&theta,&phi,&rad);

            E->trace.basicq[j][0][kk]=theta;
            E->trace.basicq[j][1][kk]=phi;
            E->trace.basicq[j][2][kk]=rad;
            E->trace.basicq[j][3][kk]=x;
            E->trace.basicq[j][4][kk]=y;
            E->trace.basicq[j][5][kk]=z;

            for (i=0; i<E->trace.number_of_extra_quantities; i++)
                E->trace.extraq[j][i][kk]=extra[i];

        }

        fprintf(E->trace.fpt,"Read %d tracers from file %s\n",numtracers,output_file);
        fflush(E->trace.fpt);

    }
    fclose(fp1);


    return;
}





/*********** CHECK SUM **************************************************/
/*                                                                      */
/* This functions checks to make sure number of tracers is preserved    */

static void check_sum(struct All_variables *E)
{

    int number, iold_number;

    number = isum_tracers(E);

    iold_number = E->trace.ilast_tracer_count;

    if (number != iold_number) {
        fprintf(E->trace.fpt,"ERROR(check_sum)-break in conservation %d %d\n",
                number,iold_number);
        fflush(E->trace.fpt);
        parallel_process_termination();
    }

    E->trace.ilast_tracer_count = number;

    return;
}

/************* ISUM TRACERS **********************************************/
/*                                                                       */
/* This function uses MPI to sum all tracers and returns number of them. */

int isum_tracers(struct All_variables *E)
{
    int imycount;
    int iallcount;
    int j;

    iallcount = 0;

    imycount = 0;
    for (j=1; j<=E->sphere.caps_per_proc; j++)
        imycount = imycount + E->trace.ntracers[j];

    MPI_Allreduce(&imycount,&iallcount,1,MPI_INT,MPI_SUM,E->parallel.world);

    return iallcount;
}



/********** CART TO SPHERE ***********************/
void cart_to_sphere(struct All_variables *E,
                    double x, double y, double z,
                    double *theta, double *phi, double *rad)
{

    double temp;
    double myatan();

    temp=x*x+y*y;

    *rad=sqrt(temp+z*z);
    *theta=atan2(sqrt(temp),z);
    *phi=myatan(y,x);


    return;
}

/********** SPHERE TO CART ***********************/
void sphere_to_cart(struct All_variables *E,
                    double theta, double phi, double rad,
                    double *x, double *y, double *z)
{

    double sint,cost,sinf,cosf;
    double temp;

    sint=sin(theta);
    cost=cos(theta);
    sinf=sin(phi);
    cosf=cos(phi);

    temp=rad*sint;

    *x=temp*cosf;
    *y=temp*sinf;
    *z=rad*cost;

    return;
}



void init_tracer_flavors(struct All_variables *E)
{
    int j, kk, number_of_tracers;
    double rad;


    /* ic_method_for_flavors == 0 (layered structure) */
    /* any tracer above z_interface is of flavor 0    */
    /* any tracer below z_interface is of flavor 1    */
    if (E->trace.ic_method_for_flavors == 0) {
        for (j=1;j<=E->sphere.caps_per_proc;j++) {

            number_of_tracers = E->trace.ntracers[j];
            for (kk=1;kk<=number_of_tracers;kk++) {
                rad = E->trace.basicq[j][2][kk];

                if (rad<=E->trace.z_interface) E->trace.extraq[j][0][kk]=1.0;
                if (rad>E->trace.z_interface) E->trace.extraq[j][0][kk]=0.0;
            }
        }
    }

    return;
}


/******************* get_neighboring_caps ************************************/
/*                                                                           */
/* Communicate with neighboring processors to get their cap boundaries,      */
/* which is later used by (E->trace.icheck_cap)()                            */
/*                                                                           */

void get_neighboring_caps(struct All_variables *E)
{
    void sphere_to_cart();

    const int ncorners = 4; /* # of top corner nodes */
    int i, j, n, d, kk, lev, idb;
    int num_ngb, neighbor_proc, tag;
    MPI_Status status[200];
    MPI_Request request[200];

    int node[ncorners];
    double xx[ncorners*3], rr[12][ncorners*3];
    int nox,noy,noz,dims;
    double x,y,z;
    double theta,phi,rad;

    dims=E->mesh.nsd;
    nox=E->lmesh.nox;
    noy=E->lmesh.noy;
    noz=E->lmesh.noz;

    node[0]=nox*noz*(noy-1)+noz;
    node[1]=noz;
    node[2]=noz*nox;
    node[3]=noz*nox*noy;

    lev = E->mesh.levmax;
    tag = 45;

    for (j=1; j<=E->sphere.caps_per_proc; j++) {

        /* loop over top corners to get their coordinates */
        n = 0;
        for (i=0; i<ncorners; i++) {
            for (d=0; d<dims; d++) {
                xx[n] = E->sx[j][d+1][node[i]];
                n++;
            }
        }

        idb = 0;
        num_ngb = E->parallel.TNUM_PASS[lev][j];
        for (kk=1; kk<=num_ngb; kk++) {
            neighbor_proc = E->parallel.PROCESSOR[lev][j].pass[kk];

            MPI_Isend(xx, n, MPI_DOUBLE, neighbor_proc,
                      tag, E->parallel.world, &request[idb]);
            idb++;

            MPI_Irecv(rr[kk], n, MPI_DOUBLE, neighbor_proc,
                      tag, E->parallel.world, &request[idb]);
            idb++;
        }

        /* Storing the current cap information */
        for (i=0; i<n; i++)
            rr[0][i] = xx[i];

        /* Wait for non-blocking calls to complete */

        MPI_Waitall(idb, request, status);

        /* Storing the received cap information
         * XXX: this part assumes:
         *      1) E->sphere.caps_per_proc==1
         *      2) E->mesh.nsd==3
         */
        for (kk=0; kk<=num_ngb; kk++) {
            for (i=1; i<=ncorners; i++) {
                theta = rr[kk][(i-1)*dims];
                phi = rr[kk][(i-1)*dims+1];
                rad = rr[kk][(i-1)*dims+2];

                sphere_to_cart(E, theta, phi, rad, &x, &y, &z);

                E->trace.xcap[kk][i] = x;
                E->trace.ycap[kk][i] = y;
                E->trace.zcap[kk][i] = z;
                E->trace.theta_cap[kk][i] = theta;
                E->trace.phi_cap[kk][i] = phi;
                E->trace.rad_cap[kk][i] = rad;
                E->trace.cos_theta[kk][i] = cos(theta);
                E->trace.sin_theta[kk][i] = sin(theta);
                E->trace.cos_phi[kk][i] = cos(phi);
                E->trace.sin_phi[kk][i] = sin(phi);
            }
        } /* end kk, number of neighbors */

        /* debugging output *
        for (kk=0; kk<=num_ngb; kk++) {
            for (i=1; i<=ncorners; i++) {
                fprintf(E->trace.fpt, "pass=%d corner=%d sx=(%g, %g, %g)\n",
                        kk, i,
                        E->trace.theta_cap[kk][i],
                        E->trace.phi_cap[kk][i],
                        E->trace.rad_cap[kk][i]);
            }
        }
        fflush(E->trace.fpt);
        /**/
    }

    return;
}


/**************** INITIALIZE TRACER ARRAYS ************************************/
/*                                                                            */
/* This function allocates memories to tracer arrays.                         */

void initialize_tracer_arrays(struct All_variables *E,
                              int j, int number_of_tracers)
{

    int kk;

    /* max_ntracers is physical size of tracer array */
    /* (initially make it 25% larger than required */

    E->trace.max_ntracers[j]=number_of_tracers+number_of_tracers/4;
    E->trace.ntracers[j]=0;

    /* make tracer arrays */

    if ((E->trace.ielement[j]=(int *) malloc(E->trace.max_ntracers[j]*sizeof(int)))==NULL) {
        fprintf(E->trace.fpt,"ERROR(make tracer array)-no memory 1a\n");
        fflush(E->trace.fpt);
        exit(10);
    }
    for (kk=1;kk<E->trace.max_ntracers[j];kk++)
        E->trace.ielement[j][kk]=-99;


    for (kk=0;kk<E->trace.number_of_basic_quantities;kk++) {
        if ((E->trace.basicq[j][kk]=(double *)malloc(E->trace.max_ntracers[j]*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"ERROR(initialize tracer arrays)-no memory 1b.%d\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
    }

    for (kk=0;kk<E->trace.number_of_extra_quantities;kk++) {
        if ((E->trace.extraq[j][kk]=(double *)malloc(E->trace.max_ntracers[j]*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"ERROR(initialize tracer arrays)-no memory 1c.%d\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
    }

    if (E->trace.nflavors > 0) {
        E->trace.ntracer_flavor[j]=(int **)malloc(E->trace.nflavors*sizeof(int*));
        for (kk=0;kk<E->trace.nflavors;kk++) {
            if ((E->trace.ntracer_flavor[j][kk]=(int *)malloc((E->lmesh.nel+1)*sizeof(int)))==NULL) {
                fprintf(E->trace.fpt,"ERROR(initialize tracer arrays)-no memory 1c.%d\n",kk);
                fflush(E->trace.fpt);
                exit(10);
            }
        }
    }


    fprintf(E->trace.fpt,"Physical size of tracer arrays (max_ntracers): %d\n",
            E->trace.max_ntracers[j]);
    fflush(E->trace.fpt);

    return;
}



/****** EXPAND TRACER ARRAYS *****************************************/

void expand_tracer_arrays(struct All_variables *E, int j)
{

    int inewsize;
    int kk;
    int icushion;

    /* expand basicq and ielement by 20% */

    icushion=100;

    inewsize=E->trace.max_ntracers[j]+E->trace.max_ntracers[j]/5+icushion;

    if ((E->trace.ielement[j]=(int *)realloc(E->trace.ielement[j],inewsize*sizeof(int)))==NULL) {
        fprintf(E->trace.fpt,"ERROR(expand tracer arrays )-no memory (ielement)\n");
        fflush(E->trace.fpt);
        exit(10);
    }

    for (kk=0;kk<=((E->trace.number_of_basic_quantities)-1);kk++) {
        if ((E->trace.basicq[j][kk]=(double *)realloc(E->trace.basicq[j][kk],inewsize*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"ERROR(expand tracer arrays )-no memory (%d)\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
    }

    for (kk=0;kk<=((E->trace.number_of_extra_quantities)-1);kk++) {
        if ((E->trace.extraq[j][kk]=(double *)realloc(E->trace.extraq[j][kk],inewsize*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"ERROR(expand tracer arrays )-no memory 78 (%d)\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
    }


    fprintf(E->trace.fpt,"Expanding physical memory of ielement, basicq, and extraq to %d from %d\n",
            inewsize,E->trace.max_ntracers[j]);

    E->trace.max_ntracers[j]=inewsize;

    return;
}




/****** REDUCE  TRACER ARRAYS *****************************************/

void reduce_tracer_arrays(struct All_variables *E)
{

    int inewsize;
    int kk;
    int iempty_space;
    int j;

    int icushion=100;

    for (j=1;j<=E->sphere.caps_per_proc;j++) {


        /* if physical size is double tracer size, reduce it */

        iempty_space=(E->trace.max_ntracers[j]-E->trace.ntracers[j]);

        if (iempty_space>(E->trace.ntracers[j]+icushion)) {


            inewsize=E->trace.ntracers[j]+E->trace.ntracers[j]/4+icushion;

            if (inewsize<1) {
                fprintf(E->trace.fpt,"Error(reduce tracer arrays)-something up (hdf3)\n");
                fflush(E->trace.fpt);
                exit(10);
            }


            if ((E->trace.ielement[j]=(int *)realloc(E->trace.ielement[j],inewsize*sizeof(int)))==NULL) {
                fprintf(E->trace.fpt,"ERROR(reduce tracer arrays )-no memory (ielement)\n");
                fflush(E->trace.fpt);
                exit(10);
            }


            for (kk=0;kk<=((E->trace.number_of_basic_quantities)-1);kk++) {
                if ((E->trace.basicq[j][kk]=(double *)realloc(E->trace.basicq[j][kk],inewsize*sizeof(double)))==NULL) {
                    fprintf(E->trace.fpt,"AKM(reduce tracer arrays )-no memory (%d)\n",kk);
                    fflush(E->trace.fpt);
                    exit(10);
                }
            }

            for (kk=0;kk<=((E->trace.number_of_extra_quantities)-1);kk++) {
                if ((E->trace.extraq[j][kk]=(double *)realloc(E->trace.extraq[j][kk],inewsize*sizeof(double)))==NULL) {
                    fprintf(E->trace.fpt,"AKM(reduce tracer arrays )-no memory 783 (%d)\n",kk);
                    fflush(E->trace.fpt);
                    exit(10);
                }
            }


            fprintf(E->trace.fpt,"Reducing physical memory of ielement, basicq, and extraq to %d from %d\n",
                    E->trace.max_ntracers[j],inewsize);

            E->trace.max_ntracers[j]=inewsize;

        } /* end if */

    } /* end j */

    return;
}


/********** PUT AWAY LATER ************************************/
/*                                             */
/* rlater has a similar structure to basicq     */
/* ilatersize is the physical memory and       */
/* ilater is the number of tracers             */

void put_away_later(struct All_variables *E, int j, int it)
{
    int kk;
    void expand_later_array();


    /* The first tracer in initiates memory allocation. */
    /* Memory is freed after parallel communications    */

    if (E->trace.ilater[j]==0) {

        E->trace.ilatersize[j]=E->trace.max_ntracers[j]/5;

        for (kk=0;kk<=((E->trace.number_of_tracer_quantities)-1);kk++) {
            if ((E->trace.rlater[j][kk]=(double *)malloc(E->trace.ilatersize[j]*sizeof(double)))==NULL) {
                fprintf(E->trace.fpt,"AKM(put_away_later)-no memory (%d)\n",kk);
                fflush(E->trace.fpt);
                exit(10);
            }
        }
    } /* end first particle initiating memory allocation */


    /* Put tracer in later array */

    E->trace.ilater[j]++;

    if (E->trace.ilater[j] >= (E->trace.ilatersize[j]-5)) expand_later_array(E,j);

    /* stack basic and extra quantities together (basic first) */

    for (kk=0;kk<=((E->trace.number_of_basic_quantities)-1);kk++)
        E->trace.rlater[j][kk][E->trace.ilater[j]]=E->trace.basicq[j][kk][it];

    for (kk=0;kk<=((E->trace.number_of_extra_quantities)-1);kk++)
        E->trace.rlater[j][E->trace.number_of_basic_quantities+kk][E->trace.ilater[j]]=E->trace.extraq[j][kk][it];


    return;
}


/****** EXPAND LATER ARRAY *****************************************/

void expand_later_array(struct All_variables *E, int j)
{

    int inewsize;
    int kk;
    int icushion;

    /* expand rlater by 20% */

    icushion=100;

    inewsize=E->trace.ilatersize[j]+E->trace.ilatersize[j]/5+icushion;

    for (kk=0;kk<=((E->trace.number_of_tracer_quantities)-1);kk++) {
        if ((E->trace.rlater[j][kk]=(double *)realloc(E->trace.rlater[j][kk],inewsize*sizeof(double)))==NULL) {
            fprintf(E->trace.fpt,"AKM(expand later array )-no memory (%d)\n",kk);
            fflush(E->trace.fpt);
            exit(10);
        }
    }


    fprintf(E->trace.fpt,"Expanding physical memory of rlater to %d from %d\n",
            inewsize,E->trace.ilatersize[j]);

    E->trace.ilatersize[j]=inewsize;

    return;
}


/***** EJECT TRACER ************************************************/

void eject_tracer(struct All_variables *E, int j, int it)
{

    int ilast_tracer;
    int kk;


    ilast_tracer=E->trace.ntracers[j];

    /* put last tracer in ejected tracer position */

    E->trace.ielement[j][it]=E->trace.ielement[j][ilast_tracer];

    for (kk=0;kk<=((E->trace.number_of_basic_quantities)-1);kk++)
        E->trace.basicq[j][kk][it]=E->trace.basicq[j][kk][ilast_tracer];

    for (kk=0;kk<=((E->trace.number_of_extra_quantities)-1);kk++)
        E->trace.extraq[j][kk][it]=E->trace.extraq[j][kk][ilast_tracer];



    E->trace.ntracers[j]--;

    return;
}



/********** ICHECK PROCESSOR SHELL *************/
/* returns -99 if rad is below current shell  */
/* returns 0 if rad is above current shell    */
/* returns 1 if rad is within current shell   */
/*                                            */
/* Shell, here, refers to processor shell     */
/*                                            */
/* shell is defined as bottom boundary up to  */
/* and not including the top boundary unless  */
/* the shell in question is the top shell     */

int icheck_processor_shell(struct All_variables *E,
                           int j, double rad)
{

    const int noz = E->lmesh.noz;
    const int nprocz = E->parallel.nprocz;
    double top_r, bottom_r;

    if (nprocz==1) return 1;

    top_r = E->sx[j][3][noz];
    bottom_r = E->sx[j][3][1];

    /* First check bottom */

    if (rad<bottom_r) return -99;


    /* Check top */

    if (rad<top_r) return 1;

    /* top processor */

    if ( (rad<=top_r) && (E->parallel.me_loc[3]==nprocz-1) ) return 1;

    /* If here, means point is above processor */
    return 0;
}


/********* ICHECK THAT PROCESSOR SHELL ********/
/*                                            */
/* Checks whether a given radius is within    */
/* a given processors radial domain.          */
/* Returns 0 if not, 1 if so.                 */
/* The domain is defined as including the bottom */
/* radius, but excluding the top radius unless   */
/* we the processor domain is the one that       */
/* is at the surface (then both boundaries are   */
/* included).                                    */

int icheck_that_processor_shell(struct All_variables *E,
                                int j, int nprocessor, double rad)
{
    int icheck_processor_shell();
    int me = E->parallel.me;

    /* nprocessor is right on top of me */
    if (nprocessor == me+1) {
        if (icheck_processor_shell(E, j, rad) == 0) return 1;
        else return 0;
    }

    /* nprocessor is right on bottom of me */
    if (nprocessor == me-1) {
        if (icheck_processor_shell(E, j, rad) == -99) return 1;
        else return 0;
    }

    /* Shouldn't be here */
    fprintf(E->trace.fpt, "Should not be here\n");
    fprintf(E->trace.fpt, "Error(check_shell) nprocessor: %d, radius: %f\n",
            nprocessor, rad);
    fflush(E->trace.fpt);
    exit(10);

    return 0;
}


