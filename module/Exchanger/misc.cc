// -*- C++ -*-
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  <LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#include <portinfo>
#include <Python.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "journal/journal.h"
#include "mpi.h"
#include "mpi/Communicator.h"
#include "BoundedBox.h"
#include "misc.h"

extern "C" {
#include "global_defs.h"
#include "element_definitions.h"
#include "citcom_init.h"
}

void commonE(All_variables*);
void initTemperatureTest(const BoundedBox&, All_variables*);
void hot_blob(const BoundedBox& bbox, All_variables* E);
void five_hot_blobs(const BoundedBox& bbox, All_variables* E);
void add_hot_blob(All_variables* E,
		  double x_center, double y_center, double z_center,
		  double radius, double baseline, double amp);
void debug_output(const All_variables* E);


// copyright

char pyExchanger_copyright__doc__[] = "";
char pyExchanger_copyright__name__[] = "copyright";

static char pyExchanger_copyright_note[] =
    "Exchanger python module: Copyright (c) 1998-2003 California Institute of Technology";


PyObject * pyExchanger_copyright(PyObject *, PyObject *)
{
    return Py_BuildValue("s", pyExchanger_copyright_note);
}

// hello

char pyExchanger_hello__doc__[] = "";
char pyExchanger_hello__name__[] = "hello";

PyObject * pyExchanger_hello(PyObject *, PyObject *)
{
    return Py_BuildValue("s", "hello");
}

//
//



// return (All_variables* E)

char pyExchanger_FinereturnE__doc__[] = "";
char pyExchanger_FinereturnE__name__[] = "FinereturnE";

PyObject * pyExchanger_FinereturnE(PyObject *, PyObject *args)
{
    PyObject *Obj;

    if (!PyArg_ParseTuple(args, "O:FinereturnE", &Obj))
        return NULL;

    mpi::Communicator * comm = (mpi::Communicator *) PyCObject_AsVoidPtr(Obj);
    MPI_Comm world = comm->handle();

    All_variables *E = citcom_init(&world);

    E->lmesh.nox = 4;
    E->lmesh.noy = 4;
    E->lmesh.noz = 3;

    E->control.theta_max = 2.1;
    E->control.theta_min = 0.9;
    E->control.fi_max = 2.1;
    E->control.fi_min = 0.9;
    E->sphere.ro = 1.9;
    E->sphere.ri = 0.9;

    commonE(E);

    PyObject *cobj = PyCObject_FromVoidPtr(E, NULL);
    return Py_BuildValue("O", cobj);
}


char pyExchanger_CoarsereturnE__doc__[] = "";
char pyExchanger_CoarsereturnE__name__[] = "CoarsereturnE";

PyObject * pyExchanger_CoarsereturnE(PyObject *, PyObject *args)
{
    PyObject *Obj;

    if (!PyArg_ParseTuple(args, "O:CoarsereturnE", &Obj))
        return NULL;

    mpi::Communicator * comm = (mpi::Communicator *) PyCObject_AsVoidPtr(Obj);
    MPI_Comm world = comm->handle();

    All_variables *E = citcom_init(&world);

    E->lmesh.nox = 4;
    E->lmesh.noy = 4;
    E->lmesh.noz = 3;

    E->control.theta_max = 3.0;
    E->control.theta_min = 0.0;
    E->control.fi_max = 3.0;
    E->control.fi_min = 0.0;
    E->sphere.ro = 2.0;
    E->sphere.ri = 0.0;

    commonE(E);

//    std::ofstream cfile("coarse.dat");
    for(int m=1;m<=E->sphere.caps_per_proc;m++)
        for(int k=1;k<=E->lmesh.noy;k++)
	    for(int j=1;j<=E->lmesh.nox;j++)
		for(int i=1;i<=E->lmesh.noz;i++)  {
		    int node = i + (j-1)*E->lmesh.noz
			     + (k-1)*E->lmesh.noz*E->lmesh.nox;

 		    E->T[m][node] = E->sx[m][1][node]
			+ E->sx[m][2][node]
			+ E->sx[m][3][node];

		    E->sphere.cap[m].V[1][node] = E->T[m][node];
		    E->sphere.cap[m].V[2][node] = 2.0*E->T[m][node];
		    E->sphere.cap[m].V[3][node] = 3.0*E->T[m][node];

//  		    cfile << "in CoarsereturnE (T, v1,v2,v3): "
// 			  <<  node << " "
// 			  << E->sx[m][1][node] << " "
// 			  << E->sx[m][2][node] << " "
// 			  << E->sx[m][3][node] << " "
// 			  << E->T[m][node] << " "
// 			  << E->sphere.cap[m].V[1][node] << " "
// 			  << E->sphere.cap[m].V[2][node] << " "
// 			  << E->sphere.cap[m].V[3][node] << " "
// 			  << std::endl;
	    }
//     cfile.close();


    PyObject *cobj = PyCObject_FromVoidPtr(E, NULL);
    return Py_BuildValue("O", cobj);
}


void commonE(All_variables *E)
{
    E->control.accuracy = 1e-6;
    E->control.tole_comp = 1e-7;

    E->parallel.nprocxy = 1;

    E->parallel.nprocx = E->parallel.nproc;
    E->parallel.nprocy = 1;
    E->parallel.nprocz = 1;

    E->parallel.me_loc[1] = E->parallel.me;
    E->parallel.me_loc[2] = 0;
    E->parallel.me_loc[3] = 0;

    E->sphere.caps_per_proc = 1;

    E->mesh.levmax = 1;
    E->mesh.levmin = 1;
    E->mesh.dof = 3;

    E->lmesh.elx = E->lmesh.nox - 1;
    E->lmesh.elz = E->lmesh.noz - 1;
    E->lmesh.ely = E->lmesh.noy - 1;

    E->lmesh.nno = E->lmesh.noz * E->lmesh.nox * E->lmesh.noy;
    E->lmesh.nel = E->lmesh.ely * E->lmesh.elx * E->lmesh.elz;
    E->lmesh.npno = E->lmesh.nel;

    int noz = E->lmesh.noz;
    int noy = E->lmesh.noy;
    int nox = E->lmesh.nox;

    E->lmesh.ELX[E->mesh.levmax] = nox-1;
    E->lmesh.ELY[E->mesh.levmax] = noy-1;
    E->lmesh.ELZ[E->mesh.levmax] = noz-1;
    E->lmesh.NOZ[E->mesh.levmax] = noz;
    E->lmesh.NOY[E->mesh.levmax] = noy;
    E->lmesh.NOX[E->mesh.levmax] = nox;
    E->lmesh.NNO[E->mesh.levmax] = nox * noz * noy;
    E->lmesh.NEL[E->mesh.levmax] = (nox-1) * (noz-1) * (noy-1);

    E->mesh.elx = E->lmesh.elx * E->parallel.nprocx;
    E->mesh.elz = E->lmesh.elz * E->parallel.nprocz;
    E->mesh.ely = E->lmesh.ely * E->parallel.nprocy;

    E->mesh.nox = E->mesh.elx + 1;
    E->mesh.noz = E->mesh.elz + 1;
    E->mesh.noy = E->mesh.ely + 1;

    E->mesh.nno = E->mesh.nox * E->mesh.noy * E->mesh.noz;
    E->mesh.nel = E->mesh.ely * E->mesh.elx * E->mesh.elz;

    for (int j=1;j<=E->sphere.caps_per_proc;j++)  {
	E->sphere.capid[j] = 1;
    }

    for (int lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {
	for (int j=1;j<=E->sphere.caps_per_proc;j++)  {
	    E->IEN[lev][j] = new IEN [E->lmesh.nel+1];
	    E->ECO[lev][j] = (struct COORD *) malloc((E->lmesh.nno+2)*sizeof(struct COORD));
	}
    }

    for (int j=1;j<=E->sphere.caps_per_proc;j++) {
	E->ien[j] = E->IEN[E->mesh.levmax][j];
	E->eco[j] = E->ECO[E->mesh.levmax][j];
    }

    for (int lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)
	for (int j=1;j<=E->sphere.caps_per_proc;j++)
	    for(int n=1; n<=E->lmesh.nno; n++) {
		E->ECO[lev][j][n].area = 1.0;
	    }

    for(int m=1;m<=E->sphere.caps_per_proc;m++) {
	E->sphere.cap[m].TB[1] = new float[E->lmesh.nno+1];
	E->sphere.cap[m].TB[2] = new float[E->lmesh.nno+1];
	E->sphere.cap[m].TB[3] = new float[E->lmesh.nno+1];
	E->sphere.cap[m].VB[1] = new float[E->lmesh.nno+1];
	E->sphere.cap[m].VB[2] = new float[E->lmesh.nno+1];
	E->sphere.cap[m].VB[3] = new float[E->lmesh.nno+1];
	E->node[m] = new unsigned int[E->lmesh.nno+1];

	for(int n=1; n<=E->lmesh.nno; n++) {
	    E->sphere.cap[m].TB[1][n] = 0;
	    E->sphere.cap[m].TB[2][n] = 0;
	    E->sphere.cap[m].TB[3][n] = 0;
	    E->sphere.cap[m].VB[1][n] = 0;
	    E->sphere.cap[m].VB[2][n] = 0;
	    E->sphere.cap[m].VB[3][n] = 0;
	    E->node[m][n] = 0;
	}
    }

    for (int lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {
	for (int j=1;j<=E->sphere.caps_per_proc;j++)  {

	    int elx = E->lmesh.ELX[lev];
	    int elz = E->lmesh.ELZ[lev];
	    int ely = E->lmesh.ELY[lev];
	    int nox = E->lmesh.NOX[lev];
	    int noz = E->lmesh.NOZ[lev];

	    for(int r=1;r<=ely;r++)
		for(int q=1;q<=elx;q++)
		    for(int p=1;p<=elz;p++)     {
			int element = (r-1)*elx*elz + (q-1)*elz  + p;
			int start = (r-1)*noz*nox + (q-1)*noz + p;
			for(int rr=1;rr<=8;rr++) {
			    E->IEN[lev][j][element].node[rr]= start
				+ offset[rr].vector[0]
				+ offset[rr].vector[1]*noz
				+ offset[rr].vector[2]*noz*nox;

// 			    std::cout << "  el = " << element
// 				      << "  lnode = " << rr
// 				      << "  ien = "
// 				      << E->IEN[lev][j][element].node[rr]
// 				      << std::endl;
			}
		    }
	}     /* end for cap j */
    }     /* end loop for lev */

    const int n = E->lmesh.nno;
    for(int m=1;m<=E->sphere.caps_per_proc;m++) {
	for(int i=1; i<=E->mesh.dof; i++) {
	    // Don't forget to delete these later
	    E->sx[m][i] = new double [n+1];
	    E->sphere.cap[m].V[i] = new float [n+1];
	}
  	E->T[m] = new double [n+1];
    }

    for(int m=1;m<=E->sphere.caps_per_proc;m++)
	for(int k=1;k<=E->lmesh.noy;k++)
	    for(int j=1;j<=E->lmesh.nox;j++)
		for(int i=1;i<=E->lmesh.noz;i++)  {
		    int node = i + (j-1)*E->lmesh.noz
			     + (k-1)*E->lmesh.noz*E->lmesh.nox;
		    E->sx[m][1][node] =
			(E->control.theta_max - E->control.theta_min)
  			  / E->mesh.elx * (j-1)
			+ (E->control.theta_max - E->control.theta_min)
			  * E->parallel.me_loc[1] / E->parallel.nprocx
			+ E->control.theta_min;
		    E->sx[m][2][node] =
			(E->control.fi_max - E->control.fi_min)
			  / E->mesh.ely * (k-1)
			+ (E->control.fi_max - E->control.fi_min)
			  * E->parallel.me_loc[2] / E->parallel.nprocy
			+ E->control.fi_min;
		    E->sx[m][3][node] =
			(E->sphere.ro -  E->sphere.ri)
			  / E->mesh.elz * (i-1)
			+ (E->sphere.ro - E->sphere.ri)
			  * E->parallel.me_loc[3] / E->parallel.nprocz
			+  E->sphere.ri;

//  		    std::cout <<  node << " "
//  			      << E->sx[m][1][node] << " "
//  			      << E->sx[m][2][node] << " "
//  			      << E->sx[m][3][node] << " "
//  			      << std::endl;
		}

    return;
}


char pyExchanger_initTemperatureTest__doc__[] = "";
char pyExchanger_initTemperatureTest__name__[] = "initTemperatureTest";

PyObject * pyExchanger_initTemperatureTest(PyObject *, PyObject *args)
{
    PyObject *obj1, *obj2;

    if (!PyArg_ParseTuple(args, "OO:initTemperatureTest", &obj1, &obj2))
	return NULL;

    BoundedBox* bbox = static_cast<BoundedBox*>(PyCObject_AsVoidPtr(obj1));
    All_variables* E = static_cast<All_variables*>(PyCObject_AsVoidPtr(obj2));

    initTemperatureTest(*bbox, E);

    Py_INCREF(Py_None);
    return Py_None;
}



void initTemperatureTest(const BoundedBox& bbox, All_variables* E)
{
    journal::debug_t debug("Exchanger");
    debug << journal::loc(__HERE__)
          << "in initTemperatureTest" << journal::end;

    //hot_blob(bbox, E);
    five_hot_blobs(bbox, E);

    debug_output(E);
}


void hot_blob(const BoundedBox& bbox, All_variables* E)
{
    // put a hot blob in the center of fine grid mesh and T=0 elsewhere

    for(int m=1;m<=E->sphere.caps_per_proc;m++)
	for(int i=1; i<E->lmesh.nno; ++i)
	    E->T[m][i] = 0;

    const double theta_min = bbox[0][0];
    const double theta_max = bbox[1][0];
    const double fi_min = bbox[0][1];
    const double fi_max = bbox[1][1];
    const double ri = bbox[0][2];
    const double ro = bbox[1][2];

    // radius of the blob is one third of the smallest dimension
    double d = std::min(std::min(theta_max - theta_min,
				 fi_max - fi_min),
                        ro - ri) / 3;

    // center of fine grid mesh
    double theta_center = 0.5 * (theta_max + theta_min);
    double fi_center = 0.5 * (fi_max + fi_min);
    double r_center = 0.5 * (ro + ri);

    double x_center = r_center * sin(fi_center) * cos(theta_center);
    double y_center = r_center * sin(fi_center) * sin(theta_center);
    double z_center = r_center * cos(fi_center);

    // compute temperature field according to nodal coordinate
    add_hot_blob(E, x_center, y_center, z_center, d, 0.5, 0.5);
}


void five_hot_blobs(const BoundedBox& bbox, All_variables* E)
{
    // put a hot blob in the center of fine grid mesh and T=0 elsewhere
    // also put 4 hot blobs around bbox in coarse mesh

    for(int m=1;m<=E->sphere.caps_per_proc;m++)
	for(int i=1; i<E->lmesh.nno; ++i)
	    E->T[m][i] = 0;

    const double theta_min = bbox[0][0];
    const double theta_max = bbox[1][0];
    const double fi_min = bbox[0][1];
    const double fi_max = bbox[1][1];
    const double ri = bbox[0][2];
    const double ro = bbox[1][2];

    // radius of blobs is one third of the smallest dimension
    double d = std::min(std::min(theta_max - theta_min,
				 fi_max - fi_min),
			ro - ri) / 3;

    // center of hot blob is in the center of bbx
    double theta_center = 0.5 * (theta_max + theta_min);
    double fi_center = 0.5 * (fi_max + fi_min);
    double r_center = 0.5 * (ro + ri);
    double x_center = r_center * sin(fi_center) * cos(theta_center);
    double y_center = r_center * sin(fi_center) * sin(theta_center);
    double z_center = r_center * cos(fi_center);
    add_hot_blob(E, x_center, y_center, z_center, d, 0.5, 0.5);

    // center of hot blob is outside bbx
    theta_center = theta_max + 0.4 * (theta_max - theta_min);
    x_center = r_center * sin(fi_center) * cos(theta_center);
    y_center = r_center * sin(fi_center) * sin(theta_center);
    z_center = r_center * cos(fi_center);
    add_hot_blob(E, x_center, y_center, z_center, d, 0.5, 0.5);

    theta_center = theta_min - 0.4 * (theta_max - theta_min);
    x_center = r_center * sin(fi_center) * cos(theta_center);
    y_center = r_center * sin(fi_center) * sin(theta_center);
    z_center = r_center * cos(fi_center);
    add_hot_blob(E, x_center, y_center, z_center, d, 0.5, 0.5);

    theta_center = 0.5 * (theta_max + theta_min);
    fi_center = fi_max + 0.4 * (fi_max - fi_min);
    x_center = r_center * sin(fi_center) * cos(theta_center);
    y_center = r_center * sin(fi_center) * sin(theta_center);
    z_center = r_center * cos(fi_center);
    add_hot_blob(E, x_center, y_center, z_center, d, 0.5, 0.5);

    fi_center = fi_min - 0.4 * (fi_max - fi_min);
    x_center = r_center * sin(fi_center) * cos(theta_center);
    y_center = r_center * sin(fi_center) * sin(theta_center);
    z_center = r_center * cos(fi_center);
    add_hot_blob(E, x_center, y_center, z_center, d, 0.5, 0.5);

}


void add_hot_blob(All_variables* E,
		  double x_center, double y_center, double z_center,
		  double radius, double baseline, double amp)
{
    // compute temperature field according to nodal coordinate
    for(int m=1;m<=E->sphere.caps_per_proc;m++)
        for(int k=1;k<=E->lmesh.noy;k++)
            for(int j=1;j<=E->lmesh.nox;j++)
                for(int i=1;i<=E->lmesh.noz;i++)  {
                    int node = i + (j-1)*E->lmesh.noz
                             + (k-1)*E->lmesh.noz*E->lmesh.nox;

                    double theta = E->sx[m][1][node];
                    double fi = E->sx[m][2][node];
                    double r = E->sx[m][3][node];

		    double x = r * sin(fi) * cos(theta);
                    double y = r * sin(fi) * sin(theta);
                    double z = r * cos(fi);

                    double distance = sqrt((x - x_center)*(x - x_center) +
                                           (y - y_center)*(y - y_center) +
                                           (z - z_center)*(z - z_center));

                    if (distance < radius)
                        E->T[m][node] += baseline
			              + amp * cos(distance/radius * M_PI);

                }
}


void debug_output(const All_variables* E)
{
    journal::debug_t debugInitT("initTemperature");
    debugInitT << journal::loc(__HERE__);

    for(int m=1;m<=E->sphere.caps_per_proc;m++)
        for(int k=1;k<=E->lmesh.noy;k++)
            for(int j=1;j<=E->lmesh.nox;j++)
                for(int i=1;i<=E->lmesh.noz;i++)  {
                    int node = i + (j-1)*E->lmesh.noz
                             + (k-1)*E->lmesh.noz*E->lmesh.nox;

                    double theta = E->sx[m][1][node];
                    double fi = E->sx[m][2][node];
                    double r = E->sx[m][3][node];

		    debugInitT << "(theta,fi,r,T) = "
			       << theta << "  "
			       << fi << "  "
			       << r << "  "
			       << E->T[m][node] << journal::newline;
                }
    debugInitT << journal::end;
}


// version
// $Id: misc.cc,v 1.24 2004/02/05 19:46:36 tan2 Exp $

// End of file
