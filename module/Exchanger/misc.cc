// -*- C++ -*-
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  <LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#include <portinfo>
#include <Python.h>
#include <iostream>
#include <fstream>

extern "C" {
#include "global_defs.h"
#include "element_definitions.h"
}

#include "misc.h"
void commonE(All_variables*);


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

PyObject * pyExchanger_FinereturnE(PyObject *, PyObject *)
{
    All_variables *E = new All_variables;

    E->parallel.nproc = 1;

    E->mesh.nox = 4;
    E->mesh.noy = 4;
    E->mesh.noz = 3;

    E->control.theta_max = 2.1;
    E->control.theta_min = 0.9;
    E->control.fi_max = 2.1;
    E->control.fi_min = 0.9;
    E->sphere.ro = 2.0;
    E->sphere.ri = 0.9;

    commonE(E);

    PyObject *cobj = PyCObject_FromVoidPtr(E, NULL);
    return Py_BuildValue("O", cobj);
}


char pyExchanger_CoarsereturnE__doc__[] = "";
char pyExchanger_CoarsereturnE__name__[] = "CoarsereturnE";

PyObject * pyExchanger_CoarsereturnE(PyObject *, PyObject *)
{
    All_variables *E = new All_variables;

    E->parallel.nproc = 12;

    E->mesh.nox = 4;
    E->mesh.noy = 4;
    E->mesh.noz = 3;

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


void commonE(All_variables *E) {
    E->parallel.me = 1;
    E->parallel.me_loc[1] = 0;
    E->parallel.me_loc[2] = 0;
    E->parallel.me_loc[3] = 0;
    E->parallel.nprocx = 1;
    E->parallel.nprocy = 1;
    E->parallel.nprocz = 1;

    E->sphere.caps_per_proc = 1;

    E->mesh.levmax = 1;
    E->mesh.levmin = 1;

    E->mesh.dof = 3;

    E->mesh.elx = E->mesh.nox - 1;
    E->mesh.ely = E->mesh.noy - 1;
    E->mesh.elz = E->mesh.noz - 1;

    E->lmesh.elx = E->mesh.elx/E->parallel.nprocx;
    E->lmesh.elz = E->mesh.elz/E->parallel.nprocz;
    E->lmesh.ely = E->mesh.ely/E->parallel.nprocy;
    E->lmesh.nox = E->lmesh.elx + 1;
    E->lmesh.noz = E->lmesh.elz + 1;
    E->lmesh.noy = E->lmesh.ely + 1;

    E->lmesh.nno = E->lmesh.noz*E->lmesh.nox*E->lmesh.noy;
    E->lmesh.nel = E->lmesh.ely*E->lmesh.elx*E->lmesh.elz;
    E->lmesh.npno = E->lmesh.nel;

    int noz = E->lmesh.noz;
    int noy = E->mesh.noy;
    int nox = E->mesh.nox;

    E->lmesh.ELX[E->mesh.levmax] = nox-1;
    E->lmesh.ELY[E->mesh.levmax] = noy-1;
    E->lmesh.ELZ[E->mesh.levmax] = noz-1;
    E->lmesh.NOZ[E->mesh.levmax] = noz;
    E->lmesh.NOY[E->mesh.levmax] = noy;
    E->lmesh.NOX[E->mesh.levmax] = nox;
    E->lmesh.NNO[E->mesh.levmax] = nox * noz * noy;
    E->lmesh.NEL[E->mesh.levmax] = (nox-1) * (noz-1) * (noy-1);

    for (int j=1;j<=E->sphere.caps_per_proc;j++)  {
	E->sphere.capid[j] = 1;
    }

    for (int lev=E->mesh.levmax;lev>=E->mesh.levmin;lev--)  {
	for (int j=1;j<=E->sphere.caps_per_proc;j++)  {
	    E->IEN[lev][j] = new IEN [E->lmesh.nel+1];
	}
    }

    for(int m=1;m<=E->sphere.caps_per_proc;m++) {
	E->sphere.cap[m].VB[1] = new float[E->lmesh.nno+1];
	E->sphere.cap[m].VB[2] = new float[E->lmesh.nno+1];
	E->sphere.cap[m].VB[3] = new float[E->lmesh.nno+1];
	E->node[m] = new unsigned int[E->lmesh.nno+1];

	for(int n=1; n<=E->lmesh.nno; n++) {
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
			(E->control.theta_max - E->control.theta_min)/(E->lmesh.nox-1)*(j-1) + E->control.theta_min;
		    E->sx[m][2][node] =
			(E->control.fi_max - E->control.fi_min)/(E->lmesh.noy-1)*(k-1) + E->control.fi_min;
		    E->sx[m][3][node] =
			(E->sphere.ro -  E->sphere.ri)/(E->lmesh.noz-1)*(i-1) +  E->sphere.ri;

//  		    std::cout <<  node << " "
//  			      << E->sx[m][1][node] << " "
//  			      << E->sx[m][2][node] << " "
//  			      << E->sx[m][3][node] << " "
//  			      << std::endl;
		}

    return;
}

// version
// $Id: misc.cc,v 1.18 2003/09/29 22:31:35 puru Exp $

// End of file
