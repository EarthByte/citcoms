// -*- C++ -*-
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//  <LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#include <portinfo>
#include <iostream>
#include <fstream>
#include <stdio.h>

#include "Array2D.h"
#include "Array2D.cc"
#include "Boundary.h"
#include "Mapping.h"
#include "CoarseGridExchanger.h"
#include "global_defs.h"


CoarseGridExchanger::CoarseGridExchanger(const MPI_Comm comm,
					 const MPI_Comm intercomm,
					 const int leader,
					 const int localLeader,
					 const int remoteLeader,
					 const All_variables *E):
    Exchanger(comm, intercomm, leader, localLeader, remoteLeader, E)
{
    std::cout << "in CoarseGridExchanger::CoarseGridExchanger" << std::endl;
}

CoarseGridExchanger::~CoarseGridExchanger() {
    std::cout << "in CoarseGridExchanger::~CoarseGridExchanger" << std::endl;
}


void CoarseGridExchanger::gather() {
    std::cout << "in CoarseGridExchanger::gather" << std::endl;

    interpretate();

    if (rank == leader) {
	const int dim = 3;
	int size = boundary->size();

	auto_array_ptr<double> tmp(new double[dim*size]);

	int nproc;
	MPI_Comm_size(comm, &nproc);
	for (int i=0; i<nproc; i++) {
	    Array2D<dim> recV(size);
	    Array2D<3>* V = localV.get();

	    if (i != leader) {
		recV.receive(comm, i);
		V = &recV;
	    }

	    for (int n=0; n<size; n++)
		if (cgmapping->bid2proc(n) == i)
		    for (int d=0; d<dim; d++)
			tmp[n*dim+d] = (*V)(d,n);
	}

	outgoingV = Velo(new Array2D<dim>(tmp, size));
	//outgoingV->print("outgoingV");
    }
    else {
	localV->send(comm, leader);
    }
}


void CoarseGridExchanger::distribute() {
    std::cout << "in CoarseGridExchanger::distribute" << std::endl;
}



void CoarseGridExchanger::interpretate() {
    std::cout << "in CoarseGridExchanger::interpretate" << std::endl;
    // interpolate velocity field to boundary nodes
    const int dim = 3;
    const int size = boundary->size();
    auto_array_ptr<double> tmp(new double[dim*size]);

    for(int i=0; i<size; i++) {
	int n1 = cgmapping->bid2elem(i);
	for(int d=0; d<dim; d++)
	    tmp[i*dim+d] = 0;

	if(n1 != 0) {
	    for(int mm=1; mm<=E->sphere.caps_per_proc; mm++)
		for(int k=0; k<8; k++) {
		    int node = E->IEN[E->mesh.levmax][mm][n1].node[k+1];
		    for(int d=0; d<dim; d++)
			tmp[i*dim+d] += cgmapping->shape(i*8+k) * E->sphere.cap[mm].V[d+1][node];
		}
	}
    }
    localV = Velo(new Array2D<dim>(tmp, size));
    //localV->print("localV");
}


void CoarseGridExchanger::mapBoundary() {
    std::cout << "in CoarseGridExchanger::mapBoundary" << std::endl;

    createMapping();
    createDataArrays();
}


void CoarseGridExchanger::createMapping() {
    cgmapping = new CoarseGridMapping(boundary, E, comm, rank, leader);
    mapping = cgmapping;
}


void CoarseGridExchanger::receiveBoundary() {
    std::cout << "in CoarseGridExchanger::receiveBoundary"
	      << "  rank = " << rank
	      << "  leader = "<< localLeader
	      << "  sender = "<< remoteLeader << std::endl;
    int size;

    if (rank == leader) {
	int tag = 0;
	MPI_Status status;

 	MPI_Recv(&size, 1, MPI_INT,
 		 remoteLeader, tag, intercomm, &status);

	boundary = new Boundary(size);
	boundary->receive(intercomm, remoteLeader);
    }

    // Broadcast info received by localLeader to the other procs
    // in the Coarse communicator.
    MPI_Bcast(&size, 1, MPI_INT, leader, comm);

    if (rank != leader)
	boundary = new Boundary(size);

    boundary->broadcast(comm, leader);
}


void CoarseGridExchanger::interpolateTemperature() {
  std::cout << "in CoarseGridExchanger::interpolateTemperature" << std::endl;

  int n1,n2,node;
  for(int i=0;i<boundary->size();i++) {
      n1 = cgmapping->bid2elem(i);
      n2 = cgmapping->bid2proc(i);

      //outgoing.T[i] = 0;
      if(n1!=0) {
	for(int mm=1;mm<=E->sphere.caps_per_proc;mm++)
	  for(int k=0; k< 8 ;k++)
	    {
	      node=E->IEN[E->mesh.levmax][mm][n1].node[k+1];
	      //outgoing.T[i]+=boundary->shape[k]*E->T[mm][node];
	    }
      }
    }
}

// version
// $Id: CoarseGridExchanger.cc,v 1.29 2003/10/11 00:38:46 tan2 Exp $

// End of file
