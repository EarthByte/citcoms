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
#include <stdlib.h>

#include "global_defs.h"
#include "Boundary.h"
#include "FineGridExchanger.h"

using namespace std;


FineGridExchanger::FineGridExchanger(MPI_Comm communicator,
				     MPI_Comm icomm,
				     int local,
				     int remote,
				     const All_variables *e):
    Exchanger(communicator, icomm, local, remote, e)
{
    std::cout << "in FineGridExchanger::FineGridExchanger" << std::endl;
}



FineGridExchanger::~FineGridExchanger() {
    std::cout << "in FineGridExchanger::~FineGridExchanger" << std::endl;
}


void FineGridExchanger::interpretate() {
    std::cout << "in FineGridExchanger::interpretate" << std::endl;
    return;
}


void FineGridExchanger::createBoundary() {
    std::cout << "in FineGridExchanger::createBoundary" << std::endl;

    if (rank == localLeader) {
      // Face nodes + Edge nodes + vertex nodes
	const int size = 2*((E->mesh.nox-2)*(E->mesh.noy-2)+(E->mesh.noy-2)*(E->mesh.noz-2)+(E->mesh.noz-2)*(E->mesh.nox-2))+4*(E->mesh.nox+E->mesh.noy+E->mesh.noz-6)+8;

	boundary = new Boundary(size);

	// initialize...
	boundary->init(E);

	//boundary->printX();
    }
}


void FineGridExchanger::sendBoundary() {
    std::cout << "in FineGridExchanger::sendBoundary"
	      << "  rank = " << rank
	      << "  leader = "<< localLeader
	      << "  receiver = "<< remoteLeader << std::endl;

    if (rank == localLeader) {
	int tag = 0;
	int itmp = boundary->size;
	MPI_Send(&itmp, 1, MPI_INT,
		 remoteLeader, tag, intercomm);

	boundary->send(intercomm, remoteLeader);
    }

    return;
}


void FineGridExchanger::mapBoundary() {
    std::cout << "in FineGridExchanger::mapBoundary" << std::endl;
    boundary->mapFineGrid(E);

    return;
}


// void FineGridExchanger::interpolate() {
  
//   double finex[3],crsex[24],xi[3],shape[8];
//   int n,node;
  
//   std::cout << "in CoarseGridExchanger::interpolate" << std::endl;
  
//   for(int i=0;i<boundary->size;i++)
//     {
//       n=boundary->bid2crseelem[i];
//       for(int j=0; j < boundary->dim ;j++)
// 	{
// 	  crsex[j]=boundary->X[j][i];
// 	  for(int k=0; k< 8 ;k++)
// 	    {
// 	      node=E->IEN[E->mesh.levmax][1][n].node[k+1];
// 	      finex[k*3+j]=E->X[E->mesh.levmax][1][j+1][node];
// 	    }
// 	}
//       std::cout << "n = " << n << " | "
// 		<< crsex[0] << " " << crsex[1] << " " << crsex[2] << " | "
// 		<< finex[0] << " " << finex[3] << " "
// 		<< finex[1] << " " << finex[10] << " "
// 		<< finex[2] << " " << finex[14] << " "
// 		<< std::endl;	
//       xi[0]=(finex[0]-crsex[0])/(finex[3]-finex[0]);
//       xi[1]=(finex[1]-crsex[1])/(finex[10]-finex[1]);
//       xi[2]=(finex[2]-crsex[2])/(finex[14]-finex[2]);
// //       std::cout << n << " " << "xi[0] = " << xi[0] << " "	
// // 		<< "xi[1] = " << xi[1] << " "	
// // 		<< "xi[2] = " << xi[2] << " "	
// // 		<< std::endl;
//       shape[0]=(1.-xi[0])*(1.-xi[1])*(1.-xi[2]);
//       shape[1]=xi[0]*(1.-xi[1])*(1.-xi[2]);
//       shape[2]=xi[0]*xi[1]*(1.-xi[2]);
//       shape[3]=(1.-xi[0])*xi[1]*(1.-xi[2]);
      
//       shape[4]=(1.-xi[0])*(1.-xi[1])*xi[2];
//       shape[5]=xi[0]*(1.-xi[1])*xi[2];
//       shape[6]=xi[0]*xi[1]*xi[2];
//       shape[7]=(1.-xi[0])*xi[1]*xi[2];
      
//       outgoing.T[i]=outgoing.v[0][i]=outgoing.v[1][i]=outgoing.v[2][i];
      
//       for(int k=0; k< 8 ;k++)
// 	{
// 	  node=E->IEN[E->mesh.levmax][1][n].node[k+1];
// 	  // 	    std::cout << "node = " << node << " "
// 	  // 		      << "k = " << k << " "
// 	  // 		      << "shape = " << shape[k] << " "
// 	  // 		      << "T = " << E->T[E->sphere.caps_per_proc][node]<<" "
// 	  // 		      << "v1 = " << E->V[E->sphere.caps_per_proc][1][node]<<" "
// 	  // 		      << "v2 = " << E->V[E->sphere.caps_per_proc][2][node]<<" "
// 	  // 		      << "v3 = " << E->V[E->sphere.caps_per_proc][3][node]<<" "
// 	  // 		      << std::endl;
// 	  outgoing.T[i]+=shape[k]*E->T[E->sphere.caps_per_proc][node];
// 	  outgoing.v[0][i]+=shape[k]*E->V[E->sphere.caps_per_proc][1][node];
// 	  outgoing.v[1][i]+=shape[k]*E->V[E->sphere.caps_per_proc][2][node];
// 	  outgoing.v[2][i]+=shape[k]*E->V[E->sphere.caps_per_proc][3][node];
// 	}
//     }
  
//   // Test
//   //     std::cout << "in CoarseGridExchanger::interpolated fields" << std::endl;
//   //     for(int i=0;i<boundary->size;i++)
//   //       {
//   // 	std::cout << i << " " << outgoing.T[i] << std::endl;
//   //       }
  
//   return;
// }

// version
// $Id: FineGridExchanger.cc,v 1.16 2003/09/27 00:27:35 tan2 Exp $

// End of file
