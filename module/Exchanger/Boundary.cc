// -*- C++ -*-
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//  <LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#include <iostream>
#include "global_defs.h"
#include "Boundary.h"
#include "ExchangerClass.h"

//using std::auto_ptr;

Boundary::Boundary(const int n) : size(n){
    std::cout << "in Boundary::Boundary  size = " << size << std::endl;

    connectivity = new int[size];
    for(int i=0; i<dim; i++)
	X[i] = new double[size];

    bid2proc = new int[size];
    bid2gid = new int[size];

    // use auto_ptr for exception-proof
    //Boundary(n, auto_ptr<int>(new int[size]),
    //     auto_ptr<double>(new double[size]),
    //     auto_ptr<double>(new double[size]),
    //     auto_ptr<double>(new double[size]));
}


// // because initialization of X[dim] involves a loop,
// // here we assume dim=3 and use a private constructor

// Boundary::Boundary(const int n,
// 		   auto_ptr<int> c,
// 		   auto_ptr<double> x0,
// 		   auto_ptr<double> x1,
// 		   auto_ptr<double> x2)
//     :    size(n), connectivity_(c),
//         X_[0](x0), X_[1](x1), X_[2](x2)
// {
//     std::cout << "in Boundary::Boundary  size = " << size << std::endl;
//     assert(dim == 3);

//     // setup traditional pointer for convenience
//     connectivity = connectivity_.get();
//     for(int i=0; i<dim; i++)
// 	X[i] = X_[i].get();

// }



Boundary::~Boundary() {
    std::cout << "in Boundary::~Boundary" << std::endl;

    delete [] connectivity;
    for(int i=0; i<dim; i++)
	delete [] X[i];

    delete [] bid2proc;
    delete [] bid2gid;
};


void Boundary::init(const All_variables *E) {
    int nodes,node1,node2,nodest;
    int *nid;

    nodest = E->lmesh.nox * E->lmesh.noy * E->lmesh.noz;
    nid = new int[nodest];
    for(int i=0;i<nodest;i++)nid[i]=0;

    nodes=0;
        // test
    for(int j=0; j<size; j++)
        connectivity[j] = j;

        //  for two YOZ planes

    if (E->parallel.me_loc[1]==0 || E->parallel.me_loc[1]==E->parallel.nprocx-1)        for (int m=1;m<=E->sphere.caps_per_proc;m++)
            for(int j=1;j<=E->lmesh.noy;j++)
                for(int i=1;i<=E->lmesh.noz;i++)  {
                    node1 = i + (j-1)*E->lmesh.noz*E->lmesh.nox;
                    node2 = node1 + (E->lmesh.nox-1)*E->lmesh.noz;

                    if ((E->parallel.me_loc[1]==0) && (!nid[node1-1]))  {
                        for(int k=0;k<dim;k++)X[k][nodes]=E->X[E->mesh.levmax][1][k+1][node1];
                        nid[node1-1]++;
                        nodes++;
                    }
                    if ((E->parallel.me_loc[1]==E->parallel.nprocx-1) && (!nid[node2-1])) {
                        for(int k=0;k<dim;k++)X[k][nodes]=E->X[E->mesh.levmax][1][k+1][node2];
                        nid[node2-1]++;
                        nodes++;
                    }
                }

        //  for two XOZ planes

    if (E->parallel.me_loc[2]==0 || E->parallel.me_loc[2]==E->parallel.nprocy-1)        for (int m=1;m<=E->sphere.caps_per_proc;m++)
            for(int j=1;j<=E->lmesh.nox;j++)
                for(int i=1;i<=E->lmesh.noz;i++)  {
                    node1 = i + (j-1)*E->lmesh.noz;
                    node2 = node1 + (E->lmesh.noy-1)*E->lmesh.noz*E->lmesh.nox;
                    if ((E->parallel.me_loc[2]==0) && (!nid[node1-1]))  {
                        for(int k=0;k<dim;k++)X[k][nodes]=E->X[E->mesh.levmax][1][k+1][node1];
                        nid[node1-1]++;
                        nodes++;
                    }
                    if((E->parallel.me_loc[2]==E->parallel.nprocy-1)&& (!nid[node2-1]))  {
                        for(int k=0;k<dim;k++)X[k][nodes]=E->X[E->mesh.levmax][1][k+1][node2];
                        nid[node2-1]++;
                        nodes++;
                    }
                }
        //  for two XOY planes

    if (E->parallel.me_loc[3]==0 || E->parallel.me_loc[3]==E->parallel.nprocz-1)        for (int m=1;m<=E->sphere.caps_per_proc;m++)
            for(int j=1;j<=E->lmesh.noy;j++)
                for(int i=1;i<=E->lmesh.nox;i++)  {
                    node1 = 1 + (i-1)*E->lmesh.noz+(j-1)*E->lmesh.nox*E->lmesh.noz;
                    node2 = node1 + E->lmesh.noz-1;

                    if ((E->parallel.me_loc[3]==0 ) && (!nid[node1-1])) {
                        for(int k=0;k<dim;k++)X[k][nodes]=E->X[E->mesh.levmax][1][k+1][node1];
                        nid[node1-1]++;
                        nodes++;
                    }
                    if ((E->parallel.me_loc[3]==E->parallel.nprocz-1) &&(!nid[node2-1])) {
                        for(int k=0;k<dim;k++)X[k][nodes]=E->X[E->mesh.levmax][1][k+1][node2];
                        nid[node2-1]++;
                        nodes++;
                    }
                }
    if(nodes != size) std::cout << " nodes != size ";

        //    for(int i=0; i<dim; i++)
        //      for(int j=0; j<size; j++) {
        //          X[i][j] = i+j;
        //      }
    delete nid;
}



void Boundary::mapFineGrid(const All_variables *E, int localLeader) {

    for(int i=0; i<size; i++)
        bid2proc[i]=localLeader;

    int n=0;
    for(int m=1;m<=E->sphere.caps_per_proc;m++)
	for(int k=1;k<=E->lmesh.noy;k++)
	    for(int j=1;j<=E->lmesh.nox;j++)
		for(int i=1;i<=E->lmesh.noz;i++)  {
		    int node = i + (j-1)*E->lmesh.noz
                        + (k-1)*E->lmesh.noz*E->lmesh.nox;

		    if((k==1)||(k==E->lmesh.noy)||(j==1)||(j==E->lmesh.nox)||(i==1)||(i==E->lmesh.noz))
		    {
                        bid2gid[n]=node;
                        n++;
		    }
		}
    if(n != size) std::cout << " nodes != size ";
    printBid2gid();
}



void Boundary::printConnectivity() const {
    int *c = connectivity;
    for(int j=0; j<size; j++)
	std::cout << "  C:  " << j << ":  " << c[j] << std::endl;
}


void Boundary::printX() const {
    for(int j=0; j<size; j++) {
	std::cout << "  X:  " << j << ":  ";
	    for(int i=0; i<dim; i++)
		std::cout << X[i][j] << " ";
	std::cout << std::endl;
    }
}


void Boundary::printBid2gid() const {
    int *c = bid2gid;
    for(int j=0; j<size; j++)
	std::cout << "  B:  " << j << ":  " << c[j] << std::endl;
}



// version
// $Id: Boundary.cc,v 1.7 2003/09/18 16:12:15 puru Exp $

// End of file
