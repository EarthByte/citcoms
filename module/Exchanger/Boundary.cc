// -*- C++ -*-
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  <LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#include <portinfo>
#include <vector>
#include "global_defs.h"
#include "journal/journal.h"
#include "Boundary.h"


Boundary::Boundary() :
    Exchanger::Boundary()
{}


Boundary::Boundary(const All_variables* E,
		   bool excludeTop,
		   bool excludeBottom) :
    Exchanger::Boundary(),
    bnode_(E->lmesh.nno+1, -1)
{
    journal::debug_t debug("CitcomS-Exchanger");
    debug << journal::loc(__HERE__) << journal::end;

    // boundary = all - interior
    int maxNodes = E->lmesh.nno - (E->lmesh.nox-2)
	                         * (E->lmesh.noy-2)
	                         * (E->lmesh.noz-2);
    X_.reserve(maxNodes);
    nodeID_.reserve(maxNodes);
    normal_.reserve(maxNodes);

    initX(E, excludeTop, excludeBottom);
    initBBox(E);
    bbox_.print("CitcomS-Boundary-BBox");

    X_.shrink();
    X_.print("CitcomS-Boundary-X");
    nodeID_.shrink();
    nodeID_.print("CitcomS-Boundary-nodeID");
    normal_.shrink();
    normal_.print("CitcomS-Boundary-normal");
}


Boundary::~Boundary()
{}


// private functions

void Boundary::initBBox(const All_variables *E)
{
    bbox_ = tightBBox();

    bbox_[0][2] = E->sphere.ri;
    bbox_[1][2] = E->sphere.ro;
}


void Boundary::initX(const All_variables* E,
		     bool excludeTop, bool excludeBottom)
{
    const int itop = E->lmesh.noz;

    if(E->parallel.me_loc[3] == E->parallel.nprocz - 1) {
	if(!excludeTop) {

	    for(int k=1; k<=E->lmesh.noy; k++)
		for(int j=1; j<=E->lmesh.nox; j++) {
		    std::vector<int> normalFlag(Exchanger::DIM,0);
		    normalFlag[2] = 1;
		    checkSidewalls(E, j, k, normalFlag);
		    int node = ijk2node(E, itop, j, k);
		    appendNode(E, node, normalFlag);
		}
	}
    } 
    else
	addSidewalls(E, itop, 0);

    
    const int ibottom = 1;

    if(E->parallel.me_loc[3] == 0)
	if(!excludeBottom) {
	    
	    for(int k=1; k<=E->lmesh.noy; k++)
		for(int j=1; j<=E->lmesh.nox; j++) {
		    std::vector<int> normalFlag(Exchanger::DIM,0);
		    normalFlag[2] = -1;
		    checkSidewalls(E, j, k, normalFlag);
		    int node = ijk2node(E, ibottom, j, k);
		    appendNode(E, node, normalFlag);
		}
	} 
	else
	    addSidewalls(E, ibottom, -1);
    else
	addSidewalls(E, ibottom, 0);


    // add sidewall nodes which aren't at top or bottom
    for(int i=2; i<E->lmesh.noz; i++)
	addSidewalls(E, i, 0);

}


void Boundary::addSidewalls(const All_variables* E, int znode, int r_normal)
{
    const int i = znode;

    for(int k=1; k<=E->lmesh.noy; k++)
	for(int j=1; j<=E->lmesh.nox; j++) {
	    std::vector<int> normalFlag(Exchanger::DIM,0);
	    normalFlag[2] = r_normal;
	    if(checkSidewalls(E, j, k, normalFlag)) {
		int node = ijk2node(E, i, j, k);
		appendNode(E, node, normalFlag);
	    }
	}
}


bool Boundary::checkSidewalls(const All_variables* E,
			      int j, int k, std::vector<int>& normalFlag)
{
    bool isBoundary = false;

    if((E->parallel.me_loc[1] == 0) && (j == 1)) {
	isBoundary |= true;
	normalFlag[0] = -1;
    }

    if((E->parallel.me_loc[1] == E->parallel.nprocx - 1)
       && (j == E->lmesh.nox)) {
	isBoundary |= true;
	normalFlag[0] = 1;
    }

    if((E->parallel.me_loc[2] == 0) && (k == 1)) {
	isBoundary |= true;
	normalFlag[1] = -1;
    }

    if((E->parallel.me_loc[2] == E->parallel.nprocy - 1)
       && (k == E->lmesh.noy)) {
	isBoundary |= true;
	normalFlag[1] = 1;
    }

    return isBoundary;
}


int Boundary::ijk2node(const All_variables* E, int i, int j, int k)
{
    return i + (j-1)*E->lmesh.noz + (k-1)*E->lmesh.noz*E->lmesh.nox;
}


void Boundary::appendNode(const All_variables* E,
			  int node, const std::vector<int>& normalFlag)
{
    const int m = 1;

    std::vector<double> x(Exchanger::DIM);
    for(int d=0; d<Exchanger::DIM; d++)
	x[d] = E->sx[m][d+1][node];

    bnode_[node] = X_.size();
    X_.push_back(x);
    nodeID_.push_back(node);
    normal_.push_back(normalFlag);
}


// version
// $Id: Boundary.cc,v 1.59 2005/02/02 21:49:17 tan2 Exp $

// End of file
