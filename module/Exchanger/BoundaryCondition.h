// -*- C++ -*-
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  <LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#if !defined(pyCitcom_BoundaryCondition_h)
#define pyCitcom_BoundaryCondition_h

#include "AreaWeightedNormal.h"
#include "Array2D.h"
#include "DIM.h"

struct All_variables;
class Boundary;
class Sink;
class Source;


class BoundaryConditionSink {
    MPI_Comm comm;
    All_variables* E;
    const Boundary& boundary;
    const Sink& sink;
    AreaWeightedNormal awnormal;
    Array2D<double,DIM> vbc;
    Array2D<double,DIM> old_vbc;
    Array2D<double,1> tbc;
    Array2D<double,1> old_tbc;
    double fge_t;
    double cge_t;

public:
    BoundaryConditionSink(const MPI_Comm& comm,
			  const Boundary& boundary,
			  const Sink& sink,
			  All_variables* E);
    ~BoundaryConditionSink() {};

    void recvTandV();
    void imposeBC();
    void storeTimestep(double fge_t, double cge_t);

private:
    void setVBCFlag();
    void setTBCFlag();
    void imposeConstraint();
    void imposeTBC();
    void imposeVBC();

};


class BoundaryConditionSource {
    All_variables* E;
    const Source& source;
    Array2D<double,DIM> vbc;
    Array2D<double,1> tbc;
    Array2D<double,DIM> fbc;

public:
    BoundaryConditionSource(const Source& src, All_variables* E);
    ~BoundaryConditionSource() {};

    void sendTandV();
    void sendTraction();

private:

};


#endif

// version
// $Id: BoundaryCondition.h,v 1.4 2003/11/23 18:21:48 ces74 Exp $

// End of file
