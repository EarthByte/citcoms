// -*- C++ -*-
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  <LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#if !defined(pyCitcom_VTSource_h)
#define pyCitcom_VTSource_h

#include <vector>
#include "mpi.h"
#include "utility.h"
#include "AbstractSource.h"
#include "Array2D.h"
#include "BoundedBox.h"
#include "DIM.h"

struct All_variables;
class BoundedMesh;
class VTInterpolator;


class VTSource : public AbstractSource {

public:
    VTSource(MPI_Comm comm, int sink,
	     BoundedMesh& mesh, const All_variables* E,
	     const BoundedBox& mybbox);
    virtual ~VTSource();

    virtual void interpolateTemperature(Array2D<double,1>& T) const;
    virtual void interpolateVelocity(Array2D<double,DIM>& V) const;

private:
    void recvMesh(BoundedMesh& mesh);
    void sendMeshNode() const;
    void initX(const BoundedMesh& mesh);

    virtual void interpolateForce(Array2D<double,DIM>& F) const {};
    virtual void interpolatePressure(Array2D<double,1>& P) const {};
    virtual void interpolateStress(Array2D<double,DIM>& S) const {};
    virtual void interpolateTraction(Array2D<double,DIM>& F) const {};

    // disable copy c'tor and assignment operator
    VTSource(const VTSource&);
    VTSource& operator=(const VTSource&);

};


#endif

// version
// $Id: VTSource.h,v 1.1 2004/02/24 20:20:32 tan2 Exp $

// End of file
