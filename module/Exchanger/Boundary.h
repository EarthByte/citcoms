// -*- C++ -*-
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//  <LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#if !defined(pyCitcom_Boundary_h)
#define pyCitcom_Boundary_h

#include <memory>
#include "mpi.h"

struct All_variables;


class Boundary {

public:
    static const int dim = 3;  // spatial dimension
    const int size;            // # of boundary nodes

    double theta_max, theta_min,
	   fi_max, fi_min,
	   ro, ri;   // domain bound of FG
    double *X[dim];  // coordinate
    int *bid2gid;    // bid (local id) -> ID (ie. global id in FG)
    int *bid2elem;   // bid -> elem from which fields are interpolated in CG
    int *bid2proc;   // bid -> proc. rank
    double *shape;   // shape functions for interpolation

    explicit Boundary(const int n);     // constructor, allocating memory only
    ~Boundary();

    void init(const All_variables *E);  // initialize X and domain bound
    void mapFineGrid(const All_variables *E);  // initialize bid2gid
    void mapCoarseGrid(const All_variables *E, const int lrank);  // initialize shape, bid2elem and bid2proc

    void send(const MPI_Comm comm, const int receiver) const;
    void receive(const MPI_Comm comm, const int sender);
    void broadcast(const MPI_Comm comm, const int broadcaster);

    void sendBid2proc(const MPI_Comm comm,
		      const int lrank, const int leader);
                                         // send bid2proc to leader

    void printX() const;
    void printBid2gid() const;
    void printBid2proc() const;
    void printBid2elem() const;
    void printBound() const;
    
private:
//     Boundary(const int,
// 	     std::auto_ptr<int>,
// 	     std::auto_ptr<double>,
// 	     std::auto_ptr<double>,
// 	     std::auto_ptr<double>);

//     const std::auto_ptr<double> X_[dim];
//     std::auto_ptr<int> bid2gid_;
//     std::auto_ptr<int> bid2proc_;  

    void testMapping(const All_variables *E) const;
    double Tetrahedronvolume(double  *x1, double *x2, double *x3, double *x4) const;
    double det3_sub(double  *x1, double *x2, double *x3) const;


};

#endif

// version
// $Id: Boundary.h,v 1.15 2003/09/28 00:11:03 tan2 Exp $

// End of file
