// -*- C++ -*-
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//
//  <LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~//

#if !defined(pyCitcom_Exchanger_h)
#define pyCitcom_Exchanger_h

#include "mpi.h"


class Boundary;     // declaration only
struct All_variables;

struct Data {
    static const int npass = 8;  // # of arrays to pass
    int size;                    // length of each array
    double *x, *y, *z;    // coordinates
    double *u, *v, *w;    // velocities
    double *T, *P;       // temperature and pressure
};



class Exchanger {
  
public:
    Exchanger(const MPI_Comm communicator,
	      const MPI_Comm intercomm,
	      const int localLeader,
	      const int remoteLeader,
	      const All_variables *E);
    virtual ~Exchanger();
  
    void reset_target(const MPI_Comm intercomm,
		      const int receiver);

//     virtual void send(int& size);
//     virtual void receive(const int size);
    void sendTemperature();
    void receiveTemperature();
    void sendVelocities();
    void receiveVelocities();
//     virtual void inter_sendTemperature();
//     virtual void inter_receiveTemperature();
    void inter_sendVelocities();
    void inter_receiveVelocities();
    double exchangeTimestep(const double);

    void wait();
    void nowait();

    virtual void gather() = 0;
    virtual void distribute() = 0;
    virtual void interpretate() = 0; // interpolation or extrapolation
    virtual void impose_bc() = 0;    // set bc flag

    virtual void mapBoundary() = 0;
                                     // create mapping from Boundary object
                                     // to global id array

protected:
    const MPI_Comm comm;
    const MPI_Comm intercomm;

    const int localLeader;
    const int remoteLeader;

    const All_variables *E;    // CitcomS data structure,
                               // Exchanger only modifies bc flags

    Boundary *boundary;

    Data outgoing;
    Data incoming;

        //Data inter_

    int rank;

private:
    // disable copy constructor and copy operator
    Exchanger(const Exchanger&);
    Exchanger operator=(const Exchanger&);

};



#endif

// version
// $Id: ExchangerClass.h,v 1.9 2003/09/18 17:16:44 puru Exp $

// End of file

