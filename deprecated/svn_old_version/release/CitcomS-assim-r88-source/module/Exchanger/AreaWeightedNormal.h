// -*- C++ -*-
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//<LicenseText>
//
// CitcomS.py by Eh Tan, Eun-seo Choi, and Pururav Thoutireddy.
// Copyright (C) 2002-2005, California Institute of Technology.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//</LicenseText>
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#if !defined(pyCitcomSExchanger_AreaWeightedNormal_h)
#define pyCitcomSExchanger_AreaWeightedNormal_h

#include <vector>
#include "mpi.h"
#include "Exchanger/Array2D.h"
#include "Exchanger/DIM.h"
#include "Exchanger/Sink.h"

struct All_variables;
class Boundary;


class AreaWeightedNormal {
    Exchanger::Array2D<double,Exchanger::DIM> nwght;

public:
    AreaWeightedNormal(const MPI_Comm& comm,
		       const Boundary& boundary,
		       const Exchanger::Sink& sink,
                       const All_variables* E);
    ~AreaWeightedNormal();

    typedef Exchanger::Array2D<double,Exchanger::DIM> Velo;

    void imposeConstraint(Velo& V,
			  const MPI_Comm& comm,
			  const Exchanger::Sink& sink,
			  const All_variables* E) const;

private:
    void computeWeightedNormal(const Boundary& boundary,
			       const All_variables* E);
    double computeTotalArea(const MPI_Comm& comm,
                            const Exchanger::Sink& sink) const;
    void normalize(double total_area);
    double computeOutflow(const Velo& V,
			  const MPI_Comm& comm,
			  const Exchanger::Sink& sink) const;
    inline int sign(double number) const;
    void reduceOutflow(Velo& V, double outflow,
		       const Exchanger::Sink& sink) const;

};


#endif

// version
// $Id: AreaWeightedNormal.h 15108 2009-06-02 22:56:46Z tan2 $

// End of file

