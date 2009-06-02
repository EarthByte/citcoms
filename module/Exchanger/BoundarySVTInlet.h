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

#if !defined(pyCitcomSExchanger_BoundaryVTInlet_h)
#define pyCitcomSExchanger_BoundaryVTInlet_h

#include "mpi.h"
#include "BaseSVTInlet.h"

struct All_variables;
class AreaWeightedNormal;


class BoundarySVTInlet : public BaseSVTInlet{
private:
    const MPI_Comm& comm;
    AreaWeightedNormal* awnormal;

public:
    BoundarySVTInlet(const Boundary& boundary,
                     const Exchanger::Sink& sink,
                     All_variables* E,
                     const MPI_Comm& comm);
    virtual ~BoundarySVTInlet();

    virtual void recv();

};


#endif

// version
// $Id: BoundaryVTInlet.h 2397 2005-10-04 22:37:25Z leif $

// End of file
