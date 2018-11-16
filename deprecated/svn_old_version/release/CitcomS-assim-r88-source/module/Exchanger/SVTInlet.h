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
// Role:
//     impose normal velocity and shear traction as velocity b.c.,
//     also impose temperature as temperature b.c.
//
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//

#if !defined(pyCitcomSExchanger_SVTInlet_h)
#define pyCitcomSExchanger_SVTInlet_h

#include "BaseSVTInlet.h"

struct All_variables;


class SVTInlet : public BaseSVTInlet{

public:
    SVTInlet(const Boundary& boundary,
	     const Exchanger::Sink& sink,
	     All_variables* E);

    virtual ~SVTInlet();

    virtual void recv();
};


#endif

// version
// $Id: SVTInlet.h 15114 2009-06-03 21:08:01Z tan2 $

// End of file
