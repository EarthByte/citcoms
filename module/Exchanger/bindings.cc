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

#include "config.h"
#include <Python.h>

#include "Exchanger/exchangers.h"
#include "exchangers.h"
#include "inlets_outlets.h"
#include "misc.h"          // miscellaneous methods

#include "bindings.h"

// the method table

struct PyMethodDef pyExchanger_methods[] = {

    // from misc.h

    {PyCitcomSExchanger_copyright__name__,
     PyCitcomSExchanger_copyright,
     METH_VARARGS,
     PyCitcomSExchanger_copyright__doc__},

    {PyCitcomSExchanger_FinereturnE__name__,
     PyCitcomSExchanger_FinereturnE,
     METH_VARARGS,
     PyCitcomSExchanger_FinereturnE__doc__},

    {PyCitcomSExchanger_CoarsereturnE__name__,
     PyCitcomSExchanger_CoarsereturnE,
     METH_VARARGS,
     PyCitcomSExchanger_CoarsereturnE__doc__},


    // from inlets_outlets.h

    {PyCitcomSExchanger_SVTInlet_create__name__,
     PyCitcomSExchanger_SVTInlet_create,
     METH_VARARGS,
     PyCitcomSExchanger_SVTInlet_create__doc__},

    {PyCitcomSExchanger_TInlet_create__name__,
     PyCitcomSExchanger_TInlet_create,
     METH_VARARGS,
     PyCitcomSExchanger_TInlet_create__doc__},

    {PyCitcomSExchanger_SInlet_create__name__,
     PyCitcomSExchanger_SInlet_create,
     METH_VARARGS,
     PyCitcomSExchanger_SInlet_create__doc__},

    {PyCitcomSExchanger_VTInlet_create__name__,
     PyCitcomSExchanger_VTInlet_create,
     METH_VARARGS,
     PyCitcomSExchanger_VTInlet_create__doc__},

    {PyCitcomSExchanger_SVTOutlet_create__name__,
     PyCitcomSExchanger_SVTOutlet_create,
     METH_VARARGS,
     PyCitcomSExchanger_SVTOutlet_create__doc__},

    {PyCitcomSExchanger_TOutlet_create__name__,
     PyCitcomSExchanger_TOutlet_create,
     METH_VARARGS,
     PyCitcomSExchanger_TOutlet_create__doc__},

    {PyCitcomSExchanger_VOutlet_create__name__,
     PyCitcomSExchanger_VOutlet_create,
     METH_VARARGS,
     PyCitcomSExchanger_VOutlet_create__doc__},

    {PyCitcomSExchanger_VTOutlet_create__name__,
     PyCitcomSExchanger_VTOutlet_create,
     METH_VARARGS,
     PyCitcomSExchanger_VTOutlet_create__doc__},

    // from exchangers.h

    {PyCitcomSExchanger_createBoundary__name__,
     PyCitcomSExchanger_createBoundary,
     METH_VARARGS,
     PyCitcomSExchanger_createBoundary__doc__},

    {PyCitcomSExchanger_createEmptyBoundary__name__,
     PyCitcomSExchanger_createEmptyBoundary,
     METH_VARARGS,
     PyCitcomSExchanger_createEmptyBoundary__doc__},

    {PyCitcomSExchanger_createEmptyInterior__name__,
     PyCitcomSExchanger_createEmptyInterior,
     METH_VARARGS,
     PyCitcomSExchanger_createEmptyInterior__doc__},

    {PyCitcomSExchanger_createGlobalBoundedBox__name__,
     PyCitcomSExchanger_createGlobalBoundedBox,
     METH_VARARGS,
     PyCitcomSExchanger_createGlobalBoundedBox__doc__},

    {PyCitcomSExchanger_createInterior__name__,
     PyCitcomSExchanger_createInterior,
     METH_VARARGS,
     PyCitcomSExchanger_createInterior__doc__},

    {PyCitcomSExchanger_initConvertor__name__,
     PyCitcomSExchanger_initConvertor,
     METH_VARARGS,
     PyCitcomSExchanger_initConvertor__doc__},

    {PyCitcomSExchanger_initTemperature__name__,
     PyCitcomSExchanger_initTemperature,
     METH_VARARGS,
     PyCitcomSExchanger_initTemperature__doc__},

    {PyCitcomSExchanger_modifyT__name__,
     PyCitcomSExchanger_modifyT,
     METH_VARARGS,
     PyCitcomSExchanger_modifyT__doc__},

    {PyCitcomSExchanger_CitcomSource_create__name__,
     PyCitcomSExchanger_CitcomSource_create,
     METH_VARARGS,
     PyCitcomSExchanger_CitcomSource_create__doc__},

    // from Exchanger/exchangers.h

    {PyExchanger_exchangeBoundedBox__name__,
     PyExchanger_exchangeBoundedBox,
     METH_VARARGS,
     PyExchanger_exchangeBoundedBox__doc__},

    {PyExchanger_exchangeSignal__name__,
     PyExchanger_exchangeSignal,
     METH_VARARGS,
     PyExchanger_exchangeSignal__doc__},

    {PyExchanger_exchangeTimestep__name__,
     PyExchanger_exchangeTimestep,
     METH_VARARGS,
     PyExchanger_exchangeTimestep__doc__},

    {PyExchanger_Inlet_impose__name__,
     PyExchanger_Inlet_impose,
     METH_VARARGS,
     PyExchanger_Inlet_impose__doc__},

    {PyExchanger_Inlet_recv__name__,
     PyExchanger_Inlet_recv,
     METH_VARARGS,
     PyExchanger_Inlet_recv__doc__},

    {PyExchanger_Inlet_storeTimestep__name__,
     PyExchanger_Inlet_storeTimestep,
     METH_VARARGS,
     PyExchanger_Inlet_storeTimestep__doc__},

    {PyExchanger_Outlet_send__name__,
     PyExchanger_Outlet_send,
     METH_VARARGS,
     PyExchanger_Outlet_send__doc__},

    {PyExchanger_Sink_create__name__,
     PyExchanger_Sink_create,
     METH_VARARGS,
     PyExchanger_Sink_create__doc__},

    // Sentinel
    {0, 0}
};

// version
// $Id$

// End of file
