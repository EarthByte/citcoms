# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

include local.def
TYPE=Regional

PROJECT = CitcomS
PACKAGE = $(TYPE)module
include std-pythonmodule.def

PROJ_CXX_SRCLIB = $(BLD_LIBDIR)/lib$(PROJECT)$(TYPE).$(EXT_LIB) 
PROJ_CXX_INCLUDES = ../../lib/Common ../../lib/$(TYPE)
#PROJ_CXX_INCLUDES = $(BLD_INCDIR)/$(PROJECT)/$(TYPE)

PROJ_SRCS = \
    bindings.cc \
    exceptions.cc \
    misc.cc \
    outputs.cc

# version
# $Id: Make.mm,v 1.6 2003/05/22 18:21:38 tan2 Exp $

# End of file
