# -*- Makefile -*-
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#                               Michael A.G. Aivazis
#                                  Steve Quenette
#                        California Institute of Technology
#                        (C) 1998-2003  All Rights Reserved
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# version
# $Id: Make.mm,v 1.1 2003/03/24 01:46:37 tan2 Exp $

include local.def
TYPE=Full

PROJECT = CitcomS/$(TYPE)
PACKAGE = lib/$(TYPE)
PROJ_LIB = $(BLD_LIBDIR)/libCitcomS$(TYPE).$(EXT_LIB)

PROJ_LIBRARIES = $(EXTERNAL_LIBPATH) $(EXTERNAL_LIBS) -lm
PROJ_CC_INCLUDES = ../Common

PROJ_SRCS = \
	Boundary_conditions.c \
	Geometry_cartesian.c \
	Global_operations.c \
	Initial_temperature.c \
	Output.c \
	Parallel_related.c \
	Problem_related.c \
	Process_buoyancy.c \
	Sphere_related.c \
	Version_dependent.c

#EXPORT_HEADERS = \
#	global_defs.h

all: $(PROJ_LIB)


# version
# $Id: Make.mm,v 1.1 2003/03/24 01:46:37 tan2 Exp $

#
# End of file
