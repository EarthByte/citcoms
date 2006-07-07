#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#<LicenseText>
#
# CitcomS.py by Eh Tan, Eun-seo Choi, and Pururav Thoutireddy.
# Copyright (C) 2002-2005, California Institute of Technology.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#</LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from CitcomSLib import CPU_time
from CitcomS.Components.CitcomComponent import CitcomComponent

class Sphere(CitcomComponent):



    def setup(self):
        return



    def run(self):
        start_time = CPU_time()
        self.launch()

        import mpi
        if not mpi.world().rank:
            import sys
            print >> sys.stderr, "initialization time = %f" % \
                  (CPU_time() - start_time)

	return



    def launch(self):
	raise NotImplementedError, "not implemented"
        return



    def setProperties(self):
        from CitcomSLib import Sphere_set_properties
        Sphere_set_properties(self.all_variables, self.inventory)
        return



    class Inventory(CitcomComponent.Inventory):

        import pyre.inventory


        nprocx = pyre.inventory.int("nprocx", default=1)
        nprocy = pyre.inventory.int("nprocy", default=1)
        nprocz = pyre.inventory.int("nprocz", default=1)

        coor = pyre.inventory.bool("coor", default=False)
        coor_file = pyre.inventory.str("coor_file", default="coor.dat")

        nodex = pyre.inventory.int("nodex", default=9)
        nodey = pyre.inventory.int("nodey", default=9)
        nodez = pyre.inventory.int("nodez", default=9)
        levels = pyre.inventory.int("levels", default=1)

        radius_outer = pyre.inventory.float("radius_outer", default=1.0)
        radius_inner = pyre.inventory.float("radius_inner", default=0.55)

	    # these parameters are for spherical harmonics output
	    # put them here temporalily
        ll_max = pyre.inventory.int("ll_max", default=20)
        nlong = pyre.inventory.int("nlong", default=361)
        nlati = pyre.inventory.int("nlati", default=181)
        output_ll_max = pyre.inventory.int("output_ll_max", default=20)




# version
__id__ = "$Id$"

# End of file
