#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from pyre.components.Component import Component

class Parallel(Component):


    def __init__(self):
        Component.__init__(self, "parallel", "parallel")
        return



    def setProperties(self):
        import CitcomS.Regional as Regional
	Regional.Parallel_set_properties(self.inventory)
        return



    class Inventory(Component.Inventory):


        import pyre.properties

        inventory = [

            pyre.properties.list("nproc_surf",1,
				 pyre.properties.choice([1,12])
				 ),
            pyre.properties.int("nprocx",1),
            pyre.properties.int("nprocy",1),
            pyre.properties.int("nprocz",1),
            # debugging flags?

            ]


# version
__id__ = "$Id: Parallel.py,v 1.3 2003/07/23 05:29:58 ces74 Exp $"

# End of file
