#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from pyre.components.Component import Component


class Exchanger(Component):


    def __init__(self, name, facility):
        Component.__init__(self, name, facility)

        self.module = None
        self.exchanger = None
        return



    def selectModule(self):
        import CitcomS.Exchanger
        self.module = CitcomS.Exchanger
        return



    def createExchanger(self, solver):
        raise NotImplementedError
        return



    def findBoundary(self):
        raise NotImplementedError
        return None



    def mapBoundary(self, boundary):
        raise NotImplementedError
        return



    def initTemperature(self):
        raise NotImplementedError
        return



    def NewStep(self):
        raise NotImplementedError
        return



    def applyBoundaryConditions(self):
        raise NotImplementedError
        return



    def stableTimestep(self, dt):
        raise NotImplementedError
        return



    def endTimestep(self):
        #raise NotImplementedError
        return



    class Inventory(Component.Inventory):

        import pyre.properties as prop


        inventory = [

            ]



# version
__id__ = "$Id: Exchanger.py,v 1.6 2003/09/10 04:01:53 tan2 Exp $"

# End of file
