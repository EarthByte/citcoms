#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from pyre.components.Component import Component


class Coupler(Component):


    def __init__(self, name, facility):
        Component.__init__(self, name, facility)

        self.exchanger = None
        self.boundary = None
        return



    def initialize(self, solver):
        # exchanger could be either a FineGridExchanger (FGE)
        # or a CoarseGridExchanger (CGE)
        self.exchanger = solver.inventory.exchanger

        # choose c++ exchanger module
        self.exchanger.selectModule()
        # create c++ exchanger
        self.exchanger.createExchanger(solver)
        return



    def launch(self, solver):
        exchanger = self.exchanger

        # find the common boundary
        print 'exchanging boundary'
        self.boundary = exchanger.findBoundary()

        # create mapping from boundary -> ID and rank
        exchanger.mapBounary(self.boundary)

        # send initial temperature field from CGE to FGE
        exchanger.initTemperature()
        return



    def newStep(self):
        self.exchanger.NewStep()
        return



    def applyBoundaryConditions(self):
        self.exchanger.appliBoundaryConditions()
        return



    def stableTimestep(self, dt):
        dt = self.exchanger.exchangeTimestep(dt)
        return dt



    def endTimestep(self):
        self.exchanger.endTimestep(self)
        return


# version
__id__ = "$Id: Coupler.py,v 1.2 2003/09/09 21:04:45 tan2 Exp $"

# End of file
