#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from Exchanger import Exchanger

class CoarseGridExchanger(Exchanger):


    def createExchanger(self, solver):
        self.exchanger = self.module.createCoarseGridExchanger(
                                     solver.communicator.handle(),
                                     solver.intercomm.handle(),
                                     solver.localLeader,
                                     solver.remoteLeader,
                                     solver.coarse_E
                                     )
        return


    def createDataArrays(self):
        self.module.createDataArrays(self.exchanger)
        return

    def deleteDataArrays(self):
        self.module.deleteDataArrays(self.exchanger)
        return    

    def findBoundary(self):
        # receive boundary from FGE
        self.module.receiveBoundary(self.exchanger)

        # create mapping from boundary to id array
        self.module.mapBoundary(self.exchanger)
        return

    def getBid2crseelem(self):
        self.module.getBid2crseelem(self.exchanger)
        return

    def gather(self):
        self.module.gather(self.exchanger)
        return

    def distribute(self):
        self.module.distribute(self.exchanger)
        return

    def interpolate(self):
        self.module.interpolate(self.exchanger)
        return

    def initTemperature(self):
        # send temperture field to FGE
        self.module.sendTemperature(self.exchanger)
        return

    def exchangeVelocities(self):
        self.module.sendVelocities(self.exchanger)
        return

    def NewStep(self):
        # wait until FGE catchs up
        self.module.wait(self.exchanger)

        # receive temperture field from FGE
        #self.module.receiveTemperature(self.exchanger)
        return



    def applyBoundaryConditions(self):
        self.module.receive(self.exchanger)
        self.module.distribute(self.exchanger)
        return



    def stableTimestep(self, dt):
        self.module.exchangeTimestep(self.exchanger, dt)
        return dt





    class Inventory(Exchanger.Inventory):

        import pyre.properties as prop


        inventory = [

            ]



# version
__id__ = "$Id: CoarseGridExchanger.py,v 1.10 2003/09/21 22:24:00 ces74 Exp $"

# End of file
