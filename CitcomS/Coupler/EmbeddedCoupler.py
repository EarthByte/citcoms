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

from Coupler import Coupler

class EmbeddedCoupler(Coupler):


    def __init__(self, name, facility):
        Coupler.__init__(self, name, facility)
        self.cge_t = 0
        self.fge_t = 0
        self.toApplyBC = True

        # exchanged information is non-dimensional
        self.inventory.dimensional = False
        # exchanged information is in spherical coordinate
        self.inventory.transformational = False
        return


    def initialize(self, solver):
        print self.name, 'entering initialize'
        Coupler.initialize(self, solver)

	# restart and use temperautre field of previous run?
        self.restart = solver.restart
        if self.restart:
            self.ic_initTemperature = solver.ic_initTemperature

	self.all_variables = solver.all_variables
        self.interior = range(self.numSrc)
        self.source["Intr"] = range(self.numSrc)
        self.II = range(self.numSrc)

        # the embedded solver should set its solver.bc.side_sbcs to on
        # otherwise, we have to stop
        if not solver.inventory.bc.inventory.side_sbcs:
            raise SystemExit('\n\nError: esolver.bc.side_sbcs must be on!\n\n\n')

        from ExchangerLib import initConvertor
        initConvertor(self.inventory.dimensional,
                      self.inventory.transformational,
                      self.all_variables)

        print self.name, 'leaving initialize'
        return


    def createMesh(self):
        from ExchangerLib import createGlobalBoundedBox, exchangeBoundedBox, createBoundary, createEmptyInterior
        inv = self.inventory
        self.globalBBox = createGlobalBoundedBox(self.all_variables)
        mycomm = self.communicator
        self.remoteBBox = exchangeBoundedBox(self.globalBBox,
                                             mycomm.handle(),
                                             self.sinkComm.handle(),
                                             0)
        self.boundary, self.myBBox = createBoundary(self.all_variables,
                                                    inv.excludeTop,
                                                    inv.excludeBottom)

        if inv.two_way_communication:
            for i in range(len(self.interior)):
                self.interior[i] = createEmptyInterior()

        return


    def createSourceSink(self):
        self.createSink()

        if self.inventory.two_way_communication:
            self.createSource()
        return


    def createSink(self):
        from ExchangerLib import Sink_create
        self.sink["BC"] = Sink_create(self.sinkComm.handle(),
                                      self.numSrc,
                                      self.boundary)
        return


    def createSource(self):
        from ExchangerLib import CitcomSource_create
        for i, comm, b in zip(range(self.numSrc),
                              self.srcComm,
                              self.interior):
            # sink is always in the last rank of a communicator
            sinkRank = comm.size - 1
            self.source["Intr"][i] = CitcomSource_create(comm.handle(),
                                                         sinkRank,
                                                         b,
                                                         self.myBBox,
                                                         self.all_variables)

        return


    def createBC(self):
        import Inlet
        self.BC = Inlet.SVTInlet(self.boundary,
                                 self.sink["BC"],
                                 self.all_variables)
        '''
        if self.inventory.incompressibility:
            self.BC = Inlet.BoundaryVTInlet(self.communicator,
                                            self.boundary,
                                            self.sink["BC"],
                                            self.all_variables,
                                            "VT")
            import journal
            journal.info("incompressibility").activate()
        else:
            self.BC = Inlet.SVTInlet(self.boundary,
                                    self.sink["BC"],
                                    self.all_variables)
        '''
        return


    def createII(self):
        import Outlet
        for i, src in zip(range(self.numSrc),
                          self.source["Intr"]):
            self.II[i] = Outlet.TOutlet(src,
                                        self.all_variables)
        return


    def initTemperature(self):
        if self.restart:
            # receive temperature from CGE and postprocess
            self.restartTemperature()
        else:
            from ExchangerLib import initTemperature
            initTemperature(self.globalBBox,
                            self.all_variables)
        return


    def restartTemperature(self):
        from ExchangerLib import createInterior, Sink_create
        interior, bbox = createInterior(self.remoteBBox,
                                        self.all_variables)
        sink = Sink_create(self.sinkComm.handle(),
                           self.numSrc,
                           interior)
        import Inlet
        inlet = Inlet.TInlet(interior, sink, self.all_variables)
        inlet.recv()
        inlet.impose()

        # Any modification of read-in temperature is done here
        # Note: modifyT is called after receiving unmodified T from CGE.
        # If T is modified before sending, FGE's T will lose sharp feature.
        # CGE has to call modifyT too to ensure consistent T field.
        self.modifyT(self.globalBBox)

        return


    def preVSolverRun(self):
        self.applyBoundaryConditions()
        return


    def newStep(self):
        if self.inventory.two_way_communication:
            if self.catchup:
                # send temperture field to CGE
                for ii in self.II:
                    ii.send()

        return


    def applyBoundaryConditions(self):
        if self.toApplyBC:
            self.BC.recv()

            self.toApplyBC = False

        self.BC.impose()

        # applyBC only when previous step is a catchup step
        if self.catchup:
            self.toApplyBC = True

        return


    def stableTimestep(self, dt):
        from ExchangerLib import exchangeTimestep
        if self.catchup:
            mycomm = self.communicator
            self.cge_t = exchangeTimestep(dt,
                                          mycomm.handle(),
                                          self.sinkComm.handle(),
                                          0)
            self.fge_t = 0
            self.catchup = False

        self.fge_t += dt
        old_dt = dt

        if self.fge_t >= self.cge_t:
            dt = dt - (self.fge_t - self.cge_t)
            self.fge_t = self.cge_t
            self.catchup = True
            #print "FGE: CATCHUP!"

        # store timestep for interpolating boundary velocities
        self.BC.storeTimestep(self.fge_t, self.cge_t)

        #print "%s - old dt = %g   exchanged dt = %g" % (
        #       self.__class__, old_dt, dt)
        #print "cge_t = %g  fge_t = %g" % (self.cge_t, self.fge_t)
        return dt


    def exchangeSignal(self, signal):
        from ExchangerLib import exchangeSignal
        mycomm = self.communicator
        newsgnl = exchangeSignal(signal,
                                 mycomm.handle(),
                                 self.sinkComm.handle(),
                                 0)
        return newsgnl



    class Inventory(Coupler.Inventory):

        import pyre.inventory as prop



        excludeTop = prop.bool("excludeTop", default=False)
        excludeBottom = prop.bool("excludeBottom", default=False)
        incompressibility = prop.bool("incompressibility", default=True)




# version
__id__ = "$Id$"

# End of file
