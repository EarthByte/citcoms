#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from pyre.components.SimulationController import SimulationController
import journal

class Controller(SimulationController):


    def __init__(self, name, facility="controller"):
        SimulationController.__init__(self, name, facility)

        self.step = 0
        self.clock = 0.0
        self.done = False
        self.solver = None
        return



    def initialize(self, app):
        self.solver = app.solver
        self.solver.initialize(app)
        return



    def launch(self, app):
        self.solver.launch(app)
        return



    def march(self, totalTime=0, steps=0):
        """explicit time loop"""

        # do io for 0th step
        self.save()

        while 1:

            # notify solvers we are starting a new timestep
            self.startTimestep()

            # synchronize boundary information
            self.applyBoundaryConditions()

            # compute an acceptable timestep
            dt = self.stableTimestep()

            # advance
            self.advance(dt)

            # update smulation clock and step number
            self.clock += dt
            self.step += 1

            # do io
            self.save()

            # notify solver we finished a timestep
            self.endTimestep(totalTime, steps)

            # are we done?
            if self.done:
                break

        # end of time advance loop

        # Notify solver we are done
        self.endSimulation()

        return



    def endTimestep(self, totalTime, steps):
        # are we done?
        if steps and self.step >= steps:
            self.done = True
        if totalTime and self.clock >= totalTime:
            self.done = True

        # solver can terminate time marching by returning True
        self.done = self.solver.endTimestep(self.clock, self.step, self.done)

        return



    def endSimulation(self):
        self.solver.endSimulation(self.step)
        return



    def save(self):
        step = self.step
        if not step % self.inventory.monitoringFrequency:
            self.solver.save(step)
        return


