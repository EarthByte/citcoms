#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from mpi.Application import Application
import CitcomS.Regional as Regional
import journal


class RegionalApp(Application):


    # for test
    def run(self):
	journal.info("staging").log("setup MPI")
        import mpi
        Regional.Citcom_Init(mpi.world().handle())

	self.rank = mpi.world().rank
	self.start_time = Regional.CPU_time()
	print "my rank is ", self.rank

        self.setProperties()

        mesher = self.inventory.mesher
        mesher.init(self)

        vsolver = self.inventory.vsolver
        vsolver.init(self)

        tsolver = self.inventory.tsolver
        tsolver.init(self)

        mesher.run()

        return


    def __init__(self, inputfile):
        Application.__init__(self, "citcomsregional")
	self.filename = inputfile
        self.total_time = 0
        self.cycles = 0
        self.keep_going = True
        self.Emergency_stop = False
        self.start_time = 0.0

	#test
	self.prefix = 'test'

        return

    def init(self):
        return

    #def fini(self):
	#self.total_time = Regional.CPU_time() - self.start_time
	#Regional.finalize()
	#Application.fini()
	#return


    def output(self):
        import CitcomS.Output as Output
        output_coord = Output.outputCoord(self.prefix, self.rank)
	output_coord.go()

	output_velo = Output.outputVelo(self.prefix, self.rank, self.cycles)
        output_velo.go()

	output_visc = Output.outputVisc(self.prefix, self.rank, self.cycles)
        output_visc.go()

        return



    def setProperties(self):
        self.inventory.bc.setProperties()
        self.inventory.visc.setProperties()

        return



    class Inventory(Application.Inventory):

        import pyre.facilities

        # facilities
        from CitcomS.Facilities.Mesher import Mesher
        from CitcomS.Facilities.TSolver import TSolver
	from CitcomS.Facilities.VSolver import VSolver

        # component modules
        import CitcomS.Advection_diffusion
        import CitcomS.Sphere
	import CitcomS.Stokes_solver

        # components
        from CitcomS.Components.BC import BC
        from CitcomS.Components.Const import Const
        from CitcomS.Components.IC import IC
        #from CitcomS.Components.Mesher import Mesher
	from CitcomS.Components.Parallel import Parallel
	from CitcomS.Components.Param import Param
        from CitcomS.Components.Phase import Phase
        from CitcomS.Components.Visc import Visc

        inventory = [
            Mesher("mesher", CitcomS.Sphere.regionalSphere()),
            VSolver("vsolver", CitcomS.Stokes_solver.imcompressibleNewtonian()),
            TSolver("tsolver", CitcomS.Advection_diffusion.temperature_diffadv()),

            pyre.facilities.facility("bc", default=BC()),
            pyre.facilities.facility("const", default=Const()),
            pyre.facilities.facility("ic", default=IC()),
            #pyre.facilities.facility("mesher", default=Mesher(1)),
	    pyre.facilities.facility("parallel", default=Parallel()),
            pyre.facilities.facility("param", default=Param()),
            pyre.facilities.facility("phase", default=Phase()),
            pyre.facilities.facility("visc", default=Visc()),

            ]


    # this is the true run(), but old
    def run_old(self):
	Application.run(self)

        # read in parameters
        Regional.read_instructions(self.filename)

	#if (Control.post_proccessing):
	#    Regional.post_processing()
	#    return

	# decide which stokes solver to use
	import CitcomS.Stokes_solver
	vsolver = CitcomS.Stokes_solver.imcompressibleNewtionian('imcompressible')
	vsolver.init()

	# decide which field to advect (and to diffuse)
	import CitcomS.Advection_diffusion as Advection_diffusion
	tsolver = Advection_diffusion.temperature_diffadv('temp')
	tsolver.init()


	# solve for 0th time step velocity and pressure
	vsolver.run()

	# output phase
        self.output()

	while (self.keep_going and not self.Emergency_stop):
	    self.cycles += 1
	    print 'cycles = ', self.cycles
	    #tsolver.run()
	    #vsolver.run()
	    total_time = Regional.CPU_time() - self.start_time

            #self.output()
	return


# version
__id__ = "$Id: RegionalApp.py,v 1.16 2003/07/15 21:50:30 tan2 Exp $"

# End of file
