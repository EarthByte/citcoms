#!/usr/bin/env python
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# <LicenseText>
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#

from Stokes_solver import Stokes_solver
import CitcomS.Regional as Regional


class ImcompressibleNewtonian(Stokes_solver):


    def _form_RHS(self):
	Regional.velocities_conform_bcs()
	Regional.assemble_forces()
	return


    def _form_LHS(self):
	Regional.get_system_viscosity()
	Regional.construct_stiffness_B_matrix()
	return


    def _solve(self):
	return Regional.solve_constrained_flow_iterative()




class ImcompressibleNonNewtonian(ImcompressibleNewtonian):

    # over-ride Stokes_solver.run()
    def run(self):

	self._init()

	self._form_RHS()

	while (self.count < 50) and self.sdepv_not_convergent:

	    self._form_LHS()
	    self.viscosity_misfit = self._solve()

	    Regional.general_stokes_solver_update_velo()
	    self.Udot_mag = Regional.general_stokes_solver_Unorm()
	    self.dUdot_mag = Regional.general_stokes_solver_Udotnorm()

	    Regional.general_stokes_solver_log(self.Udot_mag, self.dUdot_mag,
					       self.count)

	    self.sdepv_not_convergent = (self.dUdot_mag > self.viscosity_misfit)
	    self.count += 1

	self._fini()

	return


    def _init(self):
	# allocate and initialize memory here
	Regional.general_stokes_solver_init()

	self.Udot_mag = 0
	self.dUdot_mag = 0
	self.count = 1
	self.sdepv_not_convergent = True
	return


    def _fini(self):
	# free memory here
	Regional.general_stokes_solver_fini()
	return








# version
__id__ = "$Id: Incompressible.py,v 1.5 2003/06/26 23:14:42 tan2 Exp $"

# End of file
