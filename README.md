# CitcomS with Data Assimilation

## Installation

This version of CitcomS only requires a valid MPI distribution (no more python / pyre).  Most clusters can load an MPI distribute using modules

e.g., in Bern I use:

```module load openmpi/1.10.2-intel```

For reasons relating to compilers, the cluster and node setup, some mpi versions may be preferred for your particular cluster.  I was recommended by the Bern sys admin to use the intel compiler with openmpi.  You might also want to find out which openmpi is preferred on your cluster (and perhaps add this above, so we have complete documentation of what works and what doesn't).

Then, to install CitcomS v3.1.1 with data assimilation go into the src/ directory and execute

```./mymake.py```

[Yes OK, for the step above you do need a python distribution, but this is only to run the commands that configure and make Citcoms.  Python is not actually used for running the code]

Example input and output configuration files are provided in:

```inputeg/```

This example is a good one to try and run first:

```src/examples/Full```

You will need to setup your job submission script.  See a slurm example in jobsubmit/  Please add your own submission scripts to the same directory so we have more examples

For currently implemented features, scroll down below.

## Dan's work area follows

### Rakib's code

Sent to me by Sabin (12/09/17) CitcomSModDtopoRelease\_Rakib.zip
See diff/ directory for complete record

1. Advection\_diffusion
    - Scaled visc and adiabatic heating by Di (COMPLETE)
1. convection\_variables.h
    - blob\_profile
    - silo parameters
    - mantle\_temp\_adiabatic\_increase
1. global\_defs.h
    - Shell-output facility
1. Initial\_temperature.c
    - bunch of silo / blob related functions  
    - evidentally seeding plumes for the IC
1. Material\_properties.c
    - refstate updates (COMPLETE)
1. Output.c
    - outputs of shells (theta, phi, r, temperature, vr)
1. Output\_vtk.c
    - perhaps not related to Rakib's work
1. Viscosity\_structures.c
    - some new viscosity structures (case 112, 113, 117, 118) (COMPLETE)

### Code features implemented

#### Slab and lithosphere assimilation (grep for 'DJB SLAB')
1. ```lith_age_depth_function``` (bool)
1. ```lith_age_exponent``` (double)
1. ```lith_age_min``` (double)
1. ```lith_age_stencil_value``` (double)
1. ```slab_assim``` (bool)
1. ```slab_assim_file``` (char)
1. ```sten_temp``` as an output option (char)
    
#### Composition (grep for 'DJB COMP')
1. ```hybrid_method``` (bool)
1. increase memory for tracer arrays (icushion parameter)
1. user note: you should turn off tracer warnings using:
    ```itracer_warnings=off```
   in input cfg file

#### Viscosity structures (grep for 'DJB VISC')
1. case 20, used in Flament et al. (2013, 2014)
1. case 21, used in Flament et al. (2014), model TC8
1. case 22, used in Flament et al. (2014), model TC7
1. case 23, used in Flament et al. (2014), model TC9
1. case 24, used in Zhang et al. (2010) and Bower et al. (2013)
1. case 25, used by Flament for Extendend-Boussinesq (EBA) models
1. case 26, used by Flament for EBA models
1. case 27, used by Flament for EBA models
1. case 28, used by Flament for EBA models
1. case 29, used by Flament for EBA models *(see issue)*
1. case 112, used by Hassan presumably for plume models(?)
1. case 113, used by Hassan presumably for plume models(?)
1. case 117, used by Hassan presumably for plume models(?)
1. case 118, used by Hassan presumably for plume models(?)

#### Output time (grep for 'DJB TIME')
1. output time by Myr as well as/rather than number of time steps
    - ```record_every_Myr``` (int)
    - modifications made in main time loop in bin/Citcom.c
    - needs testing
1. exit time loop for negative ages (\<1 Ma)
    - ```exit_at_present``` (bool)
    - modifications made in main time loop in bin/Citcom.c
    - needs testing
    
#### Extended-Boussinesq modifications (grep for 'DJB EBA')
1. depth-dependent scaling for the dissipation number
   this effectively scales the adiabatic, viscous, and latent heating
   and can be useful to avoid large heating in certain radial parts
   of the domain, notably the surface when velocity bcs are imposed

### Code features NOT implemented

1. internal velocity bcs (I don't think these are used by anyone anyway)
1. outputs of heating terms, divv
1. tracer density for elements and nodes output (added by Ting, see svn r52 through r55)
1. buoyancy restart for dynamic topography (exclude buoyancy) (see svn r85)
1. composition and temperature spherical harmonics
1. reverse gravity acceleraton (added by Ting, see svn r76) for SBI
   note that this appears to have been subsequently removed in r88?

### Log of features implemented from legacy code
1. Adv\_diff -> COMPLETE
1. BC\_util -> ivels (TODO)
1. Composition\_related.c -> tracer density (TODO)
1. composition\_related.h -> tracer density (TODO)
1. Convection.c -> COMPLETE
1. Element\_calculations -> buoy restart for dyn topo (TODO)
1. Full\_boundary\_conditions -> ivels (TODO)
1. Full\_lith\_age\_read\_files -> COMPLETE
1. Full_read_input_from_file -> ivels (TODO)
1. Full\_solver.c -> COMPLETE
1. global\_defs.h -> dyn topo restarts, ivels (TODO)
1. Instructions.c -> outputs, heating terms, divv (TODO)
1. Lith\_age.c -> COMPLETE
1. output.c -> various outputs (TODO)
1. Pan\_problem.c -> dyn topo (TODO)
1. Problem\_related.c -> ivels (TODO)
1. Regional\_bcs.c -> ivels (TODO)
1. Regional\_lith\_age\_read\_files -> COMPLETE
1. Regional\_read\_input\_files -> ivels (TODO)
1. Regional\_solver -> COMPLETE
1. solver.h -> COMPLETE
1. Stokes\_flow\_incom -> divv calculation for output (TODO)
1. tracer\_defs.h -> output tracer density (TODO)
1. Viscosity\_structure.c -> COMPLETE
