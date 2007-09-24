/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *<LicenseText>
 *
 * CitcomS by Louis Moresi, Shijie Zhong, Lijie Han, Eh Tan,
 * Clint Conrad, Michael Gurnis, and Eun-seo Choi.
 * Copyright (C) 1994-2005, California Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *</LicenseText>
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <math.h>

#include "global_defs.h"
#include "material_properties.h"


void mat_prop_allocate(struct All_variables *E)
{
    int noz = E->lmesh.noz;
    int nno = E->lmesh.nno;
    int nel = E->lmesh.nel;

    /* reference profile of density */
    E->refstate.rho = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of coefficient of thermal expansion */
    E->refstate.thermal_expansivity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of heat capacity */
    E->refstate.heat_capacity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of thermal conductivity */
    /*E->refstate.thermal_conductivity = (double *) malloc((noz+1)*sizeof(double));*/

    /* reference profile of gravity */
    E->refstate.gravity = (double *) malloc((noz+1)*sizeof(double));

    /* reference profile of temperature */
    /*E->refstate.Tadi = (double *) malloc((noz+1)*sizeof(double));*/

}


void reference_state(struct All_variables *E)
{
    int noz = E->lmesh.noz;
    int nel = E->lmesh.nel;
    int i;
    double r, z, beta;

    beta = E->control.disptn_number * E->control.inv_gruneisen;

    /* All refstate variables (except Tadi) must be 1 at the surface.
     * Otherwise, the scaling of eqns in the code might not be correct. */

    /* Adams-Williamson EoS */
    for(i=1; i<=noz; i++) {
	r = E->sx[1][3][i];
	z = 1 - r;
	E->refstate.rho[i] = exp(beta*z);
	E->refstate.thermal_expansivity[i] = 1;
	E->refstate.heat_capacity[i] = 1;
	/*E->refstate.thermal_conductivity[i] = 1;*/
	E->refstate.gravity[i] = 1;
	/*E->refstate.Tadi[i] = (E->control.adiabaticT0 + E->control.surface_temp) * exp(E->control.disptn_number * z) - E->control.surface_temp;*/
    }

    if(E->parallel.me == 0) {
        fprintf(stderr, "nz  radius   depth    rho\n");
    }
    if(E->parallel.me < E->parallel.nprocz)
        for(i=1; i<=noz; i++) {
            fprintf(stderr, "%d %f %f %e\n",
                    i+E->lmesh.nzs-1, E->sx[1][3][i], 1-E->sx[1][3][i],
                    E->refstate.rho[i]);
        }

}


