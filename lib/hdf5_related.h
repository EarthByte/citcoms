/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/* In this file we define the following data structures:
 *
 *  HDF5_INFO
 *    Used for collective output of citcom data.
 *
 *  HDF5_TIME
 *    Used to define table with timing information.
 *
 *  field_t
 *    Used to store dataspace and hyperslab parameters.
 *
 *  field_class_t
 *    Used to deduce the number of components in each point of a field.
 *
 * Any required initialization steps are performed in h5output_open().
 *
 */

enum field_class_t
{
    SCALAR_FIELD = 0,
    VECTOR_FIELD = 1,
    TENSOR_FIELD = 2
};

typedef struct field_t
{
    /* field datatype (in file) */
    hid_t dtype;

    /* field dataspace (in file) */
    int rank;
    hsize_t *dims;
    hsize_t *maxdims;
    hsize_t *chunkdims;
    
    /* hyperslab selection parameters */
    hsize_t *offset;
    hsize_t *stride;
    hsize_t *count;
    hsize_t *block;

    /* number of data points in buffer */
    int n;
    float *data;

} field_t;

struct HDF5_TIME
{
    int step;
    float time;
    float time_step;
    float cpu;
    float cpu_step;
};

struct HDF5_INFO
{
    /* Keep a reference to the open hdf5 output file */
    char filename[100];
    hid_t file_id;

    /* Keep a reference to the MPI_Info object */
    MPI_Info mpi_info;

    /* Cap ID for current process */
    int cap;

    /* Keep track of how many times we call h5output() */
    int count;

    /* Data structures to use in dataset writes...
     * 
     * const_vector3d: coord
     * const_vector2d: surf_coord
     *
     * tensor3d: stress
     * vector3d: velocity
     * vector2d: surf_velocity
     *
     * scalar3d: temperature, viscosity, pressure
     * scalar2d: surf_heatflux, surf_topography
     * scalar1d: horiz_avg_temperature, horiz_rms_vz, horiz_rms_vxy
     *
     */
    
    field_t *const_vector3d;    /* shape (xdim,ydim,zdim,3) */
    field_t *const_vector2d;    /* shape (xdim,ydim,2) */
    field_t *const_scalar1d;    /* shape (zdim,) */
    
    field_t *tensor3d;          /* shape (tdim,xdim,ydim,zdim,6) */
    field_t *vector3d;          /* shape (tdim,xdim,ydim,zdim,3) */
    field_t *vector2d;          /* shape (tdim,xdim,ydim,2) */

    field_t *scalar3d;          /* shape (tdim,xdim,ydim,zdim) */
    field_t *scalar2d;          /* shape (tdim,xdim,ydim) */
    field_t *scalar1d;          /* shape (tdim,zdim) */

    /* For list of extendible (time dependent) fields */
    field_t *field[6];

    /* Actual data buffer -- shared over all fields! */
    float *data;
};
