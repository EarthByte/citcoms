#!/usr/bin/env python
#=====================================================================
#                Geodynamic Framework Scripts for 
#         Preprocessing, Data Assimilation, and Postprocessing
#
#                 AUTHORS: Mark Turner, Dan J. Bower
#                  ---------------------------------
#             (c) California Institute of Technology 2014
#                  ---------------------------------
#                        ALL RIGHTS RESERVED
#=====================================================================
#=====================================================================
# grid_maker.py
#=====================================================================
# This script is a general purpose tool to process Citcoms output data
# into one ore more GMT style .grd or .nc format data files.  
# Please see the usage() function below, and the sample configuration 
# file: /sample_data/grid_maker.cfg for more info.
#=====================================================================
#=====================================================================
import sys, string, os, shutil
from multiprocessing import Pool, Process
import numpy as np
#=====================================================================
import Core_Citcom
import Core_isoGMT
import Core_Util
from Core_Util import now

#=====================================================================
#=====================================================================
def usage():
    '''print usage message, and exit'''

    print('''usage: grid_maker.py [-e] configuration_file.cfg

Options and arguments:
  
-e : if the optional -e argument is given this script will print to standard out an example configuration control file.
   The parameter values in the example config.cfg file may need to be edited or commented out depending on intended use.

'configuration_file.cfg' : is a geodynamic framework formatted control file, with at least these entries: 

    pid_file = /path/to/a/citcoms_pid_file # the path of a citcoms pid0000.cfg file 
    time_spec = multi-value specification (single value, comma delimted list, start/stop/step trio),
    level_spec = multi-value specification (single value, comma delimted list, start/stop/step trio),

and at least one sub-section:

    [Subsection], where 'Subsection' may be any string, followed by:
    field = standard Citcom field name ('temp', 'visc', 'comp', etc. - see Core_Citcom.py for more info)
      
and where each sub-section may have one or more of the following optional entries:

     dimensional = True to also generate a dimensionalized grid 
     blockmedian_I = value to pass to GMT blockmedian -I option
     surface_I = value to pass to GMT surface -I option

See the example config.cfg file for more info.
''')
    sys.exit()
#=====================================================================
#=====================================================================
def initialise_variables(configFile=sys.argv[1],verbose=False):
    # get the .cfg file as a dictionary
    control_d = Core_Util.parse_configuration_file( configFile, False, False )
    
    print(f"{now()} Config file dictionary:")
    if verbose: Core_Util.tree_print( control_d )

    # set the pid file 
    pid_file = control_d['pid_file']
    
    # get the master dictionary and define aliases
    master_d = Core_Citcom.get_all_pid_data( pid_file, verbose=verbose )
    coor_d = master_d['coor_d']
    pid_d = master_d['pid_d']
    
    # set up working variables
    # get basic info about the model run
    datadir       = pid_d['datadir']
    datafile      = pid_d['datafile']
    start_age     = pid_d['start_age']
    output_format = pid_d['output_format']

    depth_list = coor_d['depth_km']
    nodez      = pid_d['nodez']
    nproc_surf = pid_d['nproc_surf']

    return control_d, pid_file, master_d, coor_d, pid_d, datadir, datafile, start_age, output_format, depth_list, nodez, nproc_surf

def main():
    print( now(), 'grid_maker.py')

    control_d, pid_file, master_d, coor_d, pid_d, datadir, datafile, start_age, output_format, depth_list, nodez, nproc_surf = initialise_variables()
    
    # print(master_d['time_d']['triples'][2])
    # for t in master_d['time_d']['triples']:
    #     print(int(t[1]))
    # sys.exit()

    # Set the verbose settings
    if 'verbose' in control_d:
        print(f"{now()} Verbose is set to {control_d['verbose']} in the config file")
        verbose=control_d['verbose']
    else:
        print(f"{now()} Verbose is not set in config file, using defaults (probably verbose=False)")
        verbose=False

    if len(sys.argv) > 2: # Prescribe procs from 2nd argument in command line. note, if removing this, also remove the 'if name==__main__' thing back to 2 args
        procs = int(sys.argv[2])
        print(f'Number of processes prescribed from command line: {procs}')
    elif 'procs' in control_d:
        procs = control_d['procs']
    else:
        print(f"{now()} Number of processors for parallel computing is not set. Will use 2 processors for now.)")
        procs = 1

    # Adjust verbose settings for all other functions- can be manually adjusted if you want to mix and match 
    Core_Util.verbose = verbose 
    Core_Citcom.verbose = verbose 
    Core_isoGMT.verbose = verbose

    # Set the debug settings
    if 'debug' in control_d:
        print(f"{now()} Debug is set to {control_d['debug']} in the config file")
        debug=control_d['debug']
    else:
        print(f"{now()} Debug is not set in config file, using defaults (probably debug=False)")
        debug=False

    # Double check for essential data 
    if master_d['time_d'] == None : 
        print( now() )
        print('ERROR: Required file "[CASE_NAME].time:" is missing from this model run.')
        print('       Aborting processing.')
        sys.exit(-1)

    # Check how to read and parse the time spec:
    read_time_d = True
    #read_time_d = False
    
    # Compute the timesteps to process
    if read_time_d : 
        time_spec_d = Core_Citcom.get_time_spec_dictionary(control_d['time_spec'], master_d['time_d'])
        
    else :
        time_spec_d = Core_Citcom.get_time_spec_dictionary(control_d['time_spec'])
    print ( now(), 'grid_maker.py: time_spec_d = ')
    if verbose: Core_Util.tree_print( time_spec_d )

    # levels to process 
    level_spec_d = Core_Util.get_spec_dictionary( control_d['level_spec'] )
    print ( now(), 'grid_maker.py: level_spec_d = ')
    if verbose: Core_Util.tree_print( level_spec_d )

    # Get coordinate data 
    lon = []
    lat = []
    
    # Check for existing coordinate data
    lon_file_cache = '_cache_lon_coords.txt'
    lat_file_cache = '_cache_lat_coords.txt'

    if os.path.exists( lon_file_cache ) and os.path.exists( lat_file_cache ) :
        print( now(), 'grid_maker.py: loadtxt: ', lon_file_cache )
        print( now(), 'grid_maker.py: loadtxt: ', lat_file_cache )
        lon = np.loadtxt( lon_file_cache )
        lat = np.loadtxt( lat_file_cache )
    else : 
        # gets lon, lat for one depth because these are depth-independent
        coord_file_format = control_d.get('coord_dir','') + '/%(datafile)s.coord.#' % vars()
        coord = Core_Citcom.read_citcom_surface_coor( master_d['pid_d'], coord_file_format )

        # flatten data since we don't care about specific cap numbers for the loop over depth
        coord = Core_Util.flatten_nested_structure( coord )

        # extract data from tuples and make into numpy array
        lon = [line[0] for line in coord]
        lat = [line[1] for line in coord]

        # save the lat data 
        np.savetxt( lon_file_cache, lon, fmt='%f' )
        np.savetxt( lat_file_cache, lat, fmt='%f' )

    # end of get coords
    print( now(), 'grid_maker.py: len(lon) = ', len(lon) )
    print( now(), 'grid_maker.py: len(lat) = ', len(lat) )

    #
    # Main looping, first over times, then sections, then levels
    # 

    # Variables that will be updated each loop:
    # age_Ma will be a zero padded string value used for filenames and reporting 
    # depth will be a zero padded string value used for filenames and reporting 

    # Variables to hold data for all grids created
    # grid_list is a list of tuples: (grid_filename, age_Ma) 
    grid_list = []

    print(now(), '=========================================================================')
    print(now(), 'grid_maker.py: Main looping, first over times, then sections, then levels')
    print(now(), '=========================================================================')
    
    # Loop over times
    # Only append ages to var list that match ages in the citcoms data (avoids double ups)
    # Note: this can be cleaned up a little - it's essentially been copied from below
    varList = []
    for tt, time in enumerate( time_spec_d['time_list'] ):
        print(f'tt:{tt}')
        print(time)

        if 'Ma' in time:
        
            # strip off units and make a number
            time = float( time.replace('Ma', '') )

            # determine what time steps are available for this age 
            # NOTE: 'temp' is requried to set which output files to check 
            found_d = Core_Citcom.find_available_timestep_from_age( master_d, 'temp', time )

            # print( now(), 'grid_maker.py: WARNING: Adjusting times to match available data:')
            # print( now(), '  request_age =', found_d['request_age'], '; request_timestep =', found_d['request_timestep'], '; request_runtime =', found_d['request_runtime'])
            # print( now(), '  found_age =', found_d['found_age'], '; found_timestep =', found_d['found_timestep'], '; found_runtime =', found_d['found_runtime'])

            # set variables for subsequent loops
            timestep = found_d['found_timestep']
            runtime_Myr = found_d['found_runtime']

            # convert the found age to an int
            age_Ma = int(np.around( found_d['found_age'] ) )

            # save age number for grids storage
            age_Ma_storing=age_Ma
             
            age_Ma = '%03d' % age_Ma

        else:

            time = float( time ) 
             
            # determine what time steps are available for this timestep 
            # NOTE: 'temp' is requried to set which output files to check 

            found_d = Core_Citcom.find_available_timestep_from_timestep( master_d, 'temp', time )

            print( now(), 'grid_maker.py: WARNING: Adjusting times to match available data:')
            print( now(), '  request_age =', found_d['request_age'], '; request_timestep =', found_d['request_timestep'], '; request_runtime =', found_d['request_runtime'])
            print( now(), '  found_age =', found_d['found_age'], '; found_timestep =', found_d['found_timestep'], '; found_runtime =', found_d['found_runtime'])

            # set variables for subsequent loops
            timestep = found_d['found_timestep']
            runtime_Myr = found_d['found_runtime']

            # convert the found age to an int
            age_Ma = int(np.around( found_d['found_age'] ) )
            
            # make a string and pad with zeros 
            #age_Ma = '%03d' % age_Ma
            age_Ma = str(age_Ma)

        print(found_d)


        print(age_Ma)
        print('\n')

        if int(time) == int(age_Ma): # Think about changing this so it chooses the nearest age, but still doesn't double up
            varList.append((tt,f'{time}Ma',master_d, control_d, nproc_surf, datadir, datafile, level_spec_d, lon, lat, grid_list, depth_list, nodez, debug, pid_file, verbose))
        else:
            continue

    # Multiprocessing of times with old-style batching - specify number of processes with 'procs' in cfg or as an argument on command line
    batches = divmod(len(varList), procs) # Return number of batches and the remainder
    for i in range(batches[0]+1): # Loop over the number of batches
        if i == batches[0]: # If we're at the last (or only) batch 
            processes = [None] * batches[1] # Batch size of remainder
        else:
            processes = [None] * procs # Batch size specified by procs
        print(processes)
        for j in range(len(processes)): # Send each process of batch to a child process
            k = (procs*i) + j # Total index count of all times that have been sent out across all batches
            print(f'##### {k} #####')
            processes[j] = Process(target=times_parallel, args=(varList[k],))
            processes[j].start()
        for j in range(len(processes)):
            # print(f'################ {processes}')
            processes[j].join()

def times_parallel(varsList):
    tt = varsList[0]
    time = varsList[1]
    master_d = varsList[2]
    control_d = varsList[3]
    nproc_surf = varsList[4]
    datadir = varsList[5]
    datafile = varsList[6]
    level_spec_d = varsList[7]
    lon = varsList[8]
    lat = varsList[9]
    grid_list = varsList[10]
    depth_list = varsList[11]
    nodez = varsList[12]
    debug = varsList[13]
    pid_file = varsList[14]
    verbose = varsList[15]
    print( now(), 'grid_maker.py: Processing time = ', time) 
    # print( now(), 'grid_maker.py: Processing time = ', tt) 
    if 'Ma' in time:
    
        # strip off units and make a number
        time = float( time.replace('Ma', '') )

        # determine what time steps are available for this age 
        # NOTE: 'temp' is requried to set which output files to check 
        found_d = Core_Citcom.find_available_timestep_from_age( master_d, 'temp', time )

        print( now(), 'grid_maker.py: WARNING: Adjusting times to match available data:')
        print( now(), '  request_age =', found_d['request_age'], '; request_timestep =', found_d['request_timestep'], '; request_runtime =', found_d['request_runtime'])
        print( now(), '  found_age =', found_d['found_age'], '; found_timestep =', found_d['found_timestep'], '; found_runtime =', found_d['found_runtime'])

        # set variables for subsequent loops
        timestep = found_d['found_timestep']
        runtime_Myr = found_d['found_runtime']

        # convert the found age to an int
        age_Ma = int(np.around( found_d['found_age'] ) )

        # save age number for grids storage
        age_Ma_storing=age_Ma
         
        age_Ma = '%03d' % age_Ma

    else:

        time = float( time ) 
         
        # determine what time steps are available for this timestep 
        # NOTE: 'temp' is requried to set which output files to check 

        found_d = Core_Citcom.find_available_timestep_from_timestep( master_d, 'temp', time )

        print( now(), 'grid_maker.py: WARNING: Adjusting times to match available data:')
        print( now(), '  request_age =', found_d['request_age'], '; request_timestep =', found_d['request_timestep'], '; request_runtime =', found_d['request_runtime'])
        print( now(), '  found_age =', found_d['found_age'], '; found_timestep =', found_d['found_timestep'], '; found_runtime =', found_d['found_runtime'])

        # set variables for subsequent loops
        timestep = found_d['found_timestep']
        runtime_Myr = found_d['found_runtime']

        # convert the found age to an int
        age_Ma = int(np.around( found_d['found_age'] ) )
        
        # make a string and pad with zeros 
        #age_Ma = '%03d' % age_Ma
        age_Ma = str(age_Ma)
        
    # output dir add - RC
    output_dir_age = int(np.around(found_d['found_age']))

    # report on integer age
    print( now(), '  age_Ma =', age_Ma)
    
    # empty file_data
    file_data = []
    
    # cache for the file_format
    file_format_cache = ''
    

    
    # Loop over sections (fields) 
    for ss, s in enumerate (control_d['_SECTIONS_'] ) :
        print( now(), 'grid_maker.py: Processing section = ', s) 

        # check for required parameter 'field'
        if not 'field' in control_d[s] :
           print('ERROR: Required parameter "field" missing from section.')
           print('       Skipping this section.')
           continue # to next section

        # get the field name 
        field_name = control_d[s]['field']

        # check for compound field
        field_name_req = ''
        if field_name == 'horiz_vmag':
            # save the requested name
            field_name_req = field_name
            # reset to get one component 
            field_name = 'vx'
            
        if field_name == 'Vx':
            # save the requested name
            field_name_req = field_name
            # reset to get one component
            field_name = 'vx'
            
        if field_name == 'Vy':
            # save the requested name
            field_name_req = field_name
            # reset to get one component
            field_name = 'vy'
            
        if field_name == 'Vz':
            # save the requested name
            field_name_req = field_name
            # reset to get one component
            field_name = 'vz'


        print('')
        print( now(), 'grid_maker.py: Processing: field =', field_name) 

        # set the region
        if nproc_surf == 12:
            grid_R = 'g'
            # optionally adjust the lon bounds of the grid to -180/180
            if 'shift_lon' in control_d :
                print( now(), 'grid_maker.py: grid_R set to to "d" : -180/+180/-90/90')
                grid_R = 'd'
            else :
                print( now(), 'grid_maker.py: grid_R set to to "g" : 0/360/-90/90')
        else:
            grid_R  = str(pid_d['lon_min']) + '/' + str(pid_d['lon_max']) + '/'
            grid_R += str(pid_d['lat_min']) + '/' + str(pid_d['lat_max'])

        # get the data file name specifics for this field 
        file_name_component = Core_Citcom.field_to_file_map[field_name]['file']
        print( now(), 'grid_maker.py: file_name_component = ', file_name_component )

        # get the data file column name specifics for this field 
        field_column = Core_Citcom.field_to_file_map[field_name]['column']
        print( now(), 'grid_maker.py: field_column = ', field_column )
        
    
        # create the total citcoms data filenames to read 
        file_format = ''
        
        pos_for_abspath = datadir.find ("/Data")
        
        # check for various data dirs
        if os.path.exists( datadir + '/0/') :
            
            print( now(), 'grid_maker.py: path found = ', datadir + '/0/' )
            file_format = datadir + '/#/' + datafile + '.' + file_name_component + '.#.' + str(timestep)

        elif os.path.exists( datadir + '/' ) :
            print( now(), 'grid_maker.py: path found = ', datadir + '/' )
            file_format = datadir + '/' + datafile + '.' + file_name_component + '.#.' + str(timestep)

        elif os.path.exists('data') :
            print( now(), 'grid_maker.py: path found = ', 'data' )
            file_format = './data/#/' + datafile + '.' + file_name_component + '.#.' + str(timestep)

        elif os.path.exists('Data') :
            print( now(), 'grid_maker.py: path found = ', 'Data' )
            file_format = './Data/#/' + datafile + '.' + file_name_component + '.#.' + str(timestep)
            
        
        # Added path to dynamic topography post-processing - RC
        elif os.path.exists('Age'+str(age_Ma)+'Ma') :
            print( now(), 'grid_maker.py: path found = ', datadir )
            file_format = './Age'+str(age_Ma)+'Ma/#/' + datafile + '.' + file_name_component + '.#.'+ str(timestep)
        
        # Added path to dynamic topography post-processing - RC
        elif os.path.exists('Age'+str(output_dir_age)+'Ma') :
            print( now(), 'grid_maker.py: path found = ', datadir )
            file_format = './Age'+str(output_dir_age)+'Ma/#/' + datafile + '.' + file_name_component + '.#.'+ str(timestep)
        
        # Added path to post-process with remote location - RC
        elif os.path.exists(datadir[:pos_for_abspath]) :
                
             file_format = datadir[:pos_for_abspath]+'/data/#/' + datafile + '.' + file_name_component + '.#.' + str(timestep)
        
        
        # report error 
        else :
            print( now() )
            print('ERROR: Cannot find output data.')
            print('       Skipping this section.')
            print( now(), 'grid_maker.py: file_format = ', file_format)
            continue # to next section
        
        print( now(), 'grid_maker.py: file_format = ', file_format )

        # check if this file data has already been read in 
        if not file_format == file_format_cache: 

            # read data by proc, e.g., velo, visc, comp_nd, surf, botm 
            file_data = Core_Citcom.read_proc_files_to_cap_list( master_d['pid_d'], file_format, field_name)
            # flatten data since we don't care about specific cap numbers for the loop over levels/depths
            file_data = Core_Util.flatten_nested_structure( file_data )
            print( now(), 'grid_maker.py: len(file_data) = ', len(file_data) )

            # update cache for next pass in loop over fields
            file_format_cache = file_format

        # Get the specific column for this field_name
        field_data = np.array( [ line[field_column] for line in file_data ] )
        print( now(), 'grid_maker.py:  len(field_data) = ', len(field_data) )

        # Check for compound field
        if field_name_req == 'horiz_vmag':
            
            # Get the second component data ('vy')
            #field_column = 1
            # RC note - this is now extracting Vy as read by gmt N-S component 
            field_column = 0
            # read data by proc, e.g., velo, visc, comp_nd, surf, botm 
            file_data2 = Core_Citcom.read_proc_files_to_cap_list( master_d['pid_d'], file_format, field_name)
            # flatten data since we don't care about specific cap numbers for the loop over levels/depths
            file_data2 = Core_Util.flatten_nested_structure( file_data2 )
            print( now(), 'grid_maker.py: len(file_data2) = ', len(file_data2) )
            field_data2 = np.array( [ line[field_column] for line in file_data2 ] )
            print( now(), 'grid_maker.py:  len(field_data2) = ', len(field_data) )

            # combine the data and rest the main variable
            field_data3 = np.hypot( field_data, field_data2)
            field_data = field_data3

            # put back field name to requested name
            field_name = field_name_req 
        # end if check on compound field


        print( now(), 'grid_maker.py:  len(field_data) = ', len(field_data) )
        print( now() )
       
        #
        # Loop over levels 

        levels = enumerate( level_spec_d['list'] )
        varList = []
        for level in levels:
            varList.append(level + (tt, ss, s, timestep, age_Ma, runtime_Myr, field_name, field_data, age_Ma_storing, lon, lat, grid_R, grid_list, verbose))

        # Loop over levels in parralel only if multiprocess_depths is set to true
        try:
            if control_d[s]['multiprocess_depths']:
                with Pool() as p:
                    p.map(levels_parallel, varList)
                    p.close()
                    p.join()
            # If set to anything other than true, loop over levels conventionally
            else:
                print('Looping over depths one at a time')
                for level in varList:
                    levels_parallel(level)
        except KeyError:
            # Loop over levels conventioanlly if multiprocess depths doesn't exist at all
            print('Looping over depths one at a time')
            for level in varList:
                levels_parallel(level)


def levels_parallel(varsList):
    ll = varsList[0]
    level = varsList[1]
    tt = varsList[2]
    ss = varsList[3]
    s = varsList[4]
    timestep = varsList[5]
    age_Ma = varsList[6]
    runtime_Myr = varsList[7]
    field_name = varsList[8]
    field_data = varsList[9]
    age_Ma_storing = varsList[10]
    lon = varsList[11]
    lat = varsList[12]
    grid_R = varsList[13]
    grid_list = varsList[14]
    verbose = varsList[15]
    print( now(), 'grid_maker.py: Processing level = ', level) 

    control_d, pid_file, master_d, coor_d, pid_d, datadir, datafile, start_age, output_format, depth_list, nodez, nproc_surf = initialise_variables()

    # ensure level is an int value 
    level = int(level)
    depth = int(depth_list[level])
    # pad the depth value 
    #depth = '%04d' % depth
    depth=str(depth)
    print( now(), '------------------------------------------------------------------------------')
    print( now(), 'grid_maker.py: tt,ss,ll = ', tt, ',', ss, ',', ll, ';')
    print( now(), 'grid_maker.py: summary for', s, ': timestep =', timestep, '; age =', age_Ma, '; runtime_Myr =', runtime_Myr, '; level =', level, '; depth =', depth, ' km; field_name =', field_name)
    print( now(), '------------------------------------------------------------------------------')


    if field_name.startswith('vertical_'):
        # perform a z slice for citcom data 
        field_slice = field_data[level::nodez] # FIXME : how to get a v slice 
       #xyz_filename = datafile + '-' + field_name + '-' + str(age_Ma_storing) + 'Ma-' + str(depth) + 'km.xyz'
        xyz_filename = datafile + '_' + field_name + '_t' + str(age_Ma_storing)+ '_' + str(depth) + '.xyz'


    else:
        # perform a z slice for citcom data 
        field_slice = field_data[level::nodez]
        #xyz_filename = datafile + '-' + field_name + '-' + str(timestep) + '-' + str(depth) + '.xyz'
        #xyz_filename = datafile + '-' + field_name + '-' + str(age_Ma_storing) + 'Ma-' + str(depth) + 'km.xyz'
        xyz_filename = datafile + '_' + field_name + '_t' + str(age_Ma_storing)+ '_' + str(depth) + '.xyz'

    print( now(), 'grid_maker.py: xyz_filename =', xyz_filename)
    
    if field_name == 'visc': field_slice = np.log10( field_slice )

    print( now(), 'grid_maker.py: type(field_slice) = ', type(field_slice) )
    print( now(), 'grid_maker.py:  len(field_slice) = ', len(field_slice) )
    print( now() )


    # create the xyz data
    xyz_data = np.column_stack( (lon, lat, field_slice) )
    np.savetxt( xyz_filename, xyz_data, fmt='%f %f %f' )

    #print( now(), 'grid_maker.py: type(xyz_data) = ', type(xyz_data) )
    #print( now(), 'grid_maker.py:  len(xyz_data) = ', len(xyz_data) )
    #print( now() )

    # recast the slice 
    #fs = np.array( field_slice )  
    #fs.shape = ( len(lat), len(lon) )
    #print( now(), 'grid_maker.py: type(fs) = ', type(field_slice) )
    #print( now(), 'grid_maker.py:  len(fs) = ', len(field_slice) )
    #print( now() )

    # check for a grid_R 
    if 'R' in control_d[s] :
        grid_R = control_d[s]['R']

    # create the median file 
    median_xyz_filename = xyz_filename.rstrip('xyz') + 'median.xyz'

    blockmedian_I = control_d[s].get('blockmedian_I', '0.5')
    cmd = xyz_filename + ' -I' + str(blockmedian_I) + ' -R' + grid_R

    Core_isoGMT.callgmt( 'blockmedian', cmd, '', '>', median_xyz_filename )

    # get a T value for median file 
    if not 'Ll' in control_d[s] or not 'Lu' in control_d[s]:
        T = Core_isoGMT.get_T_from_minmax( median_xyz_filename )
    else:
        dt = (control_d[s]['Lu']-control_d[s]['Ll'])/10
        T = '-T' + str(control_d[s]['Ll']) + '/'
        T += str(control_d[s]['Lu']) + '/' + str(dt)

    print( now(), 'grid_maker.py: T =', T)

    
    # create the grid
    grid_filename = xyz_filename.rstrip('xyz') + 'nc'

    surface_I = control_d[s].get('surface_I', '0.25')
    cmd = median_xyz_filename + ' -I' + str(surface_I) + ' -R' + grid_R 

    if 'Ll' in control_d[s]:
        cmd += ' -Ll' + str(control_d[s]['Ll'])
    if 'Lu' in control_d[s]:
        cmd += ' -Lu' + str(control_d[s]['Lu'])
    if 'T' in control_d[s]:
        cmd += ' -T' + str(control_d[s]['T'])

    #opt_a = 
    try:
        print(f'{now()} Trying the spherical interpolator')
        Core_isoGMT.callgmt( 'sphinterpolate', cmd, '', '', ' -G' + grid_filename )
    except:
        print(f'{now()} Spherical interpolator unsuccesful. Using gmt surface instead. This may cause some issues around the poles')
        Core_isoGMT.callgmt( 'surface', cmd, '', '', ' -G' + grid_filename )
    else:
        print(f'{now()} Spherical interpolator worked succesfully')

    # Produce a grid showing deviation from average
    if control_d[s].get('deviation'):
        ave = np.median(field_slice)
        field_slice_dev = field_slice - ave

        xyz_filename_dev = 'deviation_' + datafile + '_' + field_name + '_t' + str(age_Ma_storing)+ '_' + str(depth) + '.xyz'
        # create the xyz data
        xyz_data_dev = np.column_stack( (lon, lat, field_slice_dev) )
        np.savetxt( xyz_filename_dev, xyz_data_dev, fmt='%f %f %f' )

        # create the median file 
        median_xyz_filename_dev = xyz_filename_dev.rstrip('xyz') + 'median.xyz'

        cmd = xyz_filename_dev + ' -I' + str(blockmedian_I) + ' -R' + grid_R
        Core_isoGMT.callgmt( 'blockmedian', cmd, '', '>', median_xyz_filename_dev )
        
        # create the grid
        grid_filename_dev = xyz_filename_dev.rstrip('xyz') + 'nc'

        cmd = median_xyz_filename_dev + ' -I' + str(surface_I) + ' -R' + grid_R 
        if 'Ll' in control_d[s]:
            cmd += ' -Ll' + str(control_d[s]['Ll'])
        if 'Lu' in control_d[s]:
            cmd += ' -Lu' + str(control_d[s]['Lu'])
        if 'T' in control_d[s]:
            cmd += ' -T' + str(control_d[s]['T'])

        try:
            print(f'{now()} Trying the spherical interpolator')
            Core_isoGMT.callgmt( 'sphinterpolate', cmd, '', '', ' -G' + grid_filename_dev )
        except:
            print(f'{now()} Spherical interpolator unsuccesful. Using gmt surface instead. This may cause some issues around the poles')
            Core_isoGMT.callgmt( 'surface', cmd, '', '', ' -G' + grid_filename_dev )
        else:
            print(f'{now()} Spherical interpolator worked succesfully')

        dev_dir_name = f'{field_name}_deviation'
        dev_grid_dir=f'{datafile}/{dev_dir_name}/{age_Ma_storing}'

        os.makedirs(f'{dev_grid_dir}', exist_ok=True)
        if os.path.isfile(f'{dev_grid_dir}/{grid_filename_dev}'):
            os.remove(f'{dev_grid_dir}/{grid_filename_dev}')
        shutil.move(grid_filename_dev, f'{dev_grid_dir}')

        if not 'debug' in control_d:
            os.remove(xyz_filename_dev)
            os.remove(median_xyz_filename_dev
                )

    ### Jono- uncomment below to produce plots
    if 'debug' in control_d:
        # label the variables
        
        # −Dxname/yname/zname/scale/offset/title/remark
        cmd = grid_filename + ' -D/=/=/' + str(field_name) + '/=/=/' + str(field_name) + '/' + str(field_name)
        Core_isoGMT.callgmt( 'grdedit', cmd, '', '', '')
    
    # Dimensionalize grid   
    if control_d[s].get('dimensional'):
        print( now(), 'grid_maker.py: dimensional = ', control_d[s]['dimensional'])
        dim_grid_name = grid_filename.replace('.nc', '_dimensional.nc')
        Core_Citcom.dimensionalize_grid(pid_file, field_name, grid_filename, dim_grid_name, verbose=verbose)

        dim_dir_name = f'{field_name}_dimensional'

    # save this grid and its age in a list
    if control_d[s].get('dimensional'):
        grid_list.append( (dim_grid_name, age_Ma) )
    else: 
        grid_list.append( (grid_filename, age_Ma) )

    # Optional step to transform grid to plate frame
    if 'make_plate_frame_grid' in control_d :
        cmd = 'frame_change_pygplates.py %(age_Ma)s %(grid_filename)s %(grid_R)s' % vars()
        print(now(), 'grid_maker.py: cmd =', cmd)
        os.system(cmd)


    # Assoicate this grid with GPlates exported line data in .xy format:
    # compute age value 
    age_float = 0.0

    # time_list values for citcom data uses timesteps; get age 
    time_triple = Core_Citcom.get_time_triple_from_timestep(master_d['time_d']['triples'], timestep, verbose=verbose)
    age_float = time_triple[1]

    if 'debug' in control_d:
        # truncate to nearest int and make a string for the gplates .xy file name 
        if age_float < 0: age_float = 0.0
        xy_path = master_d['geoframe_d']['gplates_line_dir']
        #xy_filename = xy_path + '/' + 'topology_platepolygons_' + str(int(age_float)) + '.00Ma.xy'
        xy_filename = xy_path + '/' + 'topology_platepolygons_' + str(int(age_Ma)) + '.00Ma.xy'
        print( now(), 'grid_maker.py: xy_filename = ', xy_filename)


        # Make a plot of the grids
        J = 'X5/3' #'R0/6'
        #J = 'M5/3'
        if 'J' in control_d[s] :
            J = control_d[s]['J']

        C = 'polar'
        if 'C' in control_d[s] :
            C = control_d[s]['C']
    
        # citcoms 
        # plot non-dimensional grid
        Core_isoGMT.plot_grid( grid_filename, xy_filename, grid_R, T, J, C)

        # also plot dimensional grid 
        if control_d[s].get('dimensional') :
            print( now(), 'grid_maker.py: plotting dimensional = ', control_d[s]['dimensional'])
            dim_grid_name = grid_filename.replace('.nc', '_dimensional.nc')
            T = Core_isoGMT.get_T_from_grdinfo( dim_grid_name )
            Core_isoGMT.plot_grid( dim_grid_name, xy_filename, grid_R, T, J)

    # plot plate frame grid 
    if 'make_plate_frame_grid' in control_d :
        plateframe_grid_name = grid_filename.replace('.nc', '-plateframe.nc')
        xy_filename = ''
        xy_path = master_d['geoframe_d']['gplates_line_dir']
        # present day plate outlines : use '0' 
        xy_filename = xy_path + '/' + 'topology_platepolygons_0.00Ma.xy' 
        print( now(), 'grid_maker.py: xy_filename = ', xy_filename)

        T = Core_isoGMT.get_T_from_grdinfo( plateframe_grid_name )
        print( now(), 'grid_maker.py: T =', T)
        Core_isoGMT.plot_grid( plateframe_grid_name, xy_filename, grid_R, T, J)
    # end of plotting 

    # For normal (non-debug) mode, the produced grids go into neat folders
    # JONO - create field and age directories if needed. Done here
    # os.makedirs(field_name, exist_ok=True)
    grid_dir=f'{datafile}/{field_name}/{age_Ma_storing}'
    os.makedirs(grid_dir, exist_ok=True)
    
    if os.path.isfile(f'{grid_dir}/{grid_filename}'):
        os.remove(f'{grid_dir}/{grid_filename}')
    shutil.move(grid_filename, f'{grid_dir}')

    # Add dimensionalised grid to its own folder
    if control_d[s].get('dimensional'):
        dim_grid_dir=f'{datafile}/{dim_dir_name}/{age_Ma_storing}'
        os.makedirs(f'{dim_grid_dir}', exist_ok=True)

        if os.path.isfile(f'{dim_grid_dir}/{dim_grid_name}'):
            os.remove(f'{dim_grid_dir}/{dim_grid_name}')
        shutil.move(dim_grid_name, f'{dim_grid_dir}')

    if 'debug' in control_d:
        if os.path.isfile(f'{grid_dir}/{xyz_filename}'):
            os.remove(f'{grid_dir}/{xyz_filename}')
        shutil.move(xyz_filename, f'{grid_dir}')

        if os.path.isfile(f'{grid_dir}/{median_xyz_filename}'):
            os.remove(f'{grid_dir}/{median_xyz_filename}')
        shutil.move(median_xyz_filename, f'{grid_dir}')

        ps = grid_filename.rstrip('.nc') + '.ps'
        if os.path.isfile(f'{grid_dir}/{ps}'):
            os.remove(f'{grid_dir}/{ps}')
        shutil.move(ps, f'{grid_dir}')  

        png = grid_filename.rstrip('.nc') + '.png'
        if os.path.isfile(f'{grid_dir}/{png}'):
            os.remove(f'{grid_dir}/{png}')
        shutil.move(png, f'{grid_dir}')                                            

        cpt = grid_filename.rstrip('.nc') + '.cpt'
        if os.path.isfile(f'{grid_dir}/{cpt}'):
            os.remove(f'{grid_dir}/{cpt}')
        shutil.move(cpt, f'{grid_dir}')  

        if control_d[s].get('dimensional'):
            ps = dim_grid_name.rstrip('.nc') + '.ps'
            if os.path.isfile(f'{dim_grid_dir}/{ps}'):
                os.remove(f'{dim_grid_dir}/{ps}')
            shutil.move(ps, f'{dim_grid_dir}')  

            png = dim_grid_name.rstrip('.nc') + '.png'
            if os.path.isfile(f'{dim_grid_dir}/{png}'):
                os.remove(f'{dim_grid_dir}/{png}')
            shutil.move(png, f'{dim_grid_dir}')                                            

            cpt = dim_grid_name.rstrip('.nc') + '.cpt'
            if os.path.isfile(f'{dim_grid_dir}/{cpt}'):
                os.remove(f'{dim_grid_dir}/{cpt}')
            shutil.move(cpt, f'{dim_grid_dir}') 


    # remove some of the unneeded files
    if not 'debug' in control_d:
        os.remove(xyz_filename)
        os.remove(median_xyz_filename)

#=====================================================================
#=====================================================================
# SAVE This code for reference:
#                    # optionally adjust the lon bounds of the grid to -180/180
#                    #if 'shift_lon' in control_d : 
#                    #    print( now(), 'grid_maker.py: shifting values to -180/+180')
#                    #    arg = grid_filename
#                    #    opts = {'R' : 'd', 'S' : '' }
#                    #    Core_isoGMT.callgmt('grdedit', arg, opts)
#=====================================================================
#=====================================================================
def make_example_config_file( ):
    '''print to standard out an example configuration file for this script'''

    text = '''#=====================================================================
# config.cfg file for the grid_maker.py script
# This example has information on creating grids for Citcom data.
# ==============================================================================
# Set the basic model coordinate information common to both source types:

# set path to the model pid file 
pid_file = pid18637.cfg 

# For RS set path
pid = model_restart_age.cfg

# CitcomS coordinate files by processor (i.e. [datafile].coord.[proc])
# first, look in this user-specified directory for all files
coord_dir = coord

# second, look in data/%RANK

# NOTE: grid_maker.py will fail if coord files cannot be located

# Optional global settings

# If 'shift_lon' is set to True, then the grids will have data in the -180/+180 longitude range
# The default is for data in the 0/360 longitude range.
# shift_lon = True

# If 'make_plate_frame_grid' is set to True, then this script with produce additional data on the plate frame
#make_plate_frame_grid = True

# ==============================================================================
# Set the times to grid 

# Citcoms : use model time steps or reconstruction age, in one of these forms:

# single value:
#time_spec = 4500
#time_spec = 7Ma

# comma separated list:
#time_spec = 5400, 6200, 6961
#time_spec = 0Ma, 50Ma, 100Ma

# range of values: start/stop/step
# time_spec = 2000/2700/100
time_spec = 0Ma/10Ma/2Ma

# ==============================================================================
# Set the levels to grid 

# Citcoms : use int values from 0 to nodez-1, in one of these forms:

# single value:
#level_spec = 63

# comma separated list:
#level_spec = 64/0/10

# range of values: start/stop/step
#level_spec = 64/0/10

# NOTE: The level_spec settings must match the Citcoms field data types:

# Volume data fields : the level_spec may be values from 0 to nodez-1
level_spec = 63

# Surface (surf) data fields : the level_spec must be set to only nodez-1 
#level_spec = 64 ; for surf data 

# Botttom (botm) data fields : lthe evel_spec must be set to only 0 
#level_spec = 0 ; for bot data 

# ==============================================================================
# Set the fields to grid 
#
# Each field will be a separate section, delimited by brackets [Section_1], 
# each field requires the field name, e.g. 
# field = temp
# Each field may set optional arguments to set GMT parameters.
# Each field may set the optional parameter 'dimensional' to 'True',
# to produce an additional dimensionalized grid with the '.dimensional' filename component.
#
# See Core_Citcoms.field_to_file_map data struct for a list of field names.
#

# Citcoms :

[Grid_1]
field = temp
dimensional = True
#blockmedian_I = 0.5
#surface_I = 0.25
#Ll = 
#Lu = 
#T = 

#[Grid_2]
#field = surf_topography
#dimensional = True
#blockmedian_I = 0.5
#surface_I = 0.25
#Ll = 
#Lu =
#T = 
#=====================================================================
'''
    print( text )
#=====================================================================
#=====================================================================
if __name__ == "__main__":

    # print ( str(sys.version_info) ) 

    # check for script called wih no arguments
    if len(sys.argv) > 3:
        usage()
        sys.exit(-1)

    # create example config file 
    if sys.argv[1] == '-e':
        make_example_config_file()
        sys.exit(0)


    # run the main script workflow
    main()
    sys.exit(0)
#=====================================================================
#=====================================================================
