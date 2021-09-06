import os
import sys
import copy

import numpy as np

from pyaspect.specfemio.utils import make_record_headers
from pyaspect.specfemio.utils import _mk_symlink
from pyaspect.specfemio.utils import _copy_recursive_dir
from pyaspect.specfemio.utils import _join_path_fname
from pyaspect.specfemio.utils import _get_header_path
from pyaspect.specfemio.utils import station_auto_data_fname_id
from pyaspect.specfemio.write import _write_header
from pyaspect.specfemio.write import write_record
from pyaspect.specfemio.write import write_records
from pyaspect.specfemio.headers import *
from pyaspect.parfile import change_multiple_parameters_in_lines
from pyaspect.parfile import readlines
from pyaspect.parfile import writelines


MAX_SPEC_SRC = int(9999) # see SPECFEM3D_Cartesian manual

# list of directories need for every event
common_dir_struct = {'DATA': {},
                     'OUTPUT_FILES' : {'DATABASES_MPI':{}},
                     'SYN': {},
                     'FILT_SYN': {} }

# extra common dirs for fwi
common_fwi_dir_struct = {'SEM': {},
                         'OBS': {},
                         'FILT_OBS': {} }

# list of directories only needed for the primary run0001 dir
primary_dir_struct= {'INPUT_GRADIENT': {},
                    'INPUT_KERNELS': {},
                    'INPUT_MODEL': {},
                    'OUTPUT_MODEL': {},
                    'OUTPUT_SUM': {},
                    'SMOOTH': {},
                    'COMBINE': {},
                    'topo': {} }

def _make_dirs(fqdn,access_rights=0o755):
    if os.path.isdir(fqdn):
        raise OSError(f'The directory {fqdn} has already been created')
    try:
        os.makedirs(fqdn, access_rights)
    except OSError:
        print(f'Creation of the directory {fqdn} failed')
        return OSError


def _recursive_proj_dirs(dl,pdir,access_rights=0o755):

    if len(dl.keys()) == 0:
        return
    else:
        for dl_key in dl.keys():
            new_dir = os.path.join(pdir, dl_key)
            _make_dirs(new_dir,access_rights=0o755)
            _recursive_proj_dirs(dl[dl_key],new_dir)


def _make_proj_dir(proj_root_fqp,
                   proj_base_name,
                   pyutils_fqp=None,
                   script_fqp=None):

        projdir_fqp = os.path.join(proj_root_fqp, proj_base_name)
        _make_dirs(projdir_fqp)

        # create project level symlinks
        if pyutils_fqp != None:
            lname = 'pyutils'
            src = pyutils_fqp
            dst = os.path.join(projdir_fqp, lname)
            _mk_symlink(src,dst)

        if script_fqp != None:
            lname = 'scriptutils'
            src = script_fqp
            dst = os.path.join(projdir_fqp, lname)
            _mk_symlink(src,dst)

        return projdir_fqp

def _make_run_dir(irdir,
                  projdir_fqp,
                  spec_bin_fqp,
                  spec_utils_fqp,
                  par_lines,
                  dir_struct,
                  record_h):

    rdir_name = 'run' + str(irdir+1).zfill(4)
    rundir_fqp = os.path.join(projdir_fqp, rdir_name)
    _make_dirs(rundir_fqp)

    # make sim links for each event dir 
    # (related to the computational node(s) filesytem
    lname = 'bin'
    src = spec_bin_fqp
    dst = os.path.join(rundir_fqp, lname)
    _mk_symlink(src,dst)

    lname = 'utils'
    src = spec_utils_fqp
    dst = os.path.join(rundir_fqp, lname)
    _mk_symlink(src,dst)

    # make subdirectorieds for each event
    _recursive_proj_dirs(common_dir_struct,rundir_fqp)

    #write Par_files in DATA dirs
    ddir_fqp = os.path.join(rundir_fqp, 'DATA')
    out_par_fqp  = os.path.join(ddir_fqp, 'Par_file')
    writelines(out_par_fqp,par_lines)
    
    #write Headers and Record
    write_record(rundir_fqp,
                 record_h,
                 fname='record',
                 write_record_h=True,
                 write_h=False,
                 auto_name=True,
                 auto_network=True)


def setup_mesh_dir(proj_fqp, mesh_fqp, copy_mesh=False):

    mesh_dst_fqp = os.path.join(proj_fqp, 'MESH-default')
    if copy_mesh:
        _copy_recursive_dir(mesh_fqp,mesh_dst_fqp)
    else:
        _mk_symlink(mesh_fqp,mesh_dst_fqp)


def make_fwd_project_dir(proj_base_name,
                         proj_root_fqp,
                         parfile_fqp,
                         mesh_fqp,
                         spec_fqp,
                         pyutils_fqp,
                         script_fqp,
                         proj_record_h,
                         sub_proj_name=None,
                         batch_srcs=False,
                         copy_mesh=False,
                         max_event_rdirs=MAX_SPEC_SRC,
                         verbose=False):
    
    
    if not isinstance(proj_record_h,RecordHeader):
        raise TypeError('arg: \'record_h\' must be a RecordHeader type')

    if not isinstance(proj_base_name,str):
        raise TypeError('proj_base_name must be a str type')

    if not isinstance(proj_root_fqp,str):
        raise TypeError('proj_root_fqp must be a str type')

    if not isinstance(parfile_fqp,str):
        raise TypeError('parfile_fqp must be a str type')

    if not isinstance(mesh_fqp,str):
        raise TypeError('mesh_fqp must be a str type')

    if not isinstance(spec_fqp,str):
        raise TypeError('spec_fqp must be a str type')

    if not isinstance(pyutils_fqp,str):
        raise TypeError('pyutils_fqp must be a str type')

    if not isinstance(script_fqp,str):
        raise TypeError('script_fqp must be a str type')

    ########################################################################
    #
    # setup project structure parameters
    #
    ########################################################################
    
    nevents = proj_record_h.nevents
    nsrc = proj_record_h.nsrc
    
    #setup nbatch
    nbatch = 1
    if batch_srcs:
        nbatch = nsrc
            
       
    # calculate number of rundirs (events) per subproject/project
    l_nrundirs = [nevents]
    if not batch_srcs:
        l_nrundirs[0] = nevents*nsrc
        
    # calculate if project or subprojects and addjust rundirs
    if max_event_rdirs < l_nrundirs[0]:
        old_nrundirs = l_nrundirs[0]
        ndiv = l_nrundirs[0]//max_event_rdirs
        l_nrundirs[0] = max_event_rdirs
        for i in range(1,ndiv):
            l_nrundirs.append(max_event_rdirs)
        rem_rundirs = old_nrundirs%max_event_rdirs 
        if rem_rundirs != 0:
            l_nrundirs.append(rem_rundirs)
            
    # actual number of subproject dirs
    nprojdirs = len(l_nrundirs)
    
    # info
    if verbose:
        print(f'nevents: {nevents}')
        print(f'nsrc: {nsrc}')
        print(f'nbatch: {nbatch}')
        print(f'l_nrundirs: {l_nrundirs}')
        print(f'total_rundirs: {sum(l_nrundirs)}')
        print(f'nprojdirs: {nprojdirs}')

    # setup paths to specfem binarys and utils/tools
    spec_bin_fqp   = os.path.join(spec_fqp, 'bin')
    spec_utils_fqp   = os.path.join(spec_fqp, 'utils')
    
    
    # Read input Par_file stub 
    par_lines = readlines(parfile_fqp)
    
    # Setup output Par_files based on user input Par_file stub
    par_keys = ['SIMULATION_TYPE',
                'SAVE_FORWARD',
                'USE_FORCE_POINT_SOURCE',
                'MODEL',
                'SAVE_MESH_FILES',
                'USE_BINARY_FOR_SEISMOGRAMS',
                'SAVE_SEISMOGRAMS_DISPLACEMENT',
                'SAVE_SEISMOGRAMS_VELOCITY']

    #set solution type
    use_force_src = ForceSolutionHeader == proj_record_h.solution_cls
    if use_force_src:
        keys_vals_dict = dict(zip(par_keys,[1,False,use_force_src,'gll',False,True,True,False]))
    else:
        keys_vals_dict = dict(zip(par_keys,[1,False,use_force_src,'gll',False,True,True,False]))
        #keys_vals_dict = dict(zip(par_keys,[1,False,use_force_src,'gll',False,True,False,True]))
    par_lines = change_multiple_parameters_in_lines(par_lines,keys_vals_dict)
    
    
    ########################################################################
    #
    # If more than one proj_dir is needed then make a main proj_dir
    # for the sub_proj_dirs.  
    #
    ########################################################################
    
    #Setup sub_dir pars here
    if sub_proj_name == None:
        sub_proj_name = 'Sub_' + proj_base_name
        
    sub_proj_root_fqp = proj_root_fqp
    sub_pdir_name = proj_base_name
    
    
    # make the main_dir if needed
    if 1 < nprojdirs:
        sub_proj_root_fqp = _make_proj_dir(proj_root_fqp,
                                           proj_base_name)

        # make or copy MESH-default to the main project dir
        setup_mesh_dir(sub_proj_root_fqp, mesh_fqp, copy_mesh=copy_mesh)
        
        
    if verbose: print(f'sub_proj_root_fqp: {sub_proj_root_fqp}')

    accum_src = 0
    for ipdir in range(nprojdirs):
        
        # Make subdir or main project dir depending on number of proj_dirs
        sub_projdir_fqp = sub_proj_root_fqp
        
        if nprojdirs != 1:
            sub_pdir_name = sub_proj_name + '_' + str(ipdir+1).zfill(4)
            
        sub_projdir_fqp = _make_proj_dir(sub_projdir_fqp,
                                         sub_pdir_name,
                                         pyutils_fqp=pyutils_fqp,
                                         script_fqp=script_fqp)
    

        # Copy/Symlink MESH-default 
        if nprojdirs == 1:
            #no sub projects so copy or make symlink (user par)
            setup_mesh_dir(sub_projdir_fqp, mesh_fqp, copy_mesh=copy_mesh)
        else: 
            #multiple subproject: only make symlink for subproj dirs
            sym_mesh_fqp = os.path.join(sub_proj_root_fqp,'MESH-default')
            rel_mesh_fqp = os.path.relpath(sub_proj_root_fqp,sym_mesh_fqp)
            rel_mesh_fqp = os.path.join(rel_mesh_fqp,'MESH-default')
            setup_mesh_dir(sub_projdir_fqp, rel_mesh_fqp, copy_mesh=False)

        
        ####################################################################
        #
        # Make all the run dirs per proj/sub_proj dirs
        #
        ####################################################################
        for irdir in range(l_nrundirs[ipdir]):
            
            #FIXME: "ie" is not Correct! need list[ipdir][irdir] returns ie
            ie   = accum_src//nsrc
            isrc = accum_src%nsrc
            rdir_record_h = proj_record_h[ie,isrc:isrc+nbatch,:,:]
            
            rdir_record_h.solutions_df['proj_id'] = ipdir
            rdir_record_h.stations_df['proj_id']  = ipdir
            
            _make_run_dir(irdir,
                          sub_projdir_fqp,
                          spec_bin_fqp,
                          spec_utils_fqp,
                          par_lines,
                          common_dir_struct,
                          rdir_record_h)
            
            accum_src += nbatch

            ####################################################################
            #
            # Setup and set the SPECFEM station data paths and file prefix names
            #
            ####################################################################

            #setup paths for changing data_fqdn's in stations
            rdir_name = 'run' + str(irdir+1).zfill(4)
            rundir_fqp = os.path.join(sub_projdir_fqp, rdir_name)
            syn_fqp = os.path.join(rundir_fqp, 'SYN')
            #rel_syn_fqp = os.path.relpath(syn_fqp,sub_proj_root_fqp)
            rel_syn_fqp = os.path.relpath(syn_fqp,sub_projdir_fqp)
            

            # reset the pandas.Multiindex to allow eid and sid filtering
            proj_record_h.reset_midx()
            rec_df = proj_record_h.stations_df

            # loop over all stations and set data_fqdn
            ibool = (rec_df['eid'] == ie) & ((rec_df['sid'] >= isrc) & (rec_df['sid'] < isrc+nbatch))
            for ridx, rec in rec_df[ibool].iterrows():
                new_data_fqdn = os.path.join(rel_syn_fqp,station_auto_data_fname_id(rec))
                rec_df.loc[ridx,"data_fqdn"] = new_data_fqdn
                

            # set index back to origanl
            proj_record_h.set_default_midx()
            
    ###############################################
    #
    # write record header for the main project
    #
    ###############################################
    proj_record_fqp = _get_header_path(sub_proj_root_fqp,'project_record')
    if nprojdirs == 1:
        single_proj_fqp = os.path.join(sub_proj_root_fqp,proj_base_name)
        proj_record_fqp = _get_header_path(single_proj_fqp,'project_record')
    _write_header(proj_record_fqp,proj_record_h)


    '''
    # finally we copy or sym-link the MESH directory. We 
    # wait becuase if copying this can take a lot of time,
    # so it's smart to make sure everything before this 
    # executes successfully
    mesh_dst_fqp = os.path.join(proj_fqdn, 'MESH-default')
    if copy_mesh:
        _copy_recursive_dir(mesh_fqp,mesh_dst_fqp)
    else:
        _mk_symlink(mesh_fqp,mesh_dst_fqp)
    '''


