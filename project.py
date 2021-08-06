import os
import sys
import copy

import numpy as np

from pyaspect.specfemio.utils import make_records
from pyaspect.specfemio.utils import _mk_symlink
from pyaspect.specfemio.utils import _copy_recursive_dir
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


def _make_proj_dirs(fqdn,access_rights=0o755):
    if os.path.isdir(fqdn):
        raise OSError(f'The directory {fqdn} has already been created')
    try:
        os.makedirs(fqdn, access_rights)
    except OSError:
        print(f'Creation of the directory {fqdn} failed')


def _recursive_proj_dirs(dl,pdir,access_rights=0o755):

    if len(dl.keys()) == 0:
        return 
    else:
        for dl_key in dl.keys():
            new_dir = os.path.join(pdir, dl_key)
            _make_proj_dirs(new_dir,access_rights=0o755)
            _recursive_proj_dirs(dl[dl_key],new_dir)


def make_project(proj_name,
                 proj_root_path,
                 parfile_fqp,
                 mesh_fqp,
                 spec_fqp,
                 pyutils_fqp,
                 script_fqp,
                 src_list,
                 rec_list,
                 obs_rec_list=None,
                 parfile_kv_dict=None,
                 copy_mesh=False,
                 ignore_spec_max=False):

        
        # Check src_list and rec_list (most likely to have a mistake)
        if not isinstance(src_list,list):
            raise TypeError('src_list must be a list type')

        if not isinstance(rec_list,list):
            raise TypeError('rec_list must be a list type')
        
        if obs_rec_list != None:
            if not isinstance(obs_rec_list,list):
                raise TypeError('obs_rec_list must be a list type')

        batch_size = None
        nevents    = None
        
        if isinstance(src_list[0],list):
            if not isinstance(src_list[0][0],SolutionHeader):
                raise TypeError(f'elemts of src_list[:][:] must be of type={type(SolutionHeader)}')
            batch_size = len(src_list[0])
        else:
            if not isinstance(src_list[0],SolutionHeader):
                raise TypeError(f'elements of src_list[:] must be of type={type(SolutionHeader)}')
            batch_size = 1

        nevents    = len(src_list)

        # NON-Batch src 
        if batch_size == 1:
            if not isinstance(rec_list[0],list):
                raise TypeError(f'rec_list[:][:] must be a list type')
            if not isinstance(rec_list[0][0],StationHeader):
                raise TypeError(f'elemts of rec_list[:][:] must be of type={type(StationHeader)}')

            # if doing fwi 
            # TODO: better to not repeat code
            if obs_rec_list != None:
                if not isinstance(obs_rec_list[0],list):
                    raise TypeError(f'obs_rec_list[:][:] must be a list type')
                if len(obs_rec_list) != nevents:
                    raise Exception(f'number of observed data traces is inconsitant with number of receiver stations')
                #TODO: include Trace when when finished
                #if not isinstance(obs_rec_list[0,0],Trace):
                    #raise TypeError(f'elemts of obs_rec_list[0,:] must be of type={type(StationHeader)}')

        # Batch src 
        elif 1 < batch_size:
            if not isinstance(rec_list[0][0],list):
                raise TypeError(f'rec_list[:][:] must be a list type')
            if not isinstance(rec_list[0][0][0],StationHeader):
                raise TypeError(f'elemts of rec_list[:][:][:] must be of type={type(StationHeader)}')

            if obs_rec_list != None:
                if len(obs_rec_list) != nevents:
                    raise Exception(f'number of observed data traces is inconsitant with number of receiver stations')
                if not isinstance(obs_rec_list[0][0],list):
                    raise TypeError(f'obs_rec_list[:][:] must be a list type')
                #TODO: include Trace when when finished
                #if not isinstance(obs_rec_list[0,0,0],Trace):
                #    raise TypeError(f'elemts of obs_rec_list[:][:][:] must be of type={type(StationHeader)}')

        # not designed/implemented for other types
        else:
            raise Exception('structure of src_list and rec_list are not complient')


        if len(rec_list) != nevents:
            raise Exception(f'number of receiver stations is inconsitant with number of sourcs')



        if not isinstance(spec_fqp,str):
            raise TypeError('spec_fqp must be a str type')

        if not isinstance(pyutils_fqp,str):
            raise TypeError('pyutils_fqp must be a str type')
        
        if not isinstance(script_fqp,str):
            raise TypeError('script_fqp must be a str type')
        
        if not isinstance(proj_root_path,str):
            raise TypeError('proj_root_path must be a str type')


        if MAX_SPEC_SRC < nevents and not ignore_spec_max:
            excpt_str  = f'The number events exceeds MAX_SPEC_SRC.\n'
            excpt_str += f'To create more than {MAX_SPEC_SRC} directoies,'
            excpt_str += f'please set \'ignore_spec_max\' to True.'
            raise Exception(excpt_str)


        # set paths
        spec_bin_fqp   = os.path.join(spec_fqp, 'bin')
        spec_utils_fqp   = os.path.join(spec_fqp, 'utils')
        proj_fqdn = os.path.join(proj_root_path, proj_name)
        
        # set directory access
        access_rights = 0o755

        # check if forward or inversion
        fwd_only = True
        if obs_rec_list != None:
            fwd_only = False


        # create records now so that headers can be checked
        l_records = make_records(l_src=src_list,l_rec=rec_list)


        ######################################3
        #
        #create_project_dirs and write records
        #
        ######################################3
        
        # create main project dir
        err = _make_proj_dirs(proj_fqdn)
        
        # create project level symlinks 
        lname = 'pyutils'
        src = pyutils_fqp
        dst = os.path.join(proj_fqdn, lname)
        _mk_symlink(src,dst)

        lname = 'scriptutils'
        src = script_fqp
        dst = os.path.join(proj_fqdn, lname)
        _mk_symlink(src,dst)

        # read and setup Par_file
        par_keys = ['SIMULATION_TYPE','SAVE_FORWARD','MODEL','SAVE_MESH_FILES','USE_BINARY_FOR_SEISMOGRAMS']

        par_lines = readlines(parfile_fqp)

        keys_vals_dict = dict(zip(par_keys,[1,False,'gll',False,True]))
        par_lines = change_multiple_parameters_in_lines(par_lines,keys_vals_dict)

        if not fwd_only:
            keys_vals_dict = dict(zip(par_keys,[1,True,'gll',False,True]))
            par_lines = change_multiple_parameters_in_lines(par_lines,keys_vals_dict)

        #loop over number of events and create (run####) dirs
        for e in range(nevents):

            # make event dir
            
            rdir = 'run' + str(e+1).zfill(4)
            edir = os.path.join(proj_fqdn, rdir)
            err  = _make_proj_dirs(edir)
            
            ddir = os.path.join(edir, 'DATA')

            # make sim links for each event dir (related to the computational node(s) filesytem
            lname = 'bin'
            src = spec_bin_fqp
            dst = os.path.join(edir, lname)
            _mk_symlink(src,dst)

            lname = 'utils'
            src = spec_utils_fqp
            dst = os.path.join(edir, lname)
            _mk_symlink(src,dst)

            
            # make subdirectorieds for each event
            _recursive_proj_dirs(common_dir_struct,edir)

            # make extra fwi dirs if needed
            if not fwd_only:
                _recursive_proj_dirs(common_fwi_dir_struct,edir)
                
            # make sub-dirs needed only in the primary run0001 dir (used for inversion, model-updating, etc.)
            if e == 0 and not fwd_only:
                _recursive_proj_dirs(primary_dir_struct,edir)

            #write Par_files in DATA dirs
            par_fqdn = os.path.join(ddir, 'Par_file')
            writelines(par_fqdn,par_lines)
                
        #end for e in range(nevents)

        # Now that directories are in place, write the src-rec records
        # Note: the records were created before the project directories
        #       so that parameter checking on headers could be done
        #       before making directories.  We wait to write records so
        #       that we can be sure the directories have been made
        write_records(proj_fqdn,
                      l_records,
                      fname='proj_records',
                      write_record_h=True,
                      write_h=False,
                      auto_name=True,
                      auto_network=True)
            
                
        # finally we copy or sym-link the MESH directory. We 
        # wait becuase if copying this can take a lot of time,
        # so it's smart to make sure everything before this 
        # executes successfully
        mesh_dst_fqp = os.path.join(proj_fqdn, 'MESH-default')
        if copy_mesh:
            _copy_recursive_dir(mesh_fqp,mesh_dst_fqp)
        else:
            _mk_symlink(mesh_fqp,mesh_dst_fqp)

