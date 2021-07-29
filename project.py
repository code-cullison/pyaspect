import os
import sys
import copy

import numpy as np

from functools import wraps
from time import time
from tqdm import tqdm

from pyaspect.specfemio.utils make_recods_from
from pyaspect.specfemio.headers import *


MAX_SPEC_SRC = int(9999) # see SPECFEM3D_Cartesian manual

# list of directories need for every event
common_dir_struct = {'/DATA': {},
                     '/OUTPUT_FILES' : {'/DATABASES_MPI':{}},
                     '/SEM': {},
                     '/OBS': {},
                     '/SYN': {},
                     '/FILT_OBS': {},
                     '/FILT_SYN': {} }

# list of directories only needed for the primary run0001 dir
primary_dir_struct= {'/INPUT_GRADIENT': {},
                    '/INPUT_KERNELS': {},
                    '/INPUT_MODEL': {},
                    '/OUTPUT_MODEL': {},
                    '/OUTPUT_SUM': {},
                    '/SMOOTH': {},
                    '/COMBINE': {},
                    '/topo': {} }


def timer(func):
    @wraps(func)
    def wrap(*args, **kwargs):
        t_start = time()
        result = func(*args, **kwargs)
        t_end = time()
        print(f'Function \'{func.__name__}({args},{kwargs})\' executed in {(t_end-t_start):4.3f}s\n')
        return result
    return wrap


class Project:

    def __init__(self,
                 proj_root_path,
                 proj_name,
                 spec_fqp,
                 pyutils_fqp,
                 script_fqp,
                 src_list,
                 rec_list,
                 obs_rec_list=None,
                 ignore_spec_max=False):

        

        # Check src_list and rec_list (most likely to have a mistake)
        if not isinstance(src_list,list):
            raise TypeError('src_list must be a list type')

        if not isinstance(rec_list,list):
            raise TypeError('rec_list must be a list type')
        
        if obs_rec_list != None:
            if not isinstance(obs_rec_list,list):
                raise TypeError('obs_rec_list must be a list type')

        self.batch_size = None
        self.nevents    = None
        
        if isinstance(src_list[0],list):
            if not isinstance(src_list[0,0],SolutionHeader):
                raise TypeError(f'elemts of src_list[0,:] must be of type={type(SolutionHeader)')
            self.batch_size = len(src_list[0])
        else:
            if not isinstance(src_list[0],SolutionHeader):
                raise TypeError(f'elements of src_list[:] must be of type={type(SolutionHeader)')
            self.batch_size = 1

        self.nevents    = len(src_list)

        # NON-Batch src 
        if self.batch_size == 1:
            if not isinstance(rec_list[0],list):
                raise TypeError(f'rec_list[0,:] must be a list type')
            if not isinstance(rec_list[0,0],StationHeader):
                raise TypeError(f'elemts of rec_list[0,:] must be of type={type(StationHeader)')

            # if doing fwi 
            # TODO: better to not repeat code
            if obs_rec_list != None:
                if not isinstance(obs_rec_list[0],list):
                    raise TypeError(f'obs_rec_list[0,:] must be a list type')
                #TODO: include Trace when when finished
                #if not isinstance(obs_rec_list[0,0],Trace):
                    #raise TypeError(f'elemts of obs_rec_list[0,:] must be of type={type(StationHeader)')

        # Batch src 
        elif 1 < self.batch_size:
            if not isinstance(rec_list[0,0],list):
                raise TypeError(f'rec_list[0,:] must be a list type')
            if not isinstance(rec_list[0,0,0],StationHeader):
                raise TypeError(f'elemts of rec_list[0,0,:] must be of type={type(StationHeader)')

            if obs_rec_list != None:
                if not isinstance(obs_rec_list[0,0],list):
                    raise TypeError(f'obs_rec_list[0,:] must be a list type')
                #TODO: include Trace when when finished
                #if not isinstance(obs_rec_list[0,0,0],Trace):
                #    raise TypeError(f'elemts of obs_rec_list[0,0,:] must be of type={type(StationHeader)')

        # not designed/implemented for other types
        else:
            raise Exception('structure of src_list and rec_list are not complient')


        if len(rec_list) != self.nevents:
            raise Exception(f'number of receiver stations is inconsitant with number of sourcs')

        if len(obs_rec_list) != self.nevents:
            raise Exception(f'number of observed data traces is inconsitant with number of receiver stations')


        if not isinstance(spec_fqp,str):
            raise TypeError('spec_fqp must be a str type')

        if not isinstance(pyutils_fqp,str):
            raise TypeError('pyutils_fqp must be a str type')
        
        if not isinstance(script_fqp,str):
            raise TypeError('script_fqp must be a str type')
        
        if not isinstance(proj_root_path,str):
            raise TypeError('proj_root_path must be a str type')

        if os.path.isdir(proj_root_path):
            raise Exception(f'{proj_root_path} already exists')


        if MAX_SPEC_SRC < nevents and not ignore_spec_max:
            excpt_str  = f'The number events exceeds MAX_SPEC_SRC.\n'
            excpt_str += f'To create more than {MAX_SPEC_SRC} directoies,'
            excpt_str += f'please set \'ignore_spec_max\' to True.'
            raise Exception(excpt_str)


        self.src_list = src_list
        self.rec_list = rec_list
        self.obs_rec_list = obs_rec_list

        self.proj_root_path  = proj_root_path
        self.proj_name = proj_name

        self.spec_fqp       = spec_fqp
        self.spec_bin_fqp   = os.path.join(spec_fqp, '/bin')
        self.spec_utils_fqp   = os.path.join(spec_fqp, '/utils')

        self.pyutils_fqp = pyutils_fqp
        self.script_fqp  = script_fqp

        self.proj_fqdn = os.path.join(proj_root_path, proj_name)

        self.fwd_only = True
        if self.obs_rec_list != None:
            self.fwd_only = False

        self.access_rights = 0o755



    def _make_proj_dirs(self,fqdn):
        if os.path.isdir(fqdn):
            raise OSError(f'The directory {fqdn} has already been created')
        try:
            os.makedirs(fqdn, self.access_rights)
        except OSError:
            print(f'Creation of the directory {fqdn} failed')


    def _recursive_mkdir(self,dl,pdir):

        if len(dl.keys()) == 0:
            return 
        else:
            for dl_key in dl.keys():
                new_dir = os.path.join(pdir, dl_key)
                err = self._make_proj_dirs(new_dir)
                self._recursive_mkdir(dl[dl_key],new_dir)


    def _mk_symlink(self,src,dst,lname,prefix):
        if not os.path.islink(dst):
            os.symlink(src, dst)


    def _create_project_dirs(self):



        # create main project dir
        err = self._make_proj_dirs(self.proj_fqdn)
        
        # create project level symlinks 
        lname = 'pyutils'
        src = self.pyutils_fqp
        dst = self.proj_fqdn + '/' + lname
        self._mk_symlink(src,dst,lname)

        lname = 'scriptutils'
        src = self.script_fqp
        dst = self.proj_fqdn + '/' + lname
        self._mk_symlink(src,dst,lname)

        #loop over number of events and create (run####) dirs
        for e in range(self.nevents):

            # make event dir
            
            rdir = '/run' + str(e+1).zfill(4)
            edir = os.path.join(self.proj_fqdn, rdir)
            pdir = edir
            err  = self._make_proj_dirs(edir)

            # make sim links for each event dir (related to the computational node(s) filesytem
            lname = '/bin'
            src = self.spec_bin_fqp
            dst = os.path.join(edir, lname)
            self._mk_symlink(src,dst,lname)

            lname = '/utils'
            src = self.spec_utils_fqp
            dst = os.path.join(edir, lname)
            self._mk_symlink(src,dst,lname)

            
            # make subdirectorieds for each event
            self._recursive_mkdir(common_dir_struct,pdir)
                
            # make sub-dirs needed only in the primary run0001 dir (used for inversion, model-updating, etc.)
            if e == 0 and not self.fwd_only:
                self._recursive_mkdir(primary_dir_struct,pdir)


    def _write_project_files(self):
