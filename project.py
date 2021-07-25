import os
import sys
import copy

import numpy as np

from functools import wraps
from time import time
from tqdm import tqdm

from pyaspect.specfemio.headers import *

# these are the same symbols at https://realpython.com/directory-tree-generator-python/
# FIXME: to be used later for making the output cool
PIPE = "│"
ELBOW = "└──"
TEE = "├──"
PIPE_PREFIX = "│   "
SPACE_PREFIX = "    "


def timer(func):
    @wraps(func)
    def wrap(*args, **kwargs):
        t_start = time()
        result = func(*args, **kwargs)
        t_end = time()
        print(f'Function \'{func.__name__}({args},{kwargs})\' executed in {(t_end-t_start):4.3f}s\n')
        return result
    return wrap


MAX_SPEC_SRC = int(9999) # see SPECFEM3D_Cartesian manual


class Project:

    def __init__(self,
                 proj_path=None,
                 proj_name=None,
                 nevents=None,
                 spec_fqp=None,
                 pyutils_fqp=None,
                 script_fqp=None,
                 fwd_only=False):
        

        if not isinstance(spec_fqp,str):
            raise TypeError('spec_fqp must be a str type')

        if not isinstance(pyutils_fqp,str):
            raise TypeError('pyutils_fqp must be a str type')
        
        if not isinstance(script_fqp,str):
            raise TypeError('script_fqp must be a str type')
        
        if not isinstance(proj_path,str):
            raise TypeError('proj_path must be a str type')

        if os.path.isdir(proj_path):
            raise Exception(f'{proj_path} already exists')


        self.proj_path  = proj_path
        self.proj_name = proj_name
        self.nevents   = nevents

        self.spec_fqp       = spec_fqp
        self.spec_bin_fqp   = spec_fqp + '/bin'
        self.spec_utils_fqp = spec_fqp + '/utils'

        self.pyutils_fqp = pyutils_fqp
        self.script_fqp  = script_fqp

        self.proj_fqdn = proj_path + '/' + proj_name + '/'

        self.fwd_only = fwd_only

        self.access_rights = 0o755

        self.dir_struct = {}


    def _make_proj_dirs(self,fqdn):
        if os.path.isdir(fqdn):
            print ("The directory %s has already been created" % fqdn)
            return 1
        try:
            os.makedirs(fqdn, self.access_rights)
        except OSError:
            print ("Creation of the directory %s failed" % fqdn)
            print()
            return -1
        else:
            return 0


    def _recursive_mkdir(self,dl,ds,pdir,prefix):

        if len(dl.keys()) == 0:
            return 
        else:
            #print('dl.keys()',dl.keys())
            for dl_key in dl.keys():
                new_dir = pdir + dl_key
                err = self._make_proj_dirs(new_dir)
                if err != 0:
                    print('Oops. Somethig went wrong trying to create %s dir!' %new_dir)
                    assert False
                ds_key = (prefix + ': ',dl_key)
                ds[ds_key] = {}

                self._recursive_mkdir(dl[dl_key],ds[ds_key],new_dir,prefix+prefix)


    def _mk_symlink(self,ds,src,dst,lname,prefix):
        if not os.path.islink(dst):
            os.symlink(src, dst)
        key = (prefix + ': ',lname)
        ds[key] = {}


    def create_project_dirs(self):

        #print(f'in project.create_project_dirs, path:{self.proj_fqdn}')

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



        # create main project dir
        err = self._make_proj_dirs(self.proj_fqdn)
        if err == 0 or err == 1:
            if err == 1:
                return self.dir_struct
        else: 
            print('Oops. Something went wrong trying to creating Project directory!')
            assert False
        key = ('Project Root: ',self.proj_name)
        self.dir_struct[key] = {}
        okey = key

        dir_prefix = '----'
        
        # create project level symlinks 
        lname = 'pyutils'
        src = self.pyutils_fqp
        dst = self.proj_fqdn + '/' + lname
        self._mk_symlink(self.dir_struct[okey],src,dst,lname,dir_prefix)

        lname = 'scriptutils'
        src = self.script_fqp
        dst = self.proj_fqdn + '/' + lname
        self._mk_symlink(self.dir_struct[okey],src,dst,lname,dir_prefix)



        #loop over number of events and create directories
        ookey = okey
        for e in range(self.nevents):

            # calculate meta data
            dir_prefix = '----'
            rdir = '/run' + str(e+1).zfill(4)
            edir = self.proj_fqdn + rdir
            dname = edir
            

            # make event dir
            err = self._make_proj_dirs(edir)
            if err != 0:
                print('Oops. Somethig went wrong trying to create run/event dirs!')
                assert False
            key = (dir_prefix + ': ',rdir)
            self.dir_struct[ookey][key] = {}
            okey = key
            dir_prefix += dir_prefix

                    
            # make sim links for each event dir (related to the computational node(s) filesytem
            lname = 'bin'
            src = self.spec_bin_fqp
            dst = edir + '/' + lname
            self._mk_symlink(self.dir_struct[ookey][okey],src,dst,lname,dir_prefix)

            lname = 'utils'
            src = self.spec_utils_fqp
            dst = edir + '/' + lname
            self._mk_symlink(self.dir_struct[ookey][okey],src,dst,lname,dir_prefix)

            
            # make subdirectorieds for each event
            self._recursive_mkdir(common_dir_struct,self.dir_struct[ookey][okey],dname,dir_prefix)
                
            # make sub-dirs needed only in the primary run0001 dir (used for inversion, model-updating, etc.)
            if e == 0 and not self.fwd_only:
                self._recursive_mkdir(primary_dir_struct,self.dir_struct[ookey][okey],dname,dir_prefix)


        return self.proj_fqdn


    def dir_info(self,ds=None):

        if ds == None:
            keys = self.dir_struct.keys()
            for key in keys:
                print(key[0] + self.proj_fqdn)
                self.dir_info(self.dir_struct[key])

        elif len(ds.keys()) == 0:
            return

        else:
            for key in ds.keys():
                print(key[0] + key[1])
                self.dir_info(ds[key])


class MultiProject(object):
    
    def __init__(self,
                 spec_fqp=None,
                 pyutils_fqp=None,
                 script_fqp=None,
                 proj_path=None,
                 sub_proj_bname=None,
                 list_ndpp=None,
                 fwd_only=True,
                 ignore_spec_max=False):
            
        
        if not isinstance(spec_fqp,str):
            raise TypeError('spec_fqp must be a str type')

        if not isinstance(pyutils_fqp,str):
            raise TypeError('pyutils_fqp must be a str type')
        
        if not isinstance(script_fqp,str):
            raise TypeError('script_fqp must be a str type')
        
        if not isinstance(proj_path,str):
            raise TypeError('proj_path must be a str type')

        if os.path.isdir(proj_path):
            raise Exception(f'{proj_path} already exists')

        #TODO: check dir(s) can be made
        if not isinstance(sub_proj_bname,str):
            raise TypeError('proj_bpath must be a str type')
        
        if not isinstance(list_ndpp,list):
            raise TypeError('list_ndpp must be a list-like type')
        for num in list_ndpp:
            if not isinstance(num,int):
                raise TypeError('values in list_ndpp must be int types')
            if not ignore_spec_max and MAX_SPEC_SRC < num:
                ve_msg  = f'values in list_ndpp must be less than {MAX_SPEC_SRC}'
                ve_msg += ' unless <ignore_spec_max>=True'
                raise ValueError(f'values in list_ndpp must be less than {MAX_SPEC_SRC}')
            

        self.spec_fqp    = spec_fqp
        self.pyutils_fqp = pyutils_fqp
        self.script_fqp  = script_fqp  

        self.proj_path      = proj_path 
        self.sub_proj_bname = sub_proj_bname 

        self.list_ndpp  = list_ndpp     #number of event dirs per project dir

        self.fwd_only = fwd_only
        
        self.ignore_spec_max = ignore_spec_max
       
        self.num_proj_dir = len(self.list_ndpp)

        self.num_event_dir = sum(self.list_ndpp)

        self.is_proj_dirs_created = False
            
        ##############################################
        #                                            # 
        # Now we make project dirs and event dirs    #
        #                                            # 
        ##############################################


        self.proj_dict = {}
        for pdir in range(self.num_proj_dir):
            sub_proj_name = f'{self.sub_proj_bname}_{pdir+1}'
            sub_nevents   = self.list_ndpp[pdir]
            sub_proj = Project(self.proj_path,
                               sub_proj_name,
                               sub_nevents,
                               self.spec_fqp, 
                               self.pyutils_fqp,
                               self.script_fqp,
                               fwd_only=self.fwd_only)
            self.proj_dict[sub_proj_name] = sub_proj


    def create_all_sub_project_dirs(self,show_info=False):
        
        if self.is_proj_dirs_created:
            if show_info:
                for k,v in self.proj_dict.items():
                    print('Project Dir Info:')
                    v.dir_info()
            else:
                print('Sub-Project(s) already created:')

        else:
            for k,v in self.proj_dict.items():
                v.create_project_dirs()
                if show_info:
                    print('Created Sub-Project Dir: {k}')
                    v.dir_info()
            self.is_proj_dirs_created = True

    # end class
            

class BatchCMTProject(Project):
    
    def __init__(self,
                 spec_fqp=None,
                 pyutils_fqp=None,
                 script_fqp=None,
                 proj_path=None,
                 records=None,
                 sub_proj_bname=None,
                 ignore_spec_max=False):

        
        self.bmt_lol = bmt_lol


        #number of sources per event dir
        self.sbsize  = sbsize
        if not isinstance(sbsize,int):
            raise TypeError('sbsize must be an int type')
        
        if sbsize <= 0:
            raise TypeError('sbsize must be a positive integer g.t. zero')
            

        #create dir list for MultiProject super()
        list_ndpp = [int(len(mtl)) for mtl in bmt_lol]

        super(BatchCMTProject,self).__init__(spec_fqp=spec_fqp,
                                             pyutils_fqp=pyutils_fqp,
                                             script_fqp=script_fqp,
                                             proj_path=proj_path,
                                             sub_proj_bname=sub_proj_bname,
                                             list_ndpp=list_ndpp,
                                             fwd_only=True,
                                             ignore_spec_max=ignore_spec_max)

        self.total_events = self.num_event_dir*self.sbsize
            

class BatchCMTMultiProject(MultiProject):
    
    def __init__(self,
                 spec_fqp=None,
                 pyutils_fqp=None,
                 script_fqp=None,
                 proj_path=None,
                 sub_proj_bname=None,
                 bmt_lol=None,
                 sbsize=None,
                 ignore_spec_max=False):

        
        #TODO: devise a beter check for this structure?
        if len(bmt_lol[0][0]) != sbsize:
            emsg  = ' Number of CMT\'s per event directory != sbsize'
            emsg += f' num CMT\'s = {len(bmt_lol[0][0])}, sbsize = {sbsize}'
            raise Exception(emsg)
        
        self.bmt_lol = bmt_lol


        #number of sources per event dir
        self.sbsize  = sbsize
        if not isinstance(sbsize,int):
            raise TypeError('sbsize must be an int type')
        
        if sbsize <= 0:
            raise TypeError('sbsize must be a positive integer g.t. zero')
            

        #create dir list for MultiProject super()
        list_ndpp = [int(len(mtl)) for mtl in bmt_lol]

        super(BatchCMTMultiProject,self).__init__(spec_fqp=spec_fqp,
                                             pyutils_fqp=pyutils_fqp,
                                             script_fqp=script_fqp,
                                             proj_path=proj_path,
                                             sub_proj_bname=sub_proj_bname,
                                             list_ndpp=list_ndpp,
                                             fwd_only=True,
                                             ignore_spec_max=ignore_spec_max)

        self.total_events = self.num_event_dir*self.sbsize

