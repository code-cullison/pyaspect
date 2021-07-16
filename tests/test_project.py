import sys
sys.path.append("..") # added

import numpy as np

from project import *


# To be used as a testing driver only.
def main(*args, **kwargs): #just a test driver

    largs = [arg for larg in args for arg in larg]

    assert '--list_ndpp' in kwargs.keys()

    list_ndpp = [int(num) for num in kwargs['--list_ndpp'].split(',')]

    bsize = 6

    proj_fqp    ='./BAMBOO_RACCOON'
    spec_fqp    ='/specfem/bin/path'
    pyutils_fqp ='/pyaspect/pyutils/path'
    script_fqp  ='/bash_script/utils/path'

    if '--proj_fqp' in kwargs.keys():
        proj_fqp = int(kwargs['proj_fqp'])

    if '--create_proj_dirs' in kwargs.keys():
        proj_fqp = kwargs['--create_proj_dirs']

    if '--spec_fqp' in kwargs.keys():
        spec_fqp = int(kwargs['spec_fqp'])

    if '--pyutils_fqp' in kwargs.keys():
        pyutils_fqp = int(kwargs['pyutils_fqp'])

    if '--script_fqp' in kwargs.keys():
        script_fqp = int(kwargs['script_fqp'])

    f_only = True
    if '--inv' in largs:
        f_only = False

    igm = False
    if '--ignore_max' in largs:
        igm = True

    make_proj_instance = None
    test_proj          = None

    if '--batch' in kwargs.keys():

        bsize = int(kwargs['--batch'])
        list_ndpp = [[[(np.arange(6),1.0,0.0,0.0)]*bsize]*int(num) for num in kwargs['--list_ndpp'].split(',')]
        
        print()
        print(f'Constructing project: list_ndpp={list_ndpp}, bsize={bsize}')

        make_proj_instance = timer((lambda sfp,pyp,bfp,p,dmt,dec,igm: 
                                    BatchCMTProject(sfp,pyp,bfp,p,'TESTING',dmt,dec,igm)))
        del test_proj
        test_proj = make_proj_instance(spec_fqp,pyutils_fqp,script_fqp,proj_fqp,list_ndpp,bsize,igm)
        print('Number of Singular-Events:',test_proj.total_events)

    else:
        print()
        print(f'Constructing project: list_ndpp={list_ndpp}, bsize={bsize}')
        make_proj_instance = timer((lambda sfp,pyp,bfp,p,dmt,fo,igm: 
                                    MultiProject(sfp,pyp,bfp,p,'TESTING',dmt,fo,igm)))
        test_proj = make_proj_instance(spec_fqp,pyutils_fqp,script_fqp,proj_fqp,list_ndpp,f_only,igm)



    print('Number of Poject Directories:',test_proj.num_proj_dir)
    print('Number of Event Directories:',test_proj.num_event_dir)
    print()

    if '--show_dict' in largs:
        print('Poject ID Dictionary:\n',list(test_proj.proj_dict))
        print()

    if '--create_proj_dirs' in kwargs.keys():
        proj_fqp = kwargs['--create_proj_dirs']
        #del test_proj
        #test_proj = make_proj_instance(spec_fqp,pyutils_fqp,script_fqp,proj_fqp,list_ndpp,bsize,f_only,igm)
        #test_proj = make_proj_instance(spec_fqp,pyutils_fqp,script_fqp,proj_fqp,list_ndpp,bsize,igm)

        @timer
        def make_all_dirs(tpj,show_info=False):
            tpj.create_all_sub_project_dirs(show_info)

        print('Creating All Dirs (this may take a lot of time)')
        make_all_dirs(test_proj,show_info=False)


# To be used as a testing driver only.
if __name__=='__main__':
    m_args = []
    m_kwargs = {}
    for arg in sys.argv[1:]: # kwargs
        up_arg = arg.split('=')
        if len(up_arg) == 2:
            m_kwargs[up_arg[0]] = up_arg[1]
        else:
            m_args.append(up_arg)

    main(*m_args, **m_kwargs) # kwargs

