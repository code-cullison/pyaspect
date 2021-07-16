import numpy as np
import copy
from pyaspect.model.gridmod3d import gridmod3d as gm

def compress_gm3d_to_file(_gm3d,_path,_prefix,**kwargs):

    # keep track of useful variables
    main_pars = {}

    # get subsurface property arrays
    props = _gm3d.getNPArray()

    # prep for header data
    xmin = _gm3d.get_gorigin()[0]
    dx   = _gm3d.get_deltas()[0]
    nx   = _gm3d.get_npoints()[0]
    main_pars['xmin'] = xmin
    main_pars['dx']   = dx
    main_pars['nx']   = nx

    ymin = _gm3d.get_gorigin()[1]
    dy   = _gm3d.get_deltas()[1]
    ny   = _gm3d.get_npoints()[1]
    main_pars['ymin'] = ymin
    main_pars['dy']   = dy
    main_pars['ny']   = ny

    zmin = _gm3d.get_gorigin()[2]
    dz   = _gm3d.get_deltas()[2]
    nz   = _gm3d.get_npoints()[2]
    main_pars['zmin'] = zmin
    main_pars['dz']   = dz
    main_pars['nz']   = nz

    # header data
    xdata = np.array([xmin,dx,nx])
    ydata = np.array([ymin,dy,ny])
    zdata = np.array([zmin,dz,nz])

    # calculate maxDepth from header data
    maxDepth = (-zmin) + (nz-1)*dz

    # create output filename
    ofqn   = _path + _prefix + '_dx' + str(int(dx)) + '_dy' + str(int(dy)) + '_dz'
    ofqn  += str(int(dz)) + '_maxdepth' + str(int(maxDepth))
    for key in kwargs:
        ofqn  += '_' + str(key) + str(kwargs[key])
    ofqn  +=  '.npz'

    #merge the two dictionaries
    main_pars.update(kwargs)

    print('Compressing file to:',ofqn)
    np.savez_compressed(ofqn,props=props,xd=xdata,yd=ydata,zd=zdata,**main_pars)



def decompress_gm3d_from_file(_fqdn):

    print('Decompressing file from:',_fqdn)

    data = np.load(_fqdn)
    props = data['props'] #4D ndarray of subsurface model

    #header/meta data arrays
    xdata = data['xd']
    ydata = data['yd']
    zdata = data['zd']

    # create gridded model headers
    xmin = xdata[0]
    dx   = xdata[1]
    nx   = int(xdata[2])
    xmax = xmin + (nx-1)*dx #notice that this can be computed

    ymin = ydata[0]
    dy   = ydata[1]
    ny   = int(ydata[2])
    ymax = ymin + (ny-1)*dy #notice that this can be computed

    zmin = zdata[0]
    dz   = zdata[1]
    nz   = int(zdata[2])
    zmax = (-zmin) + (nz-1)*dz #notice that this can be computed

    #get any other parameters from the compression
    other_pars = {}
    for key in data:
        if key != 'props' and key != 'xd' and key != 'yd' and key != 'zd':
            other_pars[key] = data[key]


    # instantiate extracted GriddedModel3D object
    nsub_props = props.shape[0]
    axes_order = {'X':0,'Y':1,'Z':2} #this dict keeps track of axes order
    gm3d = gm(props,nsub_props,axes_order,(nx,ny,nz),(dx,dy,dz),(xmin,ymin,zmin))

    return gm3d,other_pars


