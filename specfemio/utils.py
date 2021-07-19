import os
import copy

import numpy as np



################################################################################
#
# Helper reading and writing functions: 
# SPECFEM3D STATIONS and SOLUTIONS type files
#
################################################################################

def _get_file_header_paths(fqp,fname):

    path = os.path.relpath(fqp, start=os.curdir)
    fqpname = os.path.join(path, fname)
    header_fqpname = os.path.join(path, f'pyheader.{fname.lower()}')

    return fqpname, header_fqpname



################################################################################
#
# Functions for writing SPECFEM3D SOLUTION type files
#
################################################################################


def forcesolution_2_str(fs):

    fslines_str  = f'{fs.name}\n'
    fslines_str += 'time shift:    %s\n' %(str(fs.tshift))
    fslines_str += 'f0:            %s\n' %(str(fs.f0))
    fslines_str += 'latorUTM:      %s\n' %(str(fs.lat_yc))
    fslines_str += 'longorUTM:     %s\n' %(str(fs.lon_xc))
    fslines_str += 'depth:         %s\n' %(str(fs.depth_km)) #VERIFY
    fslines_str += 'factor force source:            %s\n' %(str(fs.factor_fs))
    fslines_str += 'component dir vect source E:    %s\n' %(str(fs.comp_src_EX))
    fslines_str += 'component dir vect source N:    %s\n' %(str(fs.comp_src_NY))
    fslines_str += 'component dir vect source Z_UP: %s\n' %(str(fs.comp_src_Zup))

    return fslines_str


def cmtsolution_2_str(cmts):

    cmtlines_str  = f'{cmts.name}\n'
    cmtlines_str += 'event name:    %s\n' %(str(cmts.ename))
    cmtlines_str += 'time shift:    %s\n' %(str(cmts.tshift))
    cmtlines_str += 'half duration: %s\n' %(str(cmts.hdur))
    cmtlines_str += 'latorUTM:      %s\n' %(str(cmts.lat_yc))
    cmtlines_str += 'longorUTM:     %s\n' %(str(cmts.lon_xc))
    cmtlines_str += 'depth:         %s\n' %(str(cmts.depth_km))
    cmtlines_str += 'Mrr:           %s\n' %(str(cmts.mrr))
    cmtlines_str += 'Mtt:           %s\n' %(str(cmts.mtt))
    cmtlines_str += 'Mpp:           %s\n' %(str(cmts.mpp))
    cmtlines_str += 'Mrt:           %s\n' %(str(cmts.mrt))
    cmtlines_str += 'Mrp:           %s\n' %(str(cmts.mrp))
    cmtlines_str += 'Mtp:           %s'   %(str(cmts.mtp))

    return cmtlines_str

################################################################################
#
# Functions for writing SPECFEM3D STATION type files
#
################################################################################


def station_to_str(s):
        sname = f'{s.name}_{s.trid}_{s.gid}_{s.sid}'
        snet  = s.network
        slat  = s.lat_yc # or Y coordinate
        slon  = s.lon_xc # or X coordinate
        selev = s.elevation
        sbur  = s.burial
        return '%s %s %.2f %.2f %.2f %.2f\n' %(sname,snet,slat,slon,selev,sbur)


def station_list_to_str(l_stations):

    str_stations = ''
    for i in range(len(l_stations)):
        str_stations += station_to_str(l_stations[i])

    return str_stations


################################################################################
#
# Functions for makeing different lists and groups of stations
#
################################################################################


def make_station_half_cross_members(station=None,delta=None):
    
    if delta == 0:
        raise ValueError('delta cannot equal zero')
    
    station_members = []

    l_gid = np.array([1,2,3])

    if delta < 0:
        l_gid += 3
    
    # add x + delta
    station_xp = copy.deepcopy(station)
    station_xp.lon_xc += delta
    station_xp.gid  = l_gid[0]
    station_members.append(station_xp)
    
    # add y + delta
    station_yp = copy.deepcopy(station)
    station_yp.lat_yc += delta
    station_yp.gid  = l_gid[1]
    station_members.append(station_yp)
    
    # add z + delta
    station_zp = copy.deepcopy(station)
    station_zp.burial += delta
    station_zp.gid  = l_gid[2]
    station_members.append(station_zp)

    return station_members

def make_station_cross_members(station=None,delta=None):

    
    station_members = []

    # positive incremented coordinate group
    station_members += make_station_half_cross_members(station,delta)

    # negative incremented coordinate group
    station_members += make_station_half_cross_members(station,-delta)
    
    return station_members


def make_station_cross_group(station=None,delta=None):

    
    station_group = []

    # add "central" member (no coordinate change)
    cpy_station = copy.deepcopy(station)
    cpy_station.gid = station.gid = 0
    station_group.append(cpy_station) 
    
    # add pos/neg coordinate members
    station_group += make_station_cross_members(station,delta)

    return station_group


def make_station_half_cross_group(station=None,delta=None):

    station_group = []

    # add "central" member (no coordinate change)
    cpy_station = copy.deepcopy(station)
    cpy_station.gid = station.gid = 0
    station_group.append(cpy_station) 
    
    # add pos/neg coordinate members
    station_group += make_station_half_cross_members(station,delta)

    return station_group

    
def make_station_group_list(stations,delta,full_cross=True):
    
    group_station_list = []

    if full_cross:
        for i in range(len(stations)):
            group_station = make_station_cross_group(station=stations[i],delta=delta)
            group_station_list += group_station
    else:
        for i in range(len(stations)):
            group_station = make_station_half_cross_group(station=stations[i],delta=delta)
            group_station_list += group_station
    
    return list(set(group_station_list))

    
def make_station_cross_group_list(stations,delta):
    return make_station_group_list(stations,delta,True)

    
def make_station_half_cross_group_list(stations,delta):
    return make_station_group_list(stations,delta,False)



################################################################################
#
# Functions for getting coordinates
#
################################################################################


def prune_header_list(l_headers,key,val):

    if key not in l_headers[0].keys():
        raise KeyError(f'field: {key}, does not exist')

    pruned = []
    
    for h in l_headers:
        if h[key] != val:
            pruned.append(h)

    return pruned


def get_xyz_coords_from_header(h,depth_key=None):
        return np.array([h.lon_xc,h.lat_yc,h[depth_key]])


def get_xyz_coords_from_station(s):
    return get_xyz_coords_from_header(s,depth_key='burial')


def get_xyz_coords_from_solution(s):
    return get_xyz_coords_from_header(s,depth_key='depth')


def get_xyz_coords_from_header_list(l_headers,depth_key=None):

    xyz = np.zeros((len(l_headers),3))
    
    for i in range(len(l_headers)):
        h = l_headers[i]
        xyz[i,:] = get_xyz_coords_from_header(h,depth_key)[:]

    return xyz


def get_xyz_coords_from_station_list(l_stations):
    return get_xyz_coords_from_header_list(l_stations,depth_key='burial')


def get_xyz_coords_from_solution_list(l_stations):
    return get_xyz_coords_from_header_list(l_stations,depth_key='depth')


def get_xyz_coords_from_headers_except(l_headers,key=None,val=None,depth_key=None):

    p_headers = prune_header_list(l_headers,key,val)

    return get_xyz_coords_from_header_list(p_headers,depth_key)


def get_xyz_coords_from_station_list_except(l_stations,key=None,val=None):
    return get_xyz_coords_from_headers_except(l_stations,key=key,val=val,depth_key='burial')

