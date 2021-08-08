import os
import copy
import shutil

import numpy as np

from functools import wraps
from time import time


################################################################################
#
# Misc. Helper functions
#
################################################################################


def timer(func):
    @wraps(func)
    def wrap(*args, **kwargs):
        t_start = time()
        result = func(*args, **kwargs)
        t_end = time()
        print(f'Function \'{func.__name__}({args},{kwargs})\' executed in {(t_end-t_start):4.3f}s\n')
        return result
    return wrap


################################################################################
#
# Helper reading and writing functions: 
# SPECFEM3D STATIONS and SOLUTIONS type files
#
################################################################################

def _join_path_fname(fqp,fname):
    fqpname = os.path.join(fqp, fname)
    return fqpname

def _join_relpath_fname(fqp,start_fqp,fname):
    path = os.path.relpath(fqp, start=start_fqp)
    fqpname = os.path.join(path, fname)
    return fqpname


def _get_file_path(fqp,fname):
    return _join_path_fname(fqp, fname)


def _get_header_path(fqp,fname):
    return _join_path_fname(fqp, f'pyheader.{fname.lower()}')


def _get_file_header_paths(fqp,fname):

    fqpname = _get_file_path(fqp, fname)
    header_fqpname = _get_header_path(fqp, fname)

    return fqpname, header_fqpname


def _mk_symlink(src,dst):
    if not os.path.islink(dst):
        os.symlink(src, dst)


def _mk_relative_symlink(src_fqp,start_fqp,dst_fqp):
    if not os.path.islink(dst_fqp):
        rel_src_path = os.path.relpath(src_fqp, start_fqp)
        os.symlink(rel_src_path, dst_fqp)

def _copy_recursive_dir(src_fqp,dst_fqp):
    try:
        shutil.copytree(src_fqp,dst_fqp)
    except Exception as e:
        print(e)


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

def station_auto_name(s):
    return f't{str(s.trid).zfill(6)}g{str(s.gid).zfill(2)}'

def network_auto_name(s):
    return f's{str(s.sid).zfill(2)}'

def station_auto_data_fname_id(s):
    net_code  = network_auto_name(s)
    stat_code = station_auto_name(s)
    return f'{net_code}.{stat_code}'

def station_to_str(s,auto_name=False,auto_network=False):

    sname = s.name
    if auto_name:
        sname = station_auto_name(s)
        #sname = f't{str(s.trid).zfill(6)}g{str(s.gid).zfill(2)}'
    if 10 < len(sname):
        # the SPECFEM default is 32 characters, but in my testing I could
        # only use 10 characters -> It might be a FORTRAN complier optimization 
        # issue. The "read(IIN,'(a)',iostat) line" in read_stations is where
        # things go heywire
        raise Exception('trace name is too long -> less than 32 chars requied')

    snet  = s.network
    if auto_network:
        snet  = network_auto_name(s)
        #snet  = f's{str(s.sid).zfill(2)}'

    slat  = s.lat_yc # or Y coordinate
    slon  = s.lon_xc # or X coordinate
    selev = s.elevation
    sbur  = s.depth

    return '%s %s %.2f %.2f %.2f %.2f\n' %(sname,snet,slat,slon,selev,sbur)


def station_list_to_str(l_stations,auto_name=False,auto_network=False):

    str_stations = ''
    for i in range(len(l_stations)):
        str_stations += station_to_str(l_stations[i],
                                       auto_name=auto_name,
                                       auto_network=auto_network)

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
    station_zp.depth += delta
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
    cpy_station['delta'] = delta
    station_group.append(cpy_station) 
    
    # add pos/neg coordinate members
    station_group += make_station_cross_members(cpy_station,delta)

    return station_group


def make_station_half_cross_group(station=None,delta=None):

    station_group = []

    # add "central" member (no coordinate change)
    cpy_station = copy.deepcopy(station)
    cpy_station.gid = station.gid = 0
    cpy_station['delta'] = delta
    station_group.append(cpy_station) 
    
    # add pos/neg coordinate members
    station_group += make_station_half_cross_members(cpy_station,delta)

    return station_group

    
def make_grouped_station_headers(stations,delta,full_cross=True):
    
    group_station_list = []

    if full_cross:
        for i in range(len(stations)):
            group_station = make_station_cross_group(station=stations[i],delta=delta)
            #group_station_list += group_station
            group_station_list.append(group_station)
    else:
        for i in range(len(stations)):
            group_station = make_station_half_cross_group(station=stations[i],delta=delta)
            #group_station_list += group_station
            group_station_list.append(group_station)
    
    #return list(set(group_station_list))
    return group_station_list

    
def make_grouped_cross_station_headers(stations,delta):
    return make_grouped_station_headers(stations,delta,True)

    
def make_grouped_half_cross_station_headers(stations,delta):
    return make_grouped_station_headers(stations,delta,False)


def flatten_grouped_headers(l_group):
    if not isinstance(l_group[0], list):
        raise Exception('group list must be at least 2 dimensional')
    return list(sum(l_group, []))


def flatten_grouped_headers_unique(l_group):
    return list(set(flatten_grouped_headers(l_group)))

def is_grouped_headers(headers):
    is_grouped = False
    if isinstance(headers[0], list):
        is_grouped = True
    return is_grouped



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


def get_xyz_coords_from_header(h):
        return np.array([h.lon_xc,h.lat_yc,h.depth])


def get_xyz_coords_from_station(s):
    return get_xyz_coords_from_header(s)


def get_xyz_coords_from_solution(s):
    return get_xyz_coords_from_header(s)


def get_xyz_coords_from_header_list(l_headers,unique=False):

    xyz = np.zeros((len(l_headers),3))
    
    for i in range(len(l_headers)):
        h = l_headers[i]
        xyz[i,:] = get_xyz_coords_from_header(h)[:]

    if unique:
        xyz = np.unique(tuple(xyz),axis=0)

    return xyz

def get_unique_xyz_coords_from_header_list(l_headers):
    return get_xyz_coords_from_header_list(l_headers,unique=True)


def get_xyz_coords_from_station_list(l_stations):
    return get_xyz_coords_from_header_list(l_stations)


def get_xyz_coords_from_solution_list(l_stations):
    return get_xyz_coords_from_header_list(l_stations)


def get_unique_xyz_coords_from_station_list(l_stations):
    return get_unique_xyz_coords_from_header_list(l_stations)


def get_unique_xyz_coords_from_solution_list(l_stations):
    return get_unique_xyz_coords_from_header_list(l_stations)


def get_xyz_coords_from_headers_except(l_headers,key=None,val=None,unique=False):

    p_headers = prune_header_list(l_headers,key,val)

    return get_xyz_coords_from_header_list(p_headers,unique=unique)

def get_xyz_coords_from_station_list_except(l_stations,key=None,val=None):
    return get_xyz_coords_from_headers_except(l_stations,key=key,val=val)


def get_unique_xyz_coords_from_headers_except(l_headers,key=None,val=None):
    return get_xyz_coords_from_headers_except(l_stations,key=key,val=val,unique=True)

def get_unique_xyz_coords_from_station_list_except(l_stations,key=None,val=None):
    return get_xyz_coords_from_headers_except(l_stations,key=key,val=val,unique=True)


################################################################################
#
# Functions for getting SOLUTION groups and lists
#
################################################################################


def make_triplet_force_solution_members(src=None):

    from pyaspect.specfemio.headers import ForceSolutionHeader
    
    solution_members = []

    # create Force 001 (North/Y axis)
    force_solution_1 = src.copy()
    force_solution_1.sid = 1 
    force_solution_1.comp_src_EX  = 0
    force_solution_1.comp_src_NY  = src.comp_src_EX 
    force_solution_1.comp_src_Zup = 0
    hstr_1 = ForceSolutionHeader.create_header_name(eid=force_solution_1.eid,
                                              sid=force_solution_1.sid,
                                              date=force_solution_1.date,
                                              lat_yc=force_solution_1.lat_yc,
                                              lon_xc=force_solution_1.lon_xc,
                                              depth=force_solution_1.depth)
    force_solution_1.name = hstr_1
    solution_members.append(force_solution_1)

    # create Force 002 (Z axis)
    force_solution_2 = src.copy()
    force_solution_2.sid = 2 
    force_solution_2.comp_src_EX  = 0
    force_solution_2.comp_src_NY  = 0
    force_solution_2.comp_src_Zup = src.comp_src_EX 
    hstr_2 = ForceSolutionHeader.create_header_name(eid=force_solution_2.eid,
                                              sid=force_solution_2.sid,
                                              date=force_solution_2.date,
                                              lat_yc=force_solution_2.lat_yc,
                                              lon_xc=force_solution_2.lon_xc,
                                              depth=force_solution_2.depth)
    force_solution_2.name = hstr_2
    solution_members.append(force_solution_2)
    
    return solution_members

    
def make_triplet_force_solution_group(src=None):

    if src.sid != 0:
        raise ValueError('arg: \'src\' had non-zero solution_id (src.sid)')

    if src.comp_src_EX == 0:
        raise ValueError('arg: \'src\' has comp_src_EX=0')

    return [src] + make_triplet_force_solution_members(src)
    

def make_grouped_triplet_force_solution_headers(solutions=None):
    
    triplet_solution_list = []

    for i in range(len(solutions)):
        triplet_solution = make_triplet_force_solution_group(solutions[i])
        triplet_solution_list.append(triplet_solution)

    return triplet_solution_list


################################################################################
#
# Functions for replicating Stations by Solution list 
#
################################################################################


def make_replicated_station_headers_from_src_list(l_srcs,l_recs):

    from pyaspect.specfemio.headers import SolutionHeader
    from pyaspect.specfemio.headers import StationHeader

    if not isinstance(l_srcs[0],SolutionHeader):
        raise TypeError('l_srcs[:] elements must be of type:{type(SolutionHeader}')

    if not isinstance(l_recs[0],StationHeader):
        raise TypeError('l_recs[:] elements must be of type:{type(StationHeader}')

    l_grp_recs_by_srcs = []
    for i in range(len(l_srcs)):
        src = l_srcs[i]
        cpy_recs = copy.deepcopy(l_recs)
        for rec in cpy_recs:
            rec.sid = src.sid
            rec.eid = src.eid
        l_grp_recs_by_srcs.append(cpy_recs)

    return l_grp_recs_by_srcs


################################################################################
#
# Functions for creating Reciprocal Stations and Solution lists 
#
################################################################################

def make_grouped_reciprocal_station_headers_from_cmt_list(l_cmt,delta,full_cross=True):

    from pyaspect.specfemio.headers import StationHeader

    # Make the main stations
    l_vrecs = []
    for cmt in l_cmt:
        tr_bname = 'tr'
        new_r = StationHeader(name=tr_bname,
                              network='NL', #FIXME
                              lat_yc=cmt.lat_yc,
                              lon_xc=cmt.lon_xc,
                              elevation=0.0,
                              depth=cmt.depth,
                              trid=cmt.eid)
        l_vrecs.append(new_r)

    # Make the group cross stations for the derivatives
    return make_grouped_station_headers(stations=l_vrecs,delta=delta,full_cross=full_cross)

def make_grouped_cross_reciprocal_station_headers_from_cmt_list(l_cmt,delta):
    return make_grouped_reciprocal_station_headers_from_cmt_list(l_cmt,delta,full_cross=True)

def make_grouped_half_cross_reciprocal_station_headers_from_cmt_list(l_cmt,delta):
    return make_grouped_reciprocal_station_headers_from_cmt_list(l_cmt,delta,full_cross=False)


def make_grouped_reciprocal_force_solution_triplet_headers_from_rec_list(l_rec):

    import datetime
    from pyaspect.specfemio.headers import ForceSolutionHeader

    # First we make a single virtual source per receiver
    l_vsrcs = []
    for rec in l_rec:

        #NOTE!!!! the depth is set to the NEGATIVE (makes it positive) due to sign convention
        new_s = ForceSolutionHeader(ename=f'Event-{str(rec.trid).zfill(4)}',
                                    lat_yc=rec.lat_yc,
                                    lon_xc=rec.lon_xc,
                                    depth=rec.depth,
                                    tshift=0.0,
                                    date=datetime.datetime.now(),
                                    f0=0.0,
                                    factor_fs=1,
                                    comp_src_EX=1,
                                    comp_src_NY=0,
                                    comp_src_Zup=0,
                                    proj_id=0,
                                    eid=rec.trid,
                                    sid=0)
        l_vsrcs.append(new_s)

    # Second: Make a force "triplets" used for calculationg moment tensor derivatives
    return make_grouped_triplet_force_solution_headers(solutions=l_vsrcs)


def make_replicated_reciprocal_station_headers_from_src_triplet_list(l_vsrc,l_vrec):

    from pyaspect.specfemio.headers import ForceSolutionHeader
    from pyaspect.specfemio.headers import StationHeader

    if not isinstance(l_vsrc[0][0],ForceSolutionHeader):
        raise TypeError('l_vsrc[:][:] elements must be of type:{type(ForceSolutionHeader}')

    if not isinstance(l_vrec[0][0],StationHeader):
        raise TypeError('l_vrec[:][:] elements must be of type:{type(StationHeader}')

    l_grp_vrecs_by_vsrcs = []
    for i in range(len(l_vsrc)):
        vsgrp = l_vsrc[i]
        vrecs_per_vsrc = []
        for vsrc in vsgrp:
            cpy_vrecs = copy.deepcopy(flatten_grouped_headers(l_vrec))
            for vrec in cpy_vrecs:
                vrec.sid = vsrc.sid
                vrec.eid = vsrc.eid
            vrecs_per_vsrc.append(cpy_vrecs)
        l_grp_vrecs_by_vsrcs.append(vrecs_per_vsrc)

    return l_grp_vrecs_by_vsrcs


################################################################################
#
# Functions for creating Records and lists of Records
#
################################################################################


def make_record_headers(l_src=None,l_rec=None):

    if not isinstance(l_src,list):
        raise Exception('l_src must be a list type')
    if not isinstance(l_rec,list):
        raise Exception('l_rec must be a list type')

    if len(l_src) != len(l_rec):
        raise Exception('Lengths of l_src and l_rec must be equal')

    from pyaspect.specfemio.headers import SolutionHeader
    from pyaspect.specfemio.headers import StationHeader
    from pyaspect.specfemio.headers import RecordHeader

    l_records = []

    is_src_list = isinstance(l_src[0],list)

    if is_src_list:
        if not isinstance(l_src[0][0],SolutionHeader):
            raise TypeError('l_src[:][:] elements must be of type:{type(SolutionHeader}')
        if not isinstance(l_rec[0][0],list):
            raise Exception('Sources are grouped, but receivers appear not to be.')
        if not isinstance(l_rec[0][0][0],StationHeader):
            raise TypeError('l_rec[:][:][:] elements must be of type:{type(StationHeader}')
        
        for i in range(len(l_src)):
            sgrp = l_src[i]
            for j in range(len(sgrp)):
                s = sgrp[j]
                for r in l_rec[i][j]:
                    r.eid = s.eid
                    r.sid = s.sid

            record = RecordHeader(solutions_h=sgrp,stations_h=flatten_grouped_headers(l_rec[i]),rid=i)
            l_records.append(record)

    else:
        if not isinstance(l_src[0],SolutionHeader):
            raise TypeError('l_src[:] elements must be of type:{type(SolutionHeader}')
        if not isinstance(l_rec[0],list):
            raise Exception('l_rec is not complient with l_src')
        if not isinstance(l_rec[0][0],StationHeader):
            raise TypeError('l_rec[:][:] elements must be of type:{type(StationHeader}')

        for i in range(len(l_src)):
            s = l_src[i]
            for r in l_rec[i]: 
                r.eid = s.eid
                r.sid = s.sid

            record = RecordHeader(solutions_h=s,stations_h=l_rec[i],rid=i)
            l_records.append(record)

    return l_records


