import os
import copy
import pickle
import datetime

import numpy as np



################################################################################
#
# Helper reading and writing functions: 
# SPECFEM3D STATIONS and SOLUTIONS type files
#
################################################################################

def _read_headers(fqpame):

    f = open(fqpname, 'rb')
    headers = pickle.load(f)
    f.close()

    return headers


def _write_file(fqpname,file_str):

    f = open(fqpname, 'w')
    f.write(file_str)
    f.close()

def _write_header(fqpname,header):

    f = open(fqpname, 'wb')
    pickle.dump(header,f)
    f.close()


def _get_file_header_paths(fqp,fname):

    import os

    path = os.path.relpath(fqp, start=os.curdir)
    fqpname = os.path.join(path, fname)
    header_fqpname = os.path.join(path, f'pyheader.{fname.lower()}')

    return fqpname, header_fqpname



################################################################################
#
# Functions for writing SPECFEM3D SOLUTION type fiels
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


def write_cmtsolution(fqp,cmts,fname='CMTSOLUTION'):

    cmtlines_str = cmtsolution_2_str(cmts)

    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)

    _write_file(fqp_file,cmtlines_str)
    _write_header(fqp_header,cmts)


def write_forcesolution(fqp,fs,fname='FORCESOLUTION'):

    forcelines_str = forcesolution_2_str(fs)

    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)

    _write_file(fqp_file,forcelines_str)
    _write_header(fqp_header,fs)


def read_solution(fqp,fname):
    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)
    return  _read_headers(fqp_header)


# mostly useful if single CMTSOLUTION per event
def read_cmtsolution(fqp,fname='CMTSOLUTION'):
    return read_solution(fqp,fname)


# mostly useful if single FORCESOLUTION per event
def read_forcesolution(fqp,fname='FORCESOLUTION'):
    return read_solution(fqp,fname)


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
    for i in range(len(stations)):
        str_stations += station_to_str(l_stations[i])

    return str_stations


def write_stations(fqp,l_stations,fname='STATIONS'):

    '''
    path = os.path.relpath(fqp, start=os.curdir)
    fqpname = os.path.join(path, fname)
    header_fqpname = os.path.join(path, f'pyheader.{fname.lower()}')
    '''

    str_stations = station_list_to_str(l_stations)

    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)

    _write_file(fqp_file,str_stations)
    _write_header(fqp_header,cmts)

    '''
    f = open(fqpname, 'w')
    f.write(str_stations)
    f.close()

    f = open(header_fqpname, 'wb')
    pickle.dump(stations,f)
    f.close()
    '''


def read_stations(fqp,fname='STATIONS'):

    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)

    '''
    path = os.path.relpath(fqp, start=os.curdir)
    header_fqpname = os.path.join(path, f'pyheader.{fname.lower}')

    f = open(header_fqpname, 'rb')
    stations = pickle.load(f)
    f.close()
    '''
    
    #return stations
    return  _read_headers(fqp_header)



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
# General Header Class for SPECFEM3D files (STATION and SOLUTIONS)
#
################################################################################


class SpecHeader(dict):

    def __init__(self,name=None,lat_yc=None,lon_xc=None,eid=0,sid=0):
        super(SpecHeader,self).__init__()

        self['name']      = name
        self['lat_yc']    = lat_yc
        self['lon_xc']    = lon_xc
        self['eid']       = eid
        self['sid']       = sid
        self['comp_val']  = None


    def __getattr__(self, key):
        if key in self.keys():
            return self[key]
        else:
            raise AttributeError(f'Header: {key} not found' )

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        if key in self.keys():
            del self[key]
        else:
            raise AttributeError(f'Header: {key} not found' )

    def __str__(self):
        ostr = ''
        for key, value in self.items():
            ostr += f'{key}: {value}\n'
        return ostr

    def __hash__(self):
        return hash(self.hash_val())

    def _is_valid_operand(self, other):
        if type(self) != type(other):
            raise Exception('wrong types when comparing')
        return type(self) == type(other)

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return (self.eq_comparator() == other.eq_comparator())

    def __ne__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return (self.eq_comparator() != other.eq_comparator())

    def __lt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.comparator() < other.comparator()

    def __le__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.comparator() <= other.comparator()

    def __gt__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.comparator() > other.comparator()

    def __ge__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return self.comparator() >= other.comparator()


    def hash_val(self):
        raise NotImplementedError

    def comparator(self):
        raise NotImplementedError

    def eq_comparator(self):
        raise NotImplementedError


    @property
    def name(self):
        return self['name']

    @name.setter
    def name(self, value):
        self['name'] = value


    @property
    def lat_yc(self):
        return self['lat_yc']

    @lat_yc.setter
    def lat_yc(self, value):
        self['lat_yc'] = value


    @property
    def lon_xc(self):
        return self['lon_xc']

    @lon_xc.setter
    def lon_xc(self, value):
        self['lon_xc'] = value


    @property
    def eid(self):
        return self['eid']

    @eid.setter
    def eid(self, value):
        self['eid'] = value


    @property
    def sid(self):
        return self['sid']

    @sid.setter
    def sid(self, value):
        self['sid'] = value


    @property
    def comp_val(self):
        return self['comp_val']

    @comp_val.setter
    def comp_val(self, value):
        self['comp_val'] = value



################################################################################
#
# STATION Classes: STATION[_ADJOINT] files and headers for SPECFEM3D Cart.
#
################################################################################


class StationHeader(SpecHeader):

    def __init__(self,
                 name=None,
                 network=None,
                 lat_yc=None,
                 lon_xc=None,
                 elevation=None,
                 burial=None,
                 eid=0,
                 sid=0,
                 trid=0,
                 gid=0):

        super(StationHeader,self).__init__(name=name,lat_yc=lat_yc,lon_xc=lon_xc,eid=eid,sid=sid)

        if 32 < len(name):
            raise Exception('Station.name cannot exceed 32 characters')

        self['network']   = network
        self['elevation'] = elevation
        self['burial']    = burial
        self['trid']      = trid       # trace id
        self['gid']       = gid        # trace group id (trace decomp)
        self['data_fqdn'] = None


        self.comp_val     = self.trid



    def hash_val(self):
        return np.sqrt(self.lon_xc**2 + self.lat_yc**2 + self.burial**2)

    def comparator(self):
        return self.comp_val

    def eq_comparator(self):
        return (self.lat_yc, self.lon_xc, self.elevation, self.burial)


    @property
    def network(self):
        return self['network']

    @network.setter
    def network(self, value):
        self['network'] = value


    @property
    def elevation(self):
        return self['elevation']

    @elevation.setter
    def elevation(self, value):
        self['elevation'] = value


    @property
    def burial(self):
        return self['burial']

    @burial.setter
    def burial(self, value):
        self['burial']  = value


    @property
    def trid(self):
        return self['trid']

    @trid.setter
    def trid(self, value):
        self['trid'] = value


    @property
    def gid(self):
        return self['gid']

    @gid.setter
    def gid(self, value):
        self['gid'] = value


    @property
    def data_fqdn(self):
        return self['data_fqdn']

    @data_fqdn.setter
    def data_fqdn(self, value):
        self['data_fqdn'] = value


################################################################################
#
# Source Classes: [CMT | FORCE]SOLUTION files and headers for SPECFEM3D Cart.
#
################################################################################

class SolutionHeader(SpecHeader):

    def __init__(self,
                 name=None,
                 tshift=None,
                 lat_yc=None,
                 lon_xc=None,
                 depth=None,
                 eid=0,
                 sid=0):

        super(SolutionHeader,self).__init__(name=name,lat_yc=lat_yc,lon_xc=lon_xc,eid=eid,sid=sid)

        self['tshift']    = tshift
        self['depth']     = depth

        self.comp_val     = eid


    def hash_val(self):
        return np.sqrt(self.lon_xc**2 + self.lat_yc**2 + self.depth**2 + self.tshift**2)

    def comparator(self):
        return self.comp_val

    def eq_comparator(self):
        return (self.lat_yc, self.lon_xc, self.depth, self.tshift)


    @property
    def tshift(self):
        return self['tshift']

    @tshift.setter
    def tshift(self, value):
        self['tshift'] = value


    @property
    def depth(self):
        return self['depth']

    @depth.setter
    def depth(self, value):
        self['depth'] = value

    @property
    def depth_km(self):
        return self['depth']/1000.0


class ForceSolutionHeader(SolutionHeader):
    """
    FORCE  001
    time shift:       0.0
    f0:               0.0
    latorUTM:         591000.0
    longorUTM:        246000.0
    depth:            -100
    factor force source:           1
    component dir vect source E(x):   1
    component dir vect source N(y):   0
    component dir vect source Z_UP:0
    """

    def __init__(self,
             date=None,
             tshift=None,
             f0=None,
             lat_yc=None,
             lon_xc=None,
             depth=None,
             factor_fs=None,
             comp_src_EX=None,
             comp_src_NY=None,
             comp_src_Zup=None,
             eid=0,
             sid=0):


        hstr  = f'FORCE ' + str(sid).zfill(3)
        hstr += f' {date.year} {date.month} {date.day} {date.hour} {date.minute} {date.second}'
        hstr += f' {lat_yc} {lon_xc} {depth/1000} srcid_{eid}'

        super(ForceSolutionHeader,self).__init__(name=hstr,
                                                 tshift=tshift,
                                                 lat_yc=lat_yc,
                                                 lon_xc=lon_xc,
                                                 depth=depth,
                                                 eid=eid,
                                                 sid=sid)


        self['date']         = date
        self['f0']           = f0
        self['factor_fs']    = factor_fs
        self['comp_src_EX']  = comp_src_EX
        self['comp_src_NY']  = comp_src_NY
        self['comp_src_Zup'] = comp_src_Zup


    def hash_val(self):
        return np.sqrt(super().hash_val()**2 + self.f0**2)

    def eq_comparator(self):
        return tuple([*super().eq_comparator(), self.f0])


    @property
    def date(self):
        return self['date']

    @date.setter
    def date(self, value):
        self['date'] = value


    @property
    def f0(self):
        return self['f0']

    @f0.setter
    def f0(self, value):
        self['f0'] = value


    @property
    def factor_fs(self):
        return self['factor_fs']

    @factor_fs.setter
    def factor_fs(self, value):
        self['factor_fs'] = value


    @property
    def comp_src_EX(self):
        return self['comp_src_EX']

    @comp_src_EX.setter
    def comp_src_EX(self, value):
        self['comp_src_EX'] = value


    @property
    def comp_src_NY(self):
        return self['comp_src_NY']

    @comp_src_NY.setter
    def comp_src_NY(self, value):
        self['comp_src_NY'] = value


    @property
    def comp_src_Zup(self):
        return self['comp_src_Zup']

    @comp_src_Zup.setter
    def comp_src_Zup(self, value):
        self['comp_src_Zup'] = value



class CMTSolutionHeader(SolutionHeader):
    """
    PDE  1999 01 01 00 00 00.00  23000 59000 -25000 4.2 4.2 homog_test
    event name:       test
    time shift:       0.0
    half duration:    0.0
    latorUTM:         596000.0
    longorUTM:        241000.0
    depth:            -2700
    Mrr:              0
    Mtt:              0
    Mpp:              0
    Mrt:              0
    Mrp:              0
    Mtp:              -10000000.0
    """

    def __init__(self,
                 date=None,
                 ename=None,
                 tshift=None,
                 hdur=None,
                 lat_yc=None,
                 lon_xc=None,
                 depth=None,
                 mt=None,
                 eid=0,
                 sid=0):

        from pyaspect.moment_tensor import MomentTensor

        if not isinstance(mt,MomentTensor):
            raise TypeError('mt must be of type pyaspect.moment_tensor.MomentTensor')

        hstr  = f'PDE {date.year} {date.month} {date.day} {date.hour} {date.minute} {date.second}'
        hstr += f' {lat_yc} {lon_xc} {depth/1000.0} {mt.Mw} 0 srcid_{eid}'

        super(CMTSolutionHeader,self).__init__(name=hstr,
                                               tshift=tshift,
                                               lat_yc=lat_yc,
                                               lon_xc=lon_xc,
                                               depth=depth,
                                               eid=eid,
                                               sid=sid)


        self['date']   = date
        self['ename']  = ename
        self['hdur']   = hdur
        self['strike'] = mt.strike
        self['dip']    = mt.dip
        self['rake']   = mt.rake
        self['Mw']     = mt.Mw
        self['mrr']    = mt.harvard_dcm_m6()[0]
        self['mtt']    = mt.harvard_dcm_m6()[1]
        self['mpp']    = mt.harvard_dcm_m6()[2]
        self['mrt']    = mt.harvard_dcm_m6()[3]
        self['mrp']    = mt.harvard_dcm_m6()[4]
        self['mtp']    = mt.harvard_dcm_m6()[5]


    def hash_val(self):
        hlist = self.mt
        hlist.append(self.hdur)
        hlist.append(super().hash_val())
        return np.sqrt(np.sum(np.array(hlist)))

    def eq_comparator(self):
        hlist = self.mt
        hlist.append(self.hdur)
        hlist += [*super().eq_comparator()]
        return tuple(hlist)


    @property
    def date(self):
        return self['date']

    @date.setter
    def date(self, value):
        self['date'] = value


    @property
    def ename(self):
        return self['ename']

    @ename.setter
    def ename(self, value):
        self['ename'] = value


    @property
    def hdur(self):
        return self['hdur']

    @hdur.setter
    def hdur(self, value):
        self['hdur'] = value


    @property
    def mt(self):
        return [self.mrr,self.mtt,self.mpp,self.mrt,self.mrp,self.mtp]


    @property
    def strike(self):
        return self['strike']

    @strike.setter
    def strike(self, value):
        self['strike'] = value


    @property
    def dip(self):
        return self['dip']

    @dip.setter
    def dip(self, value):
        self['dip'] = value


    @property
    def rake(self):
        return self['rake']

    @rake.setter
    def rake(self, value):
        self['rake'] = value


    @property
    def Mw(self):
        return self['Mw']

    @Mw.setter
    def Mw(self, value):
        self['Mw'] = value


    @property
    def mrr(self):
        return self['mrr']

    @mrr.setter
    def mrr(self, value):
        self['mrr'] = value


    @property
    def mtt(self):
        return self['mtt']

    @mtt.setter
    def mtt(self, value):
        self['mtt'] = value


    @property
    def mpp(self):
        return self['mpp']

    @mpp.setter
    def mpp(self, value):
        self['mpp'] = value


    @property
    def mrt(self):
        return self['mrt']

    @mrt.setter
    def mrt(self, value):
        self['mrt'] = value


    @property
    def mrp(self):
        return self['mrp']

    @mrp.setter
    def mrp(self, value):
        self['mrp'] = value


    @property
    def mtp(self):
        return self['mtp']

    @mtp.setter
    def mtp(self, value):
        self['mtp'] = value
