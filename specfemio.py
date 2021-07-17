import os
import copy
import pickle

import numpy as np


class SpecHeader(dict):

    def __init__(self,name=None,lat_yc=None,lon_xc=None):
        super(SpecHeader,self).__init__()

        self['name']      = name
        self['lat_yc']    = lat_yc
        self['lon_xc']    = lon_xc
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
        print(f'self.{self.comparator()} < other.{other.comparator()}')
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
    def comp_val(self):
        return self['comp_val']

    @comp_val.setter
    def comp_val(self, value):
        self['comp_val'] = value



class StationHeader(SpecHeader):

    def __init__(self,
                 name=None,
                 network=None,
                 lat_yc=None,
                 lon_xc=None,
                 elevation=None,
                 burial=None,
                 trid=0,
                 gid=0,
                 sid=0):

        super(StationHeader,self).__init__(name=name,lat_yc=lat_yc,lon_xc=lon_xc)

        if 32 < len(name):
            raise Exception('Station.name cannot exceed 32 characters')

        self['network']   = network
        self['elevation'] = elevation
        self['burial']    = burial
        self['trid']      = trid
        self['gid']       = gid
        self['sid']       = sid
        self['data_fqdn'] = None


        self.comp_val     = self.trid



    def hash_val(self):
        #print('calling_hash')
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
    def sid(self):
        return self['sid']

    @sid.setter
    def sid(self, value):
        self['sid'] = value


    @property
    def data_fqdn(self):
        return self['data_fqdn']

    @data_fqdn.setter
    def data_fqdn(self, value):
        self['data_fqdn'] = value


class SolutionHeader(SpecHeader):
    
    def __init__(self,
                 header=None,
                 tshift=None,
                 lat_yc=None,
                 lon_xc=None,
                 depth=None,
                 sid=0,
                 gid=0):

        super(SolutionHeader,self).__init__(name=header,lat_yc=lat_yc,lon_xc=lon_xc)

        self['tshift']    = tshift
        self['depth']     = depth
        self['sid']       = sid
        self['gid']       = gid
        
        self.comp_val     = sid
        
        
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
    def sid(self):
        return self['sid']

    @sid.setter
    def sid(self, value):
        self['sid'] = value


    @property
    def gid(self):
        return self['gid']

    @gid.setter
    def gid(self, value):
        self['gid'] = value
                     
        
class ForceSolutionHeader(SolutionHeader):
    """
    FORCE  001 
    time shift:       0.0 
    f0:               0.0 
    latorUTM:         591000.0
    longorUTM:        246000.0
    depth:            -100
    """

    def __init__(self, 
             header=None,
             tshift=None,
             f0=None,
             lat_yc=None,
             lon_xc=None,
             depth=None,
             sid=0,
             gid=0):

        super(ForceSolutionHeader,self).__init__(header=header,
                                                 tshift=tshift,
                                                 lat_yc=lat_yc,
                                                 lon_xc=lon_xc,
                                                 depth=depth,
                                                 sid=sid,
                                                 gid=gid)


        self['f0'] = f0
        
        
    def hash_val(self):
        return np.sqrt(super().hash_val()**2 + self.f0**2)
        
    def eq_comparator(self):
        return tuple([*super().eq_comparator(), self.f0])
        

    @property
    def f0(self):
        return self['f0']

    @f0.setter
    def f0(self, value):
        self['f0'] = value
        


        
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
                 header=None,
                 ename=None,
                 tshift=None,
                 hdur=None,
                 lat_yc=None,
                 lon_xc=None,
                 depth=None,
                 mt=None,
                 sid=0,
                 gid=0):
                 
        if len(mt) != 6: 
            raise Exception('mt must by array-like and have only 6 values (upper-triangle)')

        super(CMTSolutionHeader,self).__init__(header=header,
                                               tshift=tshift,
                                               lat_yc=lat_yc,
                                               lon_xc=lon_xc,
                                               depth=depth,
                                               sid=sid,
                                               gid=gid)

            
        self['ename'] = ename
        self['hdur']  = hdur
        self['mrr']   = mt[0]
        self['mtt']   = mt[1]
        self['mpp']   = mt[2]
        self['mrt']   = mt[3]
        self['mrp']   = mt[4]
        self['mtp']   = mt[5]
        
        
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



def make_station_half_cross_members(station=None,delta=None):
    
    if delta == 0:
        raise ValueError('delta cannot equal zero')
    
    station_members = []

    l_gid = np.array([1,2,3])
    if delta < 0:
        l_grid += 3
    
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
    station_members.append(make_station_half_cross_members(station,delta))

    # negative incremented coordinate group
    station_members.append(make_station_half_cross_members(station,-delta))
    
    return station_group


def make_station_cross_group(station=None,delta=None):

    
    station_group = []

    # add "central" member (no coordinate change)
    cpy_station = copy.deepcopy(station)
    cpy_station.gid = station.gid = 0
    station_group.append(cpy_station) 
    
    # add pos/neg coordinate members
    station_group.append(make_station_cross_members(station,delta))

    return station_group


def make_station_half_cross_group(station=None,delta=None):

    station_group = []

    # add "central" member (no coordinate change)
    cpy_station = copy.deepcopy(station)
    cpy_station.gid = station.gid = 0
    station_group.append(cpy_station) 
    
    # add pos/neg coordinate members
    station_group.append(make_station_half_cross_members(station,delta))

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
