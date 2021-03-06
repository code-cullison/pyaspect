import copy
import importlib

import numpy as np
import pandas as pd

from pyaspect.moment_tensor import MomentTensor

################################################################################
#
# General Header Class 
#
################################################################################

class Header(dict):

    def __init__(self,name=None):
        super(Header,self).__init__()

        self['name']      = name
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


    def copy(self):
        return copy.deepcopy(self)


    @property
    def name(self):
        return self['name']

    @name.setter
    def name(self, value):
        self['name'] = value


    @property
    def comp_val(self):
        return self['comp_val']

    @comp_val.setter
    def comp_val(self, value):
        self['comp_val'] = value


################################################################################
#
# General Coordinate Headers for SPECFEM3D files (STATION and SOLUTIONS)
#
################################################################################

class CoordHeader(Header):

    def __init__(self,name=None,lat_yc=None,lon_xc=None,depth=None):
        super(CoordHeader,self).__init__(name=name)

        self['lat_yc']    = lat_yc
        self['lon_xc']    = lon_xc
        self['depth']     = depth
        self['comp_val']  = self.depth


    def hash_val(self):
        return np.sqrt(self.lon_xc**2 + self.lat_yc**2 + self.depth**2)

    def comparator(self):
        return self.comp_val

    def eq_comparator(self):
        return (self.lat_yc, self.lon_xc, self.depth)


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
    def depth(self):
        return self['depth']

    @depth.setter
    def depth(self, value):
        self['depth'] = value



################################################################################
#
# Headers for SPECFEM3D files (STATION and SOLUTIONS)
#
################################################################################


class SpecHeader(CoordHeader):

    def __init__(self,name=None,lat_yc=None,lon_xc=None,depth=None,proj_id=0,eid=0,sid=0):
        super(SpecHeader,self).__init__(name=name,lat_yc=lat_yc,lon_xc=lon_xc,depth=depth)

        self['proj_id']   = proj_id
        self['eid']       = eid
        self['sid']       = sid


    def hash_val(self):
        sup_hv = super().hash_val()
        return np.sqrt(self.proj_id**2 + self.eid**2 + self.sid**2 + sup_hv**2)
        #return np.sqrt(self.proj_id**2 + self.eid**2 + self.sid**2 + self.lon_xc**2 + self.lat_yc**2 + self.depth**2)


    def eq_comparator(self):
        return tuple(list(super().eq_comparator()) + [self.proj_id,self.eid,self.sid])


    @property
    def proj_id(self):
        return self['proj_id']

    @proj_id.setter
    def proj_id(self, value):
        self['proj_id'] = value


    @property
    def proj_id(self):
        return self['proj_id']

    @proj_id.setter
    def proj_id(self, value):
        self['proj_id'] = value


    @property
    def sid(self):
        return self['sid']

    @sid.setter
    def sid(self, value):
        self['sid'] = value



################################################################################
#
# STATION Classes: STATION[_ADJOINT] files and headers for SPECFEM3D Cart.
#
################################################################################


class StationHeader(SpecHeader):

    def __init__(self,
                 name=None,
                 lat_yc=None,
                 lon_xc=None,
                 depth=None,
                 elevation=None,
                 network=None,
                 proj_id=0,
                 eid=0,
                 sid=0,
                 trid=0,
                 gid=0):

        super(StationHeader,self).__init__(name=name,lat_yc=lat_yc,lon_xc=lon_xc,depth=depth,proj_id=proj_id,eid=eid,sid=sid)

        if 32 < len(name):
            raise Exception('Station.name cannot exceed 32 characters')

        self['network']   = network
        self['elevation'] = elevation
        self['trid']      = trid       # trace id
        self['gid']       = gid        # trace group id (trace decomp)
        self['data_fqdn'] = None


        self.comp_val     = self.trid



    def hash_val(self):
        sup_hv = super().hash_val()
        return np.sqrt(self.elevation**2 + self.trid**2 + self.gid**2 + sup_hv**2)

    def eq_comparator(self):
        return tuple(list(super().eq_comparator()) + [self.elevation])


    @classmethod
    def from_dict(cls,h_dict):
        '''
        Alternative constructor StationHeader

        This constructor takes a :py:dict: object, 

        The ``h_dict`` argument must have 'key: value' pairs for:

        * ``name``
        * ``lat_yc``
        * ``lon_xc``
        * ``depth``
        * ``elevation``
        * ``network``

        If the following 'key: value' pairs are not specfied, then 
        they will be set equal to 0:

        * ``proj_id``
        * ``eid``
        * ``sid``
        * ``trid``
        * ``gid``

        All other 'key: value' pairs will be added

        '''

        if not isinstance(h_dict,dict):
            raise TypeError(f'arg: \'h_dict\' must be of type dict')

        required_keys = ['name','lat_yc','lon_xc','depth','elevation','network']
        
        # check and get required keys
        args_dict = {}
        for rkey in required_keys:
            if rkey not in h_dict:
                raise Exception(f'missing key: \'{rkey}\'')
            else:
                args_dict[rkey] = h_dict[rkey]
                del h_dict[rkey]

        # make StationHeader
        new_header = StationHeader(**args_dict)

        # insert additional key: value pairs
        for key in h_dict:
            new_header[key] = h_dict[key]

        return new_header


    @classmethod
    def from_series(cls,h_series):
        '''
        Alternative constructor StationHeader

        This constructor takes a :py:dict: object, 

        The ``h_series`` argument must have 'column names' for:

        * ``name``
        * ``lat_yc``
        * ``lon_xc``
        * ``depth``
        * ``elevation``
        * ``network``

        If the following 'column names' are not specfied, then 
        their values will be set equal to 0:

        * ``proj_id``
        * ``eid``
        * ``sid``
        * ``trid``
        * ``gid``

        All other 'column names' will also be added

        '''

        if not isinstance(h_series,pd.Series):
            raise TypeError(f'arg: \'h_series\' must be of type panda.Series')

        return StationHeader.from_dict(h_series.to_dict())


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
                 lat_yc=None,
                 lon_xc=None,
                 depth=None,
                 tshift=None,
                 date=None,
                 ename=None,
                 proj_id=0,
                 eid=0,
                 sid=0):

        super(SolutionHeader,self).__init__(name=name,lat_yc=lat_yc,lon_xc=lon_xc,depth=depth,proj_id=proj_id,eid=eid,sid=sid)

        self['tshift'] = tshift
        self['date']   = date
        self['ename']  = ename

        self.comp_val     = eid


    def hash_val(self):
        sup_hv = super().hash_val()
        return np.sqrt(self.tshift**2 + sup_hv**2)

    def eq_comparator(self):
        return tuple(list(super().eq_comparator()) + [self.tshift])


    @property
    def tshift(self):
        return self['tshift']

    @tshift.setter
    def tshift(self, value):
        self['tshift'] = value


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
    component dir vect source Z_UP:   0
    """

    def __init__(self,
             ename=None,
             lat_yc=None,
             lon_xc=None,
             depth=None,
             tshift=None,
             date=None,
             f0=None,
             factor_fs=None,
             comp_src_EX=None, #FIXME: maybe comp_force_EX?
             comp_src_NY=None,
             comp_src_Zup=None,
             proj_id=0,
             eid=0,
             sid=0):

        #FIXME: make header name (kind of ugley way to do this)
        hstr = ForceSolutionHeader.create_header_name(eid=eid,
                                                      sid=sid,
                                                      date=date,
                                                      lat_yc=lat_yc,
                                                      lon_xc=lon_xc,
                                                      depth=depth)

        super(ForceSolutionHeader,self).__init__(name=hstr,
                                                 lat_yc=lat_yc,
                                                 lon_xc=lon_xc,
                                                 depth=depth,
                                                 tshift=tshift,
                                                 date=date,
                                                 ename=ename,
                                                 proj_id=proj_id,
                                                 eid=eid,
                                                 sid=sid)


        self['f0']           = f0
        self['factor_fs']    = factor_fs
        self['comp_src_EX']  = comp_src_EX
        self['comp_src_NY']  = comp_src_NY
        self['comp_src_Zup'] = comp_src_Zup


    def hash_val(self):
        return np.sqrt(super().hash_val()**2 + self.f0**2)

    def eq_comparator(self):
        return tuple([*super().eq_comparator(), self.f0])


    @staticmethod
    def create_header_name(eid=None,sid=None,date=None,lat_yc=None,lon_xc=None,depth=None):
        hstr  = f'FORCE ' + str(sid).zfill(3)
        hstr += f' {date.year} {date.month} {date.day} {date.hour} {date.minute} {date.second}'
        hstr += f' {lat_yc} {lon_xc} {depth/1000} srcid_{eid}'
        return hstr


    @classmethod
    def from_dict(cls,h_dict):
        '''
        Alternative constructor ForceSolutionHeader

        This constructor takes a :py:dict: object, 

        The ``h_dict`` argument must have 'key: value' pairs for:

        * ``ename``
        * ``lat_yc``
        * ``lon_xc``
        * ``depth``
        * ``tshift``
        * ``date``
        * ``f0``
        * ``factor_fs``
        * ``comp_src_EX``
        * ``comp_src_NY``
        * ``comp_src_Zup``

        If the following 'key: value' pairs are not specfied, then 
        they will be set equal to 0:

        * ``proj_id``
        * ``eid``
        * ``sid``

        All other 'key: value' pairs will be added

        '''
        if not isinstance(h_dict,dict):
            raise TypeError(f'arg: \'h_dict\' must be of type dict')

        required_keys = ['ename',
                         'lat_yc','lon_xc','depth',
                         'tshift',
                         'date',
                         'f0',
                         'factor_fs',
                         'comp_src_EX','comp_src_NY','comp_src_Zup']
        
        # check and get required keys
        args_dict = {}
        for rkey in required_keys:
            if rkey not in h_dict:
                raise Exception(f'missing key: \'{rkey}\'')
            else:
                args_dict[rkey] = h_dict[rkey]
                del h_dict[rkey]

        # make ForceSolutionHeader
        new_header = ForceSolutionHeader(**args_dict)

        # insert additional key: value pairs
        for key in h_dict:
            new_header[key] = h_dict[key]

        return new_header


    @classmethod
    def from_series(cls,h_series):
        '''
        Alternative constructor ForceSolutionHeader

        This constructor takes a :py:pandas:Series: object, 

        The ``h_series`` argument must have 'column names' for:

        * ``ename``
        * ``lat_yc``
        * ``lon_xc``
        * ``depth``
        * ``tshift``
        * ``date``
        * ``f0``
        * ``factor_fs``
        * ``comp_src_EX``
        * ``comp_src_NY``
        * ``comp_src_Zup``

        If the following 'column names' are not specfied, then 
        their values will be set equal to 0:

        * ``proj_id``
        * ``eid``
        * ``sid``

        All other 'column names' will be added

        '''
        if not isinstance(h_series,pd.Series):
            raise TypeError(f'arg: \'h_series\' must be of type pandas.Series')

        return ForceSolutionHeader.from_dict(h_series.to_dict())


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
                 ename=None,
                 lat_yc=None,
                 lon_xc=None,
                 depth=None,
                 tshift=None,
                 date=None,
                 hdur=None,
                 mt=None,
                 proj_id=0,
                 eid=0,
                 sid=0):

        from pyaspect.moment_tensor import MomentTensor

        if not isinstance(mt,MomentTensor):
            raise TypeError('mt must be of type pyaspect.moment_tensor.MomentTensor')


        hstr  = f'PDE {date.year} {date.month} {date.day} {date.hour} {date.minute} {date.second}'
        hstr += f' {lat_yc} {lon_xc} {depth/1000.0} {mt.magnitude} 0 srcid_{eid}'

        super(CMTSolutionHeader,self).__init__(name=hstr,
                                               lat_yc=lat_yc,
                                               lon_xc=lon_xc,
                                               depth=depth,
                                               tshift=tshift,
                                               date=date,
                                               ename=ename,
                                               proj_id=proj_id,
                                               eid=eid,
                                               sid=sid)


        self['hdur']   = hdur
        self['strike'] = mt.strike
        self['dip']    = mt.dip
        self['rake']   = mt.rake
        self['mw']     = mt.magnitude
        self['mrr']    = mt.m6_up_south_east()[0]
        self['mtt']    = mt.m6_up_south_east()[1]
        self['mpp']    = mt.m6_up_south_east()[2]
        self['mrt']    = mt.m6_up_south_east()[3]
        self['mrp']    = mt.m6_up_south_east()[4]
        self['mtp']    = mt.m6_up_south_east()[5]
        '''
        self['mrr']    = mt.harvard_dcm_m6()[0]
        self['mtt']    = mt.harvard_dcm_m6()[1]
        self['mpp']    = mt.harvard_dcm_m6()[2]
        self['mrt']    = mt.harvard_dcm_m6()[3]
        self['mrp']    = mt.harvard_dcm_m6()[4]
        self['mtp']    = mt.harvard_dcm_m6()[5]
        '''


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


    @classmethod
    def from_dict(cls,h_dict):
        '''
        Alternative constructor CMTSolutionHeader

        This constructor takes a :py:dict: object, 

        The ``h_dict`` argument must have 'key: value' pairs for:

        * ``ename``
        * ``lat_yc``
        * ``lon_xc``
        * ``depth``
        * ``tshift``
        * ``date``
        * ``hdur``
        * ``mrr``
        * ``mtt``
        * ``mpp``
        * ``mrt``
        * ``mrp``
        * ``mtp``

        If the following 'key: value' pairs are not specfied, then 
        they will be set equal to 0:

        * ``proj_id``
        * ``eid``
        * ``sid``

        All other 'key: value' pairs will be added

        '''
        if not isinstance(h_dict,dict):
            raise TypeError(f'arg: \'h_dict\' must be of type dict')

        required_keys = ['ename',
                         'lat_yc','lon_xc','depth',
                         'tshift',
                         'date',
                         'hdur']
        
        required_mt_keys = ['mrr',
                            'mtt',
                            'mpp',
                            'mrt',
                            'mrp',
                            'mtp']
        
        # check and get required keys
        args_dict = {}
        for rkey in required_keys:
            if rkey not in h_dict:
                raise Exception(f'missing key: \'{rkey}\'')
            else:
                args_dict[rkey] = h_dict[rkey]
                del h_dict[rkey]

        # check and get required moment-tensor keys
        mt_args_dict = {}
        for rkey in required_mt_keys:
            if rkey not in h_dict:
                raise Exception(f'missing mt-key: \'{rkey}\'')
            else:
                mt_args_dict[rkey] = h_dict[rkey]
                del h_dict[rkey]

        # add moment tensor
        '''
        args_dict['mt'] = MomentTensor(mw=mt_args_dict['mw'],
                                       strike=mt_args_dict['strike'],
                                       dip=mt_args_dict['dip'],
                                       rake=mt_args_dict['rake'])
        '''

        mrr = mt_args_dict['mrr']
        mtt = mt_args_dict['mtt']
        mpp = mt_args_dict['mpp']
        mrt = mt_args_dict['mrt']
        mrp = mt_args_dict['mrp']
        mtp = mt_args_dict['mtp']
        h_matrix = np.array([[mrr,mrt,mrp],[mrt,mtt,mtp],[mrp,mtp,mpp]])
        args_dict['mt'] = MomentTensor(m_up_south_east=h_matrix)

        # make CMTSolutionHeader
        new_header = CMTSolutionHeader(**args_dict)

        # insert additional key: value pairs
        for key in h_dict:
            new_header[key] = h_dict[key]

        return new_header
    

    @classmethod
    def from_series(cls,h_series):
        '''
        Alternative constructor CMTSolutionHeader

        This constructor takes a :py:pandas:Series: object, 

        The ``h_series`` argument must have 'column names' for:

        * ``ename``
        * ``lat_yc``
        * ``lon_xc``
        * ``depth``
        * ``tshift``
        * ``date``
        * ``hdur``
        * ``strike``
        * ``dip``
        * ``rake``
        * ``mw``

        If the following 'column names' are not specfied, then 
        their values will be set equal to 0:

        * ``proj_id``
        * ``eid``
        * ``sid``

        All other 'column names' will be added

        '''
        if not isinstance(h_series,pd.Series):
            raise TypeError(f'arg: \'h_series\' must be of type pandas.Series')

        return CMTSolutionHeader.from_dict(h_series.to_dict())


    @property
    def date(self):
        return self['date']

    @date.setter
    def date(self, value):
        self['date'] = value


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


    #no setter for Mw but there is for mw
    @property
    def Mw(self):
        return self['mw']

    @property
    def mw(self):
        return self['mw']

    @mw.setter
    def mw(self, value):
        self['mw'] = value


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



################################################################################
#
# Record Header for pairing SPECFEM SOLUTIONS with STATIONS
#
################################################################################

#TODO this is the actual record. I need the header, then make record.py
class RecordHeader(Header):

    def __init__(self,name=None,solutions_h=None,stations_h=None,proj_id=0,rid=0,iter_id=0,is_reciprocal=False):
        super(RecordHeader,self).__init__(name=name)

        l_solutions_h = solutions_h
        if not isinstance(solutions_h,list):
            l_solutions_h = [solutions_h]

        check_all = False
        check_all = all(isinstance(s,SolutionHeader) for s in l_solutions_h)
        if not check_all:
            raise Exception('elements in arg: \'solutions_h\' must be of type SolutionHeader')

        check_s = l_solutions_h[0]
        check_all = all(type(s) == type(check_s) for s in l_solutions_h)
        if not check_all:
            raise Exception('all elements in arg: \'solutions_h\' must be of type {type(check_s)}')

        l_stations_h = stations_h
        if not isinstance(stations_h,list):
            l_stations_h = [stations_h]
        
        check_all = all(isinstance(s,StationHeader) for s in l_stations_h)
        if not check_all:
            raise Exception('elements in arg: \'stations_h\' must be of type StationHeader')


        self.comp_val   = rid

        self['proj_id'] = proj_id
        self['rid']     = rid
        self['iter_id'] = iter_id
        self['is_reciprocal'] = is_reciprocal

        #FIXME: see below. is this a good/safe trick?
        self._solution_mod_name = l_solutions_h[0].__module__
        self._solution_cls_name = l_solutions_h[0].__class__.__name__

        self._station_mod_name = l_stations_h[0].__module__
        self._station_cls_name = l_stations_h[0].__class__.__name__


        self['default_solu_midx'] = ('eid','sid')
        self['default_stat_midx'] = ('eid','sid','trid','gid')

        smidx = list(self['default_solu_midx'])
        rmidx = list(self['default_stat_midx'])

        self['solutions_df'] = pd.DataFrame.from_records(l_solutions_h, index=smidx)
        self['stations_df']  = pd.DataFrame.from_records(l_stations_h,  index=rmidx)

        l_nsrc_idx = []
        for idx, df in self.solutions_df.groupby(level='eid'):
            nsids = df.index.get_level_values('sid').unique()
            l_nsrc_idx.append(nsids)

        l_nrec_idx = []
        for idx, df in self.stations_df.groupby(level='eid'):
            nsids = df.index.get_level_values('sid').unique()
            l_nrec_idx.append(nsids)

        if len(l_nsrc_idx) != len(l_nrec_idx):
            raise Exception('Number of events does not match between Solutions and Stations')

        self['nevents'] = len(l_nsrc_idx)
        
        #use different num_sids so that they can be checked
        ns_solu = None
        ns_stat = None
        prev_solu = len(l_nsrc_idx[0])
        for ie in range(len(l_nsrc_idx)):
            ns_solu = len(l_nsrc_idx[ie])
            ns_stat = len(l_nrec_idx[ie])
            if ns_solu !=  ns_stat:
                raise Exception(f'For event-{ie}: number of sid\'s differs between Solutions and Stations')
            #this check enforces that "batch" src size is the same per event
            if ns_solu != prev_solu:
                raise Exception(f'Number sid\'s differs between events')
            for isrc in range(len(l_nsrc_idx[ie])):
                if l_nsrc_idx[ie][isrc] !=  l_nrec_idx[ie][isrc]:
                    raise Exception('src-id mismatch between Solutions and Stations')

        if not isinstance(ns_solu,int):
            raise Exception('class field "ns_solu" must be an int type')

        self['nsrc'] = ns_solu


        self['added_solution_header_words'] = []
        self['added_station_header_words']  = []


    def hash_val(self):
        return np.sqrt(self.proj_id**2 + self.rid**2 + self.iter_id**2)

    def comparator(self):
        return self.comp_val

    def eq_comparator(self,other):
        are_solutions_eq = self.solutions_df.equals(other.solutions_df)
        are_stations_eq  = self.stations_df.equals(other.stations_df)
        are_ids_eq       = (self.proj_id == other.proj_id and 
                            self.rid == other.rid and 
                            self.iter_id == other.iter_id)

        return are_solutions_eq and are_stations_eq and are_ids_eq

    def __eq__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        req = self.eq_comparator(other)
        return (self.eq_comparator(other))

    def __ne__(self, other):
        if not self._is_valid_operand(other):
            return NotImplemented
        return (not self.eq_comparator(other))


    def __str__(self):
        out_str  = f'Solution Header(s):\n{self.solutions_df}\n\n'
        out_str += f'Station Header(s):\n {self.stations_df}'
        return out_str

    def __repr__(self):
        out_str  = f'Solution Header(s):\n{self.solutions_df.__repr__()}\n\n'
        out_str += f'Station Header(s):\n {self.stations_df.__repr__()}'
        return out_str


    # Helper function turns a sequence of boolean indices into
    # a single boolean index for a dataframe
    def _slice_to_bool(self,l_slice):
        if not isinstance(l_slice,list):
            raise Exception('arugment must be a list type')
        if len(l_slice) == 0:
            raise Exception('arugment must be a list of at least len=1')
        if len(l_slice) == 1:
            return l_slice[0]
        else:
            ibool = l_slice[0]
            for i in range(1,len(l_slice)):
                ibool = ibool & l_slice[i]
            return ibool


    # Helper function needed for make dataframe slices that return
    # another dataframe instead of a Series even if only one row is
    # the result of the slice
    def _convert_slice(self,kslice,smax):
        start = kslice.start
        stop = kslice.stop
        step = kslice.step
        if kslice.start == None:
            start = 0
        if kslice.stop  == None:
            stop = smax
        if kslice.step == None:
            step = 1
        return slice(start,stop,step)

    
    # Helper function for slicing the dataframes
    # Assumes that solutions_df has 2D multiindex: order=[eid,sid]
    # Assumes that stations_df has 4D multiindex: order=[eid,sid,trid,gid]
    def _get_df_slice_index(self,kslice,df,is_stations=False):

        _df = df.copy()
        nevents = _df.index.get_level_values('eid').nunique()
        nsrcs   = _df.index.get_level_values('sid').nunique()
        ntrs    = None
        ngids   = None
        if is_stations:
            ntrs    = _df.index.get_level_values('trid').nunique()
            ngids   = _df.index.get_level_values('gid').nunique()

        _df.reset_index(inplace=True)

        l_slice = []
        if isinstance(kslice,tuple):
            if isinstance(kslice[0],int):
                if kslice[0] < 0:
                    raise KeyError(kslice[0])
                l_slice.append((_df['eid'] == kslice[0]))
            else:
                s = self._convert_slice(kslice[0],nevents)
                l_slice.append( ((_df['eid'] >= s.start) &
                                 (_df['eid'] < s.stop)   &
                                 (_df['eid']%s.step == 0)) )
            if isinstance(kslice[1],int):
                if kslice[1] < 0:
                    raise KeyError(kslice[1])
                l_slice.append((_df['sid'] == kslice[1]))
            else:
                s = self._convert_slice(kslice[1],nsrcs)
                l_slice.append( ((_df['sid'] >= s.start) &
                                 (_df['sid'] < s.stop)   &
                                 (_df['sid']%s.step == 0)) )
            if 3 <= len(kslice):
                if isinstance(kslice[2],int):
                    if kslice[2] < 0:
                        raise KeyError(kslice[2])
                    l_slice.append((_df['trid'] == kslice[2]))
                else:
                    s = self._convert_slice(kslice[2],ntrs)
                    l_slice.append( ((_df['trid'] >= s.start) &
                                     (_df['trid'] < s.stop)   &
                                     (_df['trid']%s.step == 0)) )
            if 4 == len(kslice):
                if isinstance(kslice[3],int):
                    if kslice[3] < 0:
                        raise KeyError(kslice[3])
                    l_slice.append((_df['gid'] == kslice[3]))
                else:
                    s = self._convert_slice(kslice[3],ngids)
                    l_slice.append( ((_df['gid'] >= s.start) &
                                     (_df['gid'] < s.stop)   &
                                     (_df['gid']%s.step == 0)) )
            if 4 < len(kslice):
                raise Exception('too many indexers')
                
        elif isinstance(kslice,int):
            if kslice < 0:
                raise KeyError(kslice)
            l_slice.append((_df['eid'] == kslice))
        elif isinstance(kslice,slice):
            s = self._convert_slice(kslice,nevents)
            l_slice.append(((_df['eid'] >= s.start) &
                     (_df['eid'] < s.stop)   &
                     (_df['eid']%s.step == 0)))
        else:
            raise TypeError('wrong indexer type')
        
        return self._slice_to_bool(l_slice)


    # helper funciton for slicing the solutions_df
    def _get_solutions_df_slice_index(self,kslice):

        if isinstance(kslice,tuple):
            if len(kslice) < 2 or 4 < len(kslice):
                raise Exception('too many indexers for solutions_df')
            else:
                return self._get_df_slice_index(kslice[0:2],self.solutions_df,is_stations=False)
        else:
            return self._get_df_slice_index(kslice,self.solutions_df,is_stations=False)


    # helper funciton for slicing the stations_df
    def _get_stations_df_slice_index(self,kslice):

        if isinstance(kslice,tuple):
            if len(kslice) < 2 or 4 < len(kslice):
                raise Exception('incorrect number of indexers for stations_df indexer')
            else:
                return self._get_df_slice_index(kslice,self.stations_df,is_stations=True)
        else:
            return self._get_df_slice_index(kslice,self.stations_df,is_stations=True)
            

    # This will act like a dictionary if kslice is a string.
    # Otherwise, it will act like slicing and return a NEW
    # "sliced" RecordHeader
    def __getitem__(self, kslice):

        if isinstance(kslice, str):
            return super(RecordHeader, self).__getitem__(kslice)

        # get sliced dataframes
        solu_slice_idx = self._get_solutions_df_slice_index(kslice)
        stat_slice_idx = self._get_stations_df_slice_index(kslice)

        c_solu_df = self.solutions_df.reset_index()[solu_slice_idx]
        c_stat_df = self.stations_df.reset_index()[stat_slice_idx]

        # make list of solution headers
        #HeaderCls = self._get_header_class(is_stations=False)
        HeaderCls = self.solution_cls
        '''
        print(f'IN RecordHeader: check c_solution_df type:{type(c_solu_df)}')
        print(f'IN RecordHeader: _solution_mod_name:{self._solution_mod_name}')
        print(f'IN RecordHeader: _solution_cls_name:{self._solution_cls_name}')
        print(f'IN RecordHeader: SolHeaderCls:{HeaderCls}')
        print(f'IN RecordHeader: isinstance(SolHeaderCls):{isinstance(HeaderCls,CMTSolutionHeader)}')
        if self._solution_cls_name == 'CMTSolutionHeader':
            print('***IN -- IS-Type')
            for index, row in c_solu_df.iterrows():
                print(f'index,row["mw"]:\n{index}\n{row["mw"]}')
        '''
        slice_sol_h = [HeaderCls.from_series(row) for index, row in c_solu_df.iterrows()]

        # make list of station headers
        #HeaderCls = self._get_header_class(is_stations=True)
        HeaderCls = self.station_cls
        slice_stat_h = [HeaderCls.from_series(row) for index, row in c_stat_df.iterrows()]

        #make new record which is like it had been sliced
        return RecordHeader(name=self.name,
                            solutions_h=slice_sol_h,
                            stations_h=slice_stat_h,
                            proj_id=self.proj_id,
                            rid=self.rid,
                            iter_id=self.iter_id)


    # Helper function returns a "default-index" dataframe
    # i.e. flattens the dataframe index -> multiindex-cols ->  data-cols
    def _get_reset_df(self,key=None,value=None,is_stations=True):
        
        c_df = None

        if is_stations:
            c_df = self.stations_df.copy()
        else:
            c_df = self.solutions_df.copy()

        c_df.reset_index(inplace=True)

        if key != None:
            if key not in c_df.columns:
                raise KeyError('key: {key} is not a column in stations_df')
            #FIXME: find a way to check dtype of value and those in columns
            c_df = c_df.loc[c_df[key] == value]

        return c_df


    def _get_header_class(self,is_stations=True):

        HeaderCls = None

        if is_stations:
            HeaderCls = getattr(importlib.import_module(self._station_mod_name), self._station_cls_name)
        else:
            HeaderCls = getattr(importlib.import_module(self._solution_mod_name), self._solution_cls_name)

        return HeaderCls


    def _get_list_from_df(self, key=None, value=None, is_stations=True):

        c_df = self._get_reset_df(key=key,value=value,is_stations=is_stations)

        #FIXME: Q. Is this a good trick or an ugly trick?
        #       A. It's a trick.
        HeaderCls = self._get_header_class(is_stations=is_stations)

        header_list = [HeaderCls.from_series(row) for index, row in c_df.iterrows()]

        del c_df

        return header_list


    def get_solutions_header_list(self, key=None, value=None):
        return self._get_list_from_df(key=key,value=value,is_stations=False)

    def get_stations_header_list(self, key=None, value=None):
        return self._get_list_from_df(key=key,value=value,is_stations=True)


    def _add_header_word(self, key, h_values, is_stations=True):

        if not isinstance(key,str):
            raise ValueError('arg: \'key\' must be a str type')
        if not isinstance(h_values,list):
            raise ValueError('arg: \'h_values\' must be a list type')

        if is_stations:
            if len(h_values) != len(self.stations_df.index):
                raise Exception('len(\'h_values\') must equal number of stations')
            self.added_station_header_words.append(key)
            self.stations_df[key] = h_values
        else:
            if len(h_values) != len(self.solutions_df.index):
                raise Exception('len(\'h_values\') must equal number of solutions')
            self.added_solution_header_words.append(key)
            self.solutions_df[key] = h_values

    def add_station_header_word(self, key, h_values):
        self._add_header_word(key=key,h_values=h_values,is_stations=True)

    def add_solution_header_word(self, key, h_values):
        self._add_header_word(key=key,h_values=h_values,is_stations=False)


    def get_event_nsolutions(self,ievent):
        idx = pd.IndexSlice
        ns_solu = self.solutions_df.loc[idx[:,ievent,:],:].index.get_level_values('sid').nunique()
        ns_stat = self.stations_df.loc[idx[:,ievent,:,:,:],:].index.get_level_values('sid').nunique()
        if ns_stat != ns_solu:
            raise Exception(f'For event-{ievent}: number of sid\'s differs between Solutions and Stations')

        return ns_stat

    def reset_midx(self):
        self.reset_solutions_midx()
        self.reset_stations_midx()

    def reset_solutions_midx(self):
        self.solutions_df.reset_index(inplace=True)

    def reset_stations_midx(self):
        self.stations_df.reset_index(inplace=True)

    def set_default_midx(self):
        self.set_default_solutions_midx()
        self.set_default_stations_midx()

    def set_default_solutions_midx(self):
        smidx = list(self['default_solu_midx'])
        self.solutions_df.set_index(smidx,inplace=True)

    def set_default_stations_midx(self):
        rmidx = list(self['default_stat_midx'])
        self.stations_df.set_index(rmidx,inplace=True)
        

    @property
    def proj_id(self):
        return self['proj_id']
    
    @property
    def rid(self):
        return self['rid']
    
    @property
    def iter_id(self):
        return self['iter_id']
    
    @property
    def is_reciprocal(self):
        return self['is_reciprocal']
    
    @is_reciprocal.setter
    def is_reciprocal(self, value):
        self['is_reciprocal'] = value
    
    @property
    def solution_cls(self):
        return self._get_header_class(is_stations=False)
    
    @property
    def station_cls(self):
        return self._get_header_class(is_stations=True)
    
    @property
    def solutions_df(self):
        return self['solutions_df']
    
    @property
    def stations_df(self):
        return self['stations_df']
    
    @property
    def added_solution_header_words(self):
        return self['added_solution_header_words']
    
    @property
    def added_station_header_words(self):
        return self['added_station_header_words']

    @property
    def nevents(self):
        return self['nevents']

    @property
    def nsrc(self):
        return self['nsrc']



