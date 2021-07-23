import numpy as np

################################################################################
#
# General Header Class for SPECFEM3D files (STATION and SOLUTIONS)
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
# Headers for SPECFEM3D files (STATION and SOLUTIONS)
#
################################################################################

class CoordHeader(Header):

    def __init__(self,name=None,lat_yc=None,lon_xc=None,depth=None):
        super(CoordHeader,self).__init__(self,name=name)

        self['lat_yc']    = lat_yc
        self['lon_xc']    = lon_xc
        self['depth']     = depth


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
        super(SpecHeader,self).__init__(self,name=name,lat_yc=lat_yc,lon_xc=lon_xc,depth=depth)

        self['proj_id']   = proj_id
        self['eid']       = eid
        self['sid']       = sid
        self['comp_val']  = None


    def hash_val(self):
        sup_hv = super().hash_val()
        return np.sqrt(self.proj_id**2 + self.eid**2 + self.sid**2 + sup_hv**2)
        #return np.sqrt(self.proj_id**2 + self.eid**2 + self.sid**2 + self.lon_xc**2 + self.lat_yc**2 + self.depth**2)

    # NOT NEEDED
    #def eq_comparator(self):


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
    component dir vect source Z_UP:0
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
             comp_src_EX=None,
             comp_src_NY=None,
             comp_src_Zup=None,
             proj_id=0,
             eid=0,
             sid=0):


        hstr  = f'FORCE ' + str(sid).zfill(3)
        hstr += f' {date.year} {date.month} {date.day} {date.hour} {date.minute} {date.second}'
        hstr += f' {lat_yc} {lon_xc} {depth/1000} srcid_{eid}'

        super(ForceSolutionHeader,self).__init__(name=hstr,
                                                 lat_yc=lat_yc,
                                                 lon_xc=lon_xc,
                                                 depth=depth,
                                                 tshift=tshift,
                                                 date=date,
                                                 ename=ename
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
        hstr += f' {lat_yc} {lon_xc} {depth/1000.0} {mt.Mw} 0 srcid_{eid}'

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
