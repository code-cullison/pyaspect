import numpy as np


class Header(dict):

    def __init__(self):
        super(Header,self).__init__()


    def __getattr__(self, key):
        if name in self:
            return self[key]
        else:
            raise AttributeError(f'Header: {key} not found' )

    def __setattr__(self, key, value):
        self[key] = value

    def __delattr__(self, key):
        if key in self:
            del self[key]
        else:
            raise AttributeError(f'Header: {key} not found' )

    def __hash__(self):
        return hash(self.hash_val())

    def _is_valid_operand(self, other):
        return self.keys() == other.keys()

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


class StationHeader(Header):

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

        super(StationHeader,self).__init__()

        if 32 < len(name):
            raise Exception('Station.name cannot exceed 32 characters')

        self['name']      = name
        self['network']   = network
        self['lat_yc']    = lat_yc
        self['lon_xc']    = lon_xc
        self['elevation'] = elevation
        self['burial']    = burial
        self['trid']      = trid
        self['gid']       = gid
        self['sid']       = sid

        self['comp_val']  = self.trid

        self['ids'] = {}


    def hash_val(self):
        return np.sqrt(self.lon_xc**2 + self.lat_yc**2 + self.burial**2)

    def comparator(self):
        return self.comp_val

    def eq_comparator(self):
        return (self.lat_yc, self.lon_xc, self.elevation, self.burial)


    @property
    def comp_val(self):
        return self['comp_val']

    @comp_val.setter
    def comp_val(self, key, value):
        self[key] = value


    @property
    def name(self):
        return self['name']

    @name.setter
    def name(self, value):
        self['name'] = value


    @property
    def network(self):
        return self['network']

    @network.setter
    def network(self, value):
        self['network'] = value


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
        self['burial'] = value


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
    def ids(self):
        return self['ids']
