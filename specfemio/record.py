import os
import glob
import copy
import datetime
import importlib

import numpy as np
import pandas as pd

from scipy import signal

from pyaspect.specfemio.headers import RecordHeader
from pyaspect.specfemio.read import _read_headers
from pyaspect.specfemio.headers import CMTSolutionHeader
from pyaspect.specfemio.headers import StationHeader

#TODO this is the actual record. I need the header, then make record.py
class Record(RecordHeader):

    class _TraceData(object):
        
        class _XYZ(object):

            def __init__(self,ex,ny,z):

                self.ex = ex
                self.ny = ny
                self.z  = z

            def __str__(self):
                out_str  = f'Component E/X:\n{self.ex}\n\n'
                out_str += f'Component N/Y:\n{self.ny}\n\n'
                out_str += f'Component Z:\n{self.z}'
                return out_str

            def __repr__(self):
                out_str  = f'Component E/X:\n{self.ex.__repr__()}\n\n'
                out_str += f'Component N/Y:\n{self.ny.__repr__()}\n\n'
                out_str += f'Component Z:\n{self.z.__repr__()}'
                return out_str


        def __init__(self,data_df):
            
            self.df_x = data_df['comp_EX']
            self.df_y = data_df['comp_NY']
            self.df_z = data_df['comp_Z']


        def __getitem__(self,islice):
            return self._XYZ(self.df_x.loc[islice],self.df_y.loc[islice],self.df_z.loc[islice])

        def __str__(self):
            out_str  = f'Component E/X:\n{self.df_x}\n\n'
            out_str += f'Component N/Y:\n{self.df_y}\n\n'
            out_str += f'Component Z:\n{self.df_z}'
            return out_str

        def __repr__(self):
            out_str  = f'Component E/X:\n{self.df_x.__repr__()}\n\n'
            out_str += f'Component N/Y:\n{self.df_y.__repr__()}\n\n'
            out_str += f'Component Z:\n{self.df_z.__repr__()}'
            return out_str

    @classmethod
    def from_header_and_data(cls,rheader,data_df):
        '''
        Alternative constructor StationHeader

        This constructor takes two requried arguments, a RecordHeader
        and a DataFrame with the data and indexed the same as
        the solutions and stations dataframes in the Record header

        :pyaspect:specfemio:headers:RecordHeader
        :pandas:DataFrame

        * ``rheader``
        * ``data_df``
        '''

        if not isinstance(rheader,RecordHeader):
            raise TypeError(f'arg: \'rheader\' must be of type RecordHeader')

        if not isinstance(data_df,pd.DataFrame):
            raise TypeError(f'arg: \'data_df\' must be of type pandas.DataFrame')

        return Record('None',dtype='b',rheader=rheader,data_df=data_df)
        
            

    def __init__(self,header_fqp,dtype='b',rheader=None,data_df=None):

        if not isinstance(rheader,RecordHeader):
            rheader = _read_headers(header_fqp)
        else:
            if not isinstance(data_df,pd.DataFrame):
                raise ArgumentError('data_df must be provided')

        super(Record,self).__init__(name=rheader.name,
                                    solutions_h=rheader.get_solutions_header_list(),
                                    stations_h=rheader.get_stations_header_list(),
                                    proj_id=rheader.proj_id,
                                    rid=rheader.rid,
                                    iter_id=rheader.iter_id,
                                    is_reciprocal=rheader.is_reciprocal)

        self['header_path'] = os.path.dirname(header_fqp)
        
        if not isinstance(data_df,pd.DataFrame):
            self._load_data(dtype=dtype)
        else:
            self['data_df'] = data_df
            
        
        if not self._check_dfs_are_ok():
            raise Exception('pd.DataFrame indices do not match')
            

    def __str__(self):
        out_str  = f'{super(Record, self).__str__()}\n\n'
        out_str += f'Data:\n {self.data_df}'
        return out_str
    
    def __repr__(self):
        out_str  = f'{super(Record, self).__repr__()}\n\n'
        out_str += f'Data:\n {self.data_df.__repr__()}'
        return out_str

    def __getitem__(self, kslice):
        
        if not isinstance(kslice, str):
            dslice = super(Record, self)._get_df_slice_index(kslice,self.data_df,is_stations=True)
            c_data_df = self.data_df.reset_index()[dslice]
            c_rheader = super(Record, self).__getitem__(kslice)
            return Record(rheader=c_rheader,data_df=c_data_df)
        else:
            return super(Record, self).__getitem__(kslice)
    
    def _read_specfem_bin_trace(self,fpath,dtype=np.float32):
        return np.fromfile(fpath, dtype=dtype)

    def _load_data(self,dtype='b',sl=slice(None,None,None),scale=1.0,rfunc=None):

        if dtype != 'b' and _rfunc == None:
            raise Exception('can only read binary type data for the time being')
            
        #FIXME: add read ascii
        read_func = self._read_specfem_bin_trace
        if rfunc != None:
            read_func = rfunc

        l_data = []
        for eidx, edf in self.stations_df.groupby(level='eid'):
            for sidx, sdf in edf.groupby(level='sid'):
                for tidx, tdf in sdf.groupby(level='trid'):
                    for gidx, gdf in tdf.groupby(level='gid'):
                        fp_prefix = gdf.loc[(eidx,sidx,tidx,gidx),"data_fqdn"]
                        fp = os.path.join(self.header_path,fp_prefix)
                        match_fp = fp + '.*X[XYZEN].sem*'
                        data_dict = {'eid':eidx,'sid':sidx,'trid':tidx,'gid':gidx}
                        for filepath in glob.glob(match_fp):
                            comp = filepath.split('.')[-2][-1]
                            if comp == 'X' or comp == 'E':
                                data_dict['comp_EX'] = scale*read_func(filepath)
                            elif comp == 'Y' or comp == 'N':
                                data_dict['comp_NY'] = scale*read_func(filepath)
                            elif comp == 'Z':
                                data_dict['comp_Z'] = scale*read_func(filepath)
                            else:
                                raise Exception(f'Could not find component: "{comp}"')
                                
                        l_data.append(data_dict)
                            
        self['data_df'] = pd.DataFrame.from_records(l_data, index=self['default_stat_midx'])

        
    def _check_dfs_are_ok(self):
        
        data_idx = self.data_df.index
        solu_idx = self.solutions_df.index
        stat_idx = self.stations_df.index
        stat_2_solu_idx = self.stations_df.reset_index().set_index(['eid','sid']).index

        return all(data_idx == stat_idx) and all(solu_idx == stat_2_solu_idx.unique())
        
        
    @property
    def data(self):
        return self._TraceData(self.data_df)
        
    @property
    def data_df(self):
        return self['data_df']
    
    @property
    def component_names(self):
        return ['comp_EX','comp_NY','comp_Z']



def make_recip_cmt_record_from_rgf(recip_record,rgf_cmt_data_df,mt_dict):

    solu_df = recip_record.solutions_df
    stat_df = recip_record.stations_df
    l_recip_cmtsolutions = []
    l_recip_cmtstations  = []
    proj_id = recip_record.proj_id
    for eidx, edf in rgf_cmt_data_df.groupby(level='eid'):
        eid    = eidx
        mt     = mt_dict[eid]
        print(f'mt[{eid}]:\n{mt.m6_up_south_east()}')
        for tidx, tdf in edf.groupby(level='trid'):
            date      = datetime.datetime.now()
            lon_xc    = solu_df.loc[(eidx,0),'lon_xc']
            lat_yc    = solu_df.loc[(eidx,0),'lat_yc']
            depth     = solu_df.loc[(eidx,0),'depth']
            elevation = 0.
            network   = stat_df.loc[(tidx,0,eidx,0),'network']
            stat_header = StationHeader(name=f'Reciprocal-Station:{tidx}',
                                        lat_yc=lat_yc,
                                        lon_xc=lon_xc,
                                        depth=depth,
                                        elevation=elevation,
                                        network=network,
                                        proj_id=proj_id,
                                        eid=eid,
                                        sid=0,
                                        trid=tidx,
                                        gid=0)
            l_recip_cmtstations.append(stat_header)
            cmt_lon_xc    = stat_df.loc[(tidx,0,eidx,0),'lon_xc']
            cmt_lat_yc    = stat_df.loc[(tidx,0,eidx,0),'lat_yc']
            cmt_depth     = stat_df.loc[(tidx,0,eidx,0),'depth']

        cmt_header = CMTSolutionHeader(ename=f'Reciprocal-CMT:{eid}',
                                       lat_yc=cmt_lat_yc,
                                       lon_xc=cmt_lon_xc,
                                       depth=cmt_depth,
                                       tshift=0,
                                       date=date,
                                       hdur=0,
                                       mt=mt,
                                       proj_id=proj_id,
                                       eid=eid,
                                       sid=0)

        l_recip_cmtsolutions.append(cmt_header)

    record_h = RecordHeader(name=f'Reciprocal of:{recip_record.name}',
                            solutions_h=l_recip_cmtsolutions,
                            stations_h=l_recip_cmtstations,
                            proj_id=proj_id,
                            rid=recip_record.rid,
                            iter_id=recip_record.iter_id,
                            is_reciprocal=False)
    
    return Record.from_header_and_data(record_h,rgf_cmt_data_df)



def make_cmt_data_from_rgf(rgf_df,mt_dict,force_stf,cmt_stf):

    comp_dict = {'comp_EX':0,'comp_NY':1,'comp_Z':2}

    rgf_events = list(rgf_df.index.get_level_values('eid').unique())
    #print(f'rgf_events: {rgf_events}')
    mt_events  = list(mt_dict.keys())
    #print(f'mt_events: {mt_events}')
    #print(f'all: {rgf_events == mt_events}')

    if not rgf_events == mt_events:
        raise Exception('RGF-events do not match MomentTensors-events')

    l_recip_cmt_traces = []
    for eidx, edf in rgf_df.groupby(level='eid'):

        mw = mt_dict[eidx].magnitude
        m0 = mt_dict[eidx].moment
        mt_arr = mt_dict[eidx].m6_up_south_east()

        wzz =  mt_arr[0] #mrr
        wyy =  mt_arr[1] #mtt
        wxx =  mt_arr[2] #mpp
        wyz = -mt_arr[3] #mrt
        wxz =  mt_arr[4] #mrp
        wxy = -mt_arr[5] #mtp


        for tidx, tdf in edf.groupby(level='trid'):
            d_recip_cmt = {'eid':eidx,'sid':0,'trid':tidx,'gid':0}
            for comp_key in comp_dict.keys():
                ic = comp_dict[comp_key]

                composite_trace  = wxx*1*rgf_df.loc[(eidx,tidx, 0,ic, 0),'data'] #Matrix: Mee
                composite_trace += wyy*1*rgf_df.loc[(eidx,tidx, 1,ic, 1),'data'] #Matrix: Mnn
                composite_trace += wzz*1*rgf_df.loc[(eidx,tidx, 2,ic, 2),'data'] #Matrix: Mzz

                #Matrix: M1/Mxy
                composite_trace += wxy*1*rgf_df.loc[(eidx,tidx, 1,ic, 0),'data']
                composite_trace += wxy*1*rgf_df.loc[(eidx,tidx, 0,ic, 1),'data']

                #Matrix: M2/Mxz
                composite_trace += wxz*1*rgf_df.loc[(eidx,tidx, 0,ic, 2),'data']
                composite_trace += wxz*1*rgf_df.loc[(eidx,tidx, 2,ic, 0),'data']

                #Matrix: M3/Myz
                composite_trace += wyz*1*rgf_df.loc[(eidx,tidx, 1,ic, 2),'data']
                composite_trace += wyz*1*rgf_df.loc[(eidx,tidx, 2,ic, 1),'data']


                #deconvolve and then convolved
                deconv = 1.0/force_stf[0]
                scaled_trace = deconv*np.convolve(composite_trace,cmt_stf)[:len(cmt_stf)]

                # convert back to single precision
                scaled_trace = scaled_trace.astype(np.float32)

                #add trace index to trace map dictionary
                d_recip_cmt[comp_key] = scaled_trace

            # add the mapped dictionary to a list so that it can be converted to a pd.DataFrame
            l_recip_cmt_traces.append(d_recip_cmt)

    return pd.DataFrame.from_records(l_recip_cmt_traces, index=('eid','sid','trid','gid'))



def make_rgf_data(record,fl,fh,fs):

    comp_dict = {'comp_EX':0,'comp_NY':1,'comp_Z':2}
    coord_dict = {0:'lon_xc',1:'lat_yc',2:'depth'}
    sos = signal.butter(3, [fl,fh], 'bp', fs=fs, output='sos')

    l_rgf_traces = []
    m_df = pd.merge(record.stations_df,record.data_df,on=['eid','sid','trid','gid'])
    for eidx, edf in m_df.groupby(level='eid'):
        for sidx, sdf in edf.groupby(level='sid'):
            for tidx, tdf in sdf.groupby(level='trid'):
                for comp_key in comp_dict.keys():
                    ie = tidx
                    ig = eidx
                    fi = sidx
                    for di in range(3):
                        rgf_dict = {'eid':tidx,'trid':eidx,'fid':sidx}
                        rgf_dict['cid'] = comp_dict[comp_key]
                        coord_key = coord_dict[di]
                        rgf_dict['did'] = di
                        ip1 = di+1     #coord + h
                        im1 = ip1 + 3  #coord - h
                        if di == 2:
                            tm1 = ip1
                            ip1 = im1
                            im1 = tm1
                        rgf_dict['data'] = calulate_spacial_derivative(m_df,
                                                                       eidx,
                                                                       sidx,
                                                                       tidx,
                                                                       ip1,
                                                                       im1,
                                                                       sos,
                                                                       comp_key,
                                                                       coord_key)

                        rgf_dict['component'] = comp_key
                        l_rgf_traces.append(rgf_dict)

    return pd.DataFrame.from_records(l_rgf_traces, index=('eid','trid','cid','fid','did'))


# make_rgf_data is poor choice of name, but some research code needs that
# function name
def calc_dataframe_rgf(record,fl,fh,fs):
    return make_rgf_data(record,fl,fh,fs)


def calulate_spacial_derivative(tdf,eidx,sidx,tidx,g_p1,g_m1,sos,comp_key,coord_key):
    gidx_0  = pd.IndexSlice[eidx,sidx,tidx,0]
    gidx_p1 = pd.IndexSlice[eidx,sidx,tidx,g_p1]
    gidx_m1 = pd.IndexSlice[eidx,sidx,tidx,g_m1]
    df_0    = tdf.loc[gidx_0]
    df_p1   = tdf.loc[gidx_p1]
    df_m1   = tdf.loc[gidx_m1]
    data_p1 = signal.sosfilt(sos, df_p1[comp_key].astype(np.float64))
    data_m1 = signal.sosfilt(sos, df_m1[comp_key].astype(np.float64))
    c_p1    = df_p1[coord_key]
    c_m1    = df_m1[coord_key]
    c_0     = df_0[coord_key]
    delta   = 0.5*(c_p1 - c_m1)
    h       = 2.0*np.abs(delta)
    c       = c_m1 + delta

    assert h != 0
    assert c_0-c == 0

    h_scale  = 1/h
    mt_trace = h_scale*(data_p1 - data_m1)

    return mt_trace

def calc_series_composite_recip_cmt_trace(eid,trid,mt_arr,rgf_df,force_stf,cmt_stf):

    comp_dict = {'comp_EX':0,'comp_NY':1,'comp_Z':2}

    wzz =  mt_arr[0] #mrr
    wyy =  mt_arr[1] #mtt
    wxx =  mt_arr[2] #mpp
    wyz = -mt_arr[3] #mrt
    wxz =  mt_arr[4] #mrp
    wxy = -mt_arr[5] #mtp

    cmt_trace_dict = {'eid':eid, 'trid':trid}
    for comp_key in comp_dict.keys():
        ic = comp_dict[comp_key]

        composite_trace  = wxx*rgf_df.loc[(eid,trid,0,ic, 0),'data'].copy() #Matrix: Mee
        composite_trace += wyy*rgf_df.loc[(eid,trid,1,ic, 1),'data'] #Matrix: Mnn
        composite_trace += wzz*rgf_df.loc[(eid,trid,2,ic, 2),'data'] #Matrix: Mzz

        #Matrix: M1/Mxy
        composite_trace += wxy*rgf_df.loc[(eid,trid,1,ic, 0),'data']
        composite_trace += wxy*rgf_df.loc[(eid,trid,0,ic, 1),'data']

        #Matrix: M2/Mxz
        composite_trace += wxz*rgf_df.loc[(eid,trid,0,ic, 2),'data']
        composite_trace += wxz*rgf_df.loc[(eid,trid,2,ic, 0),'data']

        #Matrix: M3/Myz
        composite_trace += wyz*rgf_df.loc[(eid,trid,1,ic, 2),'data']
        composite_trace += wyz*rgf_df.loc[(eid,trid,2,ic, 1),'data']


        #deconvolve and then convolved
        deconv = 1.0/force_stf[0]
        scaled_trace = deconv*np.convolve(composite_trace.astype(np.float64),cmt_stf.astype(np.float64))[:len(cmt_stf)]

        # convert back to single precision
        cmt_trace_dict[comp_key] = scaled_trace.astype(np.float32)

    return pd.Series(cmt_trace_dict)

def calc_dataframe_composite_recipt_cmt_traces_for_one_event(eid,mt,rgf_df,force_stf,cmt_stf):

    mt_arr = mt.m6_up_south_east()

    edf = None
    ntr = rgf_df.index.get_level_values('trid').nunique()
    for tidx in range(ntr):

        tseries = calc_series_composite_recip_cmt_trace(eid,tidx,mt_arr,rgf_df,force_stf,cmt_stf)
        edf = pd.concat([edf,tseries.to_frame().T])

    edf.set_index(['eid','trid'],inplace=True)
    return edf


def calc_dataframe_composite_recipt_cmt_for_all_events(rgf_df,mt_dict,force_stf,cmt_stf):

    comp_dict = {'comp_EX':0,'comp_NY':1,'comp_Z':2}

    rgf_events = list(rgf_df.index.get_level_values('eid').unique())
    mt_events  = list(mt_dict.keys())

    if not rgf_events == mt_events:
        raise Exception('RGF-events do not match MomentTensors-events')

    rdf = None
    ne = rgf_df.index.get_level_values('eid').nunique()
    for eidx in range(ne):

        df = calc_dataframe_composite_recipt_cmt_traces_for_one_event(eidx,mt_dict[eidx],rgf_df,force_stf,cmt_stf)
        rdf = pd.concat([rdf,df])

    return rdf
