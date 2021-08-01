import pickle
import importlib

from pyaspect.specfemio.utils import _get_file_header_paths
from pyaspect.specfemio.utils import _join_path_fname
from pyaspect.specfemio.utils import _get_file_path
from pyaspect.specfemio.utils import _get_header_path
from pyaspect.specfemio.utils import _mk_relative_symlink
from pyaspect.specfemio.utils import forcesolution_2_str
from pyaspect.specfemio.utils import cmtsolution_2_str
from pyaspect.specfemio.utils import station_to_str
from pyaspect.specfemio.utils import station_list_to_str
from pyaspect.specfemio.utils import station_auto_data_fname_id
from pyaspect.specfemio.utils import flatten_grouped_headers
from pyaspect.specfemio.utils import flatten_grouped_headers_unique
from pyaspect.specfemio.utils import is_grouped_headers

from pyaspect.specfemio.headers import StationHeader
from pyaspect.specfemio.headers import ForceSolutionHeader
from pyaspect.specfemio.headers import CMTSolutionHeader
from pyaspect.specfemio.headers import RecordHeader


################################################################################
#
# Helper writing functions: 
# SPECFEM3D STATIONS and SOLUTIONS type files
#
################################################################################


def _write_file(fqpname,file_str):

    try:
        with open(fqpname, 'w') as f:
            f.write(file_str)

    except IOError as e:
        print(e)


def _write_header(fqpname,header):

    try:
        with open(fqpname, 'wb') as f:
            pickle.dump(header,f)

    except IOError as e:
        print(e)



################################################################################
#
# Functions for writing SPECFEM3D SOLUTION type fiels
#
################################################################################


def write_solution(fqp,
                   solution,
                   fname=None,
                   postfix=None,
                   write_h=True,
                   mk_sym_link=False):

    solution_str = None

    s_fname  = fname
    sym_name = None

    if isinstance(solution,ForceSolutionHeader):
        if fname == None:
            s_fname = 'FORCESOLUTION'
        if postfix != None:
            s_fname += f'.{postfix}'
        solution_str = forcesolution_2_str(solution)
        sym_name = 'FORCESOLUTION'
    elif isinstance(solution,CMTSolutionHeader):
        if fname == None:
            s_fname = 'CMTSOLUTION'
        if postfix != None:
            s_fname += f'.{postfix}'
        solution_str = cmtsolution_2_str(solution)
        sym_name = 'CMTSOLUTION'
    else:
        raise Exception(f'solution type={type(solution)} does not exist')

    fqp_file, fqp_header = _get_file_header_paths(fqp,s_fname)

    _write_file(fqp_file,solution_str)

    if write_h:
        _write_header(fqp_header,solution)

    if mk_sym_link and s_fname != 'FORCESOLUTION' and s_fname != 'CMTSOLUTION':
        dst = _join_path_fname(fqp, sym_name)
        _mk_relative_symlink(fqp_file,fqp,dst)


def write_cmtsolution(fqp,cmts,fname='CMTSOLUTION'):
    write_solution(fqp=fqp,solution=cmts,fname=fname,write_h=write_h)


def write_forcesolution(fqp,fs,fname='FORCESOLUTION',write_h=True):
    write_solution(fqp=fqp,solution=fs,fname=fname,write_h=write_h)


def write_grouped_forcesolutions(fqp,l_fs,fname='FORCESOLUTION',write_h=True):
    for s in l_fs:
        write_forcesolution(fqp,s,fname=f'FORCESOLUTION.{s.sid}',write_h=write_h)



################################################################################
#
# Functions for writing SPECFEM3D STATION type files
#
################################################################################


def write_stations(fqp,
                   l_stations,
                   fname='STATIONS',
                   write_h=True,
                   auto_name=False,
                   auto_network=False,
                   mk_sym_link=False):

    str_stations  = None
    group_headers = None

    flat_l_stations = l_stations

    if is_grouped_headers(l_stations):
        flat_l_stations = flatten_grouped_headers(l_stations)

    str_stations = station_list_to_str(flat_l_stations,
                                       auto_name=auto_name,
                                       auto_network=auto_network)

    group_headers = flat_l_stations

    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)

    _write_file(fqp_file,str_stations)
    
    if write_h:
        _write_header(fqp_header,group_headers)

    if mk_sym_link and fname != 'STATIONS':
        dst = _join_path_fname(fqp, 'STATIONS')
        _mk_relative_symlink(fqp_file,fqp,dst)



################################################################################
#
# Functions for writing SPECFEM3D STATIONS and SOLUTIONS from Records
#
################################################################################

def write_record(proj_fqp,
                 record,
                 fname='event_record',
                 write_record_h=True,
                 write_h=False,
                 auto_name=False,
                 auto_network=False):

    # do some checking
    solution_range = set(record._get_reset_df(is_stations=False)['sid'])
    station_range = set(record._get_reset_df(is_stations=False)['sid'])

    if solution_range != station_range:
        raise Exception('The solutions.sid are not the same as stations.sid')
    
    # get solutions both for path making and for writing
    solutions = record.get_solutions_header_list()
    eid = solutions[0].eid
    rid = record.rid

    #do some more checking
    for s in solutions:
        if s.eid != eid:
            raise Exception('soution.eid does not match list of solutions')

    # make some paths
    edir_fqp   = _join_path_fname(proj_fqp,f'run{str(eid+1).zfill(4)}')
    data_fqp   = _join_path_fname(edir_fqp,'DATA')
    record_fqp = _get_header_path(data_fqp,fname)

    
    # write solution files and headers
    for i in range(len(solution_range)):
        s = solutions[i]
        if s.sid != i:
            raise Exception(f'solution.sid is not the correct value')
        mk_sym_link = False
        if i == 0: #write SOLUTION symlink
            mk_sym_link = True
        write_solution(data_fqp,s,postfix=f'sid{i}',write_h=write_h,mk_sym_link=mk_sym_link)

    # write station files and headers
    for i in range(len(station_range)):
        stations = record.get_stations_header_list(key='sid',value=i)

        # more checking
        for s in stations:
            if s.sid != i:
                raise Exception(f'station.sid is not the correct value')
            if s.eid != eid:
                raise Exception('station.eid does not match list of stations')

            #TODO: could make '/SYN' dir variable dynamic
            #      also with out knowing specfem DT, can't
            #      finish the name of the station data (trace)
            data_fname = station_auto_data_fname_id(s)
            s.data_fqdn = _join_path_fname(edir_fqp,f'/SYN/{data_fname}')
        
        mk_sym_link = False
        if i == 0: #write STATIONS symlink
            mk_sym_link = True

        s_fname = f'STATIONS.sid{i}'
        write_stations(data_fqp,
                       stations,
                       fname=s_fname,
                       write_h=write_h,
                       auto_name=auto_name,
                       auto_network=auto_network,
                       mk_sym_link=mk_sym_link)

    # write record header
    if write_record_h:
        _write_header(record_fqp,record)
   

def write_records(proj_fqp,
                  l_records,
                  fname='proj_records',
                  write_record_h=True,
                  write_h=False,
                  auto_name=False,
                  auto_network=False):

    if not isinstance(l_records,list):
        raise ArgumentError('l_records must be a list')

    if not isinstance(l_records[0],RecordHeader):
        raise TypeError('l_records elements must be type:{type(RecordHeader)}')

    if len(set([r.rid for r in l_records])) != len(l_records):
        raise Exception('some records have the same rid')

    for r in l_records:
        write_record(proj_fqp,
                     r,
                     fname=f'event_record',
                     write_record_h=write_record_h,
                     write_h=write_h,
                     auto_name=auto_name,
                     auto_network=auto_network)

    l_records_fqp = _get_header_path(proj_fqp,fname)
    _write_header(l_records_fqp,l_records)

