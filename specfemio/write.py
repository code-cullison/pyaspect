import pickle
import importlib

from pyaspect.specfemio.utils import _get_file_header_paths
from pyaspect.specfemio.utils import _join_path_fname
from pyaspect.specfemio.utils import forcesolution_2_str
from pyaspect.specfemio.utils import cmtsolution_2_str
from pyaspect.specfemio.utils import station_to_str
from pyaspect.specfemio.utils import station_list_to_str
from pyaspect.specfemio.utils import flatten_grouped_headers
from pyaspect.specfemio.utils import flatten_grouped_headers_unique
from pyaspect.specfemio.utils import is_grouped_headers

from pyaspect.specfemio.headers import ForceSolutionHeader
from pyaspect.specfemio.headers import CMTSolutionHeader


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

    '''
    f = open(fqpname, 'w')
    f.write(file_str)
    f.close()
    '''

def _write_header(fqpname,header):

    try:
        with open(fqpname, 'wb') as f:
            pickle.dump(header,f)

    except IOError as e:
        print(e)

    '''
    f = open(fqpname, 'wb')
    pickle.dump(header,f)
    f.close()
    '''


################################################################################
#
# Functions for writing SPECFEM3D SOLUTION type fiels
#
################################################################################


def write_solution(fqp,solution,fname='NONAME',write_h=True):

    solution_str = None

    s_fname = fname

    if isinstance(s,ForceSolutionHeader):
        if fname == 'NONAME':
            s_fname = 'FORCESOLUTION'
        solution_str = forcesolution_2_str(fs)
    elif isinstance(s,CMTSolutionHeader):
        if fname == 'NONAME':
            s_fname = 'CMTSOLUTION'
        solution_str = cmtsolution_2_str(cmts)
    else:
        raise Exception(f'solution type={type(solution) does not exist')

    fqp_file, fqp_header = _get_file_header_paths(fqp,s_fname)

    _write_file(fqp_file,solution_str)
    _write_header(fqp_header,solution)


def write_cmtsolution(fqp,cmts,fname='CMTSOLUTION'):

    write_solution(fqp=fqp,solution=cmts,fname=fname,write_h=write_h)

    '''
    cmtlines_str = cmtsolution_2_str(cmts)

    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)

    _write_file(fqp_file,cmtlines_str)
    _write_header(fqp_header,cmts)
    '''


def write_forcesolution(fqp,fs,fname='FORCESOLUTION',write_h=True):

    write_solution(fqp=fqp,solution=fs,fname=fname,write_h=write_h)

    '''
    forcelines_str = forcesolution_2_str(fs)

    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)

    _write_file(fqp_file,forcelines_str)

    if write_h:
        _write_header(fqp_header,fs)
    '''


def write_grouped_forcesolutions(fqp,l_fs,fname='FORCESOLUTION',write_h=True):

    for s in l_fs:
        write_forcesolution(fqp,s,fname=f'FORCESOLUTION.{s.sid}',write_h=write_h)



################################################################################
#
# Functions for writing SPECFEM3D STATION type files
#
################################################################################


def write_stations(fqp,l_stations,fname='STATIONS',write_h=True):

    str_stations  = None
    group_headers = None

    if is_grouped_headers(l_stations):
        str_stations  = station_list_to_str(flatten_grouped_headers_unique(l_stations))
        group_headers = station_list_to_str(flatten_grouped_headers(l_stations))
        #print('IN GROUPED WRITE')
    else:
        #print('IN NOT GROUPED WRITE')
        str_stations = station_list_to_str(l_stations)
        group_headers = l_stations

    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)

    _write_file(fqp_file,str_stations)
    
    if write_h:
        _write_header(fqp_header,group_headers)



################################################################################
#
# Functions for writing SPECFEM3D STATIONS and SOLUTIONS from Records
#
################################################################################

def write_record(proj_fqp,record,fname='project.records'):

    #first write solution files and headers
    solutions = record.get_solutions_header_list()
    edir_fqp  = _join_path_fname(data_fqp,f'run{str(solutions[0].eid).zfill(4)}')
    data_fqp  = _join_path_fname(edir_fqp,'DATA')

    solution_range = set(record._get_reset_df(is_stations=False)['sid'])
    station_range = set(record._get_reset_df(is_stations=False)['sid'])

    if solution_range != station_range:
        raise Exception('The solutions.sid are not the same as stations.sid')

    #for s in solutions:
    for i in range(solution_range):
        s = solutions[i]
        if s.sid != i:
            raise Exception(f'solution.sid is not the correct value')
        write_solution(data_fqp,s)

    #now write station files and headers
    for i in range(s_range):
        stations = get_stations_header_list(key='sid',value=i,is_stations=True)
        for s in stations:
            if s.sid != i:
                raise Exception(f'station.sid is not the correct value')

        s_fname = f'STATIONS.sid{i}'
        write_stations(data_fqp,stations,fname=s_fname,write_h=True):

   

