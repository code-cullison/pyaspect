import pickle

from pyaspect.specfemio.utils import _get_file_header_paths

################################################################################
#
# Helper reading functions: 
# SPECFEM3D STATIONS and SOLUTIONS type files
#
################################################################################

def _read_headers(fqpname):

    f = open(fqpname, 'rb')
    headers = pickle.load(f)
    f.close()

    return headers


################################################################################
#
# Functions for writing SPECFEM3D SOLUTION type fiels
#
################################################################################

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
# Functions for reading SPECFEM3D STATION type files
#
################################################################################


def read_stations(fqp,fname='STATIONS'):

    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)
    
    return  _read_headers(fqp_header)

