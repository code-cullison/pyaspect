import pickle

from pyaspect.specfemio.utils import _get_file_header_paths
from pyaspect.specfemio.utils import forcesolution_2_str
from pyaspect.specfemio.utils import cmtsolution_2_str
from pyaspect.specfemio.utils import station_to_str
from pyaspect.specfemio.utils import station_list_to_str


################################################################################
#
# Helper writing functions: 
# SPECFEM3D STATIONS and SOLUTIONS type files
#
################################################################################


def _write_file(fqpname,file_str):

    f = open(fqpname, 'w')
    f.write(file_str)
    f.close()

def _write_header(fqpname,header):

    f = open(fqpname, 'wb')
    pickle.dump(header,f)
    f.close()


################################################################################
#
# Functions for writing SPECFEM3D SOLUTION type fiels
#
################################################################################


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



################################################################################
#
# Functions for writing SPECFEM3D STATION type files
#
################################################################################


def write_stations(fqp,l_stations,fname='STATIONS'):

    str_stations = station_list_to_str(l_stations)

    fqp_file, fqp_header = _get_file_header_paths(fqp,fname)

    _write_file(fqp_file,str_stations)
    _write_header(fqp_header,l_stations)
