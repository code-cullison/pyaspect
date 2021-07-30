import numpy as np

PAR_FILE_DICT = {'SIMULATION_TYPE '                : 'int',
                 'NOISE_TOMOGRAPHY '               : 'int',
                 'SAVE_FORWARD '                   : 'bool',
                 'INVERSE_FWI_FULL_PROBLEM '       : 'bool',
                 'UTM_PROJECTION_ZONE '            : 'int',
                 'SUPPRESS_UTM_PROJECTION '        : 'bool',
                 'NPROC '                          : 'int' ,
                 'NSTEP '                          : 'int',
                 'DT '                             : 'float',
                 'USE_LDDRK '                      : 'bool',
                 'INCREASE_CFL_FOR_LDDRK '         : 'bool',
                 'RATIO_BY_WHICH_TO_INCREASE_IT '  : 'float',
                 'NGNOD '                          : 'int',
                 'MODEL '                          : 'model',
                 'TOMOGRAPHY_PATH '                : 'path_str',
                 'SEP_MODEL_DIRECTORY '            : 'path_str',
                 'APPROXIMATE_OCEAN_LOAD '         : 'bool',
                 'TOPOGRAPHY '                     : 'bool',
                 'ATTENUATION '                    : 'bool',
                 'ANISOTROPY '                     : 'bool',
                 'GRAVITY '                        : 'bool',
                 'ATTENUATION_f0_REFERENCE '       : 'd0',
                 'MIN_ATTENUATION_PERIOD '         : 'd0',
                 'MAX_ATTENUATION_PERIOD '         : 'd0',
                 'COMPUTE_FREQ_BAND_AUTOMATIC '    : 'bool',
                 'USE_OLSEN_ATTENUATION '          : 'bool',
                 'OLSEN_ATTENUATION_RATIO '        : 'float',
                 'PML_CONDITIONS '                 : 'bool',
                 'PML_INSTEAD_OF_FREE_SURFACE '    : 'bool',
                 'f0_FOR_PML '                     : 'float',
                 'STACEY_ABSORBING_CONDITIONS '    : 'bool',
                 'STACEY_INSTEAD_OF_FREE_SURFACE ' : 'bool',
                 'BOTTOM_FREE_SURFACE '            : 'bool',
                 'UNDO_ATTENUATION_AND_OR_PML '    : 'bool',
                 'NT_DUMP_ATTENUATION '            : 'int',
                 'CREATE_SHAKEMAP '                : 'bool',
                 'MOVIE_SURFACE '                  : 'bool',
                 'MOVIE_TYPE '                     : 'int',
                 'MOVIE_VOLUME '                   : 'bool',
                 'SAVE_DISPLACEMENT '              : 'bool',
                 'USE_HIGHRES_FOR_MOVIES '         : 'bool',
                 'NTSTEP_BETWEEN_FRAMES '          : 'int',
                 'HDUR_MOVIE '                     : 'float',
                 'SAVE_MESH_FILES '                : 'bool',
                 'LOCAL_PATH '                     : 'path_str',
                 'NTSTEP_BETWEEN_OUTPUT_INFO '     : 'int',
                 'USE_SOURCES_RECEIVERS_Z '        : 'bool',
                 'USE_FORCE_POINT_SOURCE '         : 'bool',
                 'USE_RICKER_TIME_FUNCTION '       : 'bool',
                 'USE_EXTERNAL_SOURCE_FILE '       : 'bool',
                 'PRINT_SOURCE_TIME_FUNCTION '     : 'bool',
                 'NTSTEP_BETWEEN_OUTPUT_SEISMOS '  : 'int',
                 'SAVE_SEISMOGRAMS_DISPLACEMENT '  : 'bool',
                 'SAVE_SEISMOGRAMS_VELOCITY '      : 'bool',
                 'SAVE_SEISMOGRAMS_ACCELERATION '  : 'bool',
                 'SAVE_SEISMOGRAMS_PRESSURE '      : 'bool',
                 'SAVE_SEISMOGRAMS_IN_ADJOINT_RUN ': 'bool',
                 'USE_BINARY_FOR_SEISMOGRAMS '     : 'bool',
                 'SU_FORMAT '                      : 'bool',
                 'ASDF_FORMAT '                    : 'bool',
                 'WRITE_SEISMOGRAMS_BY_MASTER '    : 'bool',
                 'SAVE_ALL_SEISMOS_IN_ONE_FILE '   : 'bool',
                 'USE_TRICK_FOR_BETTER_PRESSURE '  : 'bool',
                 'USE_SOURCE_ENCODING '            : 'bool',
                 'OUTPUT_ENERGY '                  : 'bool',
                 'NTSTEP_BETWEEN_OUTPUT_ENERGY '   : 'int',
                 'NTSTEP_BETWEEN_READ_ADJSRC '     : 'int',
                 'READ_ADJSRC_ASDF '               : 'bool',
                 'ANISOTROPIC_KL '                 : 'bool',
                 'SAVE_TRANSVERSE_KL '             : 'bool',
                 'ANISOTROPIC_VELOCITY_KL '        : 'bool',
                 'APPROXIMATE_HESS_KL '            : 'bool',
                 'SAVE_MOHO_MESH '                 : 'bool',
                 'COUPLE_WITH_INJECTION_TECHNIQUE ': 'bool',
                 'INJECTION_TECHNIQUE_TYPE '       : 'int',
                 'MESH_A_CHUNK_OF_THE_EARTH '      : 'bool',
                 'TRACTION_PATH '                  : 'path_str',
                 'FKMODEL_FILE '                   : 'FKmodel',
                 'RECIPROCITY_AND_KH_INTEGRAL '    : 'bool',
                 'NUMBER_OF_SIMULTANEOUS_RUNS '    : 'int',
                 'BROADCAST_SAME_MESH_AND_MODEL '  : 'bool',
                 'GPU_MODE '                       : 'bool',
                 'ADIOS_ENABLED '                  : 'bool',
                 'ADIOS_FOR_DATABASES '            : 'bool',
                 'ADIOS_FOR_MESH '                 : 'bool',
                 'ADIOS_FOR_FORWARD_ARRAYS '       : 'bool',
                 'ADIOS_FOR_KERNELS '              : 'bool',
                 'LTS_MODE '                       : 'bool',
                 'PARTITIONING_TYPE '              : 'int'}


PAR_INT_TYPES = (int,np.int32,np.int64,np.uint32,np.uint64)
PAR_FLOAT_TYPES = (*PAR_INT_TYPES,float,np.float32,np.float64)
#TODO: more?
PAR_MODEL_TYPES = ('aniso','external','gll','salton_trough','tomo','SEP')
PAR_FKMODEL_TYPES = ('FKmodel')

PAR_TYPE_DICT = {'bool':(bool),
                 'd0':PAR_INT_TYPES,
                 'float':PAR_FLOAT_TYPES,
                 'int':PAR_INT_TYPES,
                 'model':PAR_MODEL_TYPES,
                 'path_str':(str),
                 'FKmodel': PAR_FKMODEL_TYPES}

PAR_CASTING_DICT = {'bool': lambda x: x == '.true.',
                    'd0': lambda x: str(int(x)) + '.d0',
                    'float': lambda x: float(x),
                    'int': lambda x: str(int(x)),
                    'model': lambda x: str(x),
                    'path_str': lambda x: str(x),
                    'FKmodel': lambda x: str(x)}


writelines(out_parfile_fqp):

    try:
        with open(out_parfile_fqp, 'w') as f:
            f.writelines(lines)
    except IOError as e:
        print(e)


readlines(in_parfile_fqp):

    lines = None
    try:
        with open(in_parfile_fqp, 'r') as f:
            lines = f.readlines()
    except IOError as e:
        print(e)

    return lines

    if out_par_file == None:
        out_parfile_fqp = in_parfile_fqp


change_parameter_in_lines(lines,par_key,par_val):

    if not isinstance(par_key,str):
        raise TypeError(f'arg: \'par_key\' must be of type str')

    par_type = None
    if par_key not in PAR_FILE_DICT.keys():
        raise ValueError(f'Parameter: {par_key} is not a valid parameter')
    else:
        par_type = PAR_FILE_DICT[par_key]

        if not isinstance(par_val,PAR_TYPE_DICT[par_type]):
            raise ValueError(f'arg: \'par_val\' is not a valid type for parameter: {par_key}')

    found_match = False
    for i in range(len(lines)):
        line = lines[i]
        wlist = line.split()
        if len(wlist) == 0: # empty line
            il += 1
            continue
        if wlist[0].strip() == par_key:
            found_match = True
            val_str = line.split('=')[-1].strip()
            if len(val_str):
                raise Exception(f'error reading par_val for parameter: {par_key}')
            ival = line.find(val_str)
            new_line = line[:ival] + PAR_CASTING_DICT[par_key](par_val)
            lines[i] = new_line
            break

    if not found_match:
        raise Exception(f'error looking for parameter: {par_key}')

    return lines


change_multiple_parameters_in_lines(lines,s_keys_vals):

    if not isinstance(s_keys_vals,set):
        raise TypeError(f'arg: \'s_keys_vals\' must be of type set')

    new_lines = copy.deepcopy(lines)

    for k,v in s_keys_vals:
        new_lines = change_parameter_in_lines(new_lines,par_key,par_val)

    return new_lines


change_parameter(in_parfile_fqp,lines,par_key,par_val,out_parfile_fqp=None):

    lines = readlines(in_parfile_fqp)

    if out_parfile_fqp == None:
        out_parfile_fqp = in_parfile_fqp

    new_lines = change_parameter_in_lines(lines,par_key,par_val)

    writelines(out_parfile_fqp)
