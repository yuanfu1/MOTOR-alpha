#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
Returns values from RTTOV coefficient structure
"""
from rttov_gui_f2py import rttov_gui_get_coef_val_i0, \
    rttov_gui_get_coef_val_r0, rttov_gui_get_coef_val_c0, \
    rttov_gui_get_coef_val_i1, rttov_gui_get_coef_val_r1, \
    rttov_gui_get_size_eigenvector, \
    rttov_gui_get_eigenvector

platform_name = ['noaa', 'dmsp', 'meteosat', 'goes', 'gms',
                 'fy2', 'trmm', 'ers', 'eos', 'metop',
                 'envisat', 'msg', 'fy1', 'adeos', 'mtsat',
                 'coriolis', 'jpss', 'gifts', 'sentinel3', 'meghatr',
                 'kalpana', 'meteor', 'fy3', 'coms', 'meteor-m',
                 'gosat', 'calipso', 'dummy', 'gcom-w', 'nimbus',
                 'himawari', 'mtg', 'saral', 'metopsg', 'landsat',
                 'jason', 'gpm', 'insat1', 'insat2', 'insat3',
                 'ground', 'dscovr', 'clarreo', 'ticfire', 'aircraft',
                 'iss', 'hj1', 'g-kompsat']
instrument_name = ['hirs', 'msu', 'ssu', 'amsua', 'amsub',
                   'avhrr', 'ssmi', 'vtpr1', 'vtpr2', 'tmi',
                   'ssmis', 'airs', 'hsb', 'modis', 'atsr',
                   'mhs', 'iasi', 'amsre', 'imager', 'atms',
                   'mviri', 'seviri', 'imager', 'sounder', 'imager',
                   'vissr', 'mvisr', 'cris', 'cmis', 'viirs',
                   'windsat', 'gifts', 'ssmt1', 'ssmt2', 'saphir',
                   'madras', 'ssmisz', 'vhrr', 'imager', 'sounder',
                   'mwts', 'mwhs', 'iras', 'mwri', 'abi',
                   'mi', 'msumr', 'tansofts', 'iir', 'mwr',
                   'dummyir', 'dummymw', 'dummyhi', 'dummypo', 'scams',
                   'smmr', 'ahi', 'irs', 'altika', 'iasing',
                   'tm', 'fci', 'amsr', 'amsr2', 'vissr',
                   'slstr', 'tirs', 'amr', 'oli', 'iris',
                   'ici', 'gmi', 'mwts2', 'mwhs2', 'aster',
                   'hatpro', 'mtvzagy', 'metimage', 'mws', 'mwi',
                   'epic', 'mrir', 'si', 'mrfirs', 'mbfiri',
                   'lhr', 'ismar', 'mersi1', 'mersi2', 'ecostres',
                   'irmss', 'olci', 'thir', 'ami']


def rttov_get_coef_val_i0(name):
    """
    """
    i0, err = rttov_gui_get_coef_val_i0(name.upper())
    if (err):
        i0 = None
    return i0


def rttov_get_coef_val_r0(name):

    r0, err = rttov_gui_get_coef_val_r0(name.upper())
    if (err):
        r0 = None
    return r0


def rttov_get_coef_val_c0(name):

    c0, err = rttov_gui_get_coef_val_c0(name.upper())
    if (err):
        c0 = None
    return c0


def rttov_get_coef_val_i1(name):
    i1 = None
    m, err = rttov_gui_get_coef_val_i0('SIZE_' + name.upper())
    if (not err):
        i1, err = rttov_gui_get_coef_val_i1(name.upper(), m)
        if (err):
            i1 = None
    return i1


def rttov_get_coef_val_r1(name, i=1, j=1):
    r1 = None
    if("EIGENVECTOR" == name):
        # i est pour le numero de bande
        # j est est le numero de vecteur
        print("i, j", i, j)
        m, err = rttov_gui_get_size_eigenvector(i)
        if (not err):
            r1, err = rttov_gui_get_eigenvector(i, j, m)
            if (err):
                r1 = None
    else:
        m, err = rttov_gui_get_coef_val_i0('SIZE_' + name.upper())
        if (not err):
            r1, err = rttov_gui_get_coef_val_r1(name.upper(), m)
            if (err):
                r1 = None
    return r1


if __name__ == '__main__':
    import sys
    import os
    import numpy
    from rttov_gui_f2py import rttov_gui_load
    channels = numpy.array([2, 3], dtype=int)
    rttov_gui_coeff_dir = os.environ["RTTOV_GUI_COEFF_DIR"]
    nchannels, err = rttov_gui_load(
        channels,
        rttov_gui_coeff_dir +
        "/rttov9pred54L/rtcoef_metop_2_avhrr.dat",
        "", "", "", "")
    if err != 0:
        print("cannot load coefficients")
        sys.exit(1)
    n = rttov_get_coef_val_i0('nleVels')
    print("number of levels", n)

    n = rttov_get_coef_val_i0('NCHANNELS')
    print("number of channels", n)

    n = rttov_get_coef_val_i1('FMV_VAR')
    print("fst model variables", n)

    n = rttov_get_coef_val_i1('FF_VAL_CHN')
    print("FF_VAL_CHN", n)

    n = rttov_get_coef_val_r0('FC_SAT_HEIGHT')
    print("altitude", n)

    n = rttov_get_coef_val_r1('WAVENUMBERS')
    print("WAVENUMBERS", n)

    n = rttov_get_coef_val_r1('REF_PRESSURE')
    print("pressure levels", n)

    n = rttov_get_coef_val_r1('REF_TEMPERATURE')
    print("REF_TEMPERATURE", n)

    n = rttov_get_coef_val_r1('REF_WATERVAPOR')
    print("REF_WATERVAPOR", n)

    n = rttov_get_coef_val_r1('REF_OZONE')
    print("REF_OZONE", n)

    n = rttov_get_coef_val_r1('REF_CO2')
    print("REF_CO2", n)

    n = rttov_get_coef_val_r1('LIM_PRFL_TMAX')
    print("LIM_PRFL_TMAX", n)

    n = rttov_get_coef_val_r1('Min_Ozone')
    print("MIN_OZONE", n)

    n = rttov_get_coef_val_r1('max_Ozone')
    print("MAX_OZONE", n)

    n = rttov_get_coef_val_c0('FMV_MODEL_DEF')
    print("FMV_MODEL_DEF:", n)

    n = rttov_get_coef_val_i1('FMV_GAS_POS')
    print('FMV_GAS_POS', n)
    n = rttov_get_coef_val_i1('FMV_GAS_ID')
    print('FMV_GAS_ID ', n)

    channels = numpy.array([0], dtype=int)
    nchannels, err = rttov_gui_load(
        channels,
        (rttov_gui_coeff_dir +
         "/rttov9pred101L/rtcoef_metop_2_iasi_pcrttov_compat.H5"),
        "", "", "",
        os.environ[
            "RTTOV_GUI_COEFF_DIR"] +
        "/pc/pccoef_metop_2_iasi_landsea_trace_nlte.H5")
    if err != 0:
        print("cannot load coefficients")
        sys.exit(1)
    print("nchannels iasi", nchannels)
    n = rttov_get_coef_val_i0('FMV_PC_BANDS')
    print('FMV_PC_BANDS', n)
    n = rttov_get_coef_val_i1('FMV_PC_SETS')
    print('FMV_PC_SETS', n)
    n = rttov_get_coef_val_i0('FMV_PC_CLD')
    print('FMV_PC_CLD', n)

    n = rttov_get_coef_val_r1('EIGENVECTOR', 1, 4)
    print('EIGENVECTOR', n)
    n = rttov_get_coef_val_r1('EIGENVECTOR', 1, 400)
    print('EIGENVECTOR', n)
    n = rttov_get_coef_val_r1('NOISE_IN')
    print('INSTRUMENT NOISE', n)
    print("test OK")
