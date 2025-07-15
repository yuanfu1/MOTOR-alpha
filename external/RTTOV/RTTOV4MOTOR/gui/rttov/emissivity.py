#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from rttov.core import _V


class Emissivity(dict, _V):
    """
    The Emissivity class
    A Python dictionary which contains a copy of the RTTOV v13
    rttov_emissivity structure
    Inherits the _V class for most methods
    Method:
    - setemissivity fills a Emissivity class for a given number of channels,
        emissivities are set to 0 and logical calcemis is set to True
    """
    emissivity_list = ['EMIS_IN', 'EMIS_OUT', 'SPECULARITY', 'CALCEMIS']
    list_arr_logical = ['CALCEMIS']

    def __init__(self):
        dict.__init__(self)
        for key in self.emissivity_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def setEmissivity(self, nchannels):
        """Creates emissivity arrays for a given number of channels
        EMIS_IN and EMIS_OUT are all zeros
        """
        self['EMIS_IN'] = numpy.zeros(nchannels)
        self['EMIS_OUT'] = numpy.zeros(nchannels)
        self['SPECULARITY'] = numpy.zeros(nchannels)
        self['CALCEMIS'] = numpy.empty(nchannels, bool)
        self['CALCEMIS'][:] = True
        self['EMIS_IN_ATTRIBUTE']['LBOUND'] = 1
        self['EMIS_IN_ATTRIBUTE']['UBOUND'] = nchannels
        self['EMIS_OUT_ATTRIBUTE']['LBOUND'] = 1
        self['EMIS_OUT_ATTRIBUTE']['UBOUND'] = nchannels
        self['SPECULARITY_ATTRIBUTE']['LBOUND'] = 1
        self['SPECULARITY_ATTRIBUTE']['UBOUND'] = nchannels
        self['CALCEMIS_ATTRIBUTE']['LBOUND'] = 1
        self['CALCEMIS_ATTRIBUTE']['UBOUND'] = nchannels


if __name__ == '__main__':

    import h5py
    from rttov_gui_f2py import rttov_gui_load, rttov_gui_get_emisbrdf
    from rttov.profile import Profile, getNumberOfProfiles
    from rttov.option import Option
    from rttov.reflectance import Reflectance
    import os
    from core import rttov_gui_data_test_dir
    data_test_dir = rttov_gui_data_test_dir()
    rttov_gui_profile_dir = os.environ["RTTOV_GUI_PROFILE_DIR"]
    fileName = os.path.join(rttov_gui_profile_dir,
                            "cldaer101lev_allgas.H5")

    # Lecture du profil numero 2 et des options
    nprof = getNumberOfProfiles(fileName)
    print("number of profiles found :", nprof)

    p1 = Profile()
    o1 = Option()

    iprof = 2
    p1.loadProfileNumber(fileName, iprof)

    # changer latitude longitude date et type de surface
    p1['LATITUDE'] = 45.0
    p1['LONGITUDE'] = 5.0
    p1['DATE'] = (2012, 7, 1)
    p1['SKIN']['SURFTYPE'] = 0
    p1['AZANGLE'] = 45.0
    p1['SUNAZANGLE'] = 100.0
    p1['SUNZENANGLE'] = 30.0
    p1['ZENANGLE'] = 30.0

    # Ecriture du fichier resultat "profile.H5" avec 1 profil et 1 option
    # Output p1 to HDF5 file
    ofile = os.path.join(data_test_dir, "emis_profile.H5")
    of = h5py.File(ofile, 'w')
    p1.saveh5(of, "/PROFILES/0001/")
    o1.saveh5(of, '/OPTIONS/')
    of.close()

    # RTTOV

    # On lance la lecture des coefficients RTTOV
    channel_list = numpy.array([0], dtype=int)
    rttov_coeff_dir = os.environ["RTTOV_GUI_COEFF_DIR"]
    fcoef = rttov_coeff_dir + "/rttov7pred54L/rtcoef_msg_3_seviri.dat"
    nchannels, err = rttov_gui_load(channel_list,
                                    fcoef, "", "", "", "")

    print("nchannels", nchannels)

    e1 = Emissivity()
    r1 = Reflectance()

    e1.setEmissivity(nchannels)
    r1.setReflectance(nchannels)

    ofile = os.path.join(data_test_dir, "emis_surface.H5")
    of = h5py.File(ofile, 'w')
    e1.saveh5(of, '/EMISSIVITY/')
    r1.saveh5(of, '/REFLECTANCE/')

    of.close()
    rttov_gui_emis_dir = os.environ["RTTOV_GUI_EMISS_DIR"]

    err = rttov_gui_get_emisbrdf(rttov_gui_emis_dir,
                                 os.environ["RTTOV_GUI_PREFIX"] +
                                 "/rttov_gui_data_test/emis_profile.H5",
                                 os.environ["RTTOV_GUI_PREFIX"] +
                                 "/rttov_gui_data_test/emis_surface.H5",
                                 1, 1)

    # Read emissivity and Reflectance and display
    myfile = os.path.join(data_test_dir, "emis_surface.H5")

    f = h5py.File(myfile, 'r')
    h5 = f['/EMISSIVITY/']
    e1 = Emissivity()
    e1.loadh5(h5)
    h5 = f['/REFLECTANCE/']
    r1 = Reflectance()
    r1.loadh5(h5)
    f.close()

    # Display option
    print("\n \n e1 lu dans surface.H5 ")
    e1.display()
    print("\n")
    print("\n \n r1 lu dans surface.H5 ")
    r1.display()
    print("\n")

    # Plots emissivity
    import matplotlib
    matplotlib.use("wxagg")
    import matplotlib.pyplot as plt
    plt.plot(e1['EMIS_IN'], 'r-')
    plt.plot(r1['REFL_IN'], 'r-')
    plt.xlabel('channel number -1')

    plt.ylabel(e1['EMIS_IN_ATTRIBUTE']['COMMENT'])
    plt.title(fileName)
    plt.show()
    print("test OK")
