#! /usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
from rttov.profile import Profile
from rttov.option import Option
from rttov.chanprof import Chanprof
from rttov.emissivity import Emissivity
from rttov.reflectance import Reflectance
from rttov.misc import Misc
import numpy as np
import copy


def removeall(l, o):
    while o in l:
        l.remove(o)


class Kmatrix(object):
    """
    The Kmatrix class
    This class contains all K Matrix related information
    Its members are
    - profile: the profile considered for the K calculation
    - option: calculation options
    - kmatrix: a profile class containing the K results
    - chanprof: input channels/profiles
    - emissivity: K result for emissivity
    - reflectance: K result for reflectance
    - misc: important coef file names, list of wavenumbers...
    """

    def __init__(self):
        self.profile = Profile()
        self.option = Option()
        self.kmatrix = Profile()
        self.chanprof = Chanprof()
        self.emissivity = Emissivity()
        self.reflectance = Reflectance()
        self.misc = Misc()
        self.scaled = False
        self.scalecoef = 1.0
        self.origKmat = None

    def display(self):
        self.profile.display()
        self.option.display()
        self.kmatrix.display()
        self.chanprof.display()
        self.emissivity.display()
        self.reflectance.display()
        self.misc.display()

    def loadh5(self, h5):
        h5p = h5['PROFILES/0001']
        self.profile.loadh5(h5p)
        h5o = h5['OPTIONS']
        self.option.loadh5(h5o)
        if 'KMATRIX' in list(h5.keys()):
            h5k = h5['KMATRIX']
            self.kmatrix.loadh5(h5k)
        h5c = h5['CHANPROF']
        self.chanprof.loadh5(h5c)
        if 'EMISSIVITY_K' in list(h5.keys()):
            h5e = h5['EMISSIVITY_K']
            self.emissivity.loadh5(h5e)
        if 'REFLECTANCE_K' in list(h5.keys()):
            h5r = h5['REFLECTANCE_K']
            self.reflectance.loadh5(h5r)
        h5m = h5['MISC']
        self.misc.loadh5(h5m)

    def saveh5(self, h5, path):
        self.profile.saveh5(h5, path + '/PROFILES/0001')
        self.option.saveh5(h5, path + '/OPTIONS/')
        self.kmatrix.saveh5(h5, path + '/KMATRIX/')
        self.chanprof.saveh5(h5, path + '/CHANPROF/')
        self.emissivity.saveh5(h5, path + '/EMISSIVITY_K/')
        self.reflectance.saveh5(h5, path + '/REFLECTANCE_K/')
        self.misc.saveh5(h5, path + '/MISC/')

    def getnchannels(self):
        "gives the number of channels present in the K-Matrix"
        nchan = self.kmatrix['T'].shape[1]
        return nchan

    def getwavenumber(self, channel):
        c = int(channel)
        if np.isscalar(self.misc["WAVENUMBERS"]):
            return self.misc["WAVENUMBERS"]
        else:
            return self.misc["WAVENUMBERS"][c - 1]

    def getinstrument(self):
        return self.misc["INSTRUMENT"]

    def getsatname(self):
        return self.misc["SATELLITE"]

    def getchanprof(self, channel):
        "Extracts a jacobian for a given channel"
        "A jacobian for a channel is a profile structure"
        "Channel should start from 0"
        ip = int(channel)
        # take care dimension order is opposite as Fortran"
        chanprof = Profile()
        # a cause des predicteurs on peut avoir nchannel different
        nchan = self.kmatrix['T'].shape[1]
        if (ip < 0 or ip >= nchan):
            print(("ERROR Invalid channel index",
                   ip, " allowed is [0,", nchan - 1, "]"))
            chanprof.__init__
            return chanprof
        for key in (self.kmatrix.profile_list +
                    self.kmatrix.aerosol_list +
                    self.kmatrix.cloud_list):
            if(self.kmatrix[key] is not None):
                if(key == "DATE" or key == "TIME" or key == "ID"):
                    chanprof[key] = self.kmatrix[key]
                else:
                    if np.isscalar(self.kmatrix[key]):
                        chanprof[key] = self.kmatrix[key]
                    elif(self.kmatrix[key].ndim == 0):
                        chanprof[key] = self.kmatrix[key]
                    elif(self.kmatrix[key].ndim == 1):
                        chanprof[key] = self.kmatrix[key][ip]
                    elif(self.kmatrix[key].ndim == 2):
                        chanprof[key] = self.kmatrix[key][:, ip]
                    elif(self.kmatrix[key].ndim == 3):
                        chanprof[key] = self.kmatrix[key][:, :, ip]
                    else:
                        print("ERROR kmatrix getchanprof")
                chanprof[key + '_ATTRIBUTE'] = self.kmatrix[key + '_ATTRIBUTE']
                if ('LBOUND' in list(chanprof[key + '_ATTRIBUTE'].keys())):
                    del chanprof[key + '_ATTRIBUTE']['LBOUND']
                    del chanprof[key + '_ATTRIBUTE']['UBOUND']

        for key in self.kmatrix.s2m_list:
            if(self.kmatrix['S2M'][key] is not None):
                if(self.kmatrix['S2M'][key].ndim == 0):
                    chanprof['S2M'][key] = self.kmatrix['S2M'][key]
                elif(self.kmatrix['S2M'][key].ndim == 1):
                    chanprof['S2M'][key] = self.kmatrix['S2M'][key][ip]
                else:
                    print("ERROR kmatrix getchanprof")
                chanprof['S2M'][
                    key + '_ATTRIBUTE'] = self.kmatrix['S2M'][
                    key + '_ATTRIBUTE']
                if ('LBOUND' in list(chanprof['S2M'][
                        key + '_ATTRIBUTE'].keys())):
                    del chanprof['S2M'][key + '_ATTRIBUTE']['LBOUND']
                    del chanprof['S2M'][key + '_ATTRIBUTE']['UBOUND']

        for key in self.kmatrix.sskin_list:
            if(self.kmatrix['SKIN'][key] is not None):
                if(key == "FASTEM"):
                    chanprof['SKIN'][key] = self.kmatrix['SKIN'][key]
                else:
                    if(self.kmatrix['SKIN'][key].ndim == 0):
                        chanprof['SKIN'][key] = self.kmatrix['SKIN'][key]
                    elif(self.kmatrix['SKIN'][key].ndim == 1):
                        chanprof['SKIN'][key] = self.kmatrix['SKIN'][key][ip]
                    else:
                        print("ERROR kmatrix getchanprof")
                chanprof['SKIN'][
                    key + '_ATTRIBUTE'] = self.kmatrix['SKIN'][
                    key + '_ATTRIBUTE']
                if ('LBOUND' in list(chanprof['SKIN'][
                        key + '_ATTRIBUTE'].keys())):
                    del chanprof['SKIN'][key + '_ATTRIBUTE']['LBOUND']
                    del chanprof['SKIN'][key + '_ATTRIBUTE']['UBOUND']
        return chanprof

    def kreset(self):
        self.kscale(1.0)

    def kscale(self, coef):
        """Scale KMatrix by input profile * coef
        Temperature arrays are not modified
        This method could be called several times, scaling coefficient is
        always considered relative to the original matrix
        Could be called with coef = 1.0 even for the first call
        in that case the kMatrix is reset to the original values and units"""

        if (self.scaled):
            self.kmatrix = copy.deepcopy(self.origKmatrix)
        else:
            self.origKmatrix = copy.deepcopy(self.kmatrix)
        val = coef

        if (coef == 1.0):
            print("KMatrix reset")
            reset = True
            self.scaled = False
            return
        else:
            print("KMatrix scaling by coef ", coef)
            reset = False

        coefstr = '{}'.format(coef)
        try:
            nchan = len(self.chanprof['CHANNELS'])
        except Exception:
            nchan = 1

        for key in (self.kmatrix.profile_list +
                    self.kmatrix.aerosol_list +
                    self.kmatrix.cloud_list):
            if(self.kmatrix[key] is not None):
                if(
                        key == "DATE" or key == "TIME" or key == "ID" or
                        key == "T" or key == "GAS_UNITS"):
                    pass
                else:
                    if not isinstance(self.kmatrix[key], np.ndarray):
                        continue
                    if(self.kmatrix[key].ndim == 0):
                        pass
                    if(self.kmatrix[key].ndim == 1):
                        self.kmatrix[key][:] = self.kmatrix[key][
                            :] * self.profile[key][None] * val
                    elif(self.kmatrix[key].ndim == 2):
                        self.kmatrix[key][:, :] = self.kmatrix[key][
                            :, :] * self.profile[key][:, None] * val
                    elif(self.kmatrix[key].ndim == 3):
                        print(self.kmatrix[key].shape)
                        self.profile[key] = self.kmatrix[key][:, :, ip]
                        self.kmatrix[key][:, :, :] = self.kmatrix[key][
                            :, :, :] * self.profile[key][:, :, None] * val
                    self.kmatrix[
                        key + '_ATTRIBUTE']['UNITS'] = (
                            str(
                                self.kmatrix[key + '_ATTRIBUTE'][
                                    'UNITS']) +
                        " scaled by " + str(coefstr) + " input profile")

        for key in self.kmatrix.s2m_list:
            if(self.kmatrix['S2M'][key] is not None):
                if(self.kmatrix['S2M'][key].ndim == 0):
                    pass
                elif(key == "T"):
                    pass
                elif(self.kmatrix['S2M'][key].ndim == 1):
                    for i in range(nchan):
                        self.kmatrix['S2M'][key][i] = self.kmatrix['S2M'][
                            key][i] * self.profile['S2M'][key] * val

                    self.kmatrix['S2M'][
                        key + '_ATTRIBUTE']['UNITS'] = (str(
                            self.kmatrix['S2M'][
                                key + '_ATTRIBUTE']['UNITS']) +
                        " scaled by " + coefstr +
                        " input profile")

        for key in self.kmatrix.sskin_list:
            if(self.kmatrix['SKIN'][key] is not None):
                if(key == "FASTEM" or key == "T"):
                    pass
                else:
                    if(self.kmatrix['SKIN'][key].ndim == 0):
                        pass
                    elif(self.kmatrix['SKIN'][key].ndim == 1):
                        for i in range(nchan):
                            self.kmatrix['SKIN'][key][i] = self.kmatrix[
                                'SKIN'][
                                    key][i] * self.profile['SKIN'][key] * val
                        self.kmatrix['SKIN'][key + '_ATTRIBUTE'][
                            'UNITS'] = (str(
                                self.kmatrix['SKIN'][
                                    key + '_ATTRIBUTE']['UNITS']) +
                            " scaled by " + coefstr + " input profile")

        if (reset):
            self.scaled = False
        else:
            self.scaled = True
        self.scalecoef = coef

    def read(self, filename):
        f = h5py.File(filename, 'r')
        # get the Dataset
        h5 = f['/']
        # Load kmatrix
        self.loadh5(h5)
        print('loaded')
        # Close HDF file
        f.close()


if __name__ == '__main__':
    import os
    from core import rttov_gui_data_test_dir
    data_test_dir = rttov_gui_data_test_dir()
    filename = os.path.join(data_test_dir, "kmat.h5")

    # Open file ReadOnly
    f = h5py.File(filename, 'r')

    # get the Dataset
    h5 = f['/']

    # p1 is a rttov_hdf_mod Profile instance
    k = Kmatrix()

    # Load kmatrix
    k.loadh5(h5)
    print('loaded')
    print(list(k.kmatrix.keys()))
    print(">>>> WAVENUMBER channel 1", k.getwavenumber(1))
    # Close HDF file
    f.close()

    # Display kmatrix
    # k.display()
    # k.kmatrix.display()

    # k.profile.display()
    # k.option.display()

    # Plots temperature and gas concentrations
    nchan = len(k.chanprof['CHANNELS'])
    nlev = k.profile['NLEVELS']
    import matplotlib
    matplotlib.use("wxagg")
    import matplotlib.pyplot as plt
    plt.imshow(k.kmatrix['T'], origin='upper', interpolation='nearest')
    plt.show()

    k.kscale(0.1)
    plt.imshow(k.kmatrix['T'], origin='upper', interpolation='nearest')
    plt.show()

    p1 = Profile()
    ip = 6
    print("\n\n Get Chan Profile : ", ip)
    p1 = k.getchanprof(ip)
    p1.display()
    print("\n\n Fin Chan Profile : ", ip)

    plt.plot(p1['T'], k.profile['P'], 'ro-')
    plt.xlabel(str(p1['T_ATTRIBUTE']['COMMENT']) +
               '   (' + str(p1['T_ATTRIBUTE']['UNITS']) + ')')
    plt.ylim(1100, 0.01)
    plt.yscale('log')
    plt.ylabel(str(k.profile['P_ATTRIBUTE']['COMMENT']) +
               '  ' + str(k.profile['P_ATTRIBUTE']['UNITS']))
    plt.title(p1['ID'] + ' TEMPERATURE')
    plt.show()

    if(p1.anyCloud()):
        nlevels = k.profile['NLEVELS']
        player = (k.profile['P'][0:nlevels - 1] +
                  k.profile['P'][1:nlevels]) / 2.0
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        i = 0
        lastcle = ""
        for cle in ['CFRAC'] + p1.cloud_list:
            if(p1[cle] is None):
                continue
            lastcle = cle
            color = colors[i]
            x = p1[cle]
            plt.plot(x, player, color + 'o-')
            i += 1
            plt.xlabel(str(p1[cle + '_ATTRIBUTE']['COMMENT']) +
                       '  (' + str(p1[cle + '_ATTRIBUTE']['UNITS']) + ')')
            # plt.xscale('log')
            plt.ylim(1100, 10)
            plt.yscale('log')
            plt.ylabel(str(k.profile['P_ATTRIBUTE']['COMMENT']) +
                       '   (' + str(k.profile['P_ATTRIBUTE']['UNITS']) + ')')
            plt.title(p1['ID'] + ' CLOUDS')
            plt.show()

    if (p1.anyAerosol()):

        nlevels = k.profile['NLEVELS']
        player = (k.profile['P'][0:nlevels - 1] +
                  k.profile['P'][1:nlevels]) / 2.0
        colors = ['b', 'g', 'r', 'c', 'm', 'y']
        i = 0
        lastcle = ""
        for cle in p1.aerosol_list:
            if(p1[cle] is None):
                continue
            # print "cle ",cle
            lastcle = cle
            color = colors[i % len(colors)]
            if (i >= len(colors)):
                pattern = 'o'
            else:
                pattern = 'v'
            x = p1[cle]
            plt.plot(x, player, color + pattern + '-')
            i += 1
            plt.xlabel(str(p1[cle + '_ATTRIBUTE']['COMMENT']) +
                       '  (' + str(p1[lastcle + '_ATTRIBUTE']['UNITS']) + ')')
            # plt.xscale('log')
            plt.ylim(1100, 5)
            plt.yscale('log')
            plt.ylabel(k.profile['P_ATTRIBUTE']['COMMENT'] +
                       '   (' + k.profile['P_ATTRIBUTE']['UNITS'] + ')')
            # plt.legend()
            plt.title(p1['ID'] + ' AEROSOLS')
            plt.show()

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    i = 0
    for gas in Profile.gas_list:
        if(p1[gas] is None):
            continue
        color = colors[i]
        i += 1
        plt.plot(p1[gas], k.profile['P'], color + 'o-', label=gas)
        plt.xlabel('gas concentrations (' +
                   p1[gas + '_ATTRIBUTE']['UNITS'] + ')')
        # plt.xscale('log')
        plt.ylim(1100, 0.01)
        plt.yscale('log')
        plt.ylabel(str(k.profile['P_ATTRIBUTE']['COMMENT']) +
                   '   (' + str(k.profile['P_ATTRIBUTE']['UNITS']) + ')')
        plt.legend()
        plt.title(p1['ID'] + ' GASES')
        plt.show()

    k.kscale(0.123456)
    k.kreset()
    print("test OK")
