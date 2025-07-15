'''
Created on Feb 21, 2014

@author: pascale
'''


from rttov.profile import Profile
from rttov.kmatrix import Kmatrix
from rttov.pccomp import PCCOMP


def removeall(l, o):
    while o in l:
        l.remove(o)


class Kpcmatrix(Kmatrix):
    """
    The Kmatrix class
    This class contains all K Matrix PC related information
    Its members are
    - profile: the profile considered for the K calculation
    - option: calculation options
    - kmatrix: a profile class containing the K results
    - Kpcmatrix: a profile class containing the K PC results
    - chanprof: input channels/profiles
    - emissivity: K result for emissivity
    - reflectance: K result for reflectance
    - misc: important coef file names, list of wavenumbers...
    """

    def __init__(self):
        Kmatrix.__init__(self)
        self.kpcmatrix = Profile()
        self.pccomp = PCCOMP()
        self.pccomp_k = PCCOMP()

    def display(self):
        print("-------profile-------------")
        self.profile.display()
        print("--------option -----------")
        self.option.display()
        print("------kmatrix------------")
        self.kmatrix.display()
        print("--------kpcmatrix---------")
        self.kpcmatrix.display
        print("------chanprof------------")
        self.chanprof.display()
        print("---------emissivity---------")
        self.emissivity.display()
        print("------reflectance-----------")
        self.reflectance.display()
        print("---------misc-----------")
        self.misc.display()
        print("number of pc : ", self.nbpc)
        print("-------------pccomp--------------")
        self.pccomp.display()
        print("-------------pccomp_k------------")
        self.pccomp_k.display()

    def loadh5(self, h5):
        h5p = h5['PROFILES/0001']
        self.profile.loadh5(h5p)
        h5o = h5['OPTIONS']
        self.option.loadh5(h5o)
        if 'KMATRIX' in list(h5.keys()):
            h5k = h5['KMATRIX']
            self.kmatrix.loadh5(h5k)
        h5kpc = h5['KMATRIX_K_PC']
        self.kpcmatrix.loadh5(h5kpc)
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
        pccomp = h5['PCCOMP']
        self.pccomp.loadh5(pccomp)
        pccomp_k = h5['PCCOMP_K']
        self.pccomp_k.loadh5(pccomp_k)
        self.nbpc = pccomp["TOTAL_PCSCORES"].shape[0]

    def saveh5(self, h5, path):
        self.profile.saveh5(h5, path + '/PROFILES/0001')
        self.option.saveh5(h5, path + '/OPTIONS/')
        self.kmatrix.saveh5(h5, path + '/KMATRIX/')
        self.kpcmatrix.saveh5(h5, path + '/KMATRIX_K_PC/')
        self.chanprof.saveh5(h5, path + '/CHANPROF/')
        self.emissivity.saveh5(h5, path + '/EMISSIVITY_K/')
        self.reflectance.saveh5(h5, path + '/REFLECTANCE_K/')
        self.misc.saveh5(h5, path + '/MISC/')
        self.pccomp.saveh5(h5, path + '/PCCOMP/')
        self.pccomp_k.saveh5(h5, path + '/PCCOMP_K/')

    def getkpcprof(self, ip):
        "Extracts a jacobian for a given channel"
        "A jacobian for a channel is a profile structure"
        "Channel should start from 0"
        # take care dimension order is opposite as Fortran"
        pcprof = Profile()
        nbpc = self.nbpc
        if (ip < 0 or ip >= nbpc):
            print(("ERROR Invalid channel index", ip,
                   " allowed is [0,", nbpc - 1, "]"))
            pcprof.__init__
            return pcprof
        for key in (self.kpcmatrix.profile_list + self.kpcmatrix.aerosol_list +
                    self.kpcmatrix.cloud_list):
            if(self.kpcmatrix[key] is not None):
                if(key == "DATE" or key == "TIME" or key == "ID"):
                    pcprof[key] = self.kpcmatrix[key]
                else:
                    if(self.kpcmatrix[key].ndim == 0):
                        pcprof[key] = self.kpcmatrix[key]
                    elif(self.kpcmatrix[key].ndim == 1):
                        pcprof[key] = self.kpcmatrix[key][ip]
                    elif(self.kpcmatrix[key].ndim == 2):
                        pcprof[key] = self.kpcmatrix[key][:, ip]
                    elif(self.kpcmatrix[key].ndim == 3):
                        pcprof[key] = self.kpcmatrix[key][:, :, ip]
                    else:
                        print("ERROR Kpcmatrix getpcprof")
                pcprof[key + '_ATTRIBUTE'] = self.kpcmatrix[key + '_ATTRIBUTE']
                if ('LBOUND' in list(pcprof[key + '_ATTRIBUTE'].keys())):
                    del pcprof[key + '_ATTRIBUTE']['LBOUND']
                    del pcprof[key + '_ATTRIBUTE']['UBOUND']

        for key in self.kpcmatrix.s2m_list:
            if(self.kpcmatrix['S2M'][key] is not None):
                if(self.kpcmatrix['S2M'][key].ndim == 0):
                    pcprof['S2M'][key] = self.kpcmatrix['S2M'][key]
                elif(self.kpcmatrix['S2M'][key].ndim == 1):
                    pcprof['S2M'][key] = self.kpcmatrix['S2M'][key][ip]
                else:
                    print("ERROR Kpcmatrix getpcprof")
                pcprof['S2M'][
                    key + '_ATTRIBUTE'] = self.kpcmatrix['S2M'][
                    key + '_ATTRIBUTE']
                if ('LBOUND' in list(pcprof['S2M'][
                        key + '_ATTRIBUTE'].keys())):
                    del pcprof['S2M'][key + '_ATTRIBUTE']['LBOUND']
                    del pcprof['S2M'][key + '_ATTRIBUTE']['UBOUND']

        for key in self.kpcmatrix.sskin_list:
            if(self.kpcmatrix['SKIN'][key] is not None):
                if(key == "FASTEM"):
                    pcprof['SKIN'][key] = self.kpcmatrix['SKIN'][key]
                else:
                    if(self.kpcmatrix['SKIN'][key].ndim == 0):
                        pcprof['SKIN'][key] = self.kpcmatrix['SKIN'][key]
                    elif(self.kpcmatrix['SKIN'][key].ndim == 1):
                        pcprof['SKIN'][key] = self.kpcmatrix['SKIN'][key][ip]
                    else:
                        print("ERROR Kpcmatrix getpcprof")
                pcprof['SKIN'][
                    key + '_ATTRIBUTE'] = self.kpcmatrix['SKIN'][
                    key + '_ATTRIBUTE']
                if ('LBOUND' in list(pcprof['SKIN'][
                        key + '_ATTRIBUTE'].keys())):
                    del pcprof['SKIN'][key + '_ATTRIBUTE']['LBOUND']
                    del pcprof['SKIN'][key + '_ATTRIBUTE']['UBOUND']
        return pcprof


if __name__ == '__main__':

    import h5py
    import os
    from core import rttov_gui_data_test_dir
    data_test_dir = rttov_gui_data_test_dir()
    fileName = os.path.join(data_test_dir, "pckmat.h5")

    # Open file ReadOnly
    f = h5py.File(fileName, 'r')

    # get the Dataset
    h5 = f['/']

    # p1 is a rttov_hdf_mod Profile instance
    kpc = Kpcmatrix()

    # Load kmatrix
    kpc.loadh5(h5)
    print('loaded')

    # Close HDF file
    f.close()
    profile = Profile()
    profile = kpc.getkpcprof(10)
    kpc.display()
    print("display kpc profile number 10")
    profile.display()
    print("shape of kpcmatrix: ")
    print(kpc.kpcmatrix['T'].shape[1])
    print("test OK")
