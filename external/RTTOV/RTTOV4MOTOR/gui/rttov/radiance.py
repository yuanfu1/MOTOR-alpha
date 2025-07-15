#! /usr/bin/env python
# -*- coding: utf-8 -*-


import h5py
from rttov.core import _V
from rttov.misc import Misc

radiance_list = ['BT', 'BT_CLEAR', 'TOTAL', 'CLEAR', 'REFL', 'REFL_CLEAR']


class Radiance(dict, _V):
    """
    The Radiance class
    Radiance object is Read-Only and contains minimum information
    """

    def __init__(self):
        dict.__init__(self)
        for key in radiance_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        self.misc = None

    def loadh5(self, h5):
        self.loadh5liste(h5, radiance_list)

    def read(self, filename):
        f = h5py.File(filename, 'r')

        # get the Dataset
        h5 = f['/RADIANCE/']
        self.loadh5(h5)

        # Load
        self.loadh5(h5)

        h5 = f['/MISC/']
        self.misc = Misc()
        self.misc.loadh5(h5)
        f.close()


if __name__ == '__main__':

    # Open file ReadOnly
    from core import rttov_gui_data_test_dir
    import os
    data_test_dir = rttov_gui_data_test_dir()
    ofile = os.path.join(data_test_dir, "radr.h5")
    f = h5py.File(ofile, 'r')

    # get the Dataset
    h5 = f['/RADIANCE/']

    # p1 is a rttov_hdf_mod radiance instance
    r1 = Radiance()

    # Load
    r1.loadh5(h5)

    # Close HDF file
    f.close()

    # Display option
    r1.display()

    # Save option is not available for radiances

    # Plots radiance

    import matplotlib.pyplot as plt

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    i = 0
    for v in ['BT', 'BT_CLEAR']:
        color = colors[i]
        i += 1
        plt.plot(r1[v], color + '-', label=v)
    plt.xlabel('channel (minus one)')
    plt.ylabel('Brightness Temperatures' +
               '   (' + str(r1[v + '_ATTRIBUTE']['UNITS']) + ')')
    plt.legend()
    plt.title(ofile)
    plt.show()

    for v in ['TOTAL', 'CLEAR']:
        color = colors[i]
        i += 1
        plt.plot(r1[v], color + '-', label=v)
    plt.xlabel('channel (minus one)')
    plt.ylabel('Radiances' + '  (' + str(r1[v + '_ATTRIBUTE']['UNITS']) + ')')
    plt.yscale('log')
    plt.legend()
    plt.title(ofile)
    plt.show()
    print("test OK")
