#! /usr/bin/env python
# -*- coding: utf-8 -*-

from rttov.core import _V

transmission_list = ['TAU_TOTAL', 'TAU_LEVELS']


class Transmission(dict, _V):
    """
    The Transmission class
    transmission object is Read-Only
    """

    def __init__(self):
        dict.__init__(self)
        for key in transmission_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def loadh5(self, h5):
        self.loadh5liste(h5, transmission_list)


if __name__ == '__main__':
    import os
    import h5py
    from core import rttov_gui_data_test_dir
    data_test_dir = rttov_gui_data_test_dir()

    # Open file ReadOnly

    file = os.path.join(data_test_dir, "trns.h5")
    f = h5py.File(file, 'r')

    # get the Dataset
    h5 = f['/TRANSMISSION/']

    # p1 is a rttov_hdf_mod transmission instance
    r1 = Transmission()

    # Load
    r1.loadh5(h5)

    # Close HDF file
    f.close()

    # Display option
    r1.display()

    # Save option isnot available for transmissions

    # Plots transmission
    import matplotlib
    matplotlib.use("wxagg")
    import matplotlib.pyplot as plt
    transmission_list = ['TAU_TOTAL', 'TAU_LEVELS']
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

    plt.plot(r1['TAU_TOTAL'], 'r-', label='TAU_TOTAL')
    plt.plot(r1['TAU_LEVELS'][:, 10], 'b-', label='TAU_LEVELS lev=10')
    plt.plot(r1['TAU_LEVELS'][:, 50], 'g-', label='TAU_LEVELS lev=50')
    plt.plot(r1['TAU_LEVELS'][:, 75], 'c-', label='TAU_LEVELS lev=75')
    plt.xlabel('channel (minus one)')
    plt.ylabel('Transmission' +
               '   (' + r1['TAU_TOTAL_ATTRIBUTE']['UNITS'] + ')')
    plt.legend()
    plt.title(file)
    plt.show()
    print("test OK")
