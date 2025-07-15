#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy

from rttov.core import _V, loadDset_str, loadDset_arr, getAttributes
from rttov.core import checkAttributes


class Misc(dict, _V):
    """
    The Misc class
    This class contains miscellaneous information such as
    the coefficient file name
    but also the channel wavenumbers and different ways to name
    the satellite/instrument
    """
    misc_list_str = ['SATELLITE', 'INSTRUMENT', 'ID_COMMON_NAME',
                     'COEF_FILENAME',
                     'SCAT_COEF_FILENAME', 'CLOUD_COEF_FILENAME',
                     'PC_COEF_FILENAME']
    misc_list = misc_list_str + ['WAVENUMBERS']

    def __init__(self):
        dict.__init__(self)
        for key in self.misc_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def loadh5(self, h5):
        for cle in list(h5):

            if(cle in self.misc_list_str):
                # print cle
                self[cle] = loadDset_str(h5[cle])
                self[cle + '_ATTRIBUTE'] = getAttributes(h5[cle])
                checkAttributes(self[cle + '_ATTRIBUTE'])
            else:
                self[cle] = loadDset_arr(h5[cle])
                self[cle + '_ATTRIBUTE'] = getAttributes(h5[cle])
                checkAttributes(self[cle + '_ATTRIBUTE'])


if __name__ == '__main__':

    print("see kmatrix or radiance or transmission... for unity tests example")
