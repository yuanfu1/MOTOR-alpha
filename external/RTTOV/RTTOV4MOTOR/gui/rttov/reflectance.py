#! /usr/bin/env python
# -*- coding: utf-8 -*-

import numpy

from rttov.core import _V


class Reflectance(dict, _V):
    """
    The Reflectance class
    A Python dictionary which contains a copy of the RTTOV v13
    rttov_reflectance structure
    Inherits the _V class for most methods
    Method:
    - setreflectance fills a Reflectance class for a given number of channels,
        reflectancies are set to 0 and logical calcrefl is set to True
    """
    reflectance_list = ['REFL_IN', 'REFL_OUT',
                        'DIFFUSE_REFL_IN', 'DIFFUSE_REFL_OUT',
                        'REFL_CLOUD_TOP', 'CALCREFL']
    list_arr_logical = ['CALCREFL']

    def __init__(self):
        dict.__init__(self)
        for key in self.reflectance_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def setReflectance(self, nchannels):
        """Creates reflectance arrays for a given number of channels
        REFL_IN and REFL_OUT are all zeros
        """
        self['REFL_IN'] = numpy.zeros(nchannels)
        self['REFL_OUT'] = numpy.zeros(nchannels)
        self['DIFFUSE_REFL_IN'] = numpy.zeros(nchannels)
        self['DIFFUSE_REFL_OUT'] = numpy.zeros(nchannels)
        self['REFL_CLOUD_TOP'] = numpy.zeros(nchannels)
        self['CALCREFL'] = numpy.empty(nchannels, bool)
        self['CALCREFL'][:] = True
        self['REFL_IN_ATTRIBUTE']['LBOUND'] = 1
        self['REFL_IN_ATTRIBUTE']['UBOUND'] = nchannels
        self['REFL_OUT_ATTRIBUTE']['LBOUND'] = 1
        self['REFL_OUT_ATTRIBUTE']['UBOUND'] = nchannels
        self['DIFFUSE_REFL_IN_ATTRIBUTE']['LBOUND'] = 1
        self['DIFFUSE_REFL_IN_ATTRIBUTE']['UBOUND'] = nchannels
        self['DIFFUSE_REFL_OUT_ATTRIBUTE']['LBOUND'] = 1
        self['DIFFUSE_REFL_OUT_ATTRIBUTE']['UBOUND'] = nchannels
        self['REFL_CLOUD_TOP_ATTRIBUTE']['LBOUND'] = 1
        self['REFL_CLOUD_TOP_ATTRIBUTE']['UBOUND'] = nchannels
        self['CALCREFL_ATTRIBUTE']['LBOUND'] = 1
        self['CALCREFL_ATTRIBUTE']['UBOUND'] = nchannels


if __name__ == '__main__':

    print("see emissivity.py test")
