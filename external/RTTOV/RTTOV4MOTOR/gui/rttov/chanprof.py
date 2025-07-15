#! /usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import numpy
except ImportError:
    import sys
    sys.stderr.write('ERROR: numpy is not installed\n')
    sys.exit(1)


from rttov.core import _V


class Chanprof(dict, _V):
    """
    The Chanprof class
    A Python dictionary which contains a copy of the RTTOV v11
    rttov_chanprof structure
    Inherits the _V class for most methods
    Method:
    - setchanprof fills a chanprof class for a given number of channels,
        profiles array is set to 1
    """
    chanprof_list = ['CHANNELS', 'PROFILES']

    def __init__(self):
        dict.__init__(self)
        for key in self.chanprof_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def setChanprof(self, nchannels):
        """Creates Chanprof arrays for a given number of channels
        CHANNELS are from 1 to nchannels
        PROFILES are all ones
        """
        self['CHANNELS'] = numpy.arange(nchannels)
        self['CHANNELS'] += 1
        self['PROFILES'] = numpy.ones(nchannels)


if __name__ == '__main__':

    nchannels = 4

    chanprof = Chanprof()
    chanprof.setChanprof(nchannels)
    chanprof.display()
    print("test OK")
