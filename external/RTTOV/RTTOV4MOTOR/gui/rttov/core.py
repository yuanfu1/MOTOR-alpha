#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Core functions and classes for rttov fortran interface

Copyright:
    This software was developed within the context of
    the EUMETSAT Satellite Application Facility on
    Numerical Weather Prediction (NWP SAF), under the
    Cooperation Agreement dated 7 December 2016, between
    EUMETSAT and the Met Office, UK, by one or more partners
    within the NWP SAF. The partners in the NWP SAF are
    the Met Office, ECMWF, DWD and MeteoFrance.

    Copyright 2012, EUMETSAT, All Rights Reserved.
P. Brunel
"""

try:
    import numpy
except ImportError:
    import sys
    sys.stderr.write('ERROR: numpy is not installed\n')
    sys.exit(1)
import re


def rttov_gui_data_test_dir():
    import sys
    import os
    print(os.environ["RTTOV_GUI_PREFIX"])
    try:
        prefix = os.environ["RTTOV_GUI_PREFIX"]
    except Exception:
        print("RTTOV_GUI_PREFIX is not defined")
        sys.exit(1)
    return os.path.join(prefix, "rttov_gui_data_test")


def getAttributes(dset):
    """
    Returns all attributes of the dataset
    """
    attr = {}
    for name, value in dset.attrs.items():
        attr[name] = value
    return attr


def loadDset_arr(dset):
    """
    Loads a dataset
    returns an array except if the number of dimensions is one then
    it returns a scalar
    """
    if(dset.shape):
        val = numpy.array(dset)
        if(dset.shape[0] == 1):
            return val[0]
        else:
            return val
    else:
        return dset[()]


def loadDset_log(dset):
    """
    Loads a dataset containing a logical scalar
    """
    val = numpy.array(dset)
    return (val == 1)


def loadDset_str(dset):
    """
    Loads a dataset containing a string array
    """

    val = str(dset[()]).replace("b'", "")
    val = val.replace("'", "")
    return val


def saveAttributes(dset, attr):
    """
    Saves all attributes of the dataset dset
    """
    for name, value in attr.items():
        dset.attrs[name] = value


def checkAttributes(x):
    """
    Verifies that the COMMENT and UNITS attributes are part of the variable x
    If not it fills with defaults 'none' and 'n/a'
    """
    if ('COMMENT' not in x):
        x['COMMENT'] = 'none'
    if ('UNITS' not in x):
        x['UNITS'] = 'n/a'


def cleanpath(x):
    return re.sub('/$', '', re.sub('^/', '', re.sub('/+', '/', x, count=0)))


class _V(object):
    """
    This is a virtual class which is used by default by other public classes
    It provides methods for
    - load a class from HDF5: loadh5(self, h5)
    - save a class to an HDF5: saveh5(self, h5, path)
    - display a class content to terminal: display(self)
    """

    list_arr_logical = ['']

    def loadh5(self, h5, *listearg):
        """
        Load is driven by the file content
        Classes for which we have a list_arr_logical item have a special
        treatment
        logicals are converted to integers for HDF5 storage
        """
        if(listearg):
            liste = listearg
        else:
            liste = list(h5)
        for cle in liste:
            if(cle in self.list_arr_logical):
                val = loadDset_arr(h5[cle])
                self[cle] = numpy.where(val == 1, True, False)
            else:
                self[cle] = loadDset_arr(h5[cle])
            self[cle + '_ATTRIBUTE'] = getAttributes(h5[cle])
            checkAttributes(self[cle + '_ATTRIBUTE'])

    def loadh5liste(self, h5, liste):
        """
        Load is driven by the liste argument
        this is usefull for Transmission and Radiance
        for which we do not want to load all the members
        Classes for which we have a list_arr_logical item have a special
        treatment
        logicals are converted to integers for HDF5 storage
        """
        for cle in liste:
            if(cle in self.list_arr_logical):
                val = loadDset_arr(h5[cle])
                self[cle] = numpy.where(val == 1, True, False)
            else:
                self[cle] = loadDset_arr(h5[cle])
            self[cle + '_ATTRIBUTE'] = getAttributes(h5[cle])
            checkAttributes(self[cle + '_ATTRIBUTE'])

    def saveh5(self, h5, path):
        """
        Save a class to an HDF5 file h5 under path directory
        Attributes are saved at the same time datasets are saved
        """
        list_of_names = []
        h5.visit(list_of_names.append)
        for cle in list(self.keys()):

            if(re.search('_ATTRIBUTE', cle) or self[cle] is None):
                continue
            cpath = cleanpath(path + cle)
            if cpath in list_of_names:
                del h5[cpath]
            if(cle in self.list_arr_logical):
                val = numpy.where(self[cle], 1, 0)
                dset = h5.create_dataset(cpath, data=val)
            else:
                dset = h5.create_dataset(cpath, data=self[cle])
            saveAttributes(h5[cpath], self[cle + '_ATTRIBUTE'])

    def display(self):
        """
        Display a full class content to terminal
        Members are printed by alphabetic order
        """
        for cle in sorted(self.keys()):
            print(cle, '\t =>\t', self[cle], type(self[cle]))

    def print_value(self):
        """
        Display a full class content to terminal
        Members are printed by alphabetic order
        """
        for cle in sorted(self.keys()):
            if not (cle[-9:] == "ATTRIBUTE"):
                print(cle, '\t =>\t', self[cle])


if __name__ == '__main__':
    print("no unity test fore core functions, see other classes")
