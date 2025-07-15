'''
:Copyright: 2016, EUMETSAT, All Rights Reserved.

This software was developed within the context of
the EUMETSAT Satellite Application Facility on
Numerical Weather Prediction (NWP SAF), under the
Cooperation Agreement dated 7 December 2016, between
EUMETSAT and the Met Office, UK, by one or more partners
within the NWP SAF. The partners in the NWP SAF are
the Met Office, ECMWF, DWD and MeteoFrance.

This module defines descriptors that will be used to control the access
to the Options, Profiles and Rttov classes attributes.
'''

from __future__ import absolute_import, print_function
import os
import numpy as np
from pyrttov.rttype import wrapfloat


class GenericDescriptorRO(object):
    '''A descriptor that handles Read-Only values of any type.'''

    _doc_extra = "This attribute is Read-Only."

    def __init__(self, category, target, doc='Undocumented attribute'):
        '''
        :param str category: name of the dictionary that will store the
            attribute value
        :param str target: name of the dictionary key that will be used to store
            the attribute value
        :param str doc: optional documentation string
        '''
        self._doc_update(doc)
        self._thedict = '_{}'.format(category)
        self._target = target

    def _doc_update(self, doc):
        self.__doc__ = doc + "\n\n" + self._doc_extra

    def __get__(self, instance, objtype=None):  # @UnusedVariable
        if instance is None:
            return self
        return getattr(instance, self._thedict).get(self._target, None)

    def __set__(self, instance, value):
        raise AttributeError('This is a Read-Only attribute')


class _GenericDescriptorRW(GenericDescriptorRO):
    '''An abstract class for Read/Write descriptors.'''

    _doc_extra = "This is a Read/Write attribute."

    def _target_init(self, instance, value):
        getattr(instance, self._thedict)[self._target] = value

    def __set__(self, instance, value):
        self._target_init(instance, value)


class _GenericDescriptorRWD(_GenericDescriptorRW):
    '''An abstract class for Read/Write/Delete descriptors.'''

    _doc_extra = "This is a Read/Write attribute. " \
                 "It can be emptied using `del`."

    def __delete__(self, instance):
        if self._target in getattr(instance, self._thedict):
            del getattr(instance, self._thedict)[self._target]


class TypedDescriptorRW(_GenericDescriptorRW):
    '''A descriptor that handles Read/Write values of a given type.'''

    _doc_extra = "This is a Read/Write attribute of type {!r}."

    def __init__(self, category, target, thetype,
                 doc='Undocumented attribute'):
        '''
        :param str category: name of the dictionary that will store the
            attribute value
        :param str target: name of the dictionary key that will be used to store
            the attribute value
        :param type thetype: type of the attribute
        :param str doc: optional documentation string
        '''
        self._type = thetype
        super(TypedDescriptorRW, self).__init__(category, target, doc)

    def _doc_update(self, doc):
        self.__doc__ = doc + "\n\n" + self._doc_extra.format(self._type)

    def _target_init(self, instance, value):
        if isinstance(value, (tuple, list)):
            getattr(instance, self._thedict)[
                self._target] = self._type(* value)
        else:
            getattr(instance, self._thedict)[self._target] = self._type(value)


class FilepathDescriptorRWD(_GenericDescriptorRWD):

    """A descriptor that handles a file's path (Read/Write)."""

    _doc_extra = "This is a Read/Write attribute that contains a valid " \
                 "pathname to a file. It can be emptied using `del`."

    def _target_init(self, instance, value):
        if not os.path.isfile(value):
            raise ValueError("%s does not exists or is not a file", value)
        super(FilepathDescriptorRWD, self)._target_init(instance, value)


class DirpathDescriptorRWD(_GenericDescriptorRWD):

    """A descriptor that handles a directory's path (Read/Write/Delete)."""

    _doc_extra = "This is a Read/Write attribute that contains a valid " \
                 "pathname to a directory. It can be emptied using `del`."

    def _target_init(self, instance, value):
        if not os.path.isdir(value):
            raise ValueError("%s does not exists or is not a directory", value)
        super(DirpathDescriptorRWD, self)._target_init(instance, value)


class _GenericNumpyRW(_GenericDescriptorRW):

    '''An abstract class that handles ``numpy.ndarray`` objects.'''

    _doc_extra = "This is a Read/Write attribute that holds a " \
                 "``numpy.ndarray`` of shape [{dim:s}]. When a new array is " \
                 "assigned to this attribute it will casted to a {dtype!r} "\
                 "array."

    def __init__(self, category, target, dtype=wrapfloat, dims='?',
                 doc='Undocumented attribute'):
        '''
        :param str category: name of the dictionary that will store the
                         attribute value
        :param str target: name of the dictionary key that will be used to store
                       the attribute value
        :param str dims: comma separated list of dimensions
        :param type dtype: desired typed for the ``numpy.ndarray``
        :param str doc: optional documentation string
        '''
        self._doc_dim = dims
        self._dtype = dtype
        super(_GenericNumpyRW, self).__init__(category, target, doc)

    def _doc_update(self, doc):
        self.__doc__ = (doc + "\n\n" +
                        self._doc_extra.format(dim=self._doc_dim,
                                               dtype=self._dtype))

    def _cast2dtype(self, value):
        # Try to convert the input to a ndarray
        if not isinstance(value, np.ndarray):
            try:
                npvalue = np.array(value, dtype=self._dtype)
            except(Exception):
                raise TypeError("Should be castable in a numpy.ndarray " +
                                "with type {:s}".format(self._dtype))
        # If it's already a ndarray, check the type...
        else:
            if value.dtype is not np.dtype(self._dtype):
                try:
                    npvalue = np.array(value, dtype=self._dtype)
                except(Exception):
                    raise TypeError("Should be a numpy.ndarray of type " +
                                    "{:s}. ".format(self._dtype) +
                                    "Or something compatible")
            else:
                npvalue = value
        return npvalue

    def _target_init(self, instance, value):
        getattr(instance, self._thedict)[
            self._target] = self._cast2dtype(value)


class GenericNumpyRWD(_GenericNumpyRW, _GenericDescriptorRWD):
    '''
    A descriptor that handles Read/Write/Delete values of a nprofiles/nlevels
    array.
    '''

    _doc_extra = (_GenericNumpyRW._doc_extra +
                  ' It can be emptied using `del`.')


class VerticalProfilesRW(_GenericNumpyRW):
    '''
    A descriptor that handles Read/Write values of a nprofiles/nlevels array.
    '''

    def __init__(self, category, target, dtype=wrapfloat, deltaNlevels=0,
                 doc='Undocumented attribute'):
        '''
        :param str category: name of the dictionary that will store the
            attribute value
        :param str target: name of the dictionary key that will be used to store
            the attribute value
        :param int deltaNlevels: size of the second dimension of the
            ``numpy.ndarray`` is nlevels + deltaNlevels
        :param type dtype: desired typed for the ``numpy.ndarray``
        :param str doc: optional documentation
        '''
        doc_dim = 'nprofiles,nlevels'
        if deltaNlevels < 0: doc_dim += '-'
        if deltaNlevels > 0: doc_dim += '+'
        if deltaNlevels != 0: doc_dim += str(abs(deltaNlevels))
        self._dNlevels = deltaNlevels
        super(VerticalProfilesRW, self).__init__(
            category, target, dtype=dtype, dims=doc_dim, doc=doc)

    def _target_init(self, instance, value):
        npvalue = self._cast2dtype(value)
        # Check the array's shape
        if npvalue.shape == (instance.Nprofiles, instance.Nlevels + self._dNlevels):
            getattr(instance, self._thedict)[self._target] = npvalue
        else:
            raise ValueError("Incorrect number of profiles and/or levels")


class VerticalProfilesRWD(VerticalProfilesRW, _GenericDescriptorRWD):
    '''
    A descriptor that handles Read/Write/Delete values of a nprofiles/nlevels
    array.
    '''

    _doc_extra = (VerticalProfilesRW._doc_extra +
                  ' It can be emptied using `del`.')


class ArbitraryProfilesRW(_GenericNumpyRW):
    '''
    A descriptor that handles Read/Write values of a 2D array with a
    nprofiles first dimension and an arbitrary leading dimension.
    '''

    def __init__(self, category, target, dtype=wrapfloat, leadingdim=1,
                 doc='Undocumented attribute'):
        '''
        :param str category: name of the dictionary that will store the
            attribute value
        :param str target: name of the dictionary key that will be used to store
            the attribute value
        :param int leadingdim: size of the leading dimension of the
            ``numpy.ndarray``
        :param type dtype: desired typed for the ``numpy.ndarray``
        :param str doc: optional documentation
        '''
        self._leadingdim = leadingdim
        if leadingdim > 1:
            doc_dim = 'nprofiles,{:d}'.format(leadingdim)
        else:
            doc_dim = 'nprofiles'
        super(ArbitraryProfilesRW, self).__init__(
            category, target, dtype=dtype, dims=doc_dim, doc=doc)

    def _target_init(self, instance, value):
        npvalue = self._cast2dtype(value)
        # Check the array's shape
        if self._leadingdim > 1:
            if npvalue.shape != (instance.Nprofiles, self._leadingdim):
                raise ValueError("Incorrect dimension")
        else:
            if npvalue.shape != (instance.Nprofiles,):
                raise ValueError("Incorrect dimension")
        getattr(instance, self._thedict)[self._target] = npvalue


class ArbitraryProfilesRWD(ArbitraryProfilesRW, _GenericDescriptorRWD):
    '''
    A descriptor that handles Read/Write/Del values of a 2D array with a
    nprofiles first dimension and an arbitrary leading dimension.
    '''

    _doc_extra = (ArbitraryProfilesRW._doc_extra +
                  ' It can be emptied using `del`.')


class PickGasesRO(GenericDescriptorRO):
    '''A descriptor that handles a particular Gas of a gases array.'''

    def __init__(self, callback, gas_id, doc='Undocumented attribute'):
        '''
        :param func callback: method that will be called to return the
            appropriate gas
        :param str doc: optional documentation
        '''
        self._doc_update(doc)
        self._callback = callback
        self._gas_id = gas_id

    def __get__(self, instance, objtype=None):  # @UnusedVariable
        if instance is None:
            return self
        return self._callback(instance, self._gas_id)
