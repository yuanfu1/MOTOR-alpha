'''
:Copyright: 2016, EUMETSAT, All Rights Reserved.

This software was developed within the context of
the EUMETSAT Satellite Application Facility on
Numerical Weather Prediction (NWP SAF), under the
Cooperation Agreement dated 7 December 2016, between
EUMETSAT and the Met Office, UK, by one or more partners
within the NWP SAF. The partners in the NWP SAF are
the Met Office, ECMWF, DWD and MeteoFrance.

This module defines various decorators that ease the class creation.
'''

from __future__ import absolute_import, print_function
from copy import deepcopy
from pyrttov.descriptor import GenericDescriptorRO
from pyrttov.descriptor import TypedDescriptorRW, VerticalProfilesRWD
from pyrttov.descriptor import PickGasesRO


def property_ro(func):
    """A variant of the property decorator.

    An extra message is added to the documentation.
    """
    func.__doc__ += "\n\nThis is a Read-Only attribute."
    return property(func)


def lock_attributes(obj):
    """Impose a tight control on the instance's attribute assignment.

    It will only be possible to assign an attribute under one of these
    three conditions:

    * If it's something hidden
    * If a descriptor has been defined
    * If a property has been defined

    """
    def thesetattr(self, name, value):
        '''We do not want any wild assignments.'''
        if not (name.startswith('_') or  # hidden values are ok
                # Is a descriptor defined ?
                (hasattr(type(self), name) and
                 isinstance(getattr(type(self), name), GenericDescriptorRO)) or
                # Is a property defined ?
                (hasattr(type(self), name) and
                 isinstance(getattr(type(self), name), property))
                ):
            msg = "You are not allowed to assign the {} attribute."
            raise AttributeError(msg.format(name))
        else:
            object.__setattr__(self, name, value)
    setattr(obj, '__setattr__', thesetattr)
    return obj


def add_descriptors_opts(categories):
    '''Class decorator that automatically adds typed descriptors for each
    option key.
    '''
    def thedecorator(cls):
        # Create the descriptors
        for category in categories:
            category_d = getattr(cls, '_defaults_{}'.format(category))
            for key, item in category_d.items():
                # Naming Convention...
                if (hasattr(cls, '_naming_override') and
                        key in cls._naming_override):
                    key_remap = cls._naming_override[key]
                else:
                    key_remap = ''.join([s[0].upper() + s[1:]
                                         for s in key.split('_')])
                    if key_remap.startswith('Add'):
                        key_remap = 'Add' + \
                            key_remap[3].upper() + key_remap[4:]
                # Set the descriptor
                dstr = "{}%{} option (default: {!s}).".format(category, key,
                                                              item)
                setattr(cls, key_remap,
                        TypedDescriptorRW(category, key, type(item), doc=dstr))

        # Create the initialisation routine that should be called
        # from the constructor
        def doInitDefaults(self):
            for category in categories:
                category_d = getattr(cls, '_defaults_{}'.format(category))
                setattr(self, '_{}'.format(category), deepcopy(category_d))
        cls._initDefaults = doInitDefaults
        return cls
    return thedecorator


def add_descriptors_gases2D(gases_map):
    '''Class decorator that automatically adds descriptors for each gas.

    Data are 2D ``numpy.ndarray`` of shape [nprofiles, nlevels].
    '''
    def thedecorator(cls):
        # Create the descriptors
        for gas_item, gas_desc in gases_map.items():
            # Set the descriptor
            setattr(cls, gas_desc[0],
                    VerticalProfilesRWD(
                        'profiles', gas_item,
                        doc="{} vertical profiles.".format(gas_desc[1]))
                    )
        return cls
    return thedecorator


def add_descriptors_gasesK(gases_map):
    '''Class decorator that automatically adds descriptors for each gas of the
    gases_k array.
    '''
    def thedecorator(cls):
        # Create the descriptors
        for gas_item, gas_desc in gases_map.items():
            docstr = ("{} Jacobian. ``numpy.ndarray`` of shape " +
                      "[nprofiles, nchannels, nlevels].").format(gas_desc[1])
            # Set the descriptor
            setattr(cls, gas_desc[0] + 'K',
                    PickGasesRO(cls.getItemK, gas_item, doc=docstr))
        return cls
    return thedecorator
