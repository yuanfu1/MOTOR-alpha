'''
:Copyright: 2016, EUMETSAT, All Rights Reserved.

This software was developed within the context of
the EUMETSAT Satellite Application Facility on
Numerical Weather Prediction (NWP SAF), under the
Cooperation Agreement dated 7 December 2016, between
EUMETSAT and the Met Office, UK, by one or more partners
within the NWP SAF. The partners in the NWP SAF are
the Met Office, ECMWF, DWD and MeteoFrance.
'''

from __future__ import absolute_import, print_function
import numpy as np

from pyrttov.rttype import wrapint, wrapfloat, maxscatthydro
from pyrttov.decorator import add_descriptors_gases2D
from pyrttov.descriptor import TypedDescriptorRW
from pyrttov.descriptor import ArbitraryProfilesRW, ArbitraryProfilesRWD
from pyrttov.descriptor import VerticalProfilesRW
from pyrttov.profile import _ProfilesCommon

# Gives the mapping between gas_id / (descriptor name, full description)
ItemDescriptorNamingScatt = {'Q': ('Q', 'Water Vapor (q)'),
                             'O3': ('O3', 'Ozone (O3)'),
                             'SCATT_HYDRO_FRAC': ('HydroFrac', 'Cloud fraction'),
                             'SCATT_CLW': ('Clw', 'Cloud Liquid Water'),
                             'SCATT_CIW': ('Ciw', 'Cloud Ice Water'),
                             'SCATT_RAIN': ('Rain', 'Rain'),
                             'SCATT_SNOW': ('Snow', 'Solid precip (snow)'),
                             'SCATT_GRAUPEL': ('Graupel', 'Graupel')}

for i in range(1, maxscatthydro+1):
    ItemDescriptorNamingScatt['HYDRO'+str(i)] = ('Hydro'+str(i), 'Hydrometeor ' + str(i))

for i in range(1, maxscatthydro+1):
    ItemDescriptorNamingScatt['HYDRO_FRAC'+str(i)] = ('HydroFrac'+str(i), 'Hydrometeor cloud fraction ' + str(i))

@add_descriptors_gases2D(ItemDescriptorNamingScatt)
class ProfilesScatt(_ProfilesCommon):
    '''The ProfilesScatt class holds a batch of RTTOV-SCATT profiles.

    Two mechanisms are offered to initialise the gases array:
      * All the gases are initialised at once using the :data:`Gases`
        attribute. Then, the gas Id list have to be provided to the
        :data:`GasId` attribute.
      * Each gas/hydrometeor can be provided independently using the appropriate attribute
        (:data:`Q`, :data:`CLW`, ...). Accordingly, the :data:`Gases` and
        :data:`GasId` attributes will be automatically generated.

    Whenever the :data:`Gases` attribute is set manually, it takes precedence
    over individual gas attributes that may already be defined.

    The :data:`Zeeman` attribute has default values.

    '''

    _GASES_DESCRIPTION = ItemDescriptorNamingScatt
    _PROFILES_PRINT_LIST = ('DateTimes', 'Angles', 'SurfGeom', 'SurfType',
                            'Skin', 'S2m', 'Zeeman', 'UserCfrac')

    def __init__(self, nprofiles, nlevels):
        """
        :param int nlevels: number of vertical levels of the profiles
        :param int nprofiles: number of profiles in this batch
        """
        super(ProfilesScatt, self).__init__(nprofiles, nlevels)

    Angles = ArbitraryProfilesRW('internals', 'angles', leadingdim=2,
                                 doc="Satellite zenith, satellite azimuth.")
    SurfType = ArbitraryProfilesRW('internals', 'surftype', dtype=wrapint,
                                   doc="Surface type.")
    S2m = ArbitraryProfilesRW('internals', 's2m', leadingdim=5,
                              doc="2m p, 2m t, 2m q, 10m wind u, v.")
    Skin = ArbitraryProfilesRW('internals', 'skin', leadingdim=8,
                               doc="Skin T, salinity, " +
                               "foam_fraction, fastem_coefs(1:5).")
    UserCfrac = ArbitraryProfilesRWD('profiles', 'UserCfrac',
                                     doc="User-specified cloud fraction per profile if luser_cfrac is true.")
    P = VerticalProfilesRW('profiles', 'P', doc="Pressure vertical profiles.")
    Ph = VerticalProfilesRW('profiles', 'Ph', deltaNlevels=1,
                            doc="Pressure half-level vertical profiles.")

    def setHydroN(self, hydro, n):
        """Set profile for hydrometeor n"""
        if n < 0 or n > maxscatthydro:
            raise ValueError("Hydrometeor index must lie between 1 and "+str(maxscatthydro))
        setattr(self, 'Hydro'+str(n), hydro)

    def delUserHydroN(self, n):
        """Delete profile for user hydrometeor n"""
        if n < 0 or n > maxscatthydro:
            raise ValueError("Hydrometeor index must lie between 1 and "+str(maxscatthydro))
        delattr(self, 'Hydro'+str(n))

    def setHydroFracN(self, hydrofrac, n):
        """Set profile for hydrometeor cloud fraction n"""
        if n < 0 or n > maxscatthydro:
            raise ValueError("Hydrometeor index must lie between 1 and "+str(maxscatthydro))
        setattr(self, 'HydroFrac'+str(n), hydrofrac)

    def delUserHydroFracN(self, n):
        """Delete profile for user hydrometeor cloud fraction n"""
        if n < 0 or n > maxscatthydro:
            raise ValueError("Hydrometeor index must lie between 1 and "+str(maxscatthydro))
        delattr(self, 'HydroFrac'+str(n))

    def _actual_check(self):
        return (super(ProfilesScatt, self)._actual_check() and
                self.P is not None and
                self.Ph is not None)

    def _optional_fields_init(self):
        super(ProfilesScatt, self)._optional_fields_init()
        if self.UserCfrac is None:
            self.UserCfrac = np.zeros((self.Nprofiles,), dtype=wrapfloat)
