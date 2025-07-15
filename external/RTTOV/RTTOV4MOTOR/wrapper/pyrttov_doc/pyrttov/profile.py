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
from pprint import pformat

from pyrttov.rttype import wrapint, wrapfloat, itemIdType, gasUnitType, \
                           maxuseraer
from pyrttov.decorator import add_descriptors_gases2D
from pyrttov.decorator import property_ro, lock_attributes
from pyrttov.descriptor import TypedDescriptorRW
from pyrttov.descriptor import ArbitraryProfilesRW, ArbitraryProfilesRWD
from pyrttov.descriptor import VerticalProfilesRW, VerticalProfilesRWD


# Gives the mapping between gas_id / (descriptor name, full description)
ItemDescriptorNaming = {'Q': ('Q', 'Water Vapor (q)'),
                        'O3': ('O3', 'Ozone (O3)'),
                        'CO2': ('CO2', 'Carbon dioxide (CO2)'),
                        'N2O': ('N2O', 'Nitrous oxide (N2O)'),
                        'CO': ('CO', 'Carbon monoxide (CO)'),
                        'CH4': ('CH4', 'Methane (CH4)'),
                        'SO2': ('SO2', 'Sulphur dioxide (SO2)'),
                        'CLW': ('CLW', 'Cloud Liquid Water as absorber'),
                        'CFRAC': ('Cfrac', 'Cloud fraction'),
                        'STCO': ('Stco', 'Cloud Liquid Water Type 1 (STCO)'),
                        'STMA': ('Stma', 'Cloud Liquid Water Type 2 (STMA)'),
                        'CUCC': ('Cucc', 'Cloud Liquid Water Type 3 (CUCC)'),
                        'CUCP': ('Cucp', 'Cloud Liquid Water Type 4 (CUCP)'),
                        'CUMA': ('Cuma', 'Cloud Liquid Water Type 5 (CUMA)'),
                        'CIRR': ('Cirr', 'Ice cloud'),
                        'ICEDE': ('Icede', 'Ice cloud effective diameter'),
                        'CLWDE': ('Clwde', 'Liquid cloud effective diameter'),
                        'INSO': ('Inso', 'OPAC aerosol particle type 1 (INSO)'),
                        'WASO': ('Waso', 'OPAC aerosol particle type 2 (WASO)'),
                        'SOOT': ('Soot', 'OPAC aerosol particle type 3 (SOOT)'),
                        'SSAM': ('Ssam', 'OPAC aerosol particle type 4 (SSAM)'),
                        'SSCM': ('Sscm', 'OPAC aerosol particle type 5 (SSCM)'),
                        'MINM': ('Minm', 'OPAC aerosol particle type 6 (MINM)'),
                        'MIAM': ('Miam', 'OPAC aerosol particle type 7 (MIAM)'),
                        'MICM': ('Micm', 'OPAC aerosol particle type 8 (MICM)'),
                        'MITR': ('Mitr', 'OPAC aerosol particle type 9 (MITR)'),
                        'SUSO': ('Suso', 'OPAC aerosol particle type 10 (SUSO)'),
                        'VOLA': ('Vola', 'OPAC aerosol particle type 11 (VOLA)'),
                        'VAPO': ('Vapo', 'OPAC aerosol particle type 12 (VAPO)'),
                        'ASDU': ('Asdu', 'OPAC aerosol particle type 13 (ASDU)'),
                        'BCAR': ('Bcar', 'CAMS aerosol particle type 1 (BCAR)'),
                        'DUS1': ('Dus1', 'CAMS aerosol particle type 2 (DUS1)'),
                        'DUS2': ('Dus2', 'CAMS aerosol particle type 3 (DUS2)'),
                        'DUS3': ('Dus3', 'CAMS aerosol particle type 4 (DUS3)'),
                        'SULP': ('Sulp', 'CAMS aerosol particle type 5 (SULP)'),
                        'SSA1': ('Ssa1', 'CAMS aerosol particle type 6 (SSA1)'),
                        'SSA2': ('Ssa2', 'CAMS aerosol particle type 7 (SSA2)'),
                        'SSA3': ('Ssa3', 'CAMS aerosol particle type 8 (SSA3)'),
                        'OMAT': ('Omat', 'CAMS aerosol particle type 9 (OMAT)')}

for i in range(1, maxuseraer+1):
    ItemDescriptorNaming['AER'+str(i)] = ('Aer'+str(i), 'User aerosol particle type ' + str(i))

@lock_attributes
class _ProfilesCommon(object):
    '''Factorise some common code for the Profiles like classes.'''

    _GASES_DESCRIPTION = {}
    _MANDATORY_GASES = ('Q', )
    _PROFILES_PRINT_LIST = ('DateTimes', 'Angles', 'SurfGeom', 'SurfType',
                            'Skin', 'S2m', 'Zeeman')

    def __init__(self, nprofiles, nlevels):
        """
        :param int nlevels: number of vertical levels of the profiles
        :param int nprofiles: number of profiles in this batch
        """
        self._nlevels = nlevels
        self._nprofiles = nprofiles
        self._internals = {}
        self._profiles = {}
        self._gases_eraser()
        self._conf = {'gas_units': gasUnitType('kg_per_kg')}

    @property_ro
    def Nlevels(self):
        """Number of vertical levels in a profile."""
        return self._nlevels

    @property_ro
    def Nprofiles(self):
        """Number of profiles in this Profiles object."""
        return self._nprofiles

    GasUnits = TypedDescriptorRW('conf', 'gas_units', gasUnitType,
                                 doc='Unit used in the gas arrays.')
    DateTimes = ArbitraryProfilesRW('internals', 'datetimes', dtype=wrapint,
                                    leadingdim=6, doc="Year, month, day, " +
                                    "hour, minute, second.")
    SurfGeom = ArbitraryProfilesRW('internals', 'surfgeom', leadingdim=3,
                                   doc='Latitude, longitude, elevation.')

    Zeeman = ArbitraryProfilesRWD('internals', 'zeeman', leadingdim=2,
                                  doc="Zeeman effect scheme parameters " +
                                  "(Be, cosbk).")

    P = VerticalProfilesRWD('profiles', 'P', doc="Pressure vertical profiles.")
    T = VerticalProfilesRW('profiles', 'T', doc="Temperature vertical profiles.")

    @property_ro
    def DefaultPressureLevels(self):
        """Whether or not the default pressure levels are used."""
        return self.P is None

    def _actual_check(self):
        return (self.T is not None and
                self._gases_checker() and
                self.DateTimes is not None and
                self.Angles is not None and
                self.SurfType is not None and
                self.SurfGeom is not None and
                self.S2m is not None and
                self.Skin is not None)

    def check(self):
        '''
        Check that the profile(s) are correct.

        Checks if all mandatory fields have been provided, but does not
        perform a check upon the values (this is done within RTTOV
        itself); if simplecloud, clwscheme, icecloud or zeeman have not been set
        initialise them with default values
        '''
        complete = self._actual_check()
        if complete:
            self._optional_fields_init()
        return complete

    def _buildGasIdArray(self):
        '''Build the gas_id array from individual gas arrays'''
        self._gas_id = np.array([itemIdType(name)
                                 for name in self._GASES_DESCRIPTION.keys()
                                 if name in self._profiles],
                                dtype=wrapint)
        self._gas_id.sort()

    @property_ro
    def Ngases(self):
        """Number of gases in the gases array."""
        if self._gases_manual:
            return self._gases.shape[0]
        else:
            self._buildGasIdArray()
            return len(self._gas_id)

    def _gas_id_setter(self, ids):
        if not self._gases_manual:
            raise ValueError("The gases should be set first")
        npids = np.array(ids, dtype=np.int32)
        if npids.shape != (self.Gases.shape[0],):
            raise ValueError("Incorrect shape for gas_id; " +
                             "({:d},) expected.".format(self.Gases.shape[0]))
        self._gas_id = npids

    def _gas_id_getter(self):
        if self._gases_manual:
            return self._gas_id
        else:
            self._buildGasIdArray()
            return self._gas_id

    GasId = property(_gas_id_getter, _gas_id_setter,
                     doc="Gas Id list.\n\n" +
                         "This is a Read/Write attribute.")

    def _buildGasesArray(self):
        '''Build the gases array from individual gas arrays'''
        self._gases = np.empty((self.Ngases, self.Nprofiles, self.Nlevels),
                               dtype=wrapfloat)
        for i, gasid in enumerate(self.GasId):
            self._gases[i, ...] = self._profiles[itemIdType(gasid).name]

    def _gases_eraser(self):
        self._gases = None
        self._gas_id = np.array([], dtype=wrapint)
        self._gases_manual = False

    def _gases_setter(self, gases):
        # Check that gases is a Numpy array
        if isinstance(gases, np.ndarray):
            if gases.dtype is not wrapfloat:
                npgases = np.array(gases, dtype=wrapfloat)
            else:
                npgases = gases
        else:
            raise ValueError("A dtyped float64 numpy.ndarray is expected")
        # Check dimensions
        if not (npgases.shape[1] == self.Nprofiles and
                npgases.shape[2] == self.Nlevels):
            raise ValueError("Incorrect dimensions")
        # Reset the gas_id array if it doesn't conform
        if self._gases_manual:
            if len(self._gas_id) != npgases.shape[0]:
                self._gas_id = []
        else:
            self._gas_id = []
        # Store _gases
        self._gases_manual = True
        self._gases = npgases

    def _gases_getter(self):
        if self._gases_manual:
            return self._gases
        else:
            self._buildGasesArray()
            return self._gases

    Gases = property(_gases_getter, _gases_setter, _gases_eraser,
                     doc="Gases array.\n\n" +
                         "This is a Read/Write attribute. " +
                         "It can be emptied using `del`.")

    def _gases_checker(self):
        if self._gases_manual:
            return all([itemIdType(g) in self.GasId
                        for g in self._MANDATORY_GASES])
        else:
            return all([getattr(self, g) is not None
                        for g in self._MANDATORY_GASES])

    def _optional_fields_init(self):
        if self.Zeeman is None:
            self.Zeeman = np.zeros((self.Nprofiles, 2), dtype=wrapfloat)

    def printProfiles(self, pfunction=print):
        """Print the content of this :data:`Profiles` instance.

        To the notable exception of the gases array (see :meth:`printGases`).
        """
        pfunction("gas_units={:d}".format(self.GasUnits))
        for prof in range(self.Nprofiles):
            pstring = '#{:<5d}'.format(prof)
            pfunction(">>>>>> profile {}".format(pstring))
            for stuff in self._PROFILES_PRINT_LIST:
                if getattr(self, stuff) is not None:
                    thestuff = pformat(getattr(self, stuff)[prof, ...])
                    pfunction('{} {}: {!s}'.format(pstring, stuff, thestuff))

    def printGases(self, pfunction=print):
        """Print the gases array content on the standard output."""
        if self.T is None:
            return
        pfunction(">> Gases Ids and Names:")
        if self.DefaultPressureLevels:
            head = '{:10s} {:6s} '.format('Level', 'T')
            fmt = '{:10d} {:6.2f} ' + ('{:10.4e} ' * self.Ngases)
        else:
            head = '{:10s} {:6s} '.format('P', 'T')
            fmt = '{:10.4e} {:6.2f} ' + ('{:10.4e} ' * self.Ngases)
        for i, gasid in enumerate(self.GasId):
            gname = itemIdType(gasid).name
            head += '{:10s} '.format(gname)
            pfunction(('Gas #{:<3d}, GasId: {:<3d}, ' +
                       'GasName: {:s}').format(i, gasid, gname))
        for prof in range(self.Nprofiles):
            pfunction(">> Gases array for profile #{:d}".format(prof))
            pfunction(head)
            for lev in range(self.Nlevels):
                pfunction(fmt.format(self.P[prof, lev] if self.P is not None
                                     else lev,
                                     self.T[prof, lev],
                                     * list(self.Gases[:, prof, lev])))


@add_descriptors_gases2D(ItemDescriptorNaming)
class Profiles(_ProfilesCommon):
    '''The Profiles class holds a batch of RTTOV profiles.

    Two mechanisms are offered to initialise the gases array:
      * All the gases are initialised at once using the :data:`Gases`
        attribute. Then, the gas Id list have to be provided to the
        :data:`GasId` attribute.
      * Each gas can be provided independently using the appropriate attribute
        (:data:`Q`, :data:`CO2`, ...). Accordingly, the :data:`Gases` and
        :data:`GasId` attributes will be automatically generated.

    Whenever the :data:`Gases` attribute is set manually, it takes precedence
    over individual gas attributes that may already be defined.

    The :data:`SimpleCloud`, :data:`ClwScheme` and :data:`IceCloud` and
    :data:`Zeeman` attributes have default values.

    '''

    _GASES_DESCRIPTION = ItemDescriptorNaming
    _PROFILES_PRINT_LIST = ('DateTimes', 'Angles', 'SurfGeom', 'SurfType',
                            'Skin', 'S2m', 'SimpleCloud', 'ClwScheme',
                            'IceCloud', 'Zeeman')

    def __init__(self, nprofiles, nlevels):
        """
        :param int nlevels: number of vertical levels of the profiles
        :param int nprofiles: number of profiles in this batch
        """
        super(Profiles, self).__init__(nprofiles, nlevels)
        self._conf['mmr_cldaer'] = True
    Angles = ArbitraryProfilesRW('internals', 'angles', leadingdim=4,
                                 doc="Satellite zenith, satellite azimuth, " +
                                 "sun zenith, sun azimuth.")
    SurfType = ArbitraryProfilesRW('internals', 'surftype', dtype=wrapint,
                                   leadingdim=2,
                                   doc="Description of the surface type " +
                                   "(surftype, watertype).")
    S2m = ArbitraryProfilesRW('internals', 's2m', leadingdim=6,
                              doc="2m p, 2m t, 2m q, 10m wind u, v, wind fetch.")
    Skin = ArbitraryProfilesRW('internals', 'skin', leadingdim=9,
                               doc="Skin T, salinity, snow_fraction, " +
                               "foam_fraction, fastem_coefs(1:5).")
    MmrCldAer = TypedDescriptorRW('conf', 'mmr_cldaer', bool,
                                  doc="Unit used in the cloud/aerosol arrays " +
                                  "True: kg/kg, False: g/m^3 (cld), cm^-3 (aer).")
    ClwScheme = ArbitraryProfilesRW('internals', 'clwscheme', dtype=wrapint,
                                    leadingdim=2,
                                    doc="VIS/IR cloud liquid water scheme " +
                                    "parameters (clw_scheme, clwde_param).")
    IceCloud = ArbitraryProfilesRWD('internals', 'icecloud', dtype=wrapint,
                                    leadingdim=2,
                                    doc="VIS/IR ice cloud scheme parameters " +
                                    "parameters (ice_scheme, icede_param).")
    SimpleCloud = ArbitraryProfilesRWD('internals', 'simplecloud',
                                       leadingdim=2,
                                       doc="Simple cloud scheme parameters " +
                                       "parameters (ctp, cfraction).")

    def setUserAerN(self, aer, n):
        """Set profile for user aerosol species n"""
        if n < 0 or n > maxuseraer:
            raise ValueError("Aerosol index must lie between 1 and "+str(maxuseraer))
        setattr(self, 'Aer'+str(n), aer)

    def delUserAerN(self, n):
        """Delete profile for user aerosol species n"""
        if n < 0 or n > maxuseraer:
            raise ValueError("Aerosol index must lie between 1 and "+str(maxuseraer))
        delattr(self, 'Aer'+str(n))

    def _gases_checker(self):
        if self._gases_manual:
            return itemIdType('Q') in self.GasId
        else:
            return self.Q is not None

    def _optional_fields_init(self):
        super(Profiles, self)._optional_fields_init()
        if self.SimpleCloud is None:
            self.SimpleCloud = np.zeros((self.Nprofiles, 2), dtype=wrapfloat)
        if self.ClwScheme is None:
            self.ClwScheme = np.ones((self.Nprofiles, 2), dtype=wrapint)
        if self.IceCloud is None:
            self.IceCloud = np.ones((self.Nprofiles, 2), dtype=wrapint)
