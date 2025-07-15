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
import sys
import copy
import collections
import numpy as np

import rttov_wrapper_f2py as rtwrap

from pyrttov.option import Options, DOM_IR, DOM_VIS
from pyrttov.profile import Profiles
from pyrttov.profile import ItemDescriptorNaming as _ItemDescriptorNaming
from pyrttov.profilescatt import ProfilesScatt
from pyrttov.profilescatt import ItemDescriptorNamingScatt as _ItemDescriptorNamingScatt
from pyrttov.decorator import property_ro, lock_attributes, add_descriptors_gasesK
from pyrttov.descriptor import FilepathDescriptorRWD as _FPathDescRWD
from pyrttov.descriptor import DirpathDescriptorRWD as _DPathDescRWD
from pyrttov.descriptor import GenericDescriptorRO as _DescRO
from pyrttov.descriptor import TypedDescriptorRW as _TDescRW
from pyrttov.descriptor import GenericNumpyRWD as _ArrayRWD
from pyrttov.rttype import wrapfloat, wrapint, itemIdType, maxuseraer, maxscatthydro


# No automatic exports
__all__ = []


class RttovError(Exception):
    pass


class RttovInternalError(RttovError):
    pass


class RttovWrapperError(RttovError):
    pass


@lock_attributes
class _RttovCommons(object):

    '''Elements that are common to all Rttov like classes.'''

    _PROFILES_CLASS = Profiles

    def __init__(self):
        self._coefs_conf = {}
        self._nchannels = 0
        self._instid = -1
        self._debug = False
        self._options = Options()
        self._profiles = None
        self._bBasics = dict()
        self._bRad = dict()

    def __del__(self):
        if self._instid > 0:
            self._printVerbWrap("Deallocating this inst_id.")
            err = rtwrap.rttov_drop_inst(self._instid)
            self._errCheck(err, 'Error in rttov_drop_inst')

    def _printDbg(self, msg):
        if self._debug:
            sys.stderr.write("DEBUG: " + msg + "\n")

    def _printVerbWrap(self, msg):
        if self.Options.VerboseWrapper or self._debug:
            sys.stderr.write(msg + "\n")

    @staticmethod
    def _errCheck(err, msg, exc=RttovInternalError):
        if err != 0:
            raise exc(msg)

    def _get_options(self):
        return self._options

    def _set_options(self, obj):
        if not isinstance(obj, Options):
            raise ValueError("Wrong type for options.")
        self._options = obj

    dstr = (":class:`.Options` object where RTTOV options " +
            "can be tuned." +
            "\n\nThis is a Read/Write attribute of type :class:`.Options`.")
    Options = property(_get_options, _set_options, doc=dstr)

    @property_ro
    def CoeffsLoaded(self):
        """Whether or not :meth:`loadInst` has been called successfully."""
        return self._instid > 0

    @property_ro
    def InstId(self):
        """Return the current instrument ID."""
        return self._instid

    @property_ro
    def Nchannels(self):
        """The number of channels available to this instance."""
        return self._nchannels

    FileCoef = _FPathDescRWD('coefs_conf', 'file_coef',
                             doc="Path to the RTTOV coefficients.")

    def _pre_loadInst(self, channels):
        """Before loading the instrument: checks everything."""
        if self.FileCoef is None:
            raise RttovInternalError('The FileCoef attribute must be set')
        OptStack = ['file_coef', self.FileCoef]
        OptStack.append(self.Options.defineStrOptions())
        return (channels, ' '.join(OptStack))

    def loadInst(self, channels=None):
        '''Load instrument for a list of channels.

        At least :data:`FileCoef` must have been set previously.

        :param channels: list of channels to be loaded (all if omitted).
        :type channels: List of int or numpy.ndarray
        '''

        if self.CoeffsLoaded:
            errmsg = "The coefficients were already loaded..."
            raise RttovInternalError(errmsg)

        # Perform some additional checks and determine the option string
        channels, OptString = self._pre_loadInst(channels)
        self._printDbg(OptString)

        if channels is None:
            npchannels = np.array([0, ], dtype=wrapint)
        elif isinstance(channels, collections.Iterable):
            npchannels = np.array(channels, dtype=wrapint)
        else:
            errmsg = "Incorrect type {!s} for channels".format(type(channels))
            raise TypeError(errmsg)

        self._instid = rtwrap.rttov_load_inst(OptString, npchannels)
        if self._instid > 0:
            err, self._nchannels = rtwrap.rttov_get_coef_val_i0(self._instid,
                                                                "NCHANNELS")
            if err != 0:
                raise RttovInternalError('Error getting number of channels')
        else:
            raise RttovInternalError('Error in rttov_load_inst')

        tmpmsg = "Load successful >>>>> inst_id : {0:d}, nchannels : {1:d}."
        self._printVerbWrap(tmpmsg.format(self._instid, self._nchannels))

    @property_ro
    def CoeffsNlevels(self):
        """Number of levels of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_nlevels = rtwrap.rttov_get_coef_val_i0(self._instid,
                                                          "SIZE_REF_PRFL_P")
            self._errCheck(err, 'Error getting the number of levels')
            return l_nlevels
        else:
            return 0

    @property_ro
    def RefPressures(self):
        """Pressure levels of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_p = rtwrap.rttov_get_coef_val_r1(self._instid,
                                                    "REF_PRESSURE",
                                                    self.CoeffsNlevels)
            self._errCheck(err, 'Error getting pressures from coefficient')
            return l_p
        else:
            return None

    @property_ro
    def WaveNumbers(self):
        """Channel's central wavenumbers of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_p = rtwrap.rttov_get_coef_val_r1(self._instid,
                                                    "WAVENUMBERS",
                                                    self.Nchannels)
            self._errCheck(err, 'Error getting channel central wavenumbers')
            return l_p
        else:
            return None

    @property_ro
    def CoefBco(self):
        """Channel's bco of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_p = rtwrap.rttov_get_coef_val_r1(self._instid,
                                                    "FF_BCO",
                                                    self.Nchannels)
            self._errCheck(err, 'Error getting channel bco')
            return l_p
        else:
            return None

    @property_ro
    def CoefBcs(self):
        """Channel's bcs of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_p = rtwrap.rttov_get_coef_val_r1(self._instid,
                                                    "FF_BCS",
                                                    self.Nchannels)
            self._errCheck(err, 'Error getting channel bcs')
            return l_p
        else:
            return None

    @property_ro
    def CoefPlanckc1(self):
        """Planck c1 of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_p = rtwrap.rttov_get_coef_val_r0(self._instid,
                                                    "FC_PLANCK_C1")
            self._errCheck(err, 'Error getting channel Planck c1')
            return l_p
        else:
            return None

    @property_ro
    def CoefPlanckc2(self):
        """Planck c2 of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_p = rtwrap.rttov_get_coef_val_r0(self._instid,
                                                    "FC_PLANCK_C2")
            self._errCheck(err, 'Error getting channel Planck c2')
            return l_p
        else:
            return None

    @property_ro
    def CoefSatHeight(self):
        """Satellite height of the coefficient file."""
        if self.CoeffsLoaded:
            err, l_p = rtwrap.rttov_get_coef_val_r0(self._instid,
                                                    "FC_SAT_HEIGHT")
            self._errCheck(err, 'Error getting channel satellite height')
            return l_p
        else:
            return None

    def updateOptions(self):
        """Update RTTOV options for the currently loaded instrument."""
        if not self.CoeffsLoaded:
            raise RttovWrapperError('Coefficients not loaded')

        str_options = self.Options.defineStrOptions()
        self._printDbg('Set RTTOV options: {}'.format(str_options))
        err = rtwrap.rttov_set_options(self._instid, str_options)
        self._errCheck(err, 'Error in rttov_set_options')

    def printOptions(self):
        """Print RTTOV options for the currently loaded instrument."""
        if not self.CoeffsLoaded:
            raise RttovWrapperError('Coefficients not loaded')
        err = rtwrap.rttov_print_options(self._instid)
        self._errCheck(err, 'Error in rttov_print_options')

    def _getProfiles(self):
        if self._profiles is not None:
            return self._profiles
        else:
            raise RttovInternalError("Profiles not yet initialised.")

    def _setProfiles(self, p):
        if not isinstance(p, self._PROFILES_CLASS):
            raise TypeError("Incorrect value for Profiles")
        if not p.check():
            errmsg = "Error: some mandatory profile fields are missing"
            raise RttovWrapperError(errmsg)
        self._profiles = p

    dstr = (":class:`.Profiles` object currently associated " +
            "with this Rttov instance." +
            "\n\nThis is a Read/Write attribute of type :class:`.Profiles`.")
    Profiles = property(_getProfiles, _setProfiles, doc=dstr)

    def _pressure_fixer(self):
        '''Returns the pressure array. It may come from the coefficient file'''
        if self.Profiles.DefaultPressureLevels:
            if not self.CoeffsLoaded:
                errmsg = ("Error: instrument not loaded, " +
                          "cannot use coefficient file pressure levels")
                raise RttovWrapperError(errmsg)
            if self.CoeffsNlevels != self.Profiles.Nlevels:
                errmsg = ("Error: number of levels differs to coefficient " +
                          "file, cannot use coefficient file pressure levels")
                raise RttovWrapperError(errmsg)
            self._printVerbWrap("Using pressure levels from coefficient file.")
            p = np.empty((self.Profiles.Nprofiles, self.CoeffsNlevels),
                         dtype=wrapfloat)
            p[...] = self.RefPressures
            return p
        else:
            return self.Profiles.P

    def _prerun_bBascis_init(self, channels, variables=()):
        self._printDbg(("Allocating {!s} for {:d} profiles and " +
                        "{:d} channels").format(variables,
                                                self.Profiles.Nprofiles,
                                                len(channels)))
        # First clear _bBasics of the previous results
        self._bBasics = dict()
        for variable in variables:
            self._bBasics[variable] = np.empty((self.Profiles.Nprofiles,
                                                len(channels)), dtype=wrapfloat)

    def _prerun_commons(self, channels):
        if channels is None:
            channels = np.arange(1, self.Nchannels + 1, dtype=wrapint)
        if (not isinstance(channels, np.ndarray) or
                channels.dtype is not wrapint):
            channels = np.array(channels, dtype=wrapint)
        if not self.CoeffsLoaded:
            raise RttovWrapperError("Error: coefficients not loaded")
        if self._profiles is None:
            raise RttovWrapperError("Error: profiles are not set yet")

        # update the options which are locally stored in the class instance
        # for RTTOV
        self.updateOptions()

        # A lot of prints...
        if self._debug:
            self.Profiles.printProfiles(pfunction=self._printDbg)

        return (channels, self._pressure_fixer())


    def _storeLogical(self, outd, storeid, callback):
        err, value = callback(self._instid)
        self._errCheck(err, 'Error when getting {}'.format(storeid))
        outd[storeid] = (value != 0)

    def _1DstoreGeneric(self, outd, storeid, callback, nchannels, dtype=wrapfloat):
        outd[storeid] = np.empty((self.Profiles.Nprofiles * nchannels),
                                 dtype=dtype)
        err = callback(self._instid, outd[storeid])
        self._errCheck(err, 'Error when getting {}'.format(storeid))
        outd[storeid] = outd[storeid].reshape((self.Profiles.Nprofiles,
                                               nchannels))

    def _2DstoreGeneric(self, outd, storeid, callback, nchannels,
                        onlayers=False):
        trailind_dim = self.Profiles.Nlevels - int(onlayers)
        outd[storeid] = np.empty((self.Profiles.Nprofiles * nchannels,
                                  trailind_dim), dtype=wrapfloat)
        err = callback(self._instid, outd[storeid].transpose())
        self._errCheck(err, 'Error when getting {}'.format(storeid))
        outd[storeid] = outd[storeid].reshape((self.Profiles.Nprofiles,
                                               nchannels,
                                               trailind_dim))

    def _2DlevelsstoreGeneric(self, outd, storeid, callback, nlevels, nchannels, dtype=wrapfloat):
        outd[storeid] = np.empty((self.Profiles.Nprofiles * nchannels * nlevels),
                                 dtype=dtype)
        err = callback(self._instid, outd[storeid])
        self._errCheck(err, 'Error when getting {}'.format(storeid))
        outd[storeid] = outd[storeid].reshape((self.Profiles.Nprofiles,
                                               nchannels, nlevels))


    dstring = ("Computed pressure Jacobians from the previous run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, nlevels].")
    PK = _DescRO('bBasics', 'p_k', doc=dstring)

    dstring = ("Computed temperature Jacobians from the previous run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, nlevels].")
    TK = _DescRO('bBasics', 't_k', doc=dstring)

    dstring = ("Computed skin variables Jacobians from the previous run " +
               "[Skin T, salinity, snow_fraction, foam_fraction, " +
               "fastem_coefs(1:5)]. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, 9].")
    SkinK = _DescRO('bBasics', 'skin_k', doc=dstring)

    dstring = ("Computed 2m variables Jacobians from the previous run " +
               "[2m p, 2m t, 2m q, 10m wind u, v, wind fetch]. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, 6].")
    S2mK = _DescRO('bBasics', 's2m_k', doc=dstring)

    dstring = ("Computed gas, cloud and aerosol Jacobians from the previous " +
               "run. ``numpy.ndarray`` of shape " +
               "[ngases, nprofiles, nchannels, nlevels].")
    GasesK = _DescRO('bBasics', 'gases_k', doc=dstring)

    def getItemK(self, gas_id):
        """Computed jacobian for a given gas, cloud or aerosol.

        If the jacobian was not computed, ``None`` will be returned.

        :param gas_id: the gas_id to look for
        :type gas_id: int, str or itemIdType
        :return: computed jacobian
        :rtype: numpy.ndarray of shape [nprofiles, nchannels, nlevels].
        """
        g_ids = self._bBasics.get('gases_k_id', None)
        g = self._bBasics.get('gases_k', None)
        if g_ids is not None and g is not None:
            g_idx = np.where(g_ids == itemIdType(gas_id))
            if len(g_idx[0]):
                return g[g_idx[0][0], ...]
        # The gas was not found... (or no computation was done)
        return None


@add_descriptors_gasesK(_ItemDescriptorNaming)
class Rttov(_RttovCommons):

    '''Main wrapper class for RTTOV.'''

    def __init__(self):
        super(Rttov, self).__init__()
        del self.SurfEmisRefl
        self._opt_param = {}
        self._bTrans = dict()
        self._bRad2 = dict()

    FileSccld = _FPathDescRWD('coefs_conf', 'file_sccld',
                              doc="Path to the RTTOV cloud coefficients.")

    FileScaer = _FPathDescRWD('coefs_conf', 'file_scaer',
                              doc="Path to the RTTOV aerosol coefficients.")

    FileMfasisCld = _FPathDescRWD('coefs_conf', 'file_mfasis_cld',
                                  doc="Path to the MFASIS cloud LUT file.")

    #FileMfasisAer = _FPathDescRWD('coefs_conf', 'file_mfasis_aer',
                                  #doc="Path to the MFASIS aerosol LUT file.")

    #FilePcCoef = _FPathDescRWD('coefs_conf', 'file_pccoef',
                               #doc="Path to the RTTOV PC coefficients.")

    def updateOptions(self):
        """Update RTTOV options for the currently loaded instrument."""
        # This may prevent crash
        if not self.FileSccld and not self.Options.UserCldOptParam:
            self.Options.AddClouds = False
        if not self.FileScaer and not self.Options.UserAerOptParam:
            self.Options.AddAerosl = False
        super(Rttov, self).updateOptions()

    def _pre_loadInst(self, channels):
        """Before loading the instrument: checks everything."""
        optionsloc = copy.deepcopy(self.Options)

        if self.FileCoef is None:
            raise RttovInternalError('The FileCoef attribute must be set')
        OptStack = ['file_coef', self.FileCoef]

        if self.FileSccld is not None:
            OptStack.extend(['file_sccld', self.FileSccld])
            optionsloc.AddClouds = True

        if self.FileScaer is not None:
            OptStack.extend(['file_scaer', self.FileScaer])
            optionsloc.AddAerosl = True

        if self.FileMfasisCld is not None:
            OptStack.extend(['file_mfasis_cld', self.FileMfasisCld])

        #if self.FileMfasisAer is not None:
            #OptStack.extend(['file_mfasis_aer', self.FileMfasisAer])

        #if self.FilePcCoef is not None:
            #OptStack.extend(['file_pccoef', self.FilePcCoef])
            #optionsloc.AddPC = True

        OptStack.append(optionsloc.defineStrOptions())
        return (channels, ' '.join(OptStack))

    def _getSurfEmisRefl(self):
        if self._surfemisrefl is not None:
            return self._surfemisrefl
        else:
            raise RttovInternalError("SurfEmisRefl not yet initialised.")

    def _setSurfEmisRefl(self, a):
        if (isinstance(a, np.ndarray) and a.shape[0] == 4 and
                len(a.shape) == 3):
            if a.dtype is not np.dtype(wrapfloat) or np.isfortran(a):
                a = np.array(a, order='C', dtype=wrapfloat)
            self._surfemisrefl_manual = True
            self._surfemisrefl = a
        else:
            raise TypeError("Incorrect value type/shape for SurfEmisRefl")

    def _delSurfEmisRefl(self):
        self._surfemisrefl_manual = False
        self._surfemisrefl = None

    dstring = """Array containing input/output surface emissivity, reflectance
    (direct BRDF and diffuse reflectance) and specularity values.

    This must be a ``numpy.ndarray`` of shape [4, nprofiles, nchannels], the
    size of the last two dimensions cannot be checked here. However, it will
    eventually be checked during the :meth:`runDirect` or :meth:`runK` calls.

    The `dtype` of the ``numpy.ndarray`` must be {!r}. If the input array has
    another `dtype`, a conversion will be attempted.

    Use `del` to empty the SurfEmisRefl array.

    If this array is not set by the user (or has been emptied unsig `del`), it
    will be initialised to -1 in the :meth:`runDirect` or :meth:`runK` calls.
    """.format(wrapfloat)
    SurfEmisRefl = property(_getSurfEmisRefl, _setSurfEmisRefl,
                            _delSurfEmisRefl, doc=dstring)

    def printSurfEmisRefl(self, pfunction=print):
        """Print the :data:`SurfEmisRefl` content on the standard output."""
        fmt_e = 'Channel #{:5d} e={:6.3f} brdf={:6.3f} diff={:6.3f} spec={:6.3f}'
        try:
            emisrefl = self.SurfEmisRefl
        except RttovInternalError:
            emisrefl = None
        if emisrefl is not None:
            for prof in range(self.Profiles.Nprofiles):
                pfunction(">> SurfEmis/Refl for profile #{:d}".format(prof))
                for ich in range(emisrefl.shape[-1]):
                    pfunction(fmt_e.format(ich + 1,
                                           emisrefl[0, prof, ich],
                                           emisrefl[1, prof, ich],
                                           emisrefl[2, prof, ich],
                                           emisrefl[3, prof, ich]))
        else:
            pfunction(">> SurfEmisRefl is not (yet?) initialised.")

    def _prerun_commons(self, channels):
        (channels, fixed_p) = super(Rttov, self)._prerun_commons(channels)

        # deal with SurfEmisRefl
        surfemisrefl_shape = (4, self.Profiles.Nprofiles, len(channels))
        if self._surfemisrefl_manual:
            # Check that the dimensions are fine
            if self.SurfEmisRefl.shape != surfemisrefl_shape:
                raise RttovWrapperError("Incorrect shape for SurfEmisRefl")
        else:
            self._printVerbWrap("No surface emissivity/reflectance supplied:" +
                                " setting calcemis and calcrefl to true")
            self._surfemisrefl = - np.ones(surfemisrefl_shape, dtype=wrapfloat)
        if self._debug:
            self.printSurfEmisRefl(pfunction=self._printDbg)

        # rads and btrefl are always required
        self._prerun_bBascis_init(channels, ('btrefl', 'rads', ))

        return (channels, fixed_p)

    AerPhangle  = _ArrayRWD('opt_param', 'aer_phangle', dims="aer_nphangle",
                            doc="Aersol phase function angles.")

    AerAsb      = _ArrayRWD('opt_param', 'aer_asb', dims='3,nprofiles,nchannels,nlayers',
                            doc='Aerosol absorption coefficients, scattering coefficients ' +
                                'and bpr parameters.')

    AerLegcoef  = _ArrayRWD('opt_param', 'aer_legcoef', dims='nprofiles,nchannels,nlayers,aer_nmom+1',
                            doc='Aerosol phase function Legendre coefficients.')

    AerPha      = _ArrayRWD('opt_param', 'aer_pha', dims='nprofiles,nchannels,nlayers,aer_nphangle',
                            doc='Aerosol phase functions.')

    CldPhangle  = _ArrayRWD('opt_param', 'cld_phangle', dims="cld_nphangle",
                            doc="Cloud phase function angles.")

    CldAsb      = _ArrayRWD('opt_param', 'cld_asb', dims='3,nprofiles,nchannels,nlayers',
                            doc='Cloud absorption coefficients, scattering coefficients ' +
                                'and bpr parameters.')

    CldLegcoef  = _ArrayRWD('opt_param', 'cld_legcoef', dims='nprofiles,nchannels,nlayers,cld_nmom+1',
                            doc='Cloud phase function Legendre coefficients.')

    CldPha      = _ArrayRWD('opt_param', 'cld_pha', dims='nprofiles,nchannels,nlayers,cld_nphangle',
                            doc='Cloud phase functions.')

    def calcBpr(self, phangle, pha):
        """Calculate bpr parameter for a given phase function pha defined on angles phangle.

        :parameter phangle: Phase function angles
        :type phangle: numpy.ndarray of floats of shape [nphangle,]
        :parameter pha: Phase functions
        :type pha: numpy.ndarray of floats of shape [nphangle,]

        :rtype: float
        """
        if len(phangle) != len(pha):
            raise RttovWrapperError('Angle and phase function arrays must be the same size')
        err, bpr = rtwrap.rttov_bpr(phangle, pha, self.Options.Nthreads)
        if err == 0:
            return bpr
        else:
            raise RttovInternalError('Error calculating bpr')

    def calcLegcoef(self, phangle, pha, nmom, ngauss=0):
        """Calculate Legendre coefficients for given phase function pha defined on angles phangle.

        :parameter phangle: Phase function angles
        :type phangle: numpy.ndarray of floats of shape [nphangle,]
        :parameter pha: Phase functions
        :type pha: numpy.ndarray of floats of shape [nphangle,]
        :parameter int nmom: nmom
        :parameter int nguauss: nguauss

        :rtype: numpy.ndarray of floats of shape [nmom + 1,]
        """
        if len(phangle) != len(pha):
            raise RttovWrapperError('Angle and phase function arrays must be the same size')
        legcoef = np.empty((nmom + 1), dtype=wrapfloat)
        err = rtwrap.rttov_legcoef(phangle, pha, legcoef, ngauss)
        if err == 0:
            return legcoef
        else:
            raise RttovInternalError('Error calculating Legendre coefficients')

    def _prerun_opt_param(self, nchannels, cldaer):
        """Check Aerosol or Cloud propoerties before run.

        This method is called unconditionally whenever Aerosol or Cloud parameters
        are given by the user. Therefore, some of the arrays may be missing.
        Given the Options, this method checks that the mandatory data are properly
        defined. If necessary, it allocates temporary empty arrays to comply with the
        f2py wrapper interface.

        The 4 arrays returned by this function are ready to be used (i.e no
        transpose needed).
        """
        ang   = getattr(self, cldaer + 'Phangle')
        asb   = getattr(self, cldaer + 'Asb')
        lcoef = getattr(self, cldaer + 'Legcoef')
        pha   = getattr(self, cldaer + 'Pha')

        if cldaer.lower() == 'aer':
            dothis = self.Options.AddAerosl and self.Options.UserAerOptParam
        else:
            dothis = self.Options.AddClouds and self.Options.UserCldOptParam

        (actual_asb, actual_ang, actual_pha, actual_lcoef) = (None, None, None, None)
        # Check for mandatory inputs
        if dothis:
            if asb is None:
                raise RttovWrapperError('Absorption, scattering and bpr coefficients must be defined')
            if asb.shape != (3, self.Profiles.Nprofiles, nchannels, self.Profiles.Nlevels - 1):
                raise RttovWrapperError('Incorrect shape of {:s}Asb'.format(cldaer))
            actual_asb = asb.transpose()
            if self.Options.AddSolar:
                if ang is None or pha is None:
                    raise RttovWrapperError('Phase functions and phase function angles must be ' +
                                            'defined for solar simulations')
                if pha.shape != (self.Profiles.Nprofiles, nchannels, self.Profiles.Nlevels - 1,
                                 len(ang)):
                    raise RttovWrapperError('Incorrect shape of {:s}Pha'.format(cldaer))
                actual_ang = ang.transpose()
                actual_pha = pha.transpose()
            if (self.Options.IrScattModel == DOM_IR or
                    (self.Options.VisScattModel == DOM_VIS and self.Options.AddSolar)):
                if lcoef is None:
                    raise RttovWrapperError('Legendre coefficients must be defined for DOM simulations')
                if lcoef.shape[:-1] != (self.Profiles.Nprofiles, nchannels, self.Profiles.Nlevels - 1):
                    raise RttovWrapperError('Incorrect shape of {:s}Legcoef'.format(cldaer))
                actual_lcoef = lcoef.transpose()

        if actual_asb is None:
            actual_asb = np.empty((self.Profiles.Nlevels - 1, nchannels, self.Profiles.Nprofiles, 3),
                                  dtype=wrapfloat, order='F')
        if actual_ang is None:
            actual_ang = np.empty((1,), dtype=wrapfloat, order='F')
        if actual_pha is None:
            actual_pha = np.empty((1, self.Profiles.Nlevels - 1, nchannels, self.Profiles.Nprofiles),
                                  dtype=wrapfloat, order='F')
        if actual_lcoef is None:
            actual_lcoef = np.empty((1, self.Profiles.Nlevels - 1, nchannels, self.Profiles.Nprofiles),
                                    dtype=wrapfloat, order='F')

        return (actual_ang, actual_asb, actual_lcoef, actual_pha)

    def runDirect(self, channels=None):
        """Run the RTTOV direct model.

        :param channels: list of channels to simulate (all the channels if
                         omitted)
        :type channels: list of int or numpy.ndarray
        """
        do_opt_param = (self.Options.AddAerosl and self.Options.UserAerOptParam) or \
                       (self.Options.AddClouds and self.Options.UserCldOptParam)

        channels, p = self._prerun_commons(channels)

        if do_opt_param:
            (aer_ang,       # aerosol phase fn angles                   [aer_nphangle]
             aer_asb,       # aerosol abs, sca and bpr                  [3][nprofiles][nchannels][nlayers]
             aer_lcoef,     # aerosol phase fn Legendre coefficients    [nprofiles][nchannels][nlayers][aer_nmom+1]
             aer_pha,       # aerosol phase fns                         [nprofiles][nchannels][nlayers][aer_nphangle]
             ) = self._prerun_opt_param(len(channels), 'Aer')
            (cld_ang,       # cloud phase fn angles                     [cld_nphangle]
             cld_asb,       # cloud abs, sca and bpr                    [3][nprofiles][nchannels][nlayers]
             cld_lcoef,     # cloud phase fn Legendre coefficients      [rofiles][nchannels][nlayers][cld_nmom+1]
             cld_pha,       # cloud phase fns                           [nprofiles][nchannels][nlayers][cld_nphangle]
             ) = self._prerun_opt_param(len(channels), 'Cld')
            err = rtwrap.rttov_visir_scatt_call_direct(
                self._instid,
                channels,
                self.Profiles.DateTimes.transpose(),    # profile dates/times                                     [nprofiles][6]
                self.Profiles.Angles.transpose(),       # satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
                self.Profiles.SurfGeom.transpose(),     # lat, lon, elevation                                     [nprofiles][3]
                self.Profiles.SurfType.transpose(),     # surftype, watertype                                     [nprofiles][2]
                self.Profiles.Skin.transpose(),         # skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
                self.Profiles.S2m.transpose(),          # 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
                self.Profiles.ClwScheme.transpose(),    # clw_scheme, clwde_param                                 [nprofiles][2]
                self.Profiles.IceCloud.transpose(),     # ice_scheme, idg                                         [nprofiles][2]
                p.transpose(),                          # pressure                                                [nprofiles][nlevels]
                self.Profiles.T.transpose(),            # temperature                                             [nprofiles][nlevels]
                self.Profiles.GasUnits,                 # units for gas profiles
                self.Profiles.MmrCldAer,                # units for cloud/aerosol profiles
                self.Profiles.GasId.transpose(),        # gas ID list                                             [ngases]
                self.Profiles.Gases.transpose(),        # gas profiles                                            [ngases][nprofiles][nlevels]
                aer_ang, aer_asb, aer_lcoef, aer_pha,   # Aerosol properties (see above)
                cld_ang, cld_asb, cld_lcoef, cld_pha,   # Cloud properties (see above)
                self._surfemisrefl.transpose(),         # input/output surface emis, refl, spec                   [4][nprofiles][nchannels]
                self._bBasics['btrefl'].transpose(),    # output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
                self._bBasics['rads'].transpose(),      # output radiances                                        [nprofiles][nchannels]
                )
            self._errCheck(err, 'Error in rttov_visir_scatt_call_direct')
        else:
            err = rtwrap.rttov_call_direct(
                self._instid,
                channels,
                self.Profiles.DateTimes.transpose(),    # profile dates/times                                     [nprofiles][6]
                self.Profiles.Angles.transpose(),       # satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
                self.Profiles.SurfGeom.transpose(),     # lat, lon, elevation                                     [nprofiles][3]
                self.Profiles.SurfType.transpose(),     # surftype, watertype                                     [nprofiles][2]
                self.Profiles.Skin.transpose(),         # skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
                self.Profiles.S2m.transpose(),          # 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
                self.Profiles.SimpleCloud.transpose(),  # ctp, cfraction                                          [nprofiles][2]
                self.Profiles.ClwScheme.transpose(),    # clw_scheme, clwde_param                                 [nprofiles][2]
                self.Profiles.IceCloud.transpose(),     # ice_scheme, idg                                         [nprofiles][2]
                self.Profiles.Zeeman.transpose(),       # Be, cosbk                                               [nprofiles][2]
                p.transpose(),                          # pressure                                                [nprofiles][nlevels]
                self.Profiles.T.transpose(),            # temperature                                             [nprofiles][nlevels]
                self.Profiles.GasUnits,                 # units for gas profiles
                self.Profiles.MmrCldAer,                # units for cloud/aerosol profiles
                self.Profiles.GasId.transpose(),        # gas ID list                                             [ngases]
                self.Profiles.Gases.transpose(),        # gas profiles                                            [ngases][nprofiles][nlevels]
                self._surfemisrefl.transpose(),         # input/output surface emis, refl, spec                   [4][nprofiles][nchannels]
                self._bBasics['btrefl'].transpose(),    # output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
                self._bBasics['rads'].transpose(),      # output radiances                                        [nprofiles][nchannels]
                )
            self._errCheck(err, 'Error in rttov_call_direct')

        self._doStoreTrans(len(channels))
        self._doStoreRad(len(channels))
        if not do_opt_param:
            self._doStoreRad2(len(channels))

    def runK(self, channels=None):
        """Run the RTTOV K model.

        :param channels: list of channels to simulate (all the channels if
                         omitted)
        :type channels: list of int or numpy.ndarray
        """
        do_opt_param = (self.Options.AddAerosl and self.Options.UserAerOptParam) or \
                       (self.Options.AddClouds and self.Options.UserCldOptParam)

        channels, p = self._prerun_commons(channels)

        self._printDbg("Allocating the many _k arrays")

        def ini_k_like(tab):
            theshape = list(tab.shape)
            # An extra dimension is needed for K computations
            theshape.insert(-1, len(channels))
            return np.empty(theshape, dtype=tab.dtype)

        self._bBasics['skin_k'] = ini_k_like(self.Profiles.Skin)
        self._bBasics['s2m_k'] = ini_k_like(self.Profiles.S2m)
        self._bBasics['p_k'] = ini_k_like(p)
        self._bBasics['t_k'] = ini_k_like(self.Profiles.T)
        self._bBasics['gases_k'] = ini_k_like(self.Profiles.Gases)
        self._bBasics['surfemisrefl_k'] = np.zeros_like(self._surfemisrefl)
        self._bBasics['bt_k'] = np.ones_like(self._bBasics['btrefl'])
        self._bBasics['rad_k'] = np.ones_like(self._bBasics['rads'])

        if do_opt_param:
            (aer_ang,       # aerosol phase fn angles                   [aer_nphangle]
             aer_asb,       # aerosol abs, sca and bpr                  [3][nprofiles][nchannels][nlayers]
             aer_lcoef,     # aerosol phase fn Legendre coefficients    [nprofiles][nchannels][nlayers][aer_nmom+1]
             aer_pha,       # aerosol phase fns                         [nprofiles][nchannels][nlayers][aer_nphangle]
             ) = self._prerun_opt_param(len(channels), 'Aer')
            (cld_ang,       # cloud phase fn angles                     [cld_nphangle]
             cld_asb,       # cloud abs, sca and bpr                    [3][nprofiles][nchannels][nlayers]
             cld_lcoef,     # cloud phase fn Legendre coefficients      [rofiles][nchannels][nlayers][cld_nmom+1]
             cld_pha,       # cloud phase fns                           [nprofiles][nchannels][nlayers][cld_nphangle]
             ) = self._prerun_opt_param(len(channels), 'Cld')
            err = rtwrap.rttov_visir_scatt_call_k(
                self._instid,
                channels,
                self.Profiles.DateTimes.transpose(),        # profile dates/times                                     [nprofiles][6]
                self.Profiles.Angles.transpose(),           # satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
                self.Profiles.SurfGeom.transpose(),         # lat, lon, elevation                                     [nprofiles][3]
                self.Profiles.SurfType.transpose(),         # surftype, watertype                                     [nprofiles][2]
                self.Profiles.Skin.transpose(),             # skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
                self._bBasics['skin_k'].transpose(),        # skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][nchannels][9]
                self.Profiles.S2m.transpose(),              # 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
                self._bBasics['s2m_k'].transpose(),         # 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][nchannels][6]
                self.Profiles.ClwScheme.transpose(),        # clw_scheme, clwde_param                                 [nprofiles][2]
                self.Profiles.IceCloud.transpose(),         # ice_scheme, idg                                         [nprofiles][2]
                p.transpose(),                              # pressure                                                [nprofiles][nlevels]
                self._bBasics['p_k'].transpose(),           # pressure                                                [nprofiles][nchannels][nlevels]
                self.Profiles.T.transpose(),                # temperature                                             [nprofiles][nlevels]
                self._bBasics['t_k'].transpose(),           # temperature                                             [nprofiles][nchannels][nlevels]
                self.Profiles.GasUnits,                     # units for gas profiles
                self.Profiles.MmrCldAer,                    # units for cloud/aerosol profiles
                self.Profiles.GasId.transpose(),            # gas ID list                                             [ngases]
                self.Profiles.Gases.transpose(),            # gas profiles                                            [ngases][nprofiles][nlevels]
                self._bBasics['gases_k'].transpose(),       # gas profiles                                            [ngases][nprofiles][nchannels][nlevels]
                aer_ang, aer_asb, aer_lcoef, aer_pha,       # Aerosol properties (see above)
                cld_ang, cld_asb, cld_lcoef, cld_pha,       # Cloud properties (see above)                            [nprofiles][nchannels][nlayers][cld_nphangle]
                self._surfemisrefl.transpose(),             # input/output surface emis, refl, spec                   [4][nprofiles][nchannels]
                self._bBasics['surfemisrefl_k'].transpose(),# output surface emis, refl, spec K                       [4][nprofiles][nchannels]
                self._bBasics['btrefl'].transpose(),        # output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
                self._bBasics['rads'].transpose(),          # output radiances                                        [nprofiles][nchannels]
                self._bBasics['bt_k'].transpose(),          # output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
                self._bBasics['rad_k'].transpose(),         # output radiances                                        [nprofiles][nchannels]
                )
            self._errCheck(err, 'Error in rttov_visir_scatt_call_k')
        else:
            # It's only usefull in this case...
            self._bBasics['simplecloud_k'] = ini_k_like(self.Profiles.SimpleCloud)

            err = rtwrap.rttov_call_k(
                self._instid,
                channels,
                self.Profiles.DateTimes.transpose(),        # profile dates/times                                     [nprofiles][6]
                self.Profiles.Angles.transpose(),           # satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
                self.Profiles.SurfGeom.transpose(),         # lat, lon, elevation                                     [nprofiles][3]
                self.Profiles.SurfType.transpose(),         # surftype, watertype                                     [nprofiles][2]
                self.Profiles.Skin.transpose(),             # skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
                self._bBasics['skin_k'].transpose(),        # skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][nchannels][9]
                self.Profiles.S2m.transpose(),              # 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
                self._bBasics['s2m_k'].transpose(),         # 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][nchannels][6]
                self.Profiles.SimpleCloud.transpose(),      # ctp, cfraction                                          [nprofiles][2]
                self._bBasics['simplecloud_k'].transpose(), # ctp, cfraction                                          [nprofiles][nchannels][2]
                self.Profiles.ClwScheme.transpose(),        # clw_scheme, clwde_param                                 [nprofiles][2]
                self.Profiles.IceCloud.transpose(),         # ice_scheme, idg                                         [nprofiles][2]
                self.Profiles.Zeeman.transpose(),           # Be, cosbk                                               [nprofiles][2]
                p.transpose(),                              # pressure                                                [nprofiles][nlevels]
                self._bBasics['p_k'].transpose(),           # pressure                                                [nprofiles][nchannels][nlevels]
                self.Profiles.T.transpose(),                # temperature                                             [nprofiles][nlevels]
                self._bBasics['t_k'].transpose(),           # temperature                                             [nprofiles][nchannels][nlevels]
                self.Profiles.GasUnits,                     # units for gas profiles
                self.Profiles.MmrCldAer,                    # units for cloud/aerosol profiles
                self.Profiles.GasId.transpose(),            # gas ID list                                             [ngases]
                self.Profiles.Gases.transpose(),            # gas profiles                                            [ngases][nprofiles][nlevels]
                self._bBasics['gases_k'].transpose(),       # gas profiles                                            [ngases][nprofiles][nchannels][nlevels]
                self._surfemisrefl.transpose(),             # input/output surface emis, refl, spec                   [4][nprofiles][nchannels]
                self._bBasics['surfemisrefl_k'].transpose(),# output surface emis, refl, spec K                       [4][nprofiles][nchannels]
                self._bBasics['btrefl'].transpose(),        # output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
                self._bBasics['rads'].transpose(),          # output radiances                                        [nprofiles][nchannels]
                self._bBasics['bt_k'].transpose(),          # output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
                self._bBasics['rad_k'].transpose(),         # output radiances                                        [nprofiles][nchannels]
                )
            self._errCheck(err, 'Error in rttov_call_k')

        # Save the gas_id list for a later use
        self._bBasics['gases_k_id'] = self.Profiles.GasId

        self._doStoreTrans(len(channels))
        self._doStoreRad(len(channels))
        if not do_opt_param:
            self._doStoreRad2(len(channels))

    def _doStoreTrans(self, nchannels):
        if not self.Options.StoreTrans:
            self._bTrans = dict()
            return
        self._1DstoreGeneric(self._bTrans, 'tau_total',
                             rtwrap.rttov_get_tau_total, nchannels)
        self._2DstoreGeneric(self._bTrans, 'tau_levels',
                             rtwrap.rttov_get_tau_levels, nchannels)
        self._1DstoreGeneric(self._bTrans, 'tausun_total_path1',
                             rtwrap.rttov_get_tausun_total_path1, nchannels)
        self._2DstoreGeneric(self._bTrans, 'tausun_levels_path1',
                             rtwrap.rttov_get_tausun_levels_path1, nchannels)
        self._1DstoreGeneric(self._bTrans, 'tausun_total_path2',
                             rtwrap.rttov_get_tausun_total_path2, nchannels)
        self._2DstoreGeneric(self._bTrans, 'tausun_levels_path2',
                             rtwrap.rttov_get_tausun_levels_path2, nchannels)
        self._1DstoreGeneric(self._bTrans, 'tau_total_cld',
                             rtwrap.rttov_get_tau_total_cld, nchannels)
        self._2DstoreGeneric(self._bTrans, 'tau_levels_cld',
                             rtwrap.rttov_get_tau_levels_cld, nchannels)

    def _doStoreRad(self, nchannels):
        if not self.Options.StoreRad:
            self._bRad = dict()
            return
        self._1DstoreGeneric(self._bRad, 'rad_clear',
                             rtwrap.rttov_get_rad_clear, nchannels)
        self._1DstoreGeneric(self._bRad, 'rad_total',
                             rtwrap.rttov_get_rad_total, nchannels)
        self._1DstoreGeneric(self._bRad, 'bt_clear',
                             rtwrap.rttov_get_bt_clear, nchannels)
        self._1DstoreGeneric(self._bRad, 'bt',
                             rtwrap.rttov_get_bt, nchannels)
        self._1DstoreGeneric(self._bRad, 'refl_clear',
                             rtwrap.rttov_get_refl_clear, nchannels)
        self._1DstoreGeneric(self._bRad, 'refl',
                             rtwrap.rttov_get_refl, nchannels)
        self._1DstoreGeneric(self._bRad, 'rad_cloudy',
                             rtwrap.rttov_get_rad_cloudy, nchannels)
        self._2DstoreGeneric(self._bRad, 'overcast',
                             rtwrap.rttov_get_overcast, nchannels,
                             onlayers=True)
        self._1DstoreGeneric(self._bRad, 'rad_quality',
                             rtwrap.rttov_get_rad_quality, nchannels,
                             dtype=wrapint)
        self._storeLogical(self._bRad, 'plane_parallel',
                             rtwrap.rttov_get_plane_parallel)
        self._2DstoreGeneric(self._bRad, 'geometric_height',
                             rtwrap.rttov_get_geometric_height, nchannels)


    def _doStoreRad2(self, nchannels):
        if not self.Options.StoreRad2:
            self._bRad2 = dict()
            return
        self._1DstoreGeneric(self._bRad2, 'rad2_up_clear',
                             rtwrap.rttov_get_rad2_upclear, nchannels)
        self._1DstoreGeneric(self._bRad2, 'rad2_dn_clear',
                             rtwrap.rttov_get_rad2_dnclear, nchannels)
        self._1DstoreGeneric(self._bRad2, 'rad2_refldn_clear',
                             rtwrap.rttov_get_rad2_refldnclear, nchannels)
        self._2DstoreGeneric(self._bRad2, 'rad2_up',
                             rtwrap.rttov_get_rad2_up, nchannels,
                             onlayers=True)
        self._2DstoreGeneric(self._bRad2, 'rad2_down',
                             rtwrap.rttov_get_rad2_down, nchannels,
                             onlayers=True)
        self._2DstoreGeneric(self._bRad2, 'rad2_surf',
                             rtwrap.rttov_get_rad2_surf, nchannels,
                             onlayers=True)

    dstring = ("Computed brightness temperatures or reflectances from the " +
               "previous run. ``numpy.ndarray`` of shape " +
               "[nprofiles, nchannels].")
    BtRefl = _DescRO('bBasics', 'btrefl', doc=dstring)

    dstring = ("Computed radiances from the previous run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels].")
    Rads = _DescRO('bBasics', 'rads', doc=dstring)

    dstring = ("Computed simple cloud Jacobians from the previous run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, 2].")
    SimpleCloudK = _DescRO('bBasics', 'simplecloud_k', doc=dstring)

    @property_ro
    def SurfEmisK(self):
        """Computed surface emissivity Jacobians from the previous run. ``numpy.ndarray`` of shape [nprofiles, nchannels]."""
        res = self._bBasics.get('surfemisrefl_k', None)
        if res is not None:
            res = res[0, ...]
        return res

    @property_ro
    def SurfReflK(self):
        """Computed surface BRDF Jacobians from the previous run. ``numpy.ndarray`` of shape [nprofiles, nchannels]."""
        res = self._bBasics.get('surfemisrefl_k', None)
        if res is not None:
            res = res[1, ...]
        return res

    @property_ro
    def SurfDiffuseReflK(self):
        """Computed surface diffuse reflectance Jacobians from the previous run. ``numpy.ndarray`` of shape [nprofiles, nchannels]."""
        res = self._bBasics.get('surfemisrefl_k', None)
        if res is not None:
            res = res[2, ...]
        return res

    @property_ro
    def SpecularityK(self):
        """Computed surface specularity Jacobians from the previous run. ``numpy.ndarray`` of shape [nprofiles, nchannels]."""
        res = self._bBasics.get('surfemisrefl_k', None)
        if res is not None:
            res = res[3, ...]
        return res

    def getUserAerNK(self, n):
        """Computed Jacobians for user aerosol species n from the previous run. ``numpy.ndarray`` of shape [nprofiles, nchannels, nlevels]."""
        if n < 0 or n > maxuseraer:
            raise ValueError("Aerosol index must lie between 1 and "+str(maxuseraer))
        return getattr(self, 'Aer'+str(n)+'K')

    # Outputs under the store_trans key
    fdstring = ('RTTOV transmission {} output ``numpy.ndarray`` of shape' +
                '[nprofiles, nchannels{}]. Requires store_trans=True.')

    dstring = fdstring.format('tau_total', '')
    TauTotal = _DescRO('bTrans', 'tau_total', doc=dstring)

    dstring = fdstring.format('tau_levels', ', nlevels')
    TauLevels = _DescRO('bTrans', 'tau_levels', doc=dstring)

    dstring = fdstring.format('tausun_total_path1', '')
    TauSunTotalPath1 = _DescRO('bTrans', 'tausun_total_path1', doc=dstring)

    dstring = fdstring.format('tausun_levels_path1', ', nlevels')
    TauSunLevelsPath1 = _DescRO('bTrans', 'tausun_levels_path1', doc=dstring)

    dstring = fdstring.format('tausun_total_path1', '')
    TauSunTotalPath2 = _DescRO('bTrans', 'tausun_total_path2', doc=dstring)

    dstring = fdstring.format('tausun_levels_path2', ', nlevels')
    TauSunLevelsPath2 = _DescRO('bTrans', 'tausun_levels_path2', doc=dstring)

    dstring = fdstring.format('tau_total_cld', '')
    TauTotalCld = _DescRO('bTrans', 'tau_total_cld', doc=dstring)

    dstring = fdstring.format('tau_levels_cld', ', nlevels')
    TauLevelsCld = _DescRO('bTrans', 'tau_levels_cld', doc=dstring)

    # Outputs under the store_rad key
    fdstring = ('RTTOV radiance {} output ``numpy.ndarray`` of shape' +
                '[nprofiles, nchannels{}]. Requires store_rad=True.')

    dstring = fdstring.format('clear', '')
    RadClear = _DescRO('bRad', 'rad_clear', doc=dstring)

    dstring = fdstring.format('total', '')
    RadTotal = _DescRO('bRad', 'rad_total', doc=dstring)

    dstring = fdstring.format('bt_clear', '')
    BtClear = _DescRO('bRad', 'bt_clear', doc=dstring)

    dstring = fdstring.format('bt', '')
    Bt = _DescRO('bRad', 'bt', doc=dstring)

    dstring = fdstring.format('refl_clear', '')
    ReflClear = _DescRO('bRad', 'refl_clear', doc=dstring)

    dstring = fdstring.format('refl', '')
    Refl = _DescRO('bRad', 'refl', doc=dstring)

    dstring = fdstring.format('cloudy', '')
    RadCloudy = _DescRO('bRad', 'rad_cloudy', doc=dstring)

    dstring = fdstring.format('overcast', ', nlayers')
    Overcast = _DescRO('bRad', 'overcast', doc=dstring)

    dstring = fdstring.format('quality', '')
    RadQuality = _DescRO('bRad', 'rad_quality', doc=dstring)

    dstring = fdstring.format('plane_parallel', '')
    PlaneParallel = _DescRO('bRad', 'plane_parallel', doc=dstring)

    dstring = fdstring.format('geometric_height', ', nlevels')
    GeometricHeight = _DescRO('bRad', 'geometric_height', doc=dstring)

    # Outputs under the store_rad2 key
    fdstring = ('RTTOV radiance2 {} output ``numpy.ndarray`` of shape' +
                '[nprofiles, nchannels{}]. Requires store_rad2=True.')

    dstring = fdstring.format('upclear', '')
    Rad2UpClear = _DescRO('bRad2', 'rad2_up_clear', doc=dstring)

    dstring = fdstring.format('dnclear', '')
    Rad2DnClear = _DescRO('bRad2', 'rad2_dn_clear', doc=dstring)

    dstring = fdstring.format('refldnclear', '')
    Rad2ReflDnClear = _DescRO('bRad2', 'rad2_refldn_clear', doc=dstring)

    dstring = fdstring.format('up', ', nlayers')
    Rad2Up = _DescRO('bRad2', 'rad2_up', doc=dstring)

    dstring = fdstring.format('down', ', nlayers')
    Rad2Down = _DescRO('bRad2', 'rad2_down', doc=dstring)

    dstring = fdstring.format('surf', ', nlayers')
    Rad2Surf = _DescRO('bRad2', 'rad2_surf', doc=dstring)


@add_descriptors_gasesK(_ItemDescriptorNamingScatt)
class RttovScatt(_RttovCommons):

    '''Main wrapper class for RTTOV-SCATT.'''

    _PROFILES_CLASS = ProfilesScatt

    def __init__(self):
        super(RttovScatt, self).__init__()
        self._bEmisTerms = dict()
        self._bZef = dict()
        del self.SurfEmis

    FileHydrotable = _FPathDescRWD('coefs_conf', 'file_hydrotable',
                                 doc="Path to the RTTOV-SCATT Hydrotable file.")

    MultiHydroFrac = _TDescRW('coefs_conf', 'multi_hydro_frac', bool, 
                       doc="Single hydro_frac profile or one per hydrometeor.")

    CalcZef = _TDescRW('coefs_conf', 'calc_zef', bool, 
                       doc="Flag to enable reflectivities calculation.")

    def _pre_loadInst(self, channels):
        """Before loading the instrument: checks everything."""
        optionsloc = copy.deepcopy(self.Options)

        if self.FileCoef is None:
            raise RttovInternalError('The FileCoef attribute must be set')
        OptStack = ['file_coef', self.FileCoef]

        if self.FileHydrotable is None:
            raise RttovInternalError('The FileHydrotable attribute must be set')
        OptStack.extend(['file_hydrotable', self.FileHydrotable])

        if channels is not None:
            self._printVerbWrap("channels attribute is ignored with RttovScatt.")
            channels = None

        OptStack.append(optionsloc.defineStrOptions())
        return (channels, ' '.join(OptStack))

    def _getSurfEmis(self):
        if self._surfemis is not None:
            return self._surfemis
        else:
            raise RttovInternalError("SurfEmis not yet initialised.")

    def _setSurfEmis(self, a):
        if (isinstance(a, np.ndarray) and len(a.shape) == 2):
            if a.dtype is not np.dtype(wrapfloat) or np.isfortran(a):
                a = np.array(a, order='C', dtype=wrapfloat)
            self._surfemis_manual = True
            self._surfemis = a
        else:
            raise TypeError("Incorrect value type/shape for SurfEmis")

    def _delSurfEmis(self):
        self._surfemis_manual = False
        self._surfemis = None

    dstring = """Array containing input/output surface emissivity values.

    This must be a ``numpy.ndarray`` of shape [nprofiles, nchannels], the
    size of the array cannot be checked here. However, it will
    eventually be checked during the :meth:`runDirect` or :meth:`runK` calls.

    The `dtype` of the ``numpy.ndarray`` must be {!r}. If the input array has
    another `dtype`, a conversion will be attempted.

    Use `del` to empty the SurfEmis array.

    If this array is not set by the user (or has been emptied unsig `del`), it
    will be initialised to -1 in the :meth:`runDirect` or :meth:`runK` calls.
    """.format(wrapfloat)
    SurfEmis = property(_getSurfEmis, _setSurfEmis,
                        _delSurfEmis, doc=dstring)

    def printSurfEmis(self, pfunction=print):
        """Print the :data:`SurfEmis` content on the standard output."""
        fmt_e = 'Channel #{:5d} e={:6.3f}'
        try:
            emis = self.SurfEmis
        except RttovInternalError:
            emis = None
        if emis is not None:
            for prof in range(self.Profiles.Nprofiles):
                pfunction(">> Surfemis for profile #{:d}".format(prof))
                for ich in range(emis.shape[-1]):
                    pfunction(fmt_e.format(ich + 1,
                                           emis[prof, ich]))
        else:
            pfunction(">> Surfemis is not (yet?) initialised.")

    def _prerun_commons(self, channels):
        (channels, fixed_p) = super(RttovScatt, self)._prerun_commons(channels)

        # deal with SurfEmis
        surfemis_shape = (self.Profiles.Nprofiles, len(channels))
        if self._surfemis_manual:
            # Check that the dimensions are fine
            if self.SurfEmis.shape != surfemis_shape:
                raise RttovWrapperError("Incorrect shape for SurfEmis")
        else:
            self._printVerbWrap("No surface emissivity supplied:" +
                                " setting calcemis to true")
            self._surfemis = - np.ones(surfemis_shape, dtype=wrapfloat)
        if self._debug:
            self.printSurfEmis(pfunction=self._printDbg)

        # bt is always required
        self._prerun_bBascis_init(channels, ('bt', ))

        return (channels, fixed_p)

    def runDirect(self, channels=None):
        """Run the RTTOV-SCATT direct model.

        :param channels: list of channels to simulate (all the channels if
                         omitted)
        :type channels: list of int or numpy.ndarray
        """
        channels, p = self._prerun_commons(channels)

        err = rtwrap.rttov_scatt_call_direct(
            self._instid,
            channels,
            self.Profiles.DateTimes.transpose(),    # profile dates/times                                     [nprofiles][6]
            self.Profiles.Angles.transpose(),       # satzen, satazi angles                                   [nprofiles][2]
            self.Profiles.SurfGeom.transpose(),     # lat, lon, elevation                                     [nprofiles][3]
            self.Profiles.SurfType.transpose(),     # surftype                                                [nprofiles]
            self.Profiles.Skin.transpose(),         # skin T, salinity, foam_frac, fastem_coefsx5             [nprofiles][8]
            self.Profiles.S2m.transpose(),          # 2m p, 2m t, 2m q, 10m wind u, v                         [nprofiles][5]
            self.Profiles.Zeeman.transpose(),       # Be, cosbk                                               [nprofiles][2]
            p.transpose(),                          # pressure                                                [nprofiles][nlevels]
            self.Profiles.T.transpose(),            # temperature                                             [nprofiles][nlevels]
            self.Profiles.GasUnits,                 # units for gas profiles
            self.Profiles.GasId.transpose(),        # gas ID list                                             [ngases]
            self.Profiles.Gases.transpose(),        # gas profiles                                            [ngases][nprofiles][nlevels]
            self.Profiles.Ph.transpose(),           # pressure half-levels                                    [nprofiles][nlevels+1]
            self.Profiles.UserCfrac,                # user cfrac 
            self.MultiHydroFrac,                    # multi_hydro_frac flag
            self.CalcZef,                           # calc_zef flag
            self._surfemis.transpose(),             # input/output surface emissivities                       [nprofiles][nchannels]
            self._bBasics['bt'].transpose(),        # output BTs                                              [nprofiles][nchannels]
            )
        self._errCheck(err, 'Error in rttov_scatt_call_direct')

        self._doStoreRad(len(channels))
        self._doStoreEmisTerms(len(channels))
        
        if (self.CalcZef):
	        self._doStoreZef(self.Profiles.Nlevels, len(channels))

    def runK(self, channels=None):
        """Run the RTTOV-SCATT K model.

        :param channels: list of channels to simulate (all the channels if
                         omitted)
        :type channels: list of int or numpy.ndarray
        """
        channels, p = self._prerun_commons(channels)

        self._printDbg("Allocating the many _k arrays")

        def ini_k_like(tab):
            theshape = list(tab.shape)
            # An extra dimension is needed for K computations
            theshape.insert(-1, len(channels))
            return np.empty(theshape, dtype=tab.dtype)

        self._bBasics['skin_k'] = ini_k_like(self.Profiles.Skin)
        self._bBasics['s2m_k'] = ini_k_like(self.Profiles.S2m)
        self._bBasics['p_k'] = ini_k_like(self.Profiles.P)
        self._bBasics['ph_k'] = ini_k_like(self.Profiles.Ph)
        self._bBasics['usercfrac_k'] = np.zeros_like(self._bBasics['bt'])
        self._bBasics['t_k'] = ini_k_like(self.Profiles.T)
        self._bBasics['gases_k'] = ini_k_like(self.Profiles.Gases)

        self._bBasics['surfemis_k'] = np.zeros_like(self._surfemis)
        if (self.CalcZef):
           self._bBasics['bt_k'] = np.zeros_like(self._bBasics['bt'])
           self._bBasics['zef_k'] = np.ones((self.Profiles.Nprofiles, len(channels), self.Profiles.Nlevels),dtype=wrapfloat)
        else:
           self._bBasics['bt_k'] = np.ones_like(self._bBasics['bt'])
           self._bBasics['zef_k'] = np.zeros((self.Profiles.Nprofiles, len(channels), self.Profiles.Nlevels),dtype=wrapfloat)

        err = rtwrap.rttov_scatt_call_k(
            self._instid,
            channels,
            self.Profiles.DateTimes.transpose(),        # profile dates/times                                     [nprofiles][6]
            self.Profiles.Angles.transpose(),           # satzen, satazi angles                                   [nprofiles][2]
            self.Profiles.SurfGeom.transpose(),         # lat, lon, elevation                                     [nprofiles][3]
            self.Profiles.SurfType.transpose(),         # surftype                                                [nprofiles]
            self.Profiles.Skin.transpose(),             # skin T, salinity, foam_frac, fastem_coefsx5             [nprofiles][8]
            self._bBasics['skin_k'].transpose(),        # skin T, salinity, foam_frac, fastem_coefsx5 K           [nprofiles][nchannels][8]
            self.Profiles.S2m.transpose(),              # 2m p, 2m t, 2m q, 10m wind u, v                         [nprofiles][5]
            self._bBasics['s2m_k'].transpose(),         # 2m p, 2m t, 2m q, 10m wind u, v K                       [nprofiles][nchannels][5]
            self.Profiles.Zeeman.transpose(),           # Be, cosbk                                               [nprofiles][2]
            p.transpose(),                              # pressure                                                [nprofiles][nlevels]
            self._bBasics['p_k'].transpose(),           # pressure K                                              [nprofiles][nchannels][nlevels]
            self.Profiles.T.transpose(),                # temperature                                             [nprofiles][nlevels]
            self._bBasics['t_k'].transpose(),           # temperature K                                           [nprofiles][nchannels][nlevels]
            self.Profiles.GasUnits,                     # units for gas profiles
            self.Profiles.GasId.transpose(),            # gas ID list                                             [ngases]
            self.Profiles.Gases.transpose(),            # gas profiles                                            [ngases][nprofiles][nlevels]
            self._bBasics['gases_k'].transpose(),       # gas profiles K                                          [ngases][nprofiles][nchannels][nlevels]
            self.Profiles.Ph.transpose(),               # pressure half-levels                                    [nprofiles][nlevels+1]
            self._bBasics['ph_k'].transpose(),          # pressure half-levels K                                  [nprofiles][nchannels][nlevels+1]
            self.Profiles.UserCfrac,                    # user cfrac                                              [nprofiles]
            self._bBasics['usercfrac_k'].transpose(),   # user cfrac K                                            [nprofiles][nchannels]
            self.MultiHydroFrac,                        # multi_hydro_frac flag
            self.CalcZef,                               # calc_zef flag
            self._surfemis.transpose(),                 # input/output surface emissivities                       [nprofiles][nchannels]
            self._bBasics['surfemis_k'].transpose(),    # input/output surface emissivities K                     [nprofiles][nchannels]
            self._bBasics['bt'].transpose(),            # output BTs                                              [nprofiles][nchannels]
            self._bBasics['bt_k'].transpose(),          # input BT perturbations                                  [nprofiles][nchannels]
            self._bBasics['zef_k'].transpose(),         # input zef perturbations                                 [nprofiles][nchannels][nlevels]
            )
        self._errCheck(err, 'Error in rttov_scatt_call_k')

        # Save the gas_id list for a later use
        self._bBasics['gases_k_id'] = self.Profiles.GasId

        self._doStoreRad(len(channels))

        if (self.CalcZef):
            self._doStoreZef(self.Profiles.Nlevels, len(channels))

    def _doStoreRad(self, nchannels):
        if not self.Options.StoreRad:
            self._bRad = dict()
            return
        self._1DstoreGeneric(self._bRad, 'bt_clear',
                             rtwrap.rttov_get_bt_clear, nchannels)
        self._1DstoreGeneric(self._bRad, 'bt',
                             rtwrap.rttov_get_bt, nchannels)
        self._1DstoreGeneric(self._bRad, 'rad_quality',
                             rtwrap.rttov_get_rad_quality, nchannels,
                             dtype=wrapint)
        self._storeLogical(self._bRad, 'plane_parallel',
                             rtwrap.rttov_get_plane_parallel)
        self._2DstoreGeneric(self._bRad, 'geometric_height',
                             rtwrap.rttov_get_geometric_height, nchannels)

    def _doStoreEmisTerms(self, nchannels):
        if not self.Options.StoreEmisTerms:
            self._bEmisTerms = dict()
            return
        self._1DstoreGeneric(self._bEmisTerms, 'emis_terms_cfrac',
                             rtwrap.rttov_get_emis_terms_cfrac, nchannels)
        self._1DstoreGeneric(self._bEmisTerms, 'emis_terms_bsfc',
                             rtwrap.rttov_get_emis_terms_bsfc, nchannels)
        self._1DstoreGeneric(self._bEmisTerms, 'emis_terms_tau_cld',
                             rtwrap.rttov_get_emis_terms_tau_cld, nchannels)
        self._1DstoreGeneric(self._bEmisTerms, 'emis_terms_up_cld',
                             rtwrap.rttov_get_emis_terms_up_cld, nchannels)
        self._1DstoreGeneric(self._bEmisTerms, 'emis_terms_down_cld',
                             rtwrap.rttov_get_emis_terms_down_cld, nchannels)
        self._1DstoreGeneric(self._bEmisTerms, 'emis_terms_tau_clr',
                             rtwrap.rttov_get_emis_terms_tau_clr, nchannels)
        self._1DstoreGeneric(self._bEmisTerms, 'emis_terms_up_clr',
                             rtwrap.rttov_get_emis_terms_up_clr, nchannels)
        self._1DstoreGeneric(self._bEmisTerms, 'emis_terms_down_clr',
                             rtwrap.rttov_get_emis_terms_down_clr, nchannels)
    
    def _doStoreZef(self, nlevels, nchannels):
        self._2DlevelsstoreGeneric(self._bZef, 'zef',
                             rtwrap.rttov_get_zef, nlevels, nchannels)
        self._2DlevelsstoreGeneric(self._bZef, 'azef',
                             rtwrap.rttov_get_azef, nlevels, nchannels)

    def getHydroNK(self, n):
        """Computed Jacobians for hydrometeor n from the previous run. ``numpy.ndarray`` of shape [nprofiles, nchannels, nlevels]."""
        if n < 0 or n > maxscatthydro:
            raise ValueError("Hydrometeor index must lie between 1 and "+str(maxscatthydro))
        return getattr(self, 'Hydro'+str(n)+'K')

    def getHydroFracNK(self, n):
        """Computed Jacobians for hydrometeor cloud fraction n from the previous run. ``numpy.ndarray`` of shape [nprofiles, nchannels, nlevels]."""
        if n < 0 or n > maxscatthydro:
            raise ValueError("Hydrometeor index must lie between 1 and "+str(maxscatthydro))
        return getattr(self, 'HydroFrac'+str(n)+'K')

    dstring = ("Computed brightness temperatures from the " +
               "previous run. ``numpy.ndarray`` of shape " +
               "[nprofiles, nchannels].")
    Bt = _DescRO('bBasics', 'bt', doc=dstring)

    dstring = ("Computed pressure half-level Jacobians from the previous run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, nlevels+1].")
    PhK = _DescRO('bBasics', 'ph_k', doc=dstring)

    dstring = ("Computed user cloud fraction Jacobians from the previous run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels].")
    UserCfracK = _DescRO('bBasics', 'usercfrac_k', doc=dstring)

    dstring = ("Computed surface emissivity Jacobians from the previous run. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels].")
    SurfEmisK = _DescRO('bBasics', 'surfemis_k', doc=dstring)

    dstring = ("Computed skin variables Jacobians from the previous run " +
               "[Skin T, salinity, foam_fraction, " +
               "fastem_coefs(1:5)]. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, 8].")
    SkinK = _DescRO('bBasics', 'skin_k', doc=dstring)

    dstring = ("Computed 2m variables Jacobians from the previous run " +
               "[2m p, 2m t, 2m q, 10m wind u, v]. " +
               "``numpy.ndarray`` of shape [nprofiles, nchannels, 5].")
    S2mK = _DescRO('bBasics', 's2m_k', doc=dstring)

    # Outputs under the store_rad key
    fdstring = ('RTTOV radiance {} output ``numpy.ndarray`` of shape' +
                '[nprofiles, nchannels{}]. Requires store_rad=True.')

    dstring = fdstring.format('bt_clear', '')
    BtClear = _DescRO('bRad', 'bt_clear', doc=dstring)

    dstring = fdstring.format('quality', '')
    RadQuality = _DescRO('bRad', 'rad_quality', doc=dstring)

    dstring = fdstring.format('plane_parallel', '')
    PlaneParallel = _DescRO('bRad', 'plane_parallel', doc=dstring)

    dstring = fdstring.format('geometric_height', ', nlevels')
    GeometricHeight = _DescRO('bRad', 'geometric_height', doc=dstring)

    # Outputs under the store_emis_terms key
    fdstring = ('RTTOV-SCATT emis retrieval terms {} output ``numpy.ndarray`` ' +
                'of shape [nchannels]. Requires store_emis_terms=True.')

    dstring = fdstring.format('cfrac')
    EmisTermsCfrac = _DescRO('bEmisTerms', 'emis_terms_cfrac', doc=dstring)

    dstring = fdstring.format('bsfc')
    EmisTermsBsfc = _DescRO('bEmisTerms', 'emis_terms_bsfc', doc=dstring)

    dstring = fdstring.format('tau_cld')
    EmisTermsTauCld = _DescRO('bEmisTerms', 'emis_terms_tau_cld', doc=dstring)

    dstring = fdstring.format('up_cld')
    EmisTermsUpCld = _DescRO('bEmisTerms', 'emis_terms_up_cld', doc=dstring)

    dstring = fdstring.format('down_cld')
    EmisTermsDownCld = _DescRO('bEmisTerms', 'emis_terms_down_cld', doc=dstring)

    dstring = fdstring.format('tau_clr')
    EmisTermsTauClr = _DescRO('bEmisTerms', 'emis_terms_tau_clr', doc=dstring)

    dstring = fdstring.format('up_clr')
    EmisTermsUpClr = _DescRO('bEmisTerms', 'emis_terms_up_clr', doc=dstring)

    dstring = fdstring.format('down_clr')
    EmisTermsDownClr = _DescRO('bEmisTerms', 'emis_terms_down_clr', doc=dstring)
    
    # Outputs under the store_zef key
    fdstring = ('RTTOV reflectivities {} output ``numpy.ndarray`` of shape' +
                '[nprofiles, nchannels, nlevels].')

    dstring = fdstring.format('zef')
    Zef = _DescRO('bZef', 'zef', doc=dstring)

    dstring = fdstring.format('azef')
    AZef = _DescRO('bZef', 'azef', doc=dstring)

@lock_attributes
class Atlas(object):
    '''The Atlas class holds data for an atlas for a single month and,
    where relevant, a single instrument.
    '''

    AtlasPath = _DPathDescRWD('atlas_conf', 'atlas_path',
                              doc="Path to the atlas files.")
    Verbose = _TDescRW('atlas_conf', 'verbose', bool,
                       doc='Verbosity flag for atlas.')
    IncLand = _TDescRW('atlas_conf', 'inc_land', bool,
                       doc='Flag to indicate atlas should be used for land surface profiles.')
    IncSeaIce = _TDescRW('atlas_conf', 'inc_seaice', bool,
                         doc='Flag to indicate atlas should be used for sea-ice surface profiles.')
    IncSea = _TDescRW('atlas_conf', 'inc_sea', bool,
                      doc='Flag to indicate atlas should be used for sea surface profiles.')

    def __init__(self, verbose=True):
        """
        :param bool verbose: set verbosity flag
        """
        self._atlas_conf = {}
        self._atlas_wrap_id = -1

        self.Verbose = verbose
        self.IncLand = True
        self.IncSeaIce = True
        self.IncSea = True

    def __del__(self):
        self.dropAtlas()

    def _printVerbWrap(self, msg):
        if self.Verbose:
            sys.stderr.write(msg + "\n")

    def isAtlasLoaded(self):
        return (self._atlas_wrap_id > 0)

    def dropAtlas(self):
        if self._atlas_wrap_id > 0:
            err = rtwrap.rttov_drop_atlas(self._atlas_wrap_id)
            self._atlas_wrap_id = -1
            self._printVerbWrap("Atlas deallocated.")

    def _GenericAtlasSetupCheck(self, atlas, month, inst, atlas_id=-1):
        if month < 1 or month > 12:
            self._printVerbWrap('Warning: Incorrect month number.')
            return False

        if self.isAtlasLoaded():
            self._printVerbWrap('Warning: this atlas object has already been initialised.')
            return False

        if inst is not None:
            if not isinstance(inst, Rttov) and \
               not (isinstance(inst, RttovScatt) and atlas[:2] == 'MW'):
                self._printVerbWrap('Warning: instrument argument is of wrong type, ' +
                                    '{} atlas not initialised'.format(atlas))
                return False
            if not inst.CoeffsLoaded:
                self._printVerbWrap('Warning: instrument must be loaded, ' +
                                    '{} atlas not initialised'.format(atlas))
                return False
        elif atlas[:2] == 'MW' and atlas_id == 2:
            self._printVerbWrap('Warning: CNRM atlas requires loaded instrument, ' +
                                '{} atlas not initialised'.format(atlas))
            return False

        if not self.AtlasPath:
            self._printVerbWrap(('Warning: AtlasPath not set. ' +
                                 '{} atlas not initialised').format(atlas))
            return False

        return True

    def _GenericAtlasSetupFinalise(self, atlas, atlas_wrap_id):
        if atlas_wrap_id < 1:
            self._printVerbWrap('Warning: {} atlas not initialised'.format(atlas))
            return False
        else:
            self._printVerbWrap('{} atlas loaded successfully'.format(atlas))
            return True

    def loadBrdfAtlas(self, month, inst=None, atlas_id=-1):
        """Load data from the BRDF atlas.

        :param int month: month (1-12 => Jan-Dec) for which atlas should be
            initialised
        :param Rttov inst: if present initialise atlas for this
            instrument, optional
        :param int atlas_id: ID for the BRDF atlas, optional
        """
        if self._GenericAtlasSetupCheck('BRDF', month, inst):
            inst_id = inst.InstId if inst is not None else -1
            self._atlas_wrap_id = rtwrap.rttov_load_brdf_atlas(self.AtlasPath,
                                                month, atlas_id, inst_id)
            return self._GenericAtlasSetupFinalise('BRDF', self._atlas_wrap_id)
        else:
            return False

    def loadIrEmisAtlas(self, month, inst=None, ang_corr=False, atlas_id=-1):
        """Load data from an IR emissivity atlas.

        :param int month: month (1-12 => Jan-Dec) for which atlas should be
            initialised
        :param Rttov inst: if present initialise atlas for this
            instrument, optional
        :param bool ang_corr: apply zenith angle correction to IR emissivities,
            optional
        :param int atlas_id: ID for IR atlas, optional
        """
        if self._GenericAtlasSetupCheck('IR emissivity', month, inst):
            inst_id = inst.InstId if inst is not None else -1
            self._atlas_wrap_id = rtwrap.rttov_load_ir_emis_atlas(self.AtlasPath,
                                                month, atlas_id, inst_id, int(ang_corr))
            return self._GenericAtlasSetupFinalise('IR emissivity', self._atlas_wrap_id)
        else:
            return False

    def loadMwEmisAtlas(self, month, inst=None, year=0, atlas_id=-1):
        """Load data from an MW emissivity atlas.

        :param int month: month (1-12 => Jan-Dec) for which atlas should be
            initialised
        :param inst: if present initialise atlas for this
            instrument, optional (required for CNRM atlas, ignored by TELSEM2)
        :type inst: :class:`.Rttov`, :class:`.RttovScatt`
        :param int year: year of atlas (CNRM atlas only, uses default if
            value less than 1970)
        :param int atlas_id: ID for MW atlas, optional
        """
        if self._GenericAtlasSetupCheck('MW emissivity', month, inst, atlas_id):
            inst_id = inst.InstId if inst is not None else -1
            self._atlas_wrap_id = rtwrap.rttov_load_mw_emis_atlas(self.AtlasPath,
                                                month, atlas_id, inst_id, year)
            return self._GenericAtlasSetupFinalise('MW emissivity', self._atlas_wrap_id)
        else:
            return False

    def getEmisBrdf(self, inst, channels=None):
        """Return emissivities or BRDFs from atlas for given :class:`.Rttov`
        or :class:`.RttovScatt` object inst.

        :param inst: instrument for which to retrieve emissivity/BRDF values
        :type inst: :class:`.Rttov`, :class:`.RttovScatt`
        :param channels: list of channels for which to return values
                         (all loaded channels in Inst if omitted)
        """
        if not self.isAtlasLoaded():
            raise RttovWrapperError("Error: atlas data not loaded")
        if not isinstance(inst, (Rttov, RttovScatt)):
            raise RttovWrapperError("Error: inst argument must be an Rttov (IR, MW, BRDF) or " +
                                    "RttovScatt (MW only) instance")
        if not inst.CoeffsLoaded:
            raise RttovWrapperError("Error: coefficients not loaded in inst")
        if inst._profiles is None:
            raise RttovWrapperError("Error: profiles not set in inst")

        profiles = inst.Profiles

        inst_id = inst.InstId
        nprofiles = profiles.Nprofiles

        if channels is None:
            channels = np.arange(1, inst.Nchannels + 1, dtype=wrapint)
        if (not isinstance(channels, np.ndarray) or
                channels.dtype is not wrapint):
            channels = np.array(channels, dtype=wrapint)
        nchannels = len(channels)

        latitude      = profiles.SurfGeom[:, 0]
        longitude     = profiles.SurfGeom[:, 1]
        zenangle      = profiles.Angles[:, 0]
        azangle       = profiles.Angles[:, 1]

        if isinstance(inst, RttovScatt):
            surftype      = profiles.SurfType[:]
            watertype     = np.ones((nprofiles,), dtype=wrapint)
            sunzenangle   = np.zeros((nprofiles,), dtype=wrapfloat)
            sunazangle    = np.zeros((nprofiles,), dtype=wrapfloat)
            snow_fraction = np.zeros((nprofiles,), dtype=wrapfloat)
        else:
            surftype      = profiles.SurfType[:, 0]
            watertype     = profiles.SurfType[:, 1]
            sunzenangle   = profiles.Angles[:, 2]
            sunazangle    = profiles.Angles[:, 3]
            snow_fraction = profiles.Skin[:, 2]

        emisbrdf = np.zeros((nprofiles, nchannels), dtype=wrapfloat)

        err = rtwrap.rttov_get_emisbrdf(self._atlas_wrap_id, latitude, longitude, surftype,
                  watertype, zenangle, azangle, sunzenangle, sunazangle, snow_fraction,
                  inst_id, channels.transpose(), emisbrdf.transpose())

        if err == 0:
            for p in range(nprofiles):
                if not ((surftype[p] == 0 and self.IncLand) or
                        (surftype[p] == 1 and self.IncSea)  or
                        (surftype[p] == 2 and self.IncSeaIce)):
                    emisbrdf[p][:] = -1.
        else:
            raise RttovInternalError("Error in rttov_get_emisbrdf")

        return emisbrdf
