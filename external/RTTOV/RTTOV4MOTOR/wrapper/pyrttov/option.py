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
import collections
from pyrttov.decorator import add_descriptors_opts, lock_attributes

(BAND1, BAND2, BAND3) = range(1, 4)
(INDEX1, INDEX2, INDEX3, INDEX4) = range(1, 5)
(TESSEM2, FASTEM1, FASTEM2, FASTEM3, FASTEM4, FASTEM5, FASTEM6) = range(0, 7)
(ISEM, IREMIS) = range(1, 3)
(JONSWAP, ELFOUHAILY) = range(1, 3)
(MAXRANDOM_OVERLAP, SIMPLE_OVERLAP) = range(1, 3)
(DOM_IR, CHOU) = range(1, 3)
(DOM_VIS, SINGLESCATT, MFASIS) = range(1, 4)
(LIEBE, ROSENKRANZ, TKC) = range(1, 4)


@lock_attributes
@add_descriptors_opts(['wrapper', 'config', 'interpolation', 'rt_all', 'rt_mw', \
                       'rt_ir', 'rt_ir_pc', 'scatt', 'dev'])
class Options(object):
    '''The Options class holds and manages the RTTOV options.'''

    # Defaults for wrapper options
    _defaults_wrapper = {
        'nthreads': 1,
        'nprofs_per_call': 1,
        'verbose_wrapper': False,
        'check_opts': True,
        'store_rad': False,
        'store_rad2': False,
        'store_trans': False,
        'store_emis_terms': False,
    }
    # Defaults for RTTOV configuration: config
    _defaults_config = {
        'apply_reg_limits': False,
        'verbose': True,
        'do_checkinput': True,
        'fix_hgpl': True,
    }
    # Defaults for RTTOV configuration: interpolation
    _defaults_interpolation = {
        'addinterp': True,
        'spacetop': True,
        'lgradp': False,
        'reg_limit_extrap': True,
        'interp_mode': 1,
    }
    # Defaults for RTTOV configuration: rt_all
    _defaults_rt_all = {
        'ozone_data': False,
        'co2_data': False,
        'n2o_data': False,
        'co_data': False,
        'ch4_data': False,
        'so2_data': False,
        'addrefrac': True,
        'plane_parallel': False,
        'switchrad': False,
        'use_t2m_opdep': True,
        'use_q2m': True,
        'do_lambertian': False,
        'lambertian_fixed_angle': True,
        'rad_down_lin_tau': True,
        'dtau_test': False,
    }
    # Defaults for RTTOV configuration: rt_mw
    _defaults_rt_mw = {
        'clw_data': False,
        'clw_scheme': ROSENKRANZ,
        'clw_cloud_top': 322.,
        'supply_foam_fraction': False,
        'fastem_version': FASTEM6,
    }
    # Defaults for RTTOV configuration: rt_ir
    _defaults_rt_ir = {
        'solar_sea_brdf_model': ELFOUHAILY,
        'ir_sea_emis_model': IREMIS,
        'addsolar': False,
        'rayleigh_max_wavelength': 2.,
        'rayleigh_min_pressure': 0.,
        'rayleigh_single_scatt': True,
        'do_nlte_correction': False,
        'addaerosl': False,
        'user_aer_opt_param': False,
        'addclouds': False,
        'user_cld_opt_param': False,
        'grid_box_avg_cloud': True,
        'cldcol_threshold': -1.,
        'cloud_overlap': MAXRANDOM_OVERLAP,
        'cc_low_cloud_top': 750.,
        'ir_scatt_model': CHOU,
        'vis_scatt_model': DOM_VIS,
        'dom_nstreams': 8,
        'dom_accuracy': 0.,
        'dom_opdep_threshold': 0.,
        'dom_rayleigh': False,
    }
    # Defaults for RTTOV configuration: rt_pc
    _defaults_rt_ir_pc = {
        'addpc': False,
        'addradrec': False,
        'ipcreg': INDEX1,
        'ipcbnd': BAND1,
        'npcscores': -1,
    }
    # Defaults for RTTOV-SCATT configuration: opts_scatt
    _defaults_scatt = {
        'lusercfrac': False,
        'cc_threshold': 0.001,
        'ice_polarisation': 1.40,
        'hydro_cfrac_tlad': True,
        'zero_hydro_tlad': False,
    }
    # Defaults for developer configuration: dev
    _defaults_dev = {
        'do_opdep_calc': True,
    }
    # Override the naming convention for some of the descriptors
    _naming_override = {
        'addpc': 'AddPC',
        'npcscores': 'NPCScores',
        'clw_data': 'CLWData',
        'clw_scheme': 'CLWScheme',
        'clw_cloud_top': 'CLWCloudTop',
        'cc_low_cloud_top': 'CCLowCloudTop',
        'co2_data': 'CO2Data',
        'ch4_data': 'CH4Data',
        'co_data': 'COData',
        'n2o_data': 'N2OData',
        'so2_data': 'SO2Data',
        'lusercfrac': 'LuserCfrac',
        'cc_threshold': 'CCThreshold',
        'hydro_cfrac_tlad': 'HydroCfracTLAD',
        'zero_hydro_tlad': 'ZeroHydroTLAD',
    }

    def __init__(self):
        # Call the function define by the decorator in order to initialise
        # the defaults
        self._initDefaults()

    def _initDefaults(self):
        '''May be overwritten by the class decorator.'''
        pass

    @staticmethod
    def _prepareValue(val):
        '''Transform option's values if needed.'''
        if isinstance(val, bool):
            return str(int(val))
        else:
            return str(val)

    def _prepareOpts(self):
        '''Prepare the option string.'''
        opts_processed = collections.OrderedDict()
        for confkey, thedict in {'opts%config': self._config,
                                 'opts%interpolation': self._interpolation,
                                 'opts%rt_all': self._rt_all,
                                 'opts%rt_mw': self._rt_mw,
                                 'opts%rt_ir': self._rt_ir,
                                 'opts%rt_ir%pc': self._rt_ir_pc,
                                 'opts_scatt': self._scatt,
                                 'opts%dev': self._dev}.items():
            for varkey in sorted(thedict.keys()):
                opts_processed[confkey + '%' + varkey] = \
                    self._prepareValue(thedict[varkey])
        for varkey in sorted(self._wrapper.keys()):
            opts_processed[varkey] = self._prepareValue(self._wrapper[varkey])
        if not self.AddClouds:
            del opts_processed['opts%rt_ir%cloud_overlap']
            del opts_processed['opts%rt_ir%cldcol_threshold']
            del opts_processed['opts%rt_ir%cc_low_cloud_top']
        if not self.AddClouds and not self.AddAerosl:
            del opts_processed['opts%rt_ir%ir_scatt_model']
            del opts_processed['opts%rt_ir%vis_scatt_model']
            del opts_processed['opts%rt_ir%dom_nstreams']
            del opts_processed['opts%rt_ir%dom_accuracy']
            del opts_processed['opts%rt_ir%dom_opdep_threshold']
        if not self.AddPC:
            del opts_processed['opts%rt_ir%pc%ipcbnd']
            del opts_processed['opts%rt_ir%pc%ipcreg']
            del opts_processed['opts%rt_ir%pc%npcscores']
        return opts_processed

    def defineStrOptions(self):
        '''Returns an option string understandable by RTTOV's wrapper.'''
        return ' '.join(['{} {!s}'.format(key, val)
                         for key, val in self._prepareOpts().items()])

    def __str__(self):
        return '\n'.join(['{:40s}  {!s}'.format(key, val)
                          for key, val in self._prepareOpts().items()])
