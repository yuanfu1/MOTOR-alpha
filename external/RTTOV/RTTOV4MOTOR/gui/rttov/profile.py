#! /usr/bin/env python
# -*- coding: utf-8 -*-
try:
    import h5py
except ImportError:
    import sys
    sys.stderr.write('ERROR: h5py is not installed\n')
    sys.exit(1)
try:
    import numpy
except ImportError:
    import sys
    sys.stderr.write('ERROR: numpy is not installed\n')
    sys.exit(1)
import re

from rttov.core import loadDset_arr, loadDset_str, \
    getAttributes, saveAttributes, checkAttributes, cleanpath
from rttov_gui_f2py import rttov_gui_aer_clim_prof

from rttov.default import *
import numpy as np

# conversion constants for gas content (used when adding gas):
q_mixratio_to_ppmv = 1.60771704e+6
o3_mixratio_to_ppmv = 6.03504e+5
co2_mixratio_to_ppmv = 6.58114e+5
co_mixratio_to_ppmv = 1.0340699e+6
n2o_mixratio_to_ppmv = 6.58090e+5
ch4_mixratio_to_ppmv = 1.80548e+6
q_mixratio_to_ppmv = 1.60771704e+6
so2_mixratio_to_ppmv = 4.52116e+05
mixratio_to_ppmv = {"Q": q_mixratio_to_ppmv,
                    "O3": o3_mixratio_to_ppmv,
                    "CO2": co2_mixratio_to_ppmv,
                    "CO": co_mixratio_to_ppmv,
                    "N2O": n2o_mixratio_to_ppmv,
                    "CH4": ch4_mixratio_to_ppmv,
                    "Q": q_mixratio_to_ppmv,
                    "SO2": so2_mixratio_to_ppmv,
                    }


def fct_mixratio_to_ppmv(value, gas):
    return value * mixratio_to_ppmv[gas]


def fct_ppmv_to_mixratio(value, gas):
    return value / mixratio_to_ppmv[gas]


def getNumberOfProfiles(fileName):
    """
    Get the number of profiles in and HDF5 file
    If the /PROFILE is found number of profiles is 1
    If the /PROFILES is found then the number of profiles
      is the number of subdirectories which names are formed of 4 digits
    """
    f = h5py.File(fileName, 'r')
    nprof = 0
    keys = list(f.keys())
    if("PROFILE" in keys):
        nprof = 1
    elif("PROFILES" in keys):
        for name, value in f['PROFILES'].items():
            if (re.search('[0-9]{4}', name)):
                nprof += 1
    f.close()
    return nprof


class Profile(dict):
    """
    The Profile class
    A Python dictionary which contains a copy of the RTTOV v13
    rttov_profile structure
    All members are present, for example if the profile does not
    contain any CO2
    the CO2 key is present and is set to the Python None
    Methods are:
    - loadh5 which returns a full profile from a given dataset
    - loadprofilenumber which returns a profile defined by its number
      from a HDF5 file containing several profiles
    - display display profile to terminal
    - saveh5 saves a profile class to a HDF5 file
    - addgas will add a new gas to the profile filled witha  default
      concentration or a user value
    - removegas remove a gas from the profile
    - removeallbutgas remove all gases except the given one
    """
    # Default values concentrations for aerosols and clouds
    # aerosols : 0: number density (cm-3), 1: kg/jg
    defaultaerdensity = {0: 1.e0, 1: 1.e-7}
    # cloud : 0: layer mean content(g/m3), 1: kg/jg
    defaultclddensity = {0: 1.e-01, 1: 1.e-5}

    # All known profile_type members
    gas_list = ['Q', 'O3', 'CO2', 'CH4', 'CO', 'N2O', 'SO2']
    profile_list = ['AEROSOLS', 'AZANGLE', 'BE', 'CFRAC', 'CFRACTION', 'CH4',
                    'CLOUD',
                    'CLW', 'CO', 'CO2', 'COSBK', 'CTP',
                    'DATE', 'ELEVATION', 'GAS_UNITS', 'ICEDE', 'CLWDE', 'ID',
                    'ICEDE_PARAM', 'ICE_SCHEME', 'CLW_SCHEME', 'LATITUDE', 'LONGITUDE',
                    'MMR_CLDAER', 'N2O', 'CLWDE_PARAM',
                    'NLAYERS', 'NLEVELS', 'O3', 'P', 'Q', 'SO2',
                    'SUNAZANGLE', 'SUNZENANGLE', 'T', 'TIME', 'ZENANGLE']
    sskin_list = ['SURFTYPE', 'WATERTYPE', 'T', 'SNOW_FRACTION',
                  'SOIL_MOISTURE',
                  'SALINITY', 'FOAM_FRACTION', 'FASTEM']
    # keep above always fastem in last position
    s2m_list = ['T', 'Q', 'O', 'P', 'U', 'V', 'WFETC']
    aerosol_list = None
    aerosol_long_list = None
    opac_aerosol_list = ['INSO', 'WASO', 'SOOT', 'SSAM', 'SSCM',
                         'MINM', 'MIAM', 'MICM', 'MITR', 'SUSO',
                         'VOLA', 'VAPO', 'ASDU']
    cams_aerosol_list = ['BCAR', 'DUS1', 'DUS2', 'DUS3', 'SULP', 'SSA1',
                         'SSA2', 'SSA3', 'OMAT']
    opac_aerosol_long_list = ['Insoluble', 'Water soluble', 'Soot',
                              'Sea salt (acc mode)', 'Sea salt (coa mode)',
                              'Mineral (nuc mode)', 'Mineral (acc mode)',
                              'Mineral (coa mode)',
                              'Mineral transported', 'Sulphated droplets',
                              'OPAC Volcanic ash',
                              'New Volcanic ash', 'Asian dust']
    cams_aerosol_long_list = [
        'Black carbon, fixed refractive index at all wavelengths',
        'Dust, bin 1, 0.03-0.55 micron, refractive index: Woodward 2001',
        'Dust, bin 2, 0.55-0.90 micron, refractive index: Woodward 2001',
        'Dust, bin 3, 0.90-20.0 micron, refractive index: Woodward 2001',
        'Ammonium sulphate',
        'Sea salt, bin 1, 0.03-0.5 micron',
        'Sea salt, bin 2, 0.50-5.0 micron',
        'Sea salt, bin 3, 5.0-20.0 micron',
        'Hydrophilic organic matter']
    cloud_list = ['STCO', 'STMA', 'CUCC', 'CUCP', 'CUMA', 'CIRR']
    cloud_long_list = ['Stratus Continental', 'Stratus Maritime',
                       'Cumulus Continental Clean',
                       'Cumulus Continental Polluted', 'Cumulus Maritime',
                       'Cirrus']
    gas_units = {-1: "ppmv over dry air", 0: "ppmv over dry air",
                 1: "kg/kg over moist air", 2: "ppmv over moist air",
                 None: "ppmv over moist air"}

    cld_units = {0: "$g/m^{3}$", 1: "kg/kg", None: "kg/kg"}
    aer_units = {0: "number density ($cm^{-3}$)", 1: "kg/kg", None: "kg/kg"}

    profSurfParamList = ['DATE', 'TIME', 'ID', 'LATITUDE', 'LONGITUDE',
                         'ELEVATION', 'AZANGLE', 'ZENANGLE', 'SUNAZANGLE',
                         'SUNZENANGLE', 'CFRACTION',
                         'BE', 'COSBK', 'CTP', 'ICEDE_PARAM', 'ICE_SCHEME',
                         'CLW_SCHEME']

    maxValue = {'LATITUDE': 90., 'LONGITUDE': 180., 'ELEVATION': 9000.,
                'AZANGLE': 180., 'ZENANGLE': 90.0, 'SUNAZANGLE': 180.,
                'SUNZENANGLE': 90., 'CFRACTION': 1.,
                'BE': 0.7, 'COSBK': 1., 'CTP': 1100., 'ICEDE_PARAM': 4,
                'ICE_SCHEME': 3, 'CLW_SCHEME': 2,
                'S2M_T': 400., 'S2M_Q': 0.60E+06, 'S2M_O': 1000.,
                'S2M_P': 1100.0, 'S2M_U': 100., 'S2M_V': 100.,
                "S2M_WFETC": 100000.,
                'SKIN_SNOW_FRACTION': 1., 'SKIN_SOIL_MOISTURE': 1.,
                "SKIN_SALINITY": 100., "SKIN_FOAM_FRACTION": 1.,
                "SKIN_SURFTYPE": 2, 'SKIN_T': 400., 'SKIN_WATERTYPE': 1
                }
    minValue = {'LATITUDE': -90., 'LONGITUDE': -180.,
                'ELEVATION': -500., 'AZANGLE': -180., 'ZENANGLE': 0.0,
                'SUNAZANGLE': -180, 'SUNZENANGLE': 0.0, 'CFRACTION': 0.0,
                'BE': 0.201, 'COSBK': 0., 'CTP': 50., 'ICEDE_PARAM': 1,
                'ICE_SCHEME': 1, 'CLW_SCHEME': 1,
                'S2M_T': 90., 'S2M_Q': 0.1E-10, 'S2M_O': 0.1E-10,
                'S2M_P': 400.0, 'S2M_U': -100., 'S2M_V': -100.,
                "S2M_WFETC": 0.,
                'SKIN_SNOW_FRACTION': 0., 'SKIN_SOIL_MOISTURE': 0.,
                "SKIN_SALINITY": 0., "SKIN_FOAM_FRACTION": 0.,
                "SKIN_SURFTYPE": 0, 'SKIN_T': 90., 'SKIN_WATERTYPE': 0
                }

    typeListe = {"ICEDE_PARAM": {1: "Ou and Liou", 2: "Wyser", 3: "Boudala",
                         4: "McFarquhar"},
                 "ICE_SCHEME": {1: "Baum ice dataset", 2: "Baran 2014",
                                3: "Baran 2018"},
                 "CLW_SCHEME": {1: "OPAC", 2: "CLW Deff"},
                 "SURFTYPE": {0: "land", 1: "sea", 2: "sea-ice"},
                 "WATERTYPE": {0: "fresh water", 1: "ocean water"}}

    notUsedParams = {'S2M_O', 'SKIN_SOIL_MOISTURE'}

    notUsedParamsShort = {'O', 'SOIL_MOISTURE'}

    def __init__(self, aerosol_kind="opac"):
        dict.__init__(self)
        self['S2M'] = {}
        self['SKIN'] = {}
        # aerosol kind by default "opac" (otherway it is "cams" )
        self.aerosol_kind = aerosol_kind
        if self.aerosol_kind == "opac":
            self.aerosol_list = self.opac_aerosol_list
            self.aerosol_long_list = self.opac_aerosol_long_list
        else:
            self.aerosol_kind = "cams"
            self.aerosol_list = self.cams_aerosol_list
            self.aerosol_long_list = self.cams_aerosol_long_list

        for key in self.profile_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        for key in self.s2m_list:
            self['S2M'][key] = None
            self['S2M'][
                key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        for key in self.sskin_list:
            self['SKIN'][key] = None
            self['SKIN'][
                key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        for key in self.opac_aerosol_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        for key in self.cams_aerosol_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        for key in self.cloud_list:
            self[key] = None
            self[key + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def checkAerosolKind(self):
        kind_opac = 0
        kind_cams = 0
        for item in self.keys():
            if item in self.opac_aerosol_list:
                if self[item] is not None:
                    kind_opac = kind_opac + 1
            if item in self.cams_aerosol_list:
                if self[item] is not None:
                    kind_cams = kind_cams + 1
        if kind_cams > 0 and kind_opac > 0:
            # incoherent profile
            print("ERROR incoherent aerosol content", kind_cams, kind_opac)
            return False
        if kind_cams > kind_opac:
            self.aerosol_kind = "cams"
            self.aerosol_list = self.aerosol_list = self.cams_aerosol_list
            self.aerosol_long_list = self.cams_aerosol_long_list
        else:
            self.aerosol_kind = "opac"
            self.aerosol_list = self.opac_aerosol_list
            self.aerosol_long_list = self.opac_aerosol_long_list
        return True

    def defineUnitsforCommentAttribute(self):
        if "GAS_UNITS" in list(self.keys()):
            self.my_gas_units = self.gas_units[self["GAS_UNITS"]]
        else:
            self.my_gas_units = self.gas_units[2]
        if "MMR_CLDAER" in list(self.keys()):
            if self["MMR_CLDAER"] is None:
                self["MMR_CLDAER"] = 1
            # need to take care of the possibility of having a ndarray
            # for self["MMR_CLDAER"] (kmatrix)
            if np.isscalar(self["MMR_CLDAER"]):
                self.my_cld_units = self.cld_units[self["MMR_CLDAER"]]
                self.my_aer_units = self.aer_units[self["MMR_CLDAER"]]
            else:
                self.my_cld_units = self.cld_units[self["MMR_CLDAER"][0]]
                self.my_aer_units = self.aer_units[self["MMR_CLDAER"][0]]
        else:
            self.my_cld_units = self.cld_units[1]
            self.my_aer_units = self.aer_units[1]

    def setDefaultAttributes(self):
        """
         Set the default values for ALL the profile variables even if
         variables are set to None
        """
        self.defineUnitsforCommentAttribute()

        self["ASDU_ATTRIBUTE"] = {
            'COMMENT': 'Asian dust', 'UNITS': self.my_aer_units}
        self["AZANGLE_ATTRIBUTE"] = {
            'COMMENT': 'Satellite azimuth angle',
            'UNITS': 'degree (0-360 deg; east=90)'}
        self["BE_ATTRIBUTE"] = {
            'COMMENT': 'Earth magnetic field strength', 'UNITS': 'Gauss'}
        self["CFRACTION_ATTRIBUTE"] = {
            'COMMENT': 'Black body cloud fraction (0 - 1) 1'
                       ' for 100% cloud cover',
            'UNITS': '0:1'}
        self["CFRAC_ATTRIBUTE"] = {
            'COMMENT': 'Cloud fraction (0 - 1) 1 for 100% cloud cover',
            'UNITS': '0:1'}
        self["CH4_ATTRIBUTE"] = {
            'COMMENT': 'Methan (CH4)', 'UNITS': self.my_gas_units}
        self["CIRR_ATTRIBUTE"] = {
            'COMMENT': 'Cirrus', 'UNITS': self.my_cld_units}
        self["CLW_ATTRIBUTE"] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        self["CO2_ATTRIBUTE"] = {
            'COMMENT': 'Carbon dioxide (CO2)', 'UNITS': self.my_gas_units}
        self["SO2_ATTRIBUTE"] = {
            'COMMENT': 'Sulfur dioxide (SO2)', 'UNITS': self.my_gas_units}
        self["COSBK_ATTRIBUTE"] = {
            'COMMENT': ('Cosine of the angle between the Earth '
                        'magnetic field and wave propagation direction'),
            'UNITS': 'n/a'}
        self["CO_ATTRIBUTE"] = {
            'COMMENT': 'Carbon monoxide (CO)', 'UNITS': self.my_gas_units}
        self["CTP_ATTRIBUTE"] = {
            'COMMENT': 'Black body cloud top pressure  (hPa)', 'UNITS': 'hPa'}
        self["CUCC_ATTRIBUTE"] = {
            'COMMENT': 'Cumulus Continental Clean', 'UNITS': self.my_cld_units}
        self["CUCP_ATTRIBUTE"] = {
            'COMMENT': 'Cumulus Continental Polluted',
            'UNITS': self.my_cld_units}
        self["CUMA_ATTRIBUTE"] = {
            'COMMENT': 'Cumulus Maritime', 'UNITS': self.my_cld_units}
        self["DATE_ATTRIBUTE"] = {'COMMENT': 'Year Month Day'}
        self["ELEVATION_ATTRIBUTE"] = {
            'COMMENT': 'Surface elevaton', 'UNITS': 'km'}
        self["GAS_UNITS_ATTRIBUTE"] = {
            'COMMENT': 'Units for gas profiles :'
                       '0 or less => ppmv over dry air, '
                       '1 => kg/kg over moist air, '
                       '2 => ppmv over moist air)', 'UNITS': 'n/a'}
        self["ICEDE_ATTRIBUTE"] = {
            'COMMENT': 'Ice crystals diameter', 'UNITS': 'microns'}
        self["CLWDE_ATTRIBUTE"] = {
            'COMMENT': 'Cloud liquid water particle', 'UNITS': 'microns'}
        self["CLWDE_PARAM_ATTRIBUTE"] = {
            'COMMENT': 'Cloud liquid water effective diameter parameterisation',
            'UNITS': 'n/a'}
        self["ICEDE_PARAM_ATTRIBUTE"] = {
            'COMMENT': ('Scheme for IWC to eff diameter,Dg 1=Ou and Liou; 2='
                        'Wyser et al.; 3=Boudala et al; 4=McFarquhar et al.'),
            'UNITS': 'n/a'}
        self["ID_ATTRIBUTE"] = {
            'COMMENT': 'User may give text ID to each profile', 'UNITS': 'n/a'}
        self["INSO_ATTRIBUTE"] = {
            'COMMENT': 'Insoluble', 'UNITS': self.my_aer_units}
        self["ICE_SCHEME_ATTRIBUTE"] = {
            'COMMENT': (
                'ice scheme,1=Baum ice database;2=Baran 2014;3=Baran 2018'),
            'UNITS': 'n/a'}
        self["CLW_SCHEME_ATTRIBUTE"] = {
            'COMMENT': (
                'cloud water scheme, 1=OPAC ;2=CLW Deff'),
            'UNITS': 'n/a'}
        self["LATITUDE_ATTRIBUTE"] = {
            'COMMENT': 'Latitude (deg)', 'UNITS': 'degree'}
        self["LONGITUDE_ATTRIBUTE"] = {
            'COMMENT': 'Longitude (deg)', 'UNITS': 'degree (0-360)'}
        self["MMR_CLDAER_ATTRIBUTE"] = {
            'COMMENT': 'Units for cloud and aerosol profiles '
                       '0: g/m^3 (cld); cm^-3 (aer), 1: kg/kg', 'UNITS': 'n/a'}
        self["MIAM_ATTRIBUTE"] = {
            'COMMENT': 'Mineral (acc mode)', 'UNITS': self.my_aer_units}
        self["MICM_ATTRIBUTE"] = {
            'COMMENT': 'Mineral (coa mode)', 'UNITS': self.my_aer_units}
        self["MINM_ATTRIBUTE"] = {
            'COMMENT': 'Mineral (nuc mode)', 'UNITS': self.my_aer_units}
        self["MITR_ATTRIBUTE"] = {
            'COMMENT': 'Mineral transported', 'UNITS': self.my_aer_units}
        self["N2O_ATTRIBUTE"] = {
            'COMMENT': 'Nitrous oxide (N2O)', 'UNITS': self.my_gas_units}
        self["NLAYERS_ATTRIBUTE"] = {
            'COMMENT': 'Number of atmospheric layers', 'UNITS': 'n/a'}
        self["NLEVELS_ATTRIBUTE"] = {
            'COMMENT': 'Number of atmospheric levels', 'UNITS': 'n/a'}
        self["O3_ATTRIBUTE"] = {
            'COMMENT': 'Ozone (O3)', 'UNITS': self.my_gas_units}
        self["P_ATTRIBUTE"] = {'COMMENT': 'Pressure (hPa)', 'UNITS': 'hPa'}
        self["Q_ATTRIBUTE"] = {
            'COMMENT': 'Water vapour', 'UNITS': self.my_gas_units}
        self["S2M"]["WFETC_ATTRIBUTE"] = {
            'COMMENT': 'Wind fetch (metres)', 'UNITS': 'm'}
        self["S2M"]["P_ATTRIBUTE"] = {
            'COMMENT': 'Surface pressure (hPa)', 'UNITS': 'hPa'}
        self["S2M"]["Q_ATTRIBUTE"] = {
            'COMMENT': 'Water vapour', 'UNITS': self.my_gas_units}
        self["S2M"]["U_ATTRIBUTE"] = {
            'COMMENT': 'U 10m wind component (m/s)', 'UNITS': 'm/s'}
        self["S2M"]["V_ATTRIBUTE"] = {
            'COMMENT': 'V 10m wind component (m/s)', 'UNITS': 'm/s'}
        self["S2M"]["T_ATTRIBUTE"] = {
            'COMMENT': 'Temperature (K)', 'UNITS': 'K'}
        self["S2M"]["O_ATTRIBUTE"] = {
            'COMMENT': 'Ozone', 'UNITS': self.my_gas_units}
        self["SKIN"]["FASTEM_ATTRIBUTE"] = {
            'COMMENT': 'Land/sea-ice surface parameters for fastem',
            'UNITS': 'n/a'}
        self["SKIN"]["WATERTYPE_ATTRIBUTE"] = {
            'COMMENT': '0=fresh water, 1=ocean water', 'UNITS': 'n/a'}
        self["SKIN"]["T_ATTRIBUTE"] = {
            'COMMENT': 'Radiative skin temperature (K)', 'UNITS': 'K'}
        self["SKIN"]["SURFTYPE_ATTRIBUTE"] = {
            'COMMENT': '0=land, 1=sea, 2=sea-ice', 'UNITS': 'n/a'}
        self["SKIN"]["SALINITY_ATTRIBUTE"] = {
            'COMMENT': 'Practical salinity unit %o - FASTEM-4/5 only',
            'UNITS': 'n/a'}
        self["SKIN"]["FOAM_FRACTION_ATTRIBUTE"] = {
            'COMMENT': 'Ocean foam coverage fraction passed to FASTEM',
            'UNITS': '0:1'}
        self["SKIN"]["SNOW_FRACTION_ATTRIBUTE"] = {
            'COMMENT': ('Surface snow coverage fraction (0-1).'
                        ' Used only by IR emissivity atlas'), 'UNITS': '0:1'}
        self["SKIN"]["SOIL_MOISTURE_ATTRIBUTE"] = {
            'COMMENT': 'soil moisture', 'UNITS': 'm^3/m^3'}
        self["SOOT_ATTRIBUTE"] = {'COMMENT': 'Soot',
                                  'UNITS': self.my_aer_units}
        self["SSAM_ATTRIBUTE"] = {
            'COMMENT': 'Sea salt (acc mode)', 'UNITS': self.my_aer_units}
        self["SSCM_ATTRIBUTE"] = {
            'COMMENT': 'Sea salt (coa mode)', 'UNITS': self.my_aer_units}
        self["STCO_ATTRIBUTE"] = {
            'COMMENT': 'Stratus Continental', 'UNITS': self.my_cld_units}
        self["STMA_ATTRIBUTE"] = {
            'COMMENT': 'Stratus Maritime', 'UNITS': self.my_cld_units}
        self["SUNAZANGLE_ATTRIBUTE"] = {
            'COMMENT': 'Sun azimuth angle',
            'UNITS': 'degree (0-360 deg; east=90)'}
        self["SUNZENANGLE_ATTRIBUTE"] = {
            'COMMENT': 'Sun zenith angle', 'UNITS': 'degree'}
        self["SUSO_ATTRIBUTE"] = {
            'COMMENT': 'Sulphated droplets', 'UNITS': self.my_aer_units}
        self["TIME_ATTRIBUTE"] = {'COMMENT': 'Hour Minute Second'}
        self["T_ATTRIBUTE"] = {'COMMENT': 'Temperature', 'UNITS': 'K'}
        self["VAPO_ATTRIBUTE"] = {
            'COMMENT': 'New Volcanic ash', 'UNITS': self.my_aer_units}
        self["VOLA_ATTRIBUTE"] = {
            'COMMENT': 'OPAC Volcanic ash', 'UNITS': self.my_aer_units}
        self["WASO_ATTRIBUTE"] = {
            'COMMENT': 'Water soluble', 'UNITS': self.my_aer_units}
        self["ZENANGLE_ATTRIBUTE"] = {
            'COMMENT': 'Satellite zenith angle', 'UNITS': 'degree'}

    def setDefaultProfileAsciiInput(self):
        """Put here default values for all mandatory profile variables
           The profile shall be already filled with minimum values like
           array of pressure, temperature
           water vapour
           Defaults values such as:
           - Pressure, Temperature, WaterVapour units
           - surface types and values (taken from bottom level)"""
        if self["ICEDE_PARAM"] is None:
            self["ICEDE_PARAM"] = 2   # Ice effective diameter parameterisation
        if self["ICE_SCHEME"] is None:
            self["ICE_SCHEME"] = 3  # Ice cloud scheme (Baran is default)
        if self["CLW_SCHEME"] is None:
            self["CLW_SCHEME"] = 1
        if self["GAS_UNITS"] is None:
            self["GAS_UNITS"] = 2  # ppmv over moist air

        if self["MMR_CLDAER"] is None:
            self["MMR_CLDAER"] = 1  # kg/kg

        self.defineUnitsforCommentAttribute()

        " Skin variables "
        if self["SKIN"]["T"] is None:
            self["SKIN"]["T"] = self["T"][-1]  # (K)
        if self["SKIN"]["SURFTYPE"] is None:
            self["SKIN"]["SURFTYPE"] = 1  # (0=Land, 1=Sea, 2=sea-ice)
        if self["SKIN"]["WATERTYPE"] is None:
            self["SKIN"]["WATERTYPE"] = 1  # (0=fresh water, 1=ocean water)
        if self["SKIN"]["SNOW_FRACTION"] is None:
            self["SKIN"]["SNOW_FRACTION"] = 0.  # [0,1]
        if self["SKIN"]["SOIL_MOISTURE"] is None:
            self["SKIN"]["SOIL_MOISTURE"] = 0.  # [0,1]
        if self["SKIN"]["SALINITY"] is None:
            self["SKIN"]["SALINITY"] = 37  # (%o)
        if self["SKIN"]["FASTEM"] is None:
            self["SKIN"]["FASTEM"] = numpy.array(
                [3., 5., 15., 0.1, 0.30000001])  # (5 parameters Land/sea-ice)
        if self["SKIN"]["FOAM_FRACTION"] is None:
            self["SKIN"]["FOAM_FRACTION"] = 0.0

        " 2m and 10m air variables "
        if self["S2M"]["T"] is None:
            self["S2M"]["T"] = self["T"][-1]  # (K)
        if self["S2M"]["Q"] is None:
            self["S2M"]["Q"] = self["Q"][-1]  # (ppmv or kg/kg)
        if self["S2M"]["P"] is None:
            self["S2M"]["P"] = self["P"][-1]  # (hPa)
        if self["S2M"]["U"] is None:
            self["S2M"]["U"] = 0  # (m/s)
        if self["S2M"]["V"] is None:
            self["S2M"]["V"] = 0  # (m/s)
        if self["S2M"]["WFETC"] is None:
            self["S2M"]["WFETC"] = 100000  # (m)

        " Simple cloud "
        if self["CTP"] is None:
            self["CTP"] = 500.0  # (hPa)
        if self["CFRACTION"] is None:
            self["CFRACTION"] = 0.0   # [0,1]    Clear sky is the default

        " Viewing geometry "
        if self["AZANGLE"] is None:
            self["AZANGLE"] = 0.  # (deg)
        if self["ELEVATION"] is None:
            self["ELEVATION"] = 0.  # (km)
        if self["SUNAZANGLE"] is None:
            self["SUNAZANGLE"] = 0.  # (deg)
        if self["SUNZENANGLE"] is None:
            self["SUNZENANGLE"] = 0.  # (deg)
        if self["ZENANGLE"] is None:
            self["ZENANGLE"] = 0.  # (deg)

        # Set latitude/longitude in between Lannion and Exeter,
        # in the "Manche"/"Channel"
        # Lannion is 48.750,  -3.470
        # Exeter is  50.726,  -3.476
        if self["LATITUDE"] is None:
            self["LATITUDE"] = 49.738  # (deg)
        if self["LONGITUDE"] is None:
            # (deg)    Lat/lon compatible with suftype==ocean
            self["LONGITUDE"] = -3.473

        " Magnetic field "
        if self["BE"] is None:
            self["BE"] = 0.3  # (Gauss)
        if self["COSBK"] is None:
            self["COSBK"] = 1.

        " Mislaneous "
        if self["ID"] is None:
            self["ID"] = "This is my profile"
        if self["DATE"] is None:
            self["DATE"] = numpy.array([2014, 0o4, 30])  # Year, Month, Day
        if self["TIME"] is None:
            self["TIME"] = numpy.array([12, 0, 0])     # Hour, Minute, Second
        self.setDefaultAttributes()
        if self["GAS_UNITS"] is None:
            self["GAS_UNITS"] = 2  # defaut case
        if self["MMR_CLDAER"] is None:
            self["MMR_CLDAER"] = 1  # defaut case

    def checkProfileAsciiInput(self):
        if not self.checkAerosolKind():
            raise IOError
        if (self["P"] is None or self["T"] is None or self["Q"] is None):
            print("At least one of P, T, Q arrays is missing")
            raise IOError
        self["NLEVELS"] = numpy.shape(self["P"])[0]
        self["NLAYERS"] = self["NLEVELS"] - 1
        for cle in self.gas_list:
            if(self[cle] is not None):
                if (numpy.shape(self[cle])[0] != self["NLEVELS"]):
                    print("Bad array size for", cle)
                    raise IOError
        for cle in self.aerosol_list:
            if(self[cle] is not None):
                if (numpy.shape(self[cle])[0] != self["NLAYERS"]):
                    print("Bad array size for", cle)
                    raise IOError
        for cle in self.cloud_list:
            if(self[cle] is not None):
                if (numpy.shape(self[cle])[0] != self["NLAYERS"]):
                    print("Bad array size for", cle)
                    raise IOError
        if self.anyCloud():
            if(self["CFRAC"] is not None):
                if (numpy.shape(self["CFRAC"])[0] != self["NLAYERS"]):
                    print("Bad array size for", cle)
                    raise IOError
            else:
                print("Cloud fraction missing ")
                raise IOError
            if(self["ICEDE"] is None):
                # ICEDE array is mandatory for cloudy profiles
                self["ICEDE"] = numpy.zeros(
                    self["NLAYERS"], dtype=self["P"].dtype)
        if "GAS_UNITS" not in self:
            self["GAS_UNITS"] = 2  # defaut case
        elif self["GAS_UNITS"] is None:
            self["GAS_UNITS"] = 2  # defaut case
        if "MMR_CLDAER" not in self:
            self["MMR_CLDAER"] = 1  # defaut case
        elif self["MMR_CLDAER"] is None:
            self["MMR_CLDAER"] = 1  # defaut case
        self.defineUnitsforCommentAttribute()

    def loadh5(self, h5):
        # Load driven by the file content

        for cle in list(h5):
            if(cle in ['S2M', 'SKIN']):

                for scle in list(h5[cle]):
                    self[cle][scle + '_ATTRIBUTE'] = {}
                    self[cle][scle] = loadDset_arr(h5[cle + '/' + scle])
                    self[cle][
                        scle + '_ATTRIBUTE'] = getAttributes(h5[
                            cle + '/' + scle])
                    checkAttributes(self[cle][scle + '_ATTRIBUTE'])
            elif(cle == 'ID'):

                self[cle] = loadDset_str(h5[cle])
                self[cle + '_ATTRIBUTE'] = getAttributes(h5[cle])
                checkAttributes(self[cle + '_ATTRIBUTE'])
            else:

                self[cle] = loadDset_arr(h5[cle])
                self[cle + '_ATTRIBUTE'] = getAttributes(h5[cle])
                checkAttributes(self[cle + '_ATTRIBUTE'])

        self.defineUnitsforCommentAttribute()

        if(self['AEROSOLS'] is not None):
            if self['AEROSOLS_ATTRIBUTE']['COMMENT'] == "CAMS aerosols":
                self.aerosol_kind = "cams"
                self.aerosol_list = self.cams_aerosol_list
                self.aerosol_long_list = self.cams_aerosol_long_list
            else:
                self.aerosol_kind = "opac"
                self.aerosol_list = self.opac_aerosol_list
                self.aerosol_long_list = self.opac_aerosol_long_list
            print("areosol kind:", self.aerosol_kind)
            print(self.aerosol_list)
            na = self['AEROSOLS'].shape[1]
            if(self['AEROSOLS'].ndim == 2):
                for i in range(na):
                    if (any(self['AEROSOLS'][:, i] != 0.)):
                        self[self.aerosol_list[i]] = self['AEROSOLS'][:, i]
                        self[self.aerosol_list[i] + '_ATTRIBUTE'] = {
                            'COMMENT': self.aerosol_long_list[i],
                            'UNITS': self.my_aer_units}
                        if ('UNITS' in list(self[
                                'AEROSOLS_ATTRIBUTE'].keys())):
                            self[self.aerosol_list[
                                i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.aerosol_long_list[i]}
                            self[self.aerosol_list[
                                i] + '_ATTRIBUTE']['UNITS'] = self[
                                'AEROSOLS_ATTRIBUTE']['UNITS']
                        else:
                            self[self.aerosol_list[i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.aerosol_long_list[i],
                                'UNITS': self.my_aer_units}
            else:
                for i in range(na):
                    if (numpy.any(self['AEROSOLS'][:, i, :] != 0.)):
                        self[self.aerosol_list[i]] = self['AEROSOLS'][:, i, :]
                        if ('UNITS' in list(self[
                                'AEROSOLS_ATTRIBUTE'].keys())):
                            self[self.aerosol_list[
                                i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.aerosol_long_list[i]}
                            self[self.aerosol_list[
                                i] + '_ATTRIBUTE']['UNITS'] = self[
                                'AEROSOLS_ATTRIBUTE']['UNITS']
                        else:
                            self[self.aerosol_list[i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.aerosol_long_list[i],
                                'UNITS': self.my_aer_units}
            self['AEROSOLS'] = None
            self['AEROSOLS_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

        if(self['CLOUD'] is not None):
            na = self['CLOUD'].shape[1]
            if(self['CLOUD'].ndim == 2):
                for i in range(na):
                    if (any(self['CLOUD'][:, i] != 0.)):
                        self[self.cloud_list[i]] = self['CLOUD'][:, i]
                        if ('UNITS' in list(self['CLOUD_ATTRIBUTE'].keys())):
                            self[self.cloud_list[
                                i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.cloud_long_list[i]}
                            self[self.cloud_list[
                                i] + '_ATTRIBUTE']['UNITS'] = self[
                                'CLOUD_ATTRIBUTE']['UNITS']
                        else:
                            self[self.cloud_list[i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.cloud_long_list[i],
                                'UNITS': 'layer mean content(g/m3)'}
            else:
                for i in range(na):
                    if (numpy.any(self['CLOUD'][:, i, :] != 0.)):
                        self[self.cloud_list[i]] = self['CLOUD'][:, i, :]
                        if ('UNITS' in list(self['CLOUD_ATTRIBUTE'].keys())):
                            self[self.cloud_list[
                                i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.cloud_long_list[i]}
                            self[self.cloud_list[
                                i] + '_ATTRIBUTE']['UNITS'] = self[
                                'CLOUD_ATTRIBUTE']['UNITS']
                        else:
                            self[self.cloud_list[i] + '_ATTRIBUTE'] = {
                                'COMMENT': self.cloud_long_list[i],
                                'UNITS': 'layer mean content(g/m3)'}
            self['CLOUD'] = None
            self['CLOUD_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
            if(self['ICEDE'] is None):
                # ICEDE array is mandatory for cloudy profiles
                self['ICEDE'] = numpy.zeros(
                    self['NLAYERS'], dtype=self['P'].dtype)
                self['ICEDE_ATTRIBUTE'] = {'COMMENT': 'Ice crystals diameter',
                                           'UNITS': 'microns', 'LBOUND': 1,
                                           'UBOUND': self['NLAYERS']}
            if(self['CLWDE'] is None):
                # CLWDE array is mandatory for cloudy profiles
                self['CLWDE'] = numpy.zeros(
                    self['NLAYERS'], dtype=self['P'].dtype)
                self['CLWDE_ATTRIBUTE'] = {'COMMENT': 'Liquid water diameter',
                                           'UNITS': 'microns', 'LBOUND': 1,
                                           'UBOUND': self['NLAYERS']}
            if self['CLWDE_PARAM'] is None:
                self['CLWDE_PARAM'] = 1
                self['CLWDE_PARAM_ATTRIBUTE'] = {
                    'COMMENT':
                    'Liquid water effective diameter parametrisation',
                    'UNITS': 'n/a', 'LBOUND': 1,
                    'UBOUND': 1}

    def controlProfile(self):
        """ control some values of the profile regarding to min et max
            values """
        for label in list(self.minValue.keys()):
            if label[:3] == "S2M":
                if self["S2M"][label[4:]] < self.minValue[label]:
                    self["S2M"][label[4:]] = self.minValue[label]
                if self["S2M"][label[4:]] > self.maxValue[label]:
                    self["S2M"][label[4:]] = self.maxValue[label]
            else:
                if label[:4] == "SKIN":
                    if self["SKIN"][label[5:]] < self.minValue[label]:
                        self["SKIN"][label[5:]] = self.minValue[label]
                    if self["SKIN"][label[5:]] > self.maxValue[label]:
                        self["SKIN"][label[5:]] = self.maxValue[label]
                else:
                    if self[label] < self.minValue[label]:
                        self[label] = self.minValue[label]
                    if self[label] > self.maxValue[label]:
                        self[label] = self.maxValue[label]

    def loadProfileNumber(self, filename, iprof):
        """
        Load profile number iprof from file
        Opens and closes the file
        """
        # opens the HDF5 file read only
        f = h5py.File(filename, 'r')
        dir = '/PROFILES/{:0>4}/'.format(iprof)

        # get the Dataset
        h5 = f[dir]

        # Load profile
        self.loadh5(h5)
        self.defineUnitsforCommentAttribute()

        # Close HDF file
        f.close()

    def loadProfileAscii(self, fname):
        """
        Load an ASCII profile from file fname
        """
        try:
            exec(compile(open(fname).read(), fname, 'exec'))
        except SyntaxError:
            print("Error Reading ASCII profile file")
            return 1
        try:
            self.checkProfileAsciiInput()
        except IOError:
            print("Error Checking ASCII profile file")
            return 1
        self.setDefaultProfileAsciiInput()

        return 0

    def display(self):
        # Display profile content to terminal
        for cle in sorted(self.keys()):
            print(cle, '=>', self[cle], type(self[cle]))

    def displaySurfaceParameters(self):
        # Display profile content to terminal
        surface_list = self.profSurfParamList

        for cle in surface_list:
            if cle not in self.notUsedParamsShort:
                print(cle, '=>', self[cle], type(self[cle]))
        for cle in self.s2m_list:
            if "S2M_" + cle not in self.notUsedParamsShort:
                print(cle, '=>', self["S2M"][cle], type(self["S2M"][cle]))
        for cle in self.sskin_list:
            if "SKIN_" + cle not in self.notUsedParamsShort:
                print(cle, '=>', self["SKIN"][cle], type(self["SKIN"][cle]))

    def saveh5(self, h5, path):
        """Save profile class to an HDF5 file h5 under path directory
        Do not save empty arrays
        Attributes are saved at the same time datasets are saved
        """
        list_of_names = []
        h5.visit(list_of_names.append)

        if self.anyAerosol():
            self['AEROSOLS'] = numpy.zeros(
                (self['NLAYERS'], len(self.aerosol_list)),
                dtype=self['P'].dtype)
            if self.aerosol_kind == "cams":
                self['AEROSOLS_ATTRIBUTE'] = {
                    'COMMENT': 'CAMS aerosols', 'UNITS': self.my_aer_units}
            else:
                self['AEROSOLS_ATTRIBUTE'] = {
                    'COMMENT': 'OPAC aerosols', 'UNITS': self.my_aer_units}
            for cle in self.aerosol_list:
                if(self[cle] is None):
                    continue
                else:
                    i = self.aerosol_list.index(cle)
                    self['AEROSOLS'][:, i] = self[cle]

        if self.anyCloud():
            self['CLOUD'] = numpy.zeros(
                (self['NLAYERS'], len(self.cloud_list)), dtype=self['P'].dtype)
            self['CLOUD_ATTRIBUTE'] = {
                'COMMENT': 'Cloud water/ice - IR only',
                'UNITS': self.my_cld_units}

            for cle in self.cloud_list:
                if(self[cle] is None):
                    continue
                else:
                    i = self.cloud_list.index(cle)
                    self['CLOUD'][:, i] = self[cle]
            if(self['ICEDE'] is None):
                # ICEDE array is mandatory for cloudy profiles
                self['ICEDE'] = numpy.zeros(
                    self['NLAYERS'], dtype=self['P'].dtype)
                self['ICEDE_ATTRIBUTE'] = {'COMMENT': 'Ice crystals diameter',
                                           'UNITS': 'microns', 'LBOUND': 1,
                                           'UBOUND': self['NLAYERS']}
            if(self['CLWDE'] is None):
                # CLWDE array is mandatory for cloudy profiles
                self['CLWDE'] = numpy.zeros(
                    self['NLAYERS'], dtype=self['P'].dtype)
                self['CLWDE_ATTRIBUTE'] = {'COMMENT': 'Liquid water diameter',
                                           'UNITS': 'microns', 'LBOUND': 1,
                                           'UBOUND': self['NLAYERS']}

        for cle in list(self.keys()):

            if(re.search('_ATTRIBUTE', cle) or self[cle] is None):
                # do not save empty members
                continue

            if(cle in ['S2M', 'SKIN']):
                # S2M and SKIN are HDF5 sub-directories
                for scle in list(self[cle].keys()):
                    if(
                            re.search('_ATTRIBUTE', scle) or
                            self[cle][scle] is None):
                        continue
                    # create the dataset and save attributes
                    cpath = cleanpath(path + '/' + cle + '/' + scle)
                    if cpath in list_of_names:
                        del h5[cpath]
                    dset = h5.create_dataset(cpath, data=self[cle][scle])
                    saveAttributes(h5[cpath], self[cle][scle + '_ATTRIBUTE'])

            elif(cle in self.aerosol_list or cle in self.cloud_list):
                # do not save each individual aerosol array but "aerosol" 2
                # dimensions array
                # HDF5 file only contains the RTTOV arrays
                continue

            else:
                # create the dataset and save attributes
                cpath = cleanpath(path + '/' + cle)
                if cpath in list_of_names:
                    del h5[cpath]
                h5.create_dataset(cpath, data=self[cle])
                saveAttributes(h5[cpath], self[cle + '_ATTRIBUTE'])

        # ici on peut supprimer le tableau a deux dimensions des aerosols et
        # nuages
        self['AEROSOLS'] = None
        self['AEROSOLS_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        self['CLOUD'] = None
        self['CLOUD_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def addGas(self, gas, *conc):
        """
        on ajoute un gas de nom gas a la bonne dimension nlevels
        et avec le meme type de données que le tableau de pression
        """
        # create numpy array filled with zeros, size is based on number of
        # pressure levels
        self[gas] = numpy.zeros(self['NLEVELS'], dtype=self['P'].dtype)
        self[gas + '_ATTRIBUTE'] = {'COMMENT': gas, 'UNITS': self.my_gas_units,
                                    'LBOUND': 1, 'UBOUND': self['NLEVELS']}
        # Fills with given concentration or default value
        if(conc):
            self[gas] += conc
        else:
            gas_conc = get_default_gas(gas, self["P"])
            if self['GAS_UNITS'] in [-1, 0, 2]:
                self[gas] += gas_conc
            else:
                self[gas] += fct_ppmv_to_mixratio(gas_conc, gas)

    def removeGas(self, gas):
        """
        remove gas nammed gas
        """
        self[gas] = None
        self[gas + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def removeAllButGas(self, gas):
        """
        remove all gases except gas
        """
        for g in self.gas_list:
            if(not g == gas):
                self[g] = None
                self[g + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def removeAllAerosol(self):
        """
        remove all  aerosols
        """
        for i in self.aerosol_list:
            self[i] = None
            self[i + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

        self['AEROSOLS'] = None
        self['AEROSOLS_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def removeAerosol(self, aerosol):
        """
        remove one aerosol
        """
        self[aerosol] = None
        self[aerosol + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}

    def addAerosol(self, aerosol, *conc):
        """
        add an aerosol of concentration conc
        """
        # create numpy array filled with zeros, size is based on number of
        # pressure levels, fill it with default value

        i = self.aerosol_list.index(aerosol)
        self[aerosol] = numpy.zeros(self['NLAYERS'], dtype=self['P'].dtype)
        self[aerosol + '_ATTRIBUTE'] = {'COMMENT': self.aerosol_long_list[i],
                                        'UNITS': self.my_aer_units,
                                        'LBOUND': 1, 'UBOUND': self['NLAYERS']}
        if(aerosol == 'VOLA' or aerosol == 'VAPO'):
            # volcanic ashes
            pmin = 100.
            pmax = 300.
        else:
            pmin = 600.
            pmax = 1000.
        for i in range(self['NLAYERS']):
            # here we use presure levels indices to fill cloud LAYERS.
            if(self['P'][i] > pmin and self['P'][i] < pmax):
                if(conc):
                    self[aerosol][i] += conc
                else:
                    self[aerosol][i] += self.defaultaerdensity[
                        self["MMR_CLDAER"]]

    def replaceByAerosolClim(self, clim):
        """
        use profile from climatology for aerosols
        with climatology :
        1  -->Continental clean
        2  -->Continental average
        3  -->Continental polluted
        4  -->Urban
        5  -->Desert
        6  -->Maritime clean
        7  -->Maritime polluted
        8  -->Maritime tropical
        9  -->Arctic
        10 -->Antarctic
        """
        # create numpy array filled with zeros, size is based on number of
        # pressure levels
        tabaer = None
        tabaer, err = rttov_gui_aer_clim_prof(self['P'], self['T'], self['Q'],
                                              self["GAS_UNITS"], self[
                                                  "MMR_CLDAER"],
                                              self['NLEVELS'],
                                              self['LATITUDE'],
                                              self['ELEVATION'],
                                              1.0,
                                              self['NLEVELS'])

        self.removeAllAerosol()
        self['AEROSOLS'] = tabaer[:, clim - 1, :]
        na = self['AEROSOLS'].shape[1]
        if self.aerosol_kind != "opac":
            return 1
        for i in range(na):
            if (any(self['AEROSOLS'][:, i] != 0.)):
                self[self.aerosol_list[i]] = self['AEROSOLS'][:, i]
                self[self.aerosol_list[i] + '_ATTRIBUTE'] = {
                    'COMMENT': self.aerosol_long_list[i],
                    'UNITS': self.my_aer_units}
        # ici on peut supprimer le tableau a deux dimensions des aerosols
        self['AEROSOLS'] = None
        self['AEROSOLS_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        return 0

    def anyAerosol(self):
        for cle in self.aerosol_list:
            if(self[cle] is None):
                continue
            else:
                return True
        return False

    def removeCloud(self, cloud):
        """
        remove a cloud
        """
        self[cloud] = None
        self[cloud + '_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a'}
        if not self.anyCloud():
            self['CFRAC'] = None
            self['CFRAC_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a',
                                       'LBOUND': 0, 'UBOUND': 0}
        if not self.anyCloud() and self['CLW'] is not None:
            self['CLW'] = None
            self['CLW_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a',
                                     'LBOUND': 0, 'UBOUND': 0}
        if not self.anyCloud() and self['ICEDE'] is not None:
            self['ICEDE'] = None
            self['ICEDE_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a',
                                       'LBOUND': 0, 'UBOUND': 0}
        if not self.anyCloud() and self['CLWDE'] is not None:
            self['CLWDE'] = None
            self['CLWDE_ATTRIBUTE'] = {'COMMENT': 'none', 'UNITS': 'n/a',
                                       'LBOUND': 0, 'UBOUND': 0}

    def addCloud(self, cloud, *conc):
        """
        add a cloud of concentration conc
        """
        # create numpy array for cloud filled with zeros,
        # size is based on number of
        # pressure levels

        indexItem = self.cloud_list.index(cloud)
        print("Profile addCloud", cloud)
        self[cloud] = numpy.zeros(self['NLAYERS'], dtype=self['P'].dtype)
        self[cloud + '_ATTRIBUTE'] = {
            'COMMENT': self.cloud_long_list[indexItem],
            'UNITS': self.my_cld_units,
            'LBOUND': 1, 'UBOUND': self['NLAYERS']}

        if(self['CFRAC'] is None):
            self['CFRAC'] = numpy.zeros(self['NLAYERS'], dtype=self['P'].dtype)
            self['CFRAC_ATTRIBUTE'] = {
                'COMMENT': 'Cloud fractional cover (0:1)',
                'UNITS': '0:1',
                'LBOUND': 1, 'UBOUND': self['NLAYERS']}

        if(cloud == 'STCO' or cloud == 'STMA'):
            # Stratus
            pmin = 850.
            pmax = 1000.
        elif(cloud == 'CUCC' or cloud == 'CUCP' or cloud == 'CUMA'):
            # Cumulus
            pmin = 600.
            pmax = 800.
        else:
            # Cirrus
            pmin = 450.
            pmax = 550.

        for i in range(self['NLAYERS']):
            # here we use pressure levels indices to fill cloud LAYERS.
            if(self['P'][i] > pmin and self['P'][i] < pmax):
                self[cloud][i] = self.defaultclddensity[self["MMR_CLDAER"]]
                if (self['CFRAC'][i] == 0.):
                    self['CFRAC'][i] = 0.20

    def anyCloud(self):
        for cle in self.cloud_list:
            if(self[cle] is None):
                continue
            else:
                return True
        return False

    def hasMfasis(self):
        if self.anyCloud():
            return True

    def to1DvarVect(self):
        """ return an 1Dvar vector """
        """ 1Dvar vector are only on 43, 51 or 54 levels"""
        """ a vector contains only T values, lnq bottom lnq vales ,
            Tsurf, lnq surf and Tskin"""
        """ temperature, and ln q are stored from top of atmosphere
            to ground """

        nlevels = len(self['T'])
        if (nlevels) not in (43, 51, 54):
            print(("WARNING :", "nlevels " +
                   str(nlevels) + "not in (43,51,54) "))
        veclen = nlevels + 29 + 3
        # 1-54=temp (K), 55-83=lnq (g/kg) (bottom 29 levels only), 84=Tsurf,
        # 85=lnq surf, 86=Tskin
        vector = numpy.zeros((veclen), dtype=float)
        for i in range(nlevels):
            vector[i] = self['T'][i]
        for i in range(29):
            vector[nlevels + i] = self['Q'][nlevels - 29 + i]
        vector[nlevels + 29] = self['S2M']['T']
        vector[nlevels + 29 + 1] = self['S2M']['Q']
        vector[nlevels + 29 + 2] = self['SKIN']['T']

        return vector

    def to1DvarBackgroundFile(self, filename):
        """ write an 1Dvar Background File """
        """ 1Dvar work only on 43, 51 or 54 levels"""
        """ must respect File description """
        """ https://nwpsaf.eu/deliverables/nwpsaf_1dvar/
            nwpsaf-mo-ud-032_NWPSAF_1DVar_Manual.html#aux """
        f = open(filename, "w")
        nlevels = len(self['T'])
        if (nlevels) not in (43, 51, 54):
            print(("WARNING :", "nlevels " +
                   str(nlevels) + "not in (43,51,54) \n"))
        f.write("Backgrounnd file generated from\n")
        f.write(self["ID"] + "\n")
        f.write(str(nlevels) + " fixed layers")
        for k in range(0, 6):
            f.write("\n")
        f.write(
            "x------------------- End of Header---------------------------x\n")
        f.write("\n")
        f.write("No. Background Profiles:            1\n")
        f.write("No. of Levels/Profile:             " + str(nlevels) + "\n")
        f.write("Unit for q:                     1  (1=ppmv, 2=kg/kg, 3=RH)\n")
        f.write(
            "x------------------------------------------------------------x\n")
        f.write("Profile #    1 Follows\n")
        f.write(
            "x------------------------------------------------------------x\n")
        for i in range(nlevels - 1, -1, -1):
            sp = ("%10.5f" % self["P"][i]).rjust(13)
            st = ("%9.5f" % self["T"][i]).rjust(13)
            sq = ("%f" % self["Q"][i]).rjust(13)
            so3 = ("%f" % self["O3"][i]).rjust(13)
            f.write(sp + st + sq + so3)
            f.write("\n")
        f.write("Surface Temperature (K):" +
                ("%9.5f" % self["S2M"]["T"]).rjust(18))
        f.write("\n")
        f.write("Surface Humidity :" +
                ("%f" % self["S2M"]["Q"]).rjust(18))
        f.write("\n")
        f.write("Skin Temperature (K):" + ("%9.5f" %
                                           self["SKIN"]["T"]).rjust(18))
        f.write("\n")
        f.write("Surface Pressure (hPa):" + ("%10.5f" %
                                             self["S2M"]["P"]).rjust(18))
        f.write("\n")
        f.write("10m U-Wind (m/s):" + ("%f" % self["S2M"]["U"]).rjust(18))
        f.write("\n")
        f.write("10m U-Wind (m/s):" + ("%f" % self["S2M"]["V"]).rjust(18))
        f.write("\n")
        f.close()


if __name__ == '__main__':
    import os
    rttov_profile_dir = os.environ["RTTOV_GUI_PROFILE_DIR"]
    fname = os.path.join(rttov_profile_dir, "cldaer101lev_allgas.H5")
    if not os.path.isfile(fname):
        print("Error", fname, "not found")

    # Open file ReadOnly
    f = h5py.File(fname, 'r')

    # get the Dataset
    h5 = f['/PROFILES/0001/']

    # p1 is a rttov_hdf_mod Profile instance
    p1 = Profile()
    print("loadh5")
    print(p1.aerosol_kind)
    # Load profile
    p1.loadh5(h5)

    # Close HDF file
    f.close()

    # Display profile
    p1.display()
    p1.displaySurfaceParameters()

    ofile = "rttov_profile_diroutput_profile.H5"
    of = h5py.File(ofile, 'w')
    p1.saveh5(of, "/PROFILES/0001/")
    of.close()

    # Plots temperature and gas concentrations
    import matplotlib
    matplotlib.use("wxagg")
    import matplotlib.pyplot as plt
    plt.plot(p1['T'], p1['P'], 'ro-')
    plt.xlabel(p1['T_ATTRIBUTE']['COMMENT'] +
               '   (' + p1['T_ATTRIBUTE']['UNITS'] + ')')
    plt.ylim(1100, 0.01)
    plt.yscale('log')
    plt.ylabel(p1['P_ATTRIBUTE']['COMMENT'] +
               '  ' + p1['P_ATTRIBUTE']['UNITS'])
    plt.title(p1['ID'] + ' TEMPERATURE')
    plt.show()

    if(p1.anyCloud()):
        nlevels = p1['NLEVELS']
        player = (p1['P'][0:nlevels - 1] + p1['P'][1:nlevels]) / 2.0
        colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
        i = 0
        for cle in p1.cloud_list:
            if(p1[cle] is None):
                continue
            else:
                color = colors[i]
                x = p1[cle]
                plt.plot(x, player, color + 'o-',
                         label=p1[cle + '_ATTRIBUTE']['COMMENT'])
                i += 1

        x = p1['CFRAC'][:]
        plt.plot(x, player, 'ko-', label='Cloud Fraction')
        plt.xlabel('Cloud Liquid water/ Cloud Ice Water')
        plt.xscale('log')
        plt.ylim(1100, 10)
        plt.yscale('log')
        plt.ylabel(p1['P_ATTRIBUTE']['COMMENT'] +
                   '   (' + p1['P_ATTRIBUTE']['UNITS'] + ')')
        plt.legend()
        plt.title(p1['ID'] + ' CLOUDS')
        plt.show()

    if (p1.anyAerosol()):

        nlevels = p1['NLEVELS']
        player = (p1['P'][0:nlevels - 1] + p1['P'][1:nlevels]) / 2.0
        colors = ['b', 'g', 'r', 'c', 'm', 'y']
        i = 0
        for cle in p1.aerosol_list:
            if(p1[cle] is None):
                continue
            else:
                color = colors[i % len(colors)]
                if (i >= len(colors)):
                    pattern = 'o'
                else:
                    pattern = 'v'
                x = p1[cle]
                plt.plot(x, player, color + pattern + '-',
                         label=p1[cle + '_ATTRIBUTE']['COMMENT'])
                i += 1
        plt.xlabel('Aerosols (%s)' % p1.my_aer_units)
        plt.xscale('log')
        plt.ylim(1100, 5)
        plt.yscale('log')
        plt.ylabel(p1['P_ATTRIBUTE']['COMMENT'] +
                   '   (' + p1['P_ATTRIBUTE']['UNITS'] + ')')
        plt.legend()
        plt.title(p1['ID'] + ' AEROSOLS')
        plt.show()

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    i = 0
    for gas in Profile.gas_list:
        if(p1[gas] is None):
            continue
        color = colors[i]
        i += 1
        plt.plot(p1[gas], p1['P'], color + 'o-', label=gas)
    plt.xlabel('gas concentrations (%s)' % p1.my_gas_units)
    plt.xscale('log')
    plt.ylim(1100, 0.01)
    plt.yscale('log')
    plt.ylabel(p1['P_ATTRIBUTE']['COMMENT'] +
               '   (' + p1['P_ATTRIBUTE']['UNITS'] + ')')
    plt.legend()
    plt.title(p1['ID'] + ' GASES')
    plt.show()
    print("test OK")
