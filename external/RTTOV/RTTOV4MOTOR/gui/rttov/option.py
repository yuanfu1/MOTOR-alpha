#! /usr/bin/env python
# -*- coding: utf-8 -*-

import h5py
import re

from rttov.core import _V, loadDset_arr, loadDset_log, getAttributes
from rttov.core import saveAttributes, checkAttributes
from collections import OrderedDict
import copy


class OptionObject(object):
    """
    OptionObject class
    """

    def __init__(self, otype=bool, value=False, min_max_values=None,
                 odict=None, comment="", status=False, theme=None,
                 units="N/A", hidden=False):
        self.otype = otype
        self.theme = theme
        self.value = value
        self.min_max_values = min_max_values
        self.odict = odict  # possibles values (dictionary)
        self.comment = comment
        self.status = status  # relevant or not
        self.units = units
        self.hidden = hidden

    def __setitem__(self, value):
        self.value = value

    def __getitem__(self):
        return self.value

    def getValue(self):
        if self.otype == dict:
            theKeys = self.odict.keys()
            if self.value in theKeys:
                return self.odict[self.value]
            else:
                return self.odict[theKeys[0]]
        else:
            return self.value

    def __deepcopy__(self):

        otype = copy.deepcopy(self.otype)
        theme = copy.deepcopy(self.theme)
        value = copy.deepcopy(self.value)

        min_max_values = copy.deepcopy(self.min_max_values)
        odict = copy.deepcopy(self.odict)  # possibles values (dictionary)
        comment = copy.deepcopy(self.comment)
        status = copy.deepcopy(self.status)  # relevant or not
        units = copy.deepcopy(self.units)
        hidden = copy.deepcopy(self.hidden)
        newOptionObject = OptionObject(otype, value, min_max_values,
                                       odict, comment, status, theme,
                                       units, hidden)
        return newOptionObject


class Option(OrderedDict, _V):
    """
    The Option class
    A Python dictionary which contains a copy of the RTTOV v12
    rttov_options structure
    The creation method fills the class with the same defaults
    as the Fortran structure
    Inherits the _V class for some default methods
    """

    optionsThemes = ['General configuration options',  # 0
                     'VIS/IR-only RT options',# 1
                     'General RT options', # 2
                     'PC-RTTOV options', # 3
                     'MW-only RT options', # 4
                     'Interpolation and vertical grid options'] # 5

    fullIPCREG = {1: "1", 2: '2', 3: '3', 4: '4'}
    fullIPCBND = {1: "full spectrum",
                  2: 'longwave subset',
                  3: 'short wave subset'}
    dictOptions = {"IPCREG": fullIPCREG,
                   'IPCBND': fullIPCBND,
                   'FASTEM_VERSION': {0: "TESSEM2",
                                      1: "1", 2: '2', 3: '3', 4: '4',
                                      5: '5', 6: '6'},
                   "CLW_SCHEME":  {1: "Liebe", 2: "Rosenkranz", 3: "TKC"},
                   'INTERP_MODE': {1: "1 Rochon/Rochon OD",
                                   2: '2 Log-Lin/Log-lin OD',
                                   3: '3 Rochon/Log-lin OD',
                                   4: '4 Rochon/Rochon WF',
                                   5: '5 Rochon/Log-lin WF'},
                   "VIS_SCATT_MODEL": {1: "DOM",
                                       2: "single scattering"},
                   "ALL_VIS_SCATT_MODEL": {1: "DOM",
                                           2: "single scattering",
                                           3: "MFASIS"},
                   "IR_SCATT_MODEL": {1: "DOM",
                                      2: "chou-scaling"},
                   "IR_SEA_EMIS_MODEL":  {1: "ISEM", 2: "IRemis"},
                   "SOLAR_SEA_BRDF_MODEL":  {1: "JONSWAP", 2: "Elfouhaily"}
                   }
    scatt_options = ["VIS_SCATT_MODEL", "IR_SCATT_MODEL",
                     "DOM_NSTREAMS",
                     "DOM_ACCURACY", "DOM_OPDEP_THRESHOLD",
                     "GRID_BOX_AVG_CLOUD", "DOM_RAYLEIGH"]

    ipcregNbPredictors = {1: {"IASI": {1: 300, 2: 400, 3: 500, 4: 600},
                              "AIRS": {1: 200, 2: 300, 3: 400}},
                          2: {"IASI": {1: 35, 2: 40, 3: 45, 4: 50}},
                          3: {"IASI": {1: 30, 2: 40, 3: 50, 4: 60}}
                          }

    def __deepcopy__(self):

        newoption = Option()
        for key in list(self.keys()):
            if type(self[key]) == OptionObject:
                newoption[key] = self[key].__deepcopy__()
            else:
                newoption[key] = copy.deepcopy(self[key])
        return newoption

    def __init__(self):
        OrderedDict.__init__(self)
        self["APPLY_REG_LIMITS"] = OptionObject(
            comment="Switch to restrict input"
            " profiles to coef training limits",
            theme=0)
        self["VERBOSE"] = OptionObject(
            comment="Switch for verbose output",
            theme=0)
        self["DO_CHECKINPUT"] = OptionObject(
            comment="Switch to apply internal profile "
            "checking",
            theme=0)
        self["FIX_HGPL"] = OptionObject(
            comment="Activate fix to match 2m p with "
            "elevation in geometry calculations",
            value=True,
            theme=0)
        self["ADDREFRAC"] = OptionObject(
            comment="Switch to enable atmospheric "
            "refraction",
            value=True,
            theme=2)
        self["SWITCHRAD"] = OptionObject(
            comment="Switch for input units in AD/K models",
            value=True,
            theme=2)
        self["USE_T2M_OPDEP"] = OptionObject(
            comment="Switch to enable use of 2m T "
            "variable in optical depth calculation",
            theme=2, value=True)
        self["USE_Q2M"] = OptionObject(
            comment="Switch to enable use of 2m q variable",
            theme=2, value=True)
        self["DO_LAMBERTIAN"] = OptionObject(
            comment="Switch for setting Lambertian "
            "reflection (IR and MW)",
            theme=2)
        self["LAMBERTIAN_FIXED_ANGLE"] = OptionObject(
            comment="Switch for Lambertian fixed or parameterised "
            "angle (IR and MW)",
            theme=2, value=True)
        self["PLANE_PARALLEL"] = OptionObject(
            comment="Switch to ignore atmospheric curvature",
            theme=2)
        self["RAD_DOWN_LIN_TAU"] = OptionObject(
            comment="Switch to select down-welling radiance "
            "source term",
            theme=2, value=True,
            hidden=False)
        self["DTAU_TEST"] = OptionObject(
            comment="Switch to apply dtau test in transmit/"
            "integrate calculations",
            theme=2,
            value=False,
            hidden=True)
        self["SOLAR_SEA_BRDF_MODEL"] = OptionObject(
            otype=dict,
            value=2,
            odict=self.dictOptions[
                "SOLAR_SEA_BRDF_MODEL"],
            comment="Solar sea BRDF model",
            theme=1)
        self["IR_SEA_EMIS_MODEL"] = OptionObject(
            otype=dict,
            value=2,
            odict=self.dictOptions[
                "IR_SEA_EMIS_MODEL"],
            comment="IR sea emissivity model",
            theme=1)
        self["ADDSOLAR"] = OptionObject(
            comment="Switch to enable solar simulations",
            theme=1)
        self["RAYLEIGH_MAX_WAVELENGTH"] = OptionObject(
            value=2.,
            otype=float,
            min_max_values=(0.00, 3.00),
            comment="Ignore Rayleigh scattering for channels at wavelengths"
            " above this (microns)",
            theme=1)
        self["RAYLEIGH_MIN_PRESSURE"] = OptionObject(
            value=0,
            otype=float,
            min_max_values=(0, 1100.),
            comment="Ignore Rayleigh scattering at pressures below this (hPa)",
            theme=1)
        self["RAYLEIGH_SINGLE_SCATT"] = OptionObject(
            value=True,
            comment=" Switch to enable Rayleigh single-scattering "
            "for VIS/NIR channels",
            theme=1)
        self["DO_NLTE_CORRECTION"] = OptionObject(
            comment="Switch to enable NLTE "
            "bias correction",
            theme=1)
        self["ADDAEROSL"] = OptionObject(
            comment="Switch to enable IR "
            "aerosol calculation",
            theme=1)
        self["ADDCLOUDS"] = OptionObject(
            comment="Switch to enable IR "
            "cloudy calculations",
            theme=1)
        self["USER_AER_OPT_PARAM"] = OptionObject(
            comment="Switch to supply aerosol "
            "optical properties explicitly "
            "per channel",
            theme=1,
            hidden=True)
        self["USER_CLD_OPT_PARAM"] = OptionObject(
            comment="Switch to supply cloud optical "
            " properties explicitly per channel",
            theme=1,
            hidden=True)
        self["CLDCOL_THRESHOLD"] = OptionObject(
            value=-1,
            otype=float,
            min_max_values=(-1, 1),
            comment="Ignore cloud columns with "
            "weights lower than this",
            theme=1)
        self["IR_SCATT_MODEL"] = OptionObject(
            value=2,
            otype=dict,
            odict=self.dictOptions[
                "IR_SCATT_MODEL"],
            comment="IR scattering model"
            "1=DOM, 2=Chou-scaling,  default 2",
            theme=1)
        self["VIS_SCATT_MODEL"] = OptionObject(
            value=1,
            otype=dict,
            odict=self.dictOptions[
                "ALL_VIS_SCATT_MODEL"],
            comment="VIS/NIR scattering model"
            "1=DOM, 2=single-scattering",
            theme=1)
        self["DOM_NSTREAMS"] = OptionObject(
            value=8,
            otype=int,
            min_max_values=(2, 128),
            comment="Number of DOM streams, must be "
            " even and not less than 2",
            theme=1)
        self["DOM_ACCURACY"] = OptionObject(
            value=0,
            otype=float,
            min_max_values=(0, 0.02),
            comment="Convergence criterion for "
            "termination of DOM azimuthal loop",
            theme=1)
        self["DOM_OPDEP_THRESHOLD"] = OptionObject(
            value=0,
            otype=float,
            min_max_values=(0, 20),
            comment="DOM ignores levels below "
            "this optical depth:"
            "10. is reasonable,"
            " not applied if <= 0",
            theme=1)
        self["DOM_RAYLEIGH"] = OptionObject(theme=1)
        self["OZONE_DATA"] = OptionObject(theme=2)
        self["CO2_DATA"] = OptionObject(theme=2)
        self["N2O_DATA"] = OptionObject(theme=2)
        self["CO_DATA"] = OptionObject(theme=2)
        self["CH4_DATA"] = OptionObject(theme=2)
        self["SO2_DATA"] = OptionObject(theme=2)
        self["CLOUD_OVERLAP"] = OptionObject(theme=1, hidden=True, value=1)
        self["CC_LOW_CLOUD_TOP"] = OptionObject(
            value=750,
            otype=float,
            min_max_values=(0, 1100),
            hidden=True,
            theme=1)
        self["GRID_BOX_AVG_CLOUD"] = OptionObject(
            theme=1,
            value=False,
            hidden=True,
            comment="Switch to supply grid-box average cloud "
                    "concentration (true) or "
                    "cloud concentration in cloudy fraction of "
                    "each layer (false)")
        self["ADDPC"] = OptionObject(theme=3,
                                     comment="Switch to enable PC-RTTOV")
        self["IPCBND"] = OptionObject(value=1,
                                      otype=dict,
                                      odict=self.fullIPCBND,
                                      theme=3,
                                      comment="PC spectral band")
        self["IPCREG"] = OptionObject(value=1,
                                      otype=dict,
                                      odict=self.fullIPCREG,
                                      theme=3,
                                      comment="PC predictor channel set")
        self["NPCSCORES"] = OptionObject(
            value=300,
            otype=int,
            min_max_values=(-1, 5000),
            theme=3,
            comment="Number of PC scores to compute")
        self["ADDRADREC"] = OptionObject(
            theme=3,
            comment="Switch for calculation "
            "of reconstructed radiances")

        self["FASTEM_VERSION"] = OptionObject(
            theme=4,
            otype=dict,
            value=6,
            odict={0: "TESSEM2", 1: "1", 2: '2', 3: '3', 4: '4',
                   5: '5', 6: '6'},
            comment="TESSEM2 model or FASTEM version (1-6)")
        self["SUPPLY_FOAM_FRACTION"] = OptionObject(
            theme=4,
            comment="Supply a foam fraction to FASTEM"
        )
        theComment = "Switch to enable input of cloud liquid water profile"
        self["CLW_DATA"] = OptionObject(theme=4,
                                        comment=theComment,
                                        hidden=True)
        self["CLW_SCHEME"] = OptionObject(
            theme=4,
            otype=dict,
            value=2,
            odict=self.dictOptions[
                "CLW_SCHEME"],
            comment="CLW scheme",
            hidden=True)
        self["CLW_CLOUD_TOP"] = OptionObject(
            theme=4,
            value=322,
            min_max_values=(0, 1100),
            otype=float,
            comment="Lower pressure limit for CLW absorption",
            hidden=True)

        self["ADDINTERP"] = OptionObject(
            theme=5,
            comment="Switch to enable RTTOV interpolator")
        self["INTERP_MODE"] = OptionObject(
            theme=5,
            otype=dict,
            value=5,
            odict={1: "1 Rochon/Rochon OD",
                   2: '2 Log-Lin/Log-lin OD',
                   3: '3 Rochon/Log-lin OD',
                   4: '4 Rochon/Rochon WF',
                   5: '5 Rochon/Log-lin WF'},
            comment="Interpolation mode (see user guide)")

        self["LGRADP"] = OptionObject(
            theme=5,
            comment="Switch to make pressure an active"
            " variable in TL/AD/K models")

        self["SPACETOP"] = OptionObject(
            theme=5,
            value=True,
            comment="Switch to assume space boundary at"
            " top-most input pressure level")
        self["REG_LIMIT_EXTRAP"] = OptionObject(
            theme=5,
            value=True,
            comment="Switch to extrapolate input profiles"
            " using regression limits")

        self.options_list_logical = []
        self.options_list = []
        self.optionsThemesList = {}
        for theme in self.optionsThemes:
            self.optionsThemesList[theme] = []
        for option, v in list(self.items()):
            # at this points cle_ATTRIBUTES does not exist
            self.options_list.append(option)
            self.optionsThemesList[self.optionsThemes[v.theme]].append(option)
            if v.otype == bool:
                self.options_list_logical.append(option)

# <type 'dict'>
        for option in list(self.keys()):
            self[
                option + '_ATTRIBUTE'] = {'COMMENT': self[option].comment,
                                          'UNITS': 'n/a'}

    def loadh5(self, h5):
        # Load is driven by the content of the file
        for cle in list(h5):
            if(cle in self.options_list_logical):
                self[cle].value = loadDset_log(h5[cle])
            else:
                # Reals, Integers
                self[cle].value = loadDset_arr(h5[cle])
                self[cle + '_ATTRIBUTE'] = getAttributes(h5[cle])
            checkAttributes(self[cle + '_ATTRIBUTE'])

    def saveh5(self, h5, path):
        """Save option class to an HDF5 file h5 under path directory
        Attributes are saved at the same time datasets are saved
        """
        for cle in list(self.keys()):

            if(re.search('_ATTRIBUTE', cle) or self[cle] is None):
                continue

            if(cle in self.options_list_logical):
                # Logicals are specific
                if(self[cle].value):
                    val = 1
                else:
                    val = 0
                dset = h5.create_dataset(path + cle, data=val)
            else:
                dset = h5.create_dataset(path + cle, data=self[cle].value)

            saveAttributes(h5[path + cle], self[cle + '_ATTRIBUTE'])

    def display(self):
        """
        Display a full class content to terminal
        Members are printed by alphabetic order
        """
        for cle in sorted(self.options_list):
            print(cle, '\t =>\t', self[cle].value, type(self[cle].value))
            print('\t =>\t', self[cle + "_ATTRIBUTE"])

    def print_value(self):
        """
        Display a full class content to terminal
        Members are printed by alphabetic order
        """
        for cle in sorted(self.keys()):
            if not (cle[-9:] == "ATTRIBUTE"):
                print(cle, '\t =>\t', self[cle].value)


if __name__ == '__main__':
    from core import rttov_gui_data_test_dir
    import os
    data_test_dir = rttov_gui_data_test_dir()
    ofile = os.path.join(data_test_dir, "options_res.h5")
    o1 = Option()

    # Display option
    o1.display()
    o2 = o1.__deepcopy__()

    # change something in o1
    o1["ADDAEROSL"].value = not(o2["ADDAEROSL"].value)
    print("Option/value by alphabetic order:")
    o1.print_value()

    # Save option o1 to HDF5
    of = h5py.File(ofile, 'w')
    o2.saveh5(of, '/OPTIONS/')

    of.close()

    # read option again from file:
    f = h5py.File(ofile, 'r')
    h5 = f['/OPTIONS/']
    o3 = Option()
    o3.loadh5(h5)

    print('o2["ADDAEROSL"] =', o2["ADDAEROSL"].getValue())
    # Close HDF file
    f.close()
    print(o3["ADDAEROSL"].getValue())
    print(o1["ADDAEROSL"].getValue())
    if (o3["ADDAEROSL"].getValue() != o1["ADDAEROSL"].getValue()):
        print("test 1 OK")
