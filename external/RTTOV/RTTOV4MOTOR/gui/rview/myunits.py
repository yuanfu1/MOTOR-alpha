# -*- coding: utf-8 -*-
'''
Created on Sep 26, 2016

@author: roquetp
'''

# matplotlib syntax for units
wavenumberLabel = 'wavenumber $(cm^{-1})$'
radianceLabel = "$mW.(cm^{-1})^{-1}.sr^{-1}.m^{-2}$"
micron = u"µm"

# Unicode syntax for status bar :
# https://fr.wikipedia.org/wiki/Exposants_et_indices_Unicode
exponentMinus1 = u'\u207B\u00B9'
exponentMinus2 = u'\u207B\u00B2'
exponentMinus3 = u'\u207B\u00B3'
radianceUnit = (u'mW.(cm' + exponentMinus1 + u')' + exponentMinus1 +
                u'.sr' + exponentMinus1 + u'.m' + exponentMinus2)
cmMinus1 = u'cm' + exponentMinus1
# unicode syntax for status bar :
# gases
# Label depends of profile(GAS_UNITS)
gasUniLabel = {-1: "ppmv over dry air",
               0: "ppmv over dry air",
               1: u'kg.kg' + exponentMinus1 + u' over moist air',
               2: "ppmv over moist air",
               None: "ppmv over moist air"

               }
gasUnitLabelForJacobian = {
    -1: "ppmv" + exponentMinus1,
    0: "ppmv" + exponentMinus1,
    1: "kg.kg" + exponentMinus1,
    2: "ppmv" + exponentMinus1,
    None: "ppmv"
}
# aerosols : 0: number density (cm-3), 1: kg/kg
# depend of MMR_CLDAER
aerUniLabel = {0: u'number density (cm' + exponentMinus3 + u')',
               1: u'kg.kg' + exponentMinus1}

# cloud : 0: layer mean content(g/m3), 1: kg/kg
cldUniLabel = {0: u"layer mean content (g.m" + exponentMinus3 + u")",
               1: u"kg.kg" + exponentMinus1}

# jacobien units unicode enconding


def makePretty(gas, strunit, gas_unit):
    # look for input perturbation unit
    if "mw" in str(strunit):
        punit = radianceUnit
    else:
        punit = "K"
    if gas == "T":
        itemUnit = "K" + exponentMinus1
    else:
        itemUnit = gasUnitLabelForJacobian[gas_unit]
    return punit + "." + itemUnit

# jacobien units for matplotlib


def makePrettyForMatplotlib(gas, strunit, gas_unit):
    # look for input perturbation unit
    if "mw" in str(strunit):
        punit = radianceLabel
    else:
        punit = "$K$"
    if gas == "T":
        itemUnit = "$K^{-1}$"
    else:
        if gas_unit == 1:
            itemUnit = "$kg.kg^{-1}$"
        else:
            itemUnit = "$ppmv^{-1}$"
    return punit + "." + itemUnit
