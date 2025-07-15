'''
Created on Sep 8, 2015

@author: pascale
'''
# RTTOV GUI Colors can be customised
#  you can modify this file as your convenience but beware to respect the
#  python syntax
#  because this file is part of RTTOV GUI source code
#  to know the list of matplotlib available colors you can execute  :
#  python rview/colors.py


def print_matplotlib_color_list():
    import matplotlib
    for name, hexacode in matplotlib.colors.cnames.items():
        print((name, hexacode))


# colors for the profile and kprofile windows
# profile colors
profileItemColors = {'Q': "red", 'O3': "green", 'CO2': "magenta",
                     'CH4': "cyan",
                     'CO': "yellow", 'N2O': "black", "T": "blue",
                     "SO2": "brown",
                     'INSO': "coral", 'WASO': "crimson", 'SOOT': "chocolate",
                     'SSAM': "yellow", 'SSCM': "orange", 'MINM': "black",
                     'MIAM': "aqua", 'MICM': "blueviolet", 'MITR': "brown",
                     'SUSO': "cadetblue", 'VOLA': "magenta", 'VAPO': "green",
                     'ASDU': "red", 'STCO': "red", 'STMA': "green",
                     'CUCC': "magenta", 'CUCP': "cyan", 'CUMA': "yellow",
                     'CIRR': "blue", 'CFRAC': "brown", "CLW": "green",
                     'ICEDE': "blue", 'CLWDE': "blueviolet",
                     'BCAR': "black", 'DUS1': "chocolate",
                     'DUS2': "grey", 'DUS3': "brown", 'SULP': "crimson",
                     'SSA1': "blue", 'SSA2': "blueviolet", 'SSA3': "aqua",
                     'OMAT': "orange"}
# profile markers
profileItemMarkers = {'Q': "_", 'O3': "_", 'CO2': "_", 'CH4': "_", 'CO': "_",
                      'N2O': "_", "SO2": "_",
                      "T": "_", 'INSO': "+", 'WASO': "+",
                      'SOOT': "+", 'SSAM': "+", 'SSCM': "+", 'MINM': "+",
                      'MIAM': "+", 'MICM': "+", 'MITR': "+",  'SUSO': "+",
                      'VOLA': "+", 'VAPO': "+", 'ASDU': "+",
                      'STCO': "*", 'STMA': "*", 'CUCC': "*", 'CUCP': "*",
                      'CUMA': "*", 'CIRR': "*", 'CFRAC': "*", "CLW": "*",
                      'ICEDE': "*", 'CLWDE': "+", 'CLW': '+',
                      'BCAR': "+", 'DUS1': "+",
                      'DUS2': "+", 'DUS3': "+", 'SULP': "+",
                      'SSA1': "+", 'SSA2': "+", 'SSA3': "+",
                      'OMAT': "+"}

# colors for K profile window (old K curves )
oldKColors = ["orange", "green", "blueviolet"]
# style for K profile window (old K curves
StyleForPreviousValues = ['-', '--', ':']


# colors for radiance window
radFrameRadTotal = 'blue'
radFrameRadClear = 'blue'
radFrameBT = 'green'
radFrameBTClear = 'green'
radFrameRefl = 'cyan'
radFrameReflClear = 'cyan'

# colors for the surface window (emissivity and surface BRDF in and out)
surfaceIn = "blue"
surfaceOut = "red"
# colors for RbtView
btViewColors = {"T-Bg": "red",
                "True": "black",
                "Bg": "blue"}

# pcView colors : colors for the PC score windows (colors for successive runs)
pcViewColors = ["blue", "red", "blueviolet", "green", "black",
                "cyan", "green", "red", "yellow", "orange", "black",
                "aqua", "crimson", "brown", "cadetblue", "magenta",
                "green", "brown", "blue"]

# colors for the KPC profile window
kpcviewColors = ["black", "blue", "cyan", "green", "red", "black", "blue",
                 "cyan", "green", "red", "yellow", "orange", "black",
                 "aqua", "BlueViolet", "brown", "CadetBlue", "magenta",
                 "green", "brown"]
kpcviewLinestyles = ['-', '-', '-', '-', '-', '--', '--',
                     '--', '--', '--', '-', '--', '--', '--', '--', '--']
kpcviewLinewidths = [2.0, 2.0, 1.5, 1.5, 1, 1, 1, 1, 1, 1, 1, 1]

# Colors for 1Dvar profileFrame
# colors for the retrieved profile (list must contain ten different colors
RttovGui1DvarStepColor = ["red", "cyan", "magenta", "violet",
                          "darkviolet", "mediumpurple", "indigo",
                          "navy", "royalblue", "green"]

# colors for True and Background Colors (T and Q)
RttovGui1DvarItemColors = {"true": "black",
                           "background": "blue"}

if __name__ == "__main__":
    print_matplotlib_color_list()
