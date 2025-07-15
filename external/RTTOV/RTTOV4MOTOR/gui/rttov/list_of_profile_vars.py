""" NWP-SAF RTTOV GUI v11
    May 2014
    List of profile variables that can be set by the user
    The variable names are the ones described in the RTTOV Users guide
    for the profile structure (Annex O), capital letters; except for:
     - clouds where 2D cloud array is replaced by 1D cloud arrays, one for
       each cloud short name (table 16 of Users guide)
     - aerosols where 2D aerosols array is replaced by 1D aerosol arrays,
       one for each aerosol short name (table 19 of Users guide)
    Units should be the ones described in RTTOV Users Guide for profile
    structure (Annex O)

    Take this file as a template, do not change the variable name,
    just fill values or comment lines
    Arrays like numpy.array([...]) should be replaced by user values
    like numpy.array([10., 20., 30.]) or commented

    Default values if some variables are not provided
    Variable attributes like Comment and Units are set internally to Profile
    class
"""

" Mandatory arrays: P (hPa), T(K), Q(ppmv) on levels"
self["P"] = numpy.array([...])
self["T"] = numpy.array([...])
self["Q"] = numpy.array([...])

"--------------------------------------------"
" Optional profile variables "
"--------------------------------------------"

" Other Gases (ppmv on levels)"
self["O3"] = numpy.array([...])
self["CO2"] = numpy.array([...])
self["CH4"] = numpy.array([...])
self["CO"] = numpy.array([...])
self["N2O"] = numpy.array([...])

" Aerosols (1/cm3 on layers)"
self["INSO"] = numpy.array([...])  # Insoluble
self["WASO"] = numpy.array([...])  # Water soluble
self["SOOT"] = numpy.array([...])  # Soot
self["SSAM"] = numpy.array([...])  # Sea salt (acc mode)
self["SSCM"] = numpy.array([...])  # Sea salt (coa mode)
self["MINM"] = numpy.array([...])  # Mineral (nuc mode)
self["MIAM"] = numpy.array([...])  # Mineral (acc mode)
self["MICM"] = numpy.array([...])  # Mineral (coa mode)
self["MITR"] = numpy.array([...])  # Mineral transported
self["SUSO"] = numpy.array([...])  # Sulphated droplets
self["VOLA"] = numpy.array([...])  # OPAC Volcanic ash
self["VAPO"] = numpy.array([...])  # New Volcanic ash
self["ASDU"] = numpy.array([...])  # Asian dust

" Clouds (g/m3 on layers)"
self["STCO"] = numpy.array([...])  # Stratus Continental
self["STMA"] = numpy.array([...])  # Stratus Maritime
self["CUCC"] = numpy.array([...])  # Cumulus Continental Clean
self["CUCP"] = numpy.array([...])  # Cumulus Continental Polluted
self["CUMA"] = numpy.array([...])  # Cumulus Maritime
self["CIRR"] = numpy.array([...])  # Cirrus
# Cloud Fraction (should be set if any cloud)
self["CFRAC"] = numpy.array([...])
self["ICEDE_PARAM"] = 2  # Ice effective diameter parameterisation
self["ICE_SCHEME"] = 2   # Ice cloud scheme
self["CLWDE_PARAM"] = 1  # Cloud liquid water effective diameter parameterisation
self["CLW_SCHEME"] = 1   # Cloud liquid water scheme

" Skin variables "
self["SKIN"]["T"] = self["T"][-1]  # (K)
self["SKIN"]["SURFTYPE"] = 1  # (0=Land, 1=Sea, 2=sea-ice)
self["SKIN"]["WATERTYPE"] = 1  # (0=fresh water, 1=ocean water)
self["SKIN"]["SNOW_FRACTION"] = 0.  # [0,1]
self["SKIN"]["SALINITY"] = 37  # (%o)
self["SKIN"]["FASTEM"] = numpy.array(
    [0., 0., 0., 0., 0.])  # (5 parameters Land/sea-ice)

" 2m and 10m air variables "
self["S2M"]["T"] = self["T"][-1]  # (K)
self["S2M"]["Q"] = self["Q"][-1]  # (ppmv)
self["S2M"]["P"] = self["P"][-1]  # (hPa)
self["S2M"]["U"] = 0  # (m/s)
self["S2M"]["V"] = 0  # (m/s)
self["S2M"]["WFETC"] = 100000  # (m)

" Simple cloud "
self["CTP"] = 500.0  # (hPa)
self["CFRACTION"] = 0.0  # [0,1]    Clear sky is the default

" Viewing geometry "
self["AZANGLE"] = 0.  # (deg)
self["ELEVATION"] = 0.  # (km)
self["LATITUDE"] = 49.738  # (deg)  Lannion is 48.750,  -3.470
self["LONGITUDE"] = -3.473  # (deg)  Exeter is  50.726,  -3.476

self["SUNAZANGLE"] = 0.  # (deg)
self["SUNZENANGLE"] = 0.  # (deg)
self["ZENANGLE"] = 0.  # (deg)

" Magnetic field "
self["BE"] = 0.3  # (Gauss)
self["COSBK"] = 1.

" Mislaneous "
self["ID"] = "This is my profile"
self["DATE"] = numpy.array([2014, 0o4, 30], dtype=int)  # Year, Month, Day
self["TIME"] = numpy.array([12, 0, 0], dtype=int)     # Hour, Minute, Second


"--------------------------------------------"
" Other profile variables that should NOT be provided"
"--------------------------------------------"

" Levels and Layers are calculated "
