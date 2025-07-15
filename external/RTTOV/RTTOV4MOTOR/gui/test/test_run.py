#! /usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module tests the rttov Fortran and prepare the files which are needed
for unitary tests for each windows
This files are kept in the rttov_gui_test_data directory
Notice that you should have the PYTHONPATH environment variable
containig the directorys that holds the rttov Python package
as well as the RTTOV_GUI_PREFIX environment vairaible with the
name of the top gui directory.
This test does:
- load US76 profile
- update latitude/logitude/surface_type in order to be over land
- load MSG3 coefficients for Solar and Thermal channels.
- retrieve surafce emissivity/reflectance for all channels
- run forward model
- display surface emissivity/reflectance
- display radiances
- display brightness temperatures
- display relectances
- drops the coefficients
- load PC coefficients
- run pc
- run k pc with addradrec = True
- run k pc with addradrec = False
"""

import h5py
import os
import sys
from shutil import copy2

from rttov import emissivity, reflectance, profile, \
    option, radiance, getcoefval

import rttov_gui_f2py

from rmodel import project

try:
    myfile = os.environ["RTTOV_GUI_PROFILE_DIR"] + "/us76_50lev_allgas.H5"
except Exception:
    print("RTTOV_GUI_PROFILE_DIR is not defined")
    sys.exit(1)
try:
    coeff_dir = os.environ["RTTOV_GUI_COEFF_DIR"]
except Exception:
    print("RTTOV_GUI_COEFF_DIR is not defined")
    sys.exit(1)
try:
    prefix = os.environ["RTTOV_GUI_PREFIX"]
except Exception:
    print("RTTOV_GUI_PREFIX is not defined")
    sys.exit(1)

res_dir = os.path.join(prefix, "rttov_gui_data_test")
try:
    wrk_dir = os.environ["RTTOV_GUI_WRK_DIR"]
    if not os.path.exists(res_dir):
        print("create", res_dir, "directory")
        os.makedirs(res_dir)
except Exception:
    print("RTTOV_GUI_WRK_DIR is not defined")
    sys.exit(1)

if not os.path.isdir(res_dir):
    print((res_dir + "does not exist"))
    sys.exit(1)

if not os.access(res_dir, os.W_OK):
    print((res_dir + "not writable"))
    sys.exit(1)

# Read US76 profile
nprof = profile.getNumberOfProfiles(myfile)
print(("number of profiles found :", nprof))

p1 = profile.Profile()
o1 = option.Option()

iprof = 6
p1.loadProfileNumber(myfile, iprof)


# ADD solar calculations
o1['ADDSOLAR'].value = True
o1['ADDINTERP'].value = True

# change latitude longitude date and surface type in order to
# put the profile over land in France
# Black body cloud 800hPa cloud fraction 0.40
p1['LATITUDE'] = 45.0
p1['LONGITUDE'] = 5.0
p1['SKIN']['SURFTYPE'] = 0
p1['AZANGLE'] = 45.0
p1['SUNAZANGLE'] = 100.0
p1['SUNZENANGLE'] = 30.0
p1['ZENANGLE'] = 30.0
p1['CFRACTION'] = 0.4
p1['CTP'] = 800.0
o1["CO_DATA"].value = False
o1["CO2_DATA"].value = False
o1["N2O_DATA"].value = False
o1["CH4_DATA"].value = False
o1["ADDSOLAR"].value = False

profileFile = os.path.join(res_dir, "profile.H5")
surfaceFile = os.path.join(res_dir, "surface.H5")
radianceFile = os.path.join(res_dir, "radr.H5")
kmatFile = os.path.join(res_dir, "kmat.H5")
transFile = os.path.join(res_dir, "trans.H5")

# write file "profile.H5" with 1 profile and an option dataset
# Output p1 to HDF5 file
ofile = profileFile
of = h5py.File(ofile, 'w')
p1.saveh5(of, "/PROFILES/0001/")
o1.saveh5(of, '/OPTIONS/')
of.close()

print("load coefficient file")
# Read RTTOV v12 coefficient
nchannels, err = rttov_gui_f2py.rttov_gui_load(
    [0],
    coeff_dir + "/rttov9pred54L/rtcoef_msg_3_seviri.dat", "", "", "", "")
print(("nchannels", nchannels))
listchan = getcoefval.rttov_get_coef_val_i1('FF_ORI_CHN')

print(listchan)

e1 = emissivity.Emissivity()
r1 = reflectance.Reflectance()
e1.setEmissivity(nchannels)
r1.setReflectance(nchannels)
ofile = surfaceFile
of = h5py.File(ofile, 'w')
e1.saveh5(of, '/EMISSIVITY/')
r1.saveh5(of, '/REFLECTANCE/')
of.close()

# Get surface emissivity/reflectance from Atlas

surfaceFile = ofile
err = rttov_gui_f2py.rttov_gui_get_emisbrdf(os.environ["RTTOV_GUI_EMISS_DIR"],
                                            profileFile,
                                            surfaceFile,
                                            1, 1)


# Read emissivity and Reflectance; print and plot
fileName = os.path.join(res_dir, "surface.H5")
copy2(os.path.join(wrk_dir, "kmat.h5"), res_dir)

f = h5py.File(fileName, 'r')
h5 = f['/EMISSIVITY/']
e1 = emissivity.Emissivity()
e1.loadh5(h5)
h5 = f['/REFLECTANCE/']
r1 = reflectance.Reflectance()
r1.loadh5(h5)
f.close()

print("\n \n e1 read from surface.H5 ")
e1.display()
print("\n")
print("\n \n r1 read from surface.H5 ")
r1.display()
print("\n")

# Now RUN RTTOV Direct model
nthreads = 2
err = rttov_gui_f2py.rttov_gui_run(profileFile,
                                   surfaceFile,
                                   radianceFile,
                                   kmatFile,
                                   transFile,
                                   "DIRECT", nthreads)
assert err == 0

copy2(os.path.join(wrk_dir, "kmat.h5"), res_dir)
copy2(os.path.join(wrk_dir, "radr.h5"), res_dir)

# Read emissivity and Reflectance; print and plot

f = h5py.File(surfaceFile, 'r')
h5 = f['/EMISSIVITY/']
e1 = emissivity.Emissivity()
e1.loadh5(h5)
h5 = f['/REFLECTANCE/']
r1 = reflectance.Reflectance()
r1.loadh5(h5)
f.close()
print("\n \n e1 read from surface.H5 after RTTOV run ")
e1.display()
print("\n")
print("\n \n r1 read from surface.H5 after RTTOV run  ")
r1.display()
print("\n")


# Read Radiance results
f = h5py.File(radianceFile, 'r')
h5 = f['/RADIANCE/']

# rad is a rttov_hdf_mod radiance instance
rad = radiance.Radiance()
rad.loadh5(h5)
f.close()
rad.display()


proj = project.Project()
profileName = os.environ["RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
err = proj.openProfile(profileName, 1)
assert err == 0
coefFile = proj.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
    "/rttov9pred101L/rtcoef_metop_2_iasi_pcrttov_compat.H5"
pcFile = proj.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
    "/pc/pccoef_metop_2_iasi_landsea_trace_nlte.H5"
proj.myCoeffs.fileName["standard"] = coefFile
proj.myCoeffs.fileName["PC"] = pcFile
proj.myOption["ADDPC"].value = True
proj.ctrlCoherence()
print(proj.myCoeffs.fileName)
err = proj.loadCoefficients()
assert err == 0
err = proj.runPCK()
assert err == 0
copy2(os.path.join(wrk_dir, "pckmat.h5"), res_dir)
proj.myOption["ADDRADREC"].value = True
err = proj.runPC()
assert err == 0
copy2(os.path.join(wrk_dir, "pc.h5"), res_dir)
copy2(os.path.join(wrk_dir, "trns.h5"), res_dir)

# test data for pcView
proj.myOption["ADDRADREC"].value = False
err = proj.runPC()
assert err == 0
copy2(os.path.join(wrk_dir, "pc.h5"), os.path.join(res_dir, "pc1.h5"))
proj.myOption["OZONE_DATA"].value = True
err = proj.runPC()
assert err == 0

copy2(os.path.join(wrk_dir, "pc.h5"), os.path.join(res_dir, "pc2.h5"))
proj.myOption["CO2_DATA"].value = True
err = proj.runPC()
assert err == 0
copy2(os.path.join(wrk_dir, "pc.h5"), os.path.join(res_dir, "pc3.h5"))

print(">>>>>>>>>>TEST OK")
