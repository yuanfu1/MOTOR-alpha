#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" defaut values for profiles """
import numpy as np
import os

defaultppmv = {'Q':    2e4, 'O3': 0.5,   'CO2': 392.0,
               'CH4': 1809.0, 'CO': 0.174, 'N2O': 325.0, 'SO2':  0.000319807}


def myinterp(val_in, press_in, press_out):
    logpin = np.log(press_in)
    logpout = np.log(press_out)
    val_out = np.interp(x=logpout, xp=logpin, fp=val_in)
    return val_out


def get_default_gas(gas, press_out):
    """ get defaut gas concentration from mean profile 83 """
    try:
        dirname = os.path.join(os.environ["RTTOV_GUI_PREFIX"], "rttov", "data")
        fname = os.path.join(dirname,
                             gas + ".txt")
        pname = os.path.join(dirname, "P.txt")
        gas_in = np.genfromtxt(fname)
        press_in = np.genfromtxt(pname)
        gas_con = myinterp(gas_in, press_in, press_out)
    except Exception:
        gas_con = defaultppmv[gas]

    return gas_con


if __name__ == '__main__':
    import matplotlib
    matplotlib.use("wxagg")
    import matplotlib.pyplot as plt
    import h5py
    try:
        rttov_profile_dir = os.environ["RTTOV_GUI_PROFILE_DIR"]
    except Exception:
        rttov_profile_dir = "../../rttov_test/profile-datasets-hdf/"

    fname = os.path.join(rttov_profile_dir, "cldaer101lev_allgas.H5")
    if not os.path.isfile(fname):
        print("Error", fname, "not found")
    h5f = h5py.File(fname, "r")
    p_in = h5f["PROFILES/0001/P"]

    so2 = get_default_gas("SO2", p_in)
    if so2 is not None:
        print(so2.shape)
        levels = [x for x in range(so2.shape[0], 0, -1)]
        y = np.array(levels)
        print(y.shape)
        plt.plot(so2, y)
        plt.title("default SO2 profile")
    else:
        print("error for so2")

    plt.show()
    fname = os.path.join(rttov_profile_dir, "cldaer50lev_co2o3.H5")
    if not os.path.isfile(fname):
        print("Error", fname, "not found")
    h5f = h5py.File(fname, "r")
    p_in = h5f["PROFILES/0001/P"]
    so2 = get_default_gas("SO2", p_in)
    if so2 is not None:
        print(so2.shape)
        levels = [x for x in range(so2.shape[0], 0, -1)]
        y = np.array(levels)
        print(y.shape)
        plt.plot(so2, y)
        plt.title("default SO2 profile")
    else:
        print("error for so2")
    plt.show()
    print("test OK")
