# -*- coding: utf-8 -*-
'''
Created on Nov 24, 2014

@author: pascale
'''

import rmodel
import numpy as np
import logging
import rttov
import copy
import matplotlib.pyplot as plt
import os
import h5py
import shutil
from r1Dvar import r1dvarObjects


q_mixration_to_ppmv = 1.60771704e+6


def vector2profile(vector, exprofile):
    """ convert vector to profile """
    """ take exprofile to begin with """
    """ a vector contains only T values, """
    """ lnq bottom lnq vales , Tsurf, lnq surf and Tskin"""
    """ temperature, and ln q are stored from top of atmosphere to ground """

    nlevels = vector.shape[0] - 29 - 3
    if nlevels != exprofile['NLEVELS']:
        logging.error("vector and profile are not of the same dimension")
        return None

    profile = copy.deepcopy(exprofile)

    for i in range(nlevels):
        profile["T"][i] = vector[i]

    for i in range(29):
        profile["Q"][nlevels - 29 + i] = vector[nlevels + i]

    profile["S2M"]["T"] = vector[nlevels + 29]
    profile["S2M"]["Q"] = vector[nlevels + 29 + 1]
    profile["SKIN"]["T"] = vector[nlevels + 29 + 2]

    return profile


class Project1dvar(object):

    def __init__(self, pTrue, pBg, satellite, instrument,
                 matrixBfile=None, channel_list=None,
                 changeTheBackground=False):
        """ initialize a 1DvarProject with a True and a Background project """
        """ if channel_list is present : we work on a restricted list """
        """ of channel """
        """ but in fact channel_list is set from the Coefficient """
        """ object of pBg """
        """ so we use this list to extract a sub-matrix from the Rmatrix """
        """ satellite must be = platform-num (ex : metop-1 , noaa-19) """
        """ if changeTheBackground=False we don not change the Background """
        """ if changeTheBackground=False at each new step the previously  """
        """ retrieved profile becomes the new background """

        self.changeTheBackground = changeTheBackground
        self.ENV = {}
        self.pTrue = pTrue
        self.pBg = pBg
        # set distinct filenames for the true and the background project
        # (we will run rttov with different files)
        self.pTrue.setFileNameMark("True")
        self.pBg.setFileNameMark("Bg")
        logging.info("Initialize a Project1dvar object")
        logging.info("are pBg coefficients loaded ? " +
                     str(self.pBg.myCoeffs.loadCoeffs))
        self.satellite = satellite
        self.instrument = instrument

        self.BgInitialProfile = copy.deepcopy(self.pBg.myProfile)
        self.BmatrixFileName = None
        self.RmatrixFilename = None

        self._SetConfig()
        self.factorB = 1.
        self.factorR = 1.
        self.MaxNoise = 1.
        self.retrievedVectors = []
        self.retrievedProfiles = []
        self.Rmatrix = Rmatrix()
        self.Bmatrix = Bmatrix()
        # read R matrix
        self.initCoeffsandReadRmatrix()
        # Bmatrix
        if matrixBfile:
            self.BmatrixFileName = matrixBfile
        if not os.path.exists(self.BmatrixFileName):
            logging.error(
                "Sorry something is wrong : I cannot find the Bmatrix :"
                " check the presence of the nwpsaf-1dvar files ")
            raise IOError
        # read the Bmatrix
        B = Bmatrix()
        B.read_matrices(self.BmatrixFileName)
        # The file contains 4 matrix but according to Peter the 4 matrix are
        # identical
        self.Bmat = B.matrices[1]
        self.step_counter = 0

    def initCoeffsandReadRmatrix(self):
        """ initialise myCoeffs for pBg project and read the R matrix """
        """ the R matrix is computed according to the list of channels """
        if not self.pBg.myCoeffs.loadCoeffs:
            logging.info("load coefficients for pBg")
            err = self.pBg.loadCoefficients()
            if err != 0:
                raise RuntimeError("Cannot load coefficient in r1dvar project")
        nbChannels = self.pBg.myCoeffs.nchannels
        logging.info("nbChannels: {}".format(nbChannels))
        # initialise the  1dvarproject channel_list used by readRmatrix
        my_chan_list = self.pBg.myCoeffs.getFF_ORI_CHN()
        # must start at 0
        self.channel_list = [my_chan_list[i] -
                             1 for i in range(len(my_chan_list))]
        # read the R matrix

        self.readRmatrix(self.satellite, self.instrument)

    def reinitCoeff(self, satellite, instrument, coeffs):
        """ must be called if we change the coefficients """
        """ (select channel or change instrument ) """
        """  satellite must be = platform-num (ex : metop-1 , noaa-19) """
        self.pBg.myCoeffs = coeffs
        self.satellite = satellite
        self.instrument = instrument

        self.initCoeffsandReadRmatrix()

    def _SetConfig(self):
        try:
            self.ENV["NWPSAF_1DVAR_HOME"] = os.environ("NWPSAF_1DVAR_HOME")
        except:
            logging.warning(
                "The environment variable NWPSAF_1DVAR_HOME is not set :"
                " no worries I can use my matrix ..")
            self.ENV["NWPSAF_1DVAR_HOME"] = os.environ[
                "RTTOV_GUI_PREFIX"] + "/r1Dvar/data"
        self.BmatrixFileName = self.ENV[
            "NWPSAF_1DVAR_HOME"] + "/Sample_Bmatrices/Bmatrix_54L"

        self.RmatrixFilenames = {
            "AIRS": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                 "AIRS_COEFFS_DIR/Rmatrix_orig"),
            "ATMS": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                 "ATMS_COEFFS_DIR/Rmatrix_orig"),
            "AMSU-A": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                   "ATOVS_COEFFS_DIR/Rmatrix_orig"),
            "AMSU-B": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                   "ATOVS_COEFFS_DIR/Rmatrix_orig"),
            "HIRS": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                 "ATOVS_COEFFS_DIR/Rmatrix_orig"),
            "MHS": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                "ATOVS_COEFFS_DIR/Rmatrix_orig"),
            "mhs": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                "ATOVS_COEFFS_DIR/Rmatrix_orig"),
            "IASI": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                 "IASI_COEFFS_DIR/Rmatrix_orig"),
            "CRIS": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                 "CrIS_COEFFS_DIR/Rmatrix_orig"),
            "SSMIS": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                  "SSMIS_COEFFS_DIR/Rmatrix_orig"),
            "CrIS": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                 "CrIS_COEFFS_DIR/Rmatrix_orig"),
            "amsua": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                  "ATOVS_COEFFS_DIR/Rmatrix_orig"),
            "amsub": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                  "ATOVS_COEFFS_DIR/Rmatrix_orig"),
            "hirs": os.path.join(self.ENV["NWPSAF_1DVAR_HOME"],
                                 "ATOVS_COEFFS_DIR/Rmatrix_orig")
        }

    def readRmatrix(self, satellite, instrument, filename=None):
        """ read the observation error matrix file and init self.Rmat """
        """ if filename is specified then it has the priority """
        """ first step : determine the file name from the instrument """
        """ (stored in rge Rmatrixfilenames dictionary"""
        """ create a Rmatrix objet which will contains all the Rmatrix """
        """ read in a file """
        """ then determine the actual Rmat from satellite and instrument """
        """ ATOVS case is particular : must retrieve HIRS, AMSU-A and MHS"""
        """ for some other files it is trivial : they contains just one """
        """ matrix """
        """  satellite must be = platform-num (ex : metop-1 , noaa-19) """
        # when this method is call self.channel_list exit and is ready (no need
        # to shift indices)
        logging.info("read R matrix for " + satellite + " " + instrument)
        isATOV = False
        if instrument in ("hirs", "amsua", "amsub", "mhs"):
            isATOV = True
        else:
            if instrument == "cris":
                instrument = "CrIS"
            else:
                instrument = instrument.upper()
        if satellite == "metop-2":
            satellite = "MetOp-A"
        else:
            if satellite == "metop-1":
                satellite = "MetOp-B"
            else:
                satellite = satellite.upper()

        logging.info("instrument :" + str(instrument) + " " +
                     str(satellite) + " is ATOV ?: " + str(isATOV))
        if filename is not None:
            self.RmatrixFilename = filename
        else:
            self.instrument = instrument
            self.satellite = satellite

            self.RmatrixFilename = self.RmatrixFilenames[instrument]
            try:
                self.RmatrixFilename = self.RmatrixFilenames[instrument]
                logging.info(" RmatrixFilename : " + self.RmatrixFilename +
                             " for " +
                             instrument + " " + satellite +
                             " is ATOV ?: " + str(isATOV))
            except (KeyError):
                raise ValueError("wrong instrument" + str(instrument))

            self.Rmatrix = Rmatrix()
            logging.info("Read the Rmatrix")
            self.Rmatrix.read(self.RmatrixFilename)
            logging.info("R matrix successfully read")

            try:
                if len(list(self.Rmatrix.Rmatrix.keys())) == 1:
                    self.Rmat = list(self.Rmatrix.Rmatrix.values())[
                        0]["matrix"]
                else:
                    if isATOV:

                        if self.instrument == "amsua":
                            self.Rmat = self.Rmatrix.Rmatrix[
                                satellite + " " + "ATOVS"]["matrix"][
                                20:35, 20:35]
                        else:
                            if (
                                    self.instrument == "amsub" or
                                    self.instrument == "MHS" or
                                    self.instrument == "mhs"):
                                self.Rmat = self.Rmatrix.Rmatrix[
                                    satellite + " " + "ATOVS"]["matrix"][
                                    35:40, 35:40]
                                print(self.Rmat)
                            else:
                                if self.instrument == "hirs":
                                    # HIRS : 1 vis channel, 19 IR channels :
                                    # take only the IR channels !
                                    # in this case we start at 1 because the
                                    # first channel is vis !!!!
                                    self.Rmat = self.Rmatrix.Rmatrix[
                                        satellite + " " + "ATOVS"]["matrix"][
                                        1:20, 1:20]
                                else:
                                    raise ValueError(
                                        "wrong satellite and instrument" +
                                        str(satellite) + " " + str(instrument))
                    else:
                        raise ValueError(
                            "wrong satellite and instrument " +
                            str(satellite) +
                            " " + str(instrument))
            except (KeyError):

                raise ValueError("wrong satellite and instrument " +
                                 str(satellite) + " " + str(instrument))
        # keep only selected channels

        if self.channel_list is None:
            logging.info("no channel list take all of the matrix")
            self.channel_list = list(range(0, self.Rmat.shape[0]))
        else:
            # this channel list already begin with 0
            foo = self.Rmat[self.channel_list, :]
            bar = foo[:, self.channel_list]
            self.Rmat = bar

    def plotRmat(self):
        plt.imshow(self.Rmat, origin="lower", interpolation="nearest")
        plt.title("Rmatrix for " + self.instrument + " " + self.satellite)
        plt.colorbar()
        plt.show()

    def setFactorB(self, factor):
        self.factorB = factor

    def setFactorR(self, factor):
        self.factorR = factor

    def setMaxNoise(self, maxNoise):
        self.MaxNoise = maxNoise

    def stepRetrieve1d(self):
        """ perform a step of retrieval
          the retrieved vector becomes the new background profile"""
        if self.Xr is not None:
            if self.changeTheBackground:
                # make a new profile object with the retrieved vector and the
                # former background profile
                newBgProfile = vector2profile(self.Xr, self.pBg.myProfile)
                self.pBg.myProfile = newBgProfile
                self.pBg.myProfile = self.retrievedProfiles[-1]
                logging.info("new profile background initialised " +
                             self.pBg.profileFileName)
                self.retrieve1d()
            else:
                # we do not change the Background
                self.retrieve1d()

    def set1dvarRequiredOptionsAndValues(self, project):
        """ set the required option and valuesfor the 1dvar algorithm """
        """ this options are CFRACTION=0 """
        """ Nadir condition """
        logging.info(
            "set 1dvar required options and values (CFRACTION=0, no gas,"
            " calculation at nadir etc) ")
        project.myOption["SWITCHRAD"].value = True
        # No gas
        project.myOption["CO2_DATA"].value = False
        project.myOption["CH4_DATA"].value = False
        project.myOption["OZONE_DATA"].value = False
        project.myOption["CO_DATA"].value = False
        project.myOption["N2O_DATA"].value = False
        # no cloud
        project.myOption["ADDCLOUDS"].value = False
        project.myProfile["CFRACTION"] = 0
        # no aerosols
        project.myOption["ADDAEROSL"].value = False
        # no solar
        project.myOption["ADDSOLAR"].value = False

        # calculation at nadir
        project.myProfile["ZENANGLE"] = 0
        project.myProfile["SUNZENANGLE"] = 0

    def assumptionOnSurfaceParameters(self, profileTrue, profileBg):
        """ restrictions the surface parameters of the Bg are applied
             on the True except for T skin and T 2m  """
        profileBg["SKIN"]['SURFTYPE'] = 1  # (0=Land, 1=Sea, 2=sea-ice)
        for item in ["LATITUDE", "LONGITUDE", "ELEVATION", "AZANGLE", "BE",
                     "COSBK"]:
            profileTrue[item] = profileBg[item]
        for item in ["Q", "O", "P", "U", "V", "WFETC"]:
            profileTrue["S2M"][item] = profileBg["S2M"][item]
        for item in ["SURFTYPE", "WATERTYPE", "SALINITY", "FASTEM",
                     "SNOW_FRACTION", "SOIL_MOISTURE"]:
            profileTrue["SKIN"][item] = profileBg["SKIN"][item]

    def retrieve1d(self):
        """ perform the basic retrieval algorithm """
        """ must put switchrad a true et cfraction a 0 """
        # load coefficients
        self.step_counter = self.step_counter + 1
        self.pTrue.myCoeffs = self.pBg.myCoeffs
        self.pTrue.myOption = self.pBg.myOption
        # we are making assumption on surface parameters some must be the same
        self.assumptionOnSurfaceParameters(
            self.pTrue.myProfile, self.pBg.myProfile)

        # at this point the coefficient must be loaded
        if not self.pBg.myCoeffs.loadCoeffs:
            logging.warning(
                "strange : coefficients would have been "
                "loaded at this point ...")
            err = self.pBg.loadCoefficients()
            if err != 0:
                raise RuntimeError("Cannot load coefficient in r1dvar project")

        if len(self.channel_list) != self.pTrue.myCoeffs.nchannels:
            logging.warning("wrong channel list : take all")
            self.channel_list = list(
                range(1, self.pTrue.myCoeffs.nchannels + 1))
            self.readRmatrix(self.satellite, self.instrument)

        nbChannels = self.pTrue.myCoeffs.nchannels
        logging.info("nbChannels: {}".format(nbChannels))

        logging.debug("coefficients file: " +
                      self.pBg.myCoeffs.fileName["standard"])
        if self.step_counter > 1:
            # test if we have a different number of channels
            if len(self.Y_Xt) == 1:
                previousNumberOfChannels = 1
            else:
                previousNumberOfChannels = self.Y_Xt.shape[0]
            notSameNumberOfChannels = (
                self.pBg.myCoeffs.nchannels != previousNumberOfChannels)
            if notSameNumberOfChannels:
                logging.info("not the same number of channels")

        # Run direct on True

        logging.info(">>>>>>>>> Run direct on True")
        # must put switchrad a true et cfraction a 0 and other things on
        # options
        self.set1dvarRequiredOptionsAndValues(self.pTrue)
        err = self.pTrue.runDirect()
        if (err != 0):
            raise RuntimeError("Run Direct failed on true")
        # put radiance to 0 in the radiance file for plotting
        self.changeRadianceFile(self.pTrue.radianceFileName)
        # 1) Y(Xt) brightness temperatures from Xt
        YXt = rttov.radiance.Radiance()
        logging.info(">>>>>>>>>>>reading BT from : ")
        YXt.read(self.pTrue.radianceFileName)
        self.Y_Xt = YXt["BT_CLEAR"]  # No clouds so BT_CLEAR=BT
        Xt = self.pTrue.myProfile.to1DvarVect()
        self.Xt = Xt

        # Run direct on Background
        logging.info(">>>>>>>> Run direct on Background")
        # must put switchrad a true et cfraction a 0 and other things on
        # options
        self.set1dvarRequiredOptionsAndValues(self.pBg)
        err = self.pBg.runDirect()
        if (err != 0):
            raise RuntimeError("Run Direct failed on background")
        # put radiance to 0 in the radiance file for plotting
        self.changeRadianceFile(self.pBg.radianceFileName)

        logging.info(">>>>>>>> Run K on Background")
        err = self.pBg.runK()
        if (err != 0):
            raise RuntimeError("Run K failed on background")
        else:
            logging.info("run K on background OK")

        # 2) Y((Xb) brightness temperarures from Xb
        YXb = rttov.radiance.Radiance()
        logging.info(">>>>>>>>>>reading BT " + self.pBg.radianceFileName)
        YXb.read(self.pBg.radianceFileName)
        Y_Xb = YXb["BT_CLEAR"]

        # 3) compute jacobian matrix K(Xb) and transpose KT(Xb)
        KXb = r1dvarObjects.R1dvarKmatrix()
        logging.info("reading kmatrix " + self.pBg.KMatrixFileName)
        KXb.read(self.pBg.KMatrixFileName)

        logging.info(">>>>>>>KmatrixVectorAndTranspose ")
        # extract the jacobian matrix from the KMatrix Object
        K_Xb = KXb.toKmatrix1dvar()
        # transpose this matrix
        KT_Xb = K_Xb.T

        # 4) apply scaling factor to background errors B --> fb x B
        # the Bmatrix has been read in init and has been put in self.Bmat
        myBmat = self.Bmat * self.factorB

        # 5) apply scaling factor to observation errors R --> fr x R
        # The R matrix has been read in init and has been put in self.Rmat

        myRmat = self.Rmat * self.factorR

        # 6) compute add random noise to Y(Xt)

#
        # Noise in BT = +/- MaxNoise degree
        logging.info("Max Noise : " + str(self.MaxNoise))
        self.N = self.MaxNoise * (np.random.rand(nbChannels) * 2 - 1)
        self.Y_Xt = self.Y_Xt + self.N

        # save the file with True BT+Noise in order to be able to plot it after
        new_name = self.pTrue.radianceFileName.replace(".h5", "_Noise.h5")
        self.changeRadianceFile(self.pTrue.radianceFileName, self.N, new_name)

        # 7) compute linear 1DVAR weights W where
        # W = B . KT . [ K .B . KT + R ](-1)
        # W =  fact1 . [ foo1 .KT  + R ](-1)
        # W =  fact1 .  [  foo2   +  R ](-1)
        # W =  fact1 .  [       foo3   ](-1)
        # W =  fact1 . Matrice_inverse
        fact1 = myBmat.dot(KT_Xb)
        foo1 = K_Xb.dot(myBmat)
        foo2 = foo1.dot(KT_Xb)
        foo3 = foo2 + myRmat
        # matrix inversion
        logging.info("matrix inversion")
        Matrice_inverse = np.linalg.inv(foo3)
        logging.info("compute the linear 1DVAR weight ")
        W = fact1.dot(Matrice_inverse)

        # 8) compute linear 1Dvar retrieved profile (Xr)
        Xb = self.pBg.myProfile.to1DvarVect()
        # we have already added N to Y_Xt earlier
        YtminusYb = self.Y_Xt - Y_Xb

        Xr = Xb + W.dot(YtminusYb)
        # if Q <0 alors Q=0
        for k in range(Xr.shape[0] - 29 - 3, (Xr.shape[0] - 3)):
            if Xr[k] <= 0:
                # if find something negative for Q : stick to the background
                # ....
                Xr[k] = Xb[k]
        k = Xr.shape[0] - 2
        if Xr[k] <= 0:
            # if find something negative for Q : stick to the background ....
            Xr[k] = Xb[k]

        self.Xr = Xr
        self.Y_Xb = Y_Xb
        self.Xb = Xb
        # append the retrieved Vector and the retrieved profile to list
        self.retrievedVectors.append(Xr)
        retrievedProfile = vector2profile(self.Xr, self.pBg.myProfile)
        self.retrievedProfiles.append(retrievedProfile)
        return Xr

    def changeRadianceFile(self, radianceFileName, Noise=None, NewName=None):
        """   change the radiance file set radiance to 0 add Noise to BT) """
        if NewName:
            shutil.copyfile(radianceFileName, NewName)
            rf = h5py.File(radianceFileName, 'r+')
        else:
            rf = h5py.File(radianceFileName, 'r+')
        group = rf["RADIANCE"]
        bt = rf["RADIANCE/BT"][:]
        bt_clear = rf["RADIANCE/BT_CLEAR"][:]
        rad_clear = rf["RADIANCE/CLEAR"]
        a = rad_clear[:]

        a_new = np.zeros_like(a)
        del rf["RADIANCE/CLEAR"]
        del rf["RADIANCE/CLOUDY"]
        del rf["RADIANCE/TOTAL"]
        del rf["RADIANCE/OVERCAST"]
        group.create_dataset("CLEAR", a_new.shape, data=a_new)
        group.create_dataset("CLOUDY", a_new.shape, data=a_new)
        group.create_dataset("TOTAL", a_new.shape, data=a_new)
        group.create_dataset("OVERCAST", a_new.shape, data=a_new)
        if Noise is not None:
            bt_new = bt + Noise
            bt_clear_new = bt_clear + Noise
            del rf["RADIANCE/BT"]
            del rf["RADIANCE/BT_CLEAR"]
            group.create_dataset("BT", bt_new.shape, data=bt_new)
            group.create_dataset("BT_CLEAR", bt_new.shape, data=bt_clear_new)
        rf.close()

    def plot(self):
        colors = {0: "red",
                  1: "orange",
                  2: "yellow",
                  3: "magenta",
                  4: "brown"
                  }
        ax = plt.subplot(2, 2, 1)
        Xb0 = self.BgInitialProfile.to1DvarVect()
        Y = self.pTrue.myProfile["P"]
        pression = Y
        print("pression")
        for k in range(0, pression.shape[0]):
            print(("%.2f" % pression[k]))
        ax.set_ylim((pression[-1] + 30, pression[0]))
        ax.set_yscale("log")
        ax.plot(self.Xt[:54], Y, color="black", label="T true")
        ax.plot(Xb0[:54], Y, color="blue", label="T background")
        for k in range(len(self.retrievedVectors)):
            ax.plot(self.retrievedVectors[k][
                    :54], Y, color="red", label="T retrieved step " + str(k))
        ax.plot(self.Xt[85], Y[-1] + 20, color="black",
                marker="*", label="T skin True")
        ax.plot(Xb0[85], Y[-1] + 20, color="blue",
                marker="*", label="T skin background")
        for k in range(len(self.retrievedVectors)):
            ax.plot(self.retrievedVectors[k][85],
                    Y[-1] + 20, color="red", marker="*",
                    label="T skin retrieved step" + str(k))
        ax.legend(prop={'size': 10}, shadow=True, fancybox=True, loc='best')
        ax = plt.subplot(2, 2, 2)

        pression = Y[:29]
        print((pression.shape))

        ax.set_ylim((pression[-1], pression[0]))
        ax.set_yscale("log")

        ax.plot(self.Xt[54:83], pression, color="black", label="Q true")
        ax.plot(Xb0[54:83], pression, color="blue", label="Q background")
        print("q true")
        print(self.Xt[54:83])
        print("Q initial")
        print(Xb0[54:83])
        for k in range(len(self.retrievedVectors)):
            ax.plot(self.retrievedVectors[k][54:83], pression, color=colors[
                    k], label="Q retrieved step " + str(k))
            print("Q retrieved k", k)
            print(self.retrievedVectors[k][54:83])
        ax.legend(prop={'size': 10}, shadow=True, fancybox=True, loc='best')

        ax = plt.subplot(2, 2, 3)
        if len(self.Y_Xb.shape) == 0:
            nbchannels = 1
        else:
            nbchannels = self.Y_Xb.shape[0]

        channels = np.arange(1, nbchannels + 1)

        print((self.Y_Xb.shape))
        print((channels.shape))
        ax.plot(channels, self.Y_Xt, color="black", label="BT True")
        ax.plot(channels, self.Y_Xb, color="blue", label="BT Background")
        ax.plot(channels, self.Y_Xt + self.N,
                color="yellow", label="BT True+Noise")
        ax.legend(prop={'size': 10}, shadow=True, fancybox=True, loc='best')
        plt.title(self.satellite + " " + self.instrument)
        plt.show()


class Bmatrix(object):
    """ Class for dealing with the background error covariance
        matrices used by the 1DVar sheme """

    def __init__(self):
        self.matrices = {}
        self.header = {}
        self.dimension = {}

    def _read(self, f, dim=86):

        data = np.zeros((dim * dim), dtype=np.float64)
        indice = 0
        while indice < dim * dim:
            lin = f.readline().split()
            for value_str in lin:
                data[indice] = np.float64(value_str)
                indice = indice + 1
        return data.reshape(dim, dim)

    def read_matrices(self, filename):
        f = open(filename)
        for j in (1, 2, 3, 4):
            self.header[j] = []
            for i in (0, 1, 3):
                k = f.readline()
                self.header[j].append(k)
            self.dimension[j] = int(k)
            self.matrices[j] = self._read(f, dim=self.dimension[j])

            # convert lnq g/kg in ppmv for Q
            dim = self.dimension[j]
            for k in range(dim - 29 - 3, dim - 3):
                for l in range(dim - 29 - 3, dim - 3):
                    self.matrices[j][k, l] = q_mixration_to_ppmv * \
                        np.exp(self.matrices[j][k, l]) * 0.001
            self.matrices[j][dim - 2, dim - 2] = q_mixration_to_ppmv * \
                np.exp(self.matrices[j][dim - 2, dim - 2]) * 0.001

        f.close()

    def plot(self, j):
        print("Bmatrice :", j)
        print(self.matrices[j])

        plt.imshow(self.matrices[j], origin="lower", interpolation="nearest")
        plt.title(self.header[j][0])
        plt.colorbar()
        plt.show()


class RmatrixBase(object):
    """ Class for dealing with mesurement error covariance matrix
        used by the 1DVar scheme """

    def __init__(self):
        pass

    def read(self, filename):
        pass


class Rmatrix(RmatrixBase):

    def read(self, filename):
        """ Rmatrix fila structure is explained  here : """
        """ https://nwpsaf.eu/deliverables/nwpsaf_1dvar/
            nwpsaf-mo-ud-032_NWPSAF_1DVar_Manual.html#aux """
        """ read the entire file : populate a dictionary named Rmatrx with :"""
        """ nchannels : number of channels"""
        """ nband : number of band """
        """ inverse : """
        """ type :"""
        """ data : the data"""
        """ finaly with the __mkMatrix method is called to make the Matrix """
        """ dictionary which will contain numpy arrays """

        self.Rmatrix = {}
        logging.debug("open " + filename)
        f = open(filename)
        cond = True
        while cond:
            li = f.readline()
            if len(li) == 0:
                cond = False
                continue
            satint = li.replace("\n", "")
            logging.debug("read" + str(satint))
            self.Rmatrix[satint] = {}
            line = f.readline()
            self.Rmatrix[satint]["type"], self.Rmatrix[satint][
                "nchannel"], self.Rmatrix[
                    satint]["nband"], self.Rmatrix[satint]["inverse"] = [
                        int(x) for x in line.split()]
            logging.debug(str(self.Rmatrix[satint]["type"]) +
                          str(self.Rmatrix[satint]["nchannel"]) +
                          str(self.Rmatrix[satint]["nband"]) +
                          str(self.Rmatrix[satint]["inverse"]))
            if self.Rmatrix[satint]["type"] == 2:
                self.Rmatrix[satint]["channel"] = []
                self.Rmatrix[satint]["data"] = {}
                while (len(self.Rmatrix[satint]["channel"]) <
                       self.Rmatrix[satint]["nchannel"]):
                    li = f.readline().split()
                    for x in li:
                        self.Rmatrix[satint]["channel"].append(int(x))
                for band in range(self.Rmatrix[satint]["nband"]):
                    self.Rmatrix[satint]["data"][band] = []
                    while (len(self.Rmatrix[satint]["data"][band]) <
                           self.Rmatrix[satint]["nchannel"]):
                        li = f.readline().split()
                        for x in li:
                            self.Rmatrix[satint]["data"][band].append(float(x))
                logging.debug("end")
            else:
                logging.debug('cannot read matrix')
                cond = False
        self._mkMatrix()

    def _mkMatrix(self):
        """ create the matrix from the data read """
        logging.debug("create the matrix")
        for satint in list(self.Rmatrix.keys()):
            logging.debug("create the matrix for satellite" + satint +
                          str(self.Rmatrix[satint]["type"]) +
                          str(self.Rmatrix[satint]["nband"]))
            self.Rmatrix[satint]["matrix"] = np.zeros(
                (self.Rmatrix[satint]["nchannel"],
                 self.Rmatrix[satint]["nchannel"]), dtype=float)
            if (
                    self.Rmatrix[satint]["type"] == 2 and
                    self.Rmatrix[satint]["nband"] == 1):
                for i in range(self.Rmatrix[satint]["nchannel"]):
                    self.Rmatrix[satint]["matrix"][
                        i, i] = self.Rmatrix[satint]["data"][0][i]

            if (
                    self.Rmatrix[satint]["type"] == 2 and
                    self.Rmatrix[satint]["nband"] != 1):
                # diagonal case
                for i in range(self.Rmatrix[satint]["nchannel"]):
                    self.Rmatrix[satint]["matrix"][
                        i, i] = self.Rmatrix[satint]["data"][0][i]

                for band in range(1, self.Rmatrix[satint]["nband"]):
                    logging.debug("band" + str(band))

                    for i in range(0, self.Rmatrix[satint]["nchannel"] - band):
                        self.Rmatrix[satint]["matrix"][
                            i + band, i] = self.Rmatrix[satint]["data"][band][
                            i]
                        self.Rmatrix[satint]["matrix"][
                            i, i + band] = self.Rmatrix[satint]["data"][band][
                            i]

    def plot_matrix(self):
        for satint in list(self.Rmatrix.keys()):
            plt.imshow(self.Rmatrix[satint]["matrix"],
                       origin="upper", interpolation="nearest")
            plt.title(satint)
            plt.colorbar()
            plt.show()


class Nwpsaf1dvarProfile(object):

    def __init__(self):
        self.data = {}

    def toRttovGuiProfile(self):
        """ convert Background to profile """

        self.profile = rmodel.project.pProfile()
        self.option = rmodel.project.pOption()
        self.profile.setDefaultAttributes()
        self.profile["NLEVELS"] = self.data["P"][::-1].shape[0]
        self.profile["LAYERS"] = self.profile["NLEVELS"] - 1
        self.profile["P"] = self.data["P"][::-1]
        self.profile["T"] = self.data["T"][::-1]
        self.profile["Q"] = self.data["Q"][::-1]
        self.profile["O3"] = self.data["O3"][::-1]
        self.profile["S2M"]["T"] = self.data["Surface Temperature (K)"]
        self.profile["SKIN"]["T"] = self.data["Skin Temperature (K)"]
        for k, v in list(self.data.items()):
            if k.startswith("Surface Humidity"):
                self.profile["S2M"]["Q"] = v
        self.profile["S2M"]["U"] = self.data["10m U-Wind (m/s)"]
        self.profile["S2M"]["V"] = self.data["10m U-Wind (m/s)"]
        self.profile["S2M"]["P"] = self.data["Surface Pressure (hPa)"]
        self.profile.setDefaultProfileAsciiInput()
        return self.profile


class RetrievedProfile(object):
    """ class for dealing with 1D var ascii Retrieved_Profile files
        from the NWPDAF 1DVAR software """
    """ read them and transform it in profile object """

    def __init__(self):
        self.retrieved = Nwpsaf1dvarProfile()
        self.background = Nwpsaf1dvarProfile()
        self.header = {}

    def read(self, filename):
        data = np.genfromtxt(filename, skip_header=3, skip_footer=8)
        self.retrieved.data["P"] = data[:, 0]
        self.retrieved.data["T"] = data[:, 1]
        self.retrieved.data["Q"] = data[:, 2]
        self.retrieved.data["O3"] = data[:, 3]
        self.background.data["P"] = data[:, 0]
        self.background.data["T"] = data[:, 4]
        self.background.data["Q"] = data[:, 5]
        self.background.data["O3"] = data[:, 5]

    def plot(self):
        colors = {0: "red",
                  1: "orange",
                  2: "yellow",
                  3: "magenta",
                  4: "brown"
                  }
        ax = plt.subplot(2, 2, 1)
        Xb = self.background
        Xr = self.retrieved

        Y = self.background.data["P"][:]
        pression = Y
        print("pression")
        for k in range(0, pression.shape[0]):
            print(("%.2f" % pression[k]))
        ax.set_ylim((pression[0] + 30, pression[-1]))
        ax.set_yscale("log")
        ax.plot(Xb.data["T"], Y, color="blue", label="T background")
        ax.plot(Xr.data["T"], Y, color="red", label="T retrieved")

        ax.legend(prop={'size': 10}, shadow=True, fancybox=True, loc='best')
        ax = plt.subplot(2, 2, 2)
        ax.set_ylim((pression[0] + 30, pression[-1]))

        ax.set_yscale("log")

        ax.plot(Xb.data["Q"], pression, color="blue", label="Q background")

        ax.plot(Xr.data["Q"], pression, color="red", label="Q retrieved  ")

        ax.legend(prop={'size': 10}, shadow=True, fancybox=True, loc='best')

        plt.show()


class Background(object):
    """ class for dealing with 1D var ascii Backround files from the
        NWPDAF 1DVAR software """
    """ read them and transform it in profile object """

    def __init__(self):
        self.prof = Nwpsaf1dvarProfile()
        self.data = {}
        self.header = {}

    def read(self, filename):
        data = np.genfromtxt(filename, skip_header=16, skip_footer=6)

        self.data["P"] = data[:, 0]
        self.data["T"] = data[:, 1]
        self.data["Q"] = data[:, 2]
        self.data["O3"] = data[:, 3]
        self.nlevels = data.shape[0]

        labels = np.genfromtxt(
            filename, usecols=0, dtype=str, skip_header=16 + data.shape[0],
            delimiter=":")
        raw_data = np.genfromtxt(
            filename, skip_header=16 + data.shape[0], delimiter=":")[:, 1:]

        self.data_surf = {label: row for label,
                          row in zip(labels, raw_data[:, 0])}
        for key in ("P", "Q", "T", "O3"):
            self.prof.data[key] = self.data[key]
        self.prof.data["Surface Temperature (K)"] = self.data_surf[
            "Surface Temperature (K)"]
        self.prof.data["Skin Temperature (K)"] = self.data_surf[
            "Skin Temperature (K)"]
        self.prof.data["10m U-Wind (m/s)"] = self.data_surf["10m U-Wind (m/s)"]
        self.prof.data["10m U-Wind (m/s)"] = self.data_surf["10m U-Wind (m/s)"]
        self.prof.data["Surface Pressure (hPa)"] = self.data_surf[
            "Surface Pressure (hPa)"]

    def toProfile(self):
        return self.prof.toRttovGuiProfile()

    def print_data(self):
        print("P")
        print((self.data["P"]))
        print("T")
        print((self.data["T"]))
        print("Q")
        print((self.data["Q"]))
        print("O3")
        print((self.data["P"]))
        print((self.data_surf))


if __name__ == '__main__':

    rmatrixfilename = "./data/IASI_COEFFS_DIR/Rmatrix_orig"

    print("rmatrixfilename")
    R = Rmatrix()
    R.read(rmatrixfilename)
    R.plot_matrix()
    print("End ok")
