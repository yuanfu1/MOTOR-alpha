# -*- coding: utf-8 -*-

# wxPython 2.8
import sys
try:
    import wx
    if wx.VERSION[1] == 8:
        sys.stderr.write(
            ' wx python version 2.8  some features of RTTOVGUI might'
            ' not work properly\n')

except ImportError:
    sys.stderr.write(
        'ERROR: wx is not installed or not compatible with python \n')
    sys.exit(1)

try:
    import datetime
except ImportError:
    sys.stderr.write('ERROR: datetime is not installed\n')
    sys.exit(1)

try:
    import h5py
except ImportError:
    sys.stderr.write('ERROR: h5py  is not installed\n')
    sys.exit(1)

import rttov
from rmodel import config
import os
import rttov_gui_f2py
import logging
import numpy
import signal
from rmodel import karchive
import copy

if wx.VERSION[1] == 8:
    from wx.lib.pubsub import Publisher as pub
else:
    if wx.VERSION[0] <= 3:
        # wxPython 2.9 or 3.0:
        # wx.lib.pubsub: Pusub now defaults to the new "kwarg" version of
        # the API.
        # In order to continue using the original "arg1" API you
        # will need to import wx.lib.pubsub.setuparg1 before importing any
        # other
        # pubsub modules.
        import wx.lib.pubsub.setuparg1
        from wx.lib.pubsub import pub
    else:
        # wxPython 4
        from pubsub import pub


class Coefficients(object):

    def __init__(self):
        self.loadCoeffs = False
        self.channel_list = None
        self.WAVENUMBERS = None
        self.fileName = {"standard": "",
                         "aerosols": "", "clouds": "", "mfasis": "", "PC": ""}
        self.nchannels = 0
        self.Filenchannels = 0

        # gas_id: cf rttov_const.F90
#        & gas_id_mixed       = 1, &
#        & gas_id_watervapour = 2, &
#        & gas_id_ozone       = 3, &
#        & gas_id_wvcont      = 4, &
#        & gas_id_co2         = 5, &
#        & gas_id_n2o         = 6, &
#        & gas_id_co          = 7, &
#        & gas_id_ch4         = 8  &
#        & gas_id_so2         = 9
        self.gas_id = {'Q': 1, 'O3': 3, 'CO2': 5,
                       "N2O": 6, "CO": 7, "CH4": 8, "SO2": 9}

    def get_FMV_MODEL_VER(self):
        """ return the version of predictors """
        version = None
        if self.loadCoeffs:
            version = rttov.getcoefval.rttov_get_coef_val_i0("FMV_MODEL_VER")
        return version

    def is_V13(self):
        """ return True if the coeffile is V13 """
        version = self.get_FMV_MODEL_VER()
        if version == 13:
            return True
        else:
            return False

    def isMW(self):
        """ return True if the sensor is a MicroWave sensor """
        if self.loadCoeffs:
            ID_SENSOR = rttov.getcoefval.rttov_get_coef_val_i0("ID_SENSOR")
            if ID_SENSOR == 2 or ID_SENSOR == 4:
                return True
            else:
                return False
        else:
            return False

    def hasSolar(self):
        if self.loadCoeffs:
            SS_VAL_CHN = rttov.getcoefval.rttov_get_coef_val_i1("SS_VAL_CHN")
            if SS_VAL_CHN.sum() > 0:
                return True
            else:
                return False
        else:
            return False

    def hasGas(self, gas):
        posgaslist = self.getFMV_GAS_POS()
        if posgaslist is not None:
            if len(posgaslist) < self.gas_id[gas]:
                return False
            if posgaslist[self.gas_id[gas] - 1] <= 0:
                return False
            else:
                return True
        else:
            return False

    def getFMV_GAS_POS(self):
        if self.loadCoeffs:
            return rttov.getcoefval.rttov_get_coef_val_i1('FMV_GAS_POS')
        else:
            return None

    def hasPC(self):
        if (self.loadCoeffs) and self.fileName["PC"] != "":
            return True
        else:
            return False

    def hasClouds(self):
        if (self.loadCoeffs) and self.fileName["clouds"] != "":
            return True
        else:
            return False

    def hasMfasis(self):
        if (self.loadCoeffs) and self.fileName["mfasis"] != "":
            return True
        else:
            return False

    def hasAerosols(self):
        if (self.loadCoeffs) and self.fileName["aerosols"] != "":
            return True
        else:
            return False

    def isPCClouds(self):
        if self.loadCoeffs and self.fileName["PC"] != "":
            cloud = rttov.getcoefval.rttov_get_coef_val_i0('FMV_PC_CLD')
            if cloud != 0 and self.hasClouds():
                return True
            else:
                return False
        else:
            return False

    def getFMV_PC_BANDS(self):
        if self.loadCoeffs:
            if self.hasPC():
                return rttov.getcoefval.rttov_get_coef_val_i0('FMV_PC_BANDS')
        else:
            return None

    def getMaxIPCREG(self, band):
        if self.loadCoeffs and self.hasPC():
            FMV_PC_SET = self.getFMV_PC_SETS()
            if band in range(1, len(FMV_PC_SET) + 1):
                maxIPCREG = FMV_PC_SET[band - 1]
            else:
                maxIPCREG = 1
        else:
            maxIPCREG = 1
        return maxIPCREG

    def getFMV_PC_SETS(self):
        if self.loadCoeffs:
            if self.hasPC():
                return rttov.getcoefval.rttov_get_coef_val_i1('FMV_PC_SETS')
        else:
            return None

    def getFF_ORI_CHN(self):
        """ return FF_ORI_CHN: the list of selected channels """
        if self.loadCoeffs:
            mylist = rttov.getcoefval.rttov_get_coef_val_i1('FF_ORI_CHN')
            return mylist
        else:
            return None

    def isHighResolution(self):
        """ return True if the sensor is a Hyper spectral sensor """
        if self.loadCoeffs:
            if rttov.getcoefval.rttov_get_coef_val_i0('ID_SENSOR') == 3:
                return True
            else:
                return False
        else:
            return False

    def hasNLTE(self):
        """ return True if the sensor can do an NLTE correction """

        if self.loadCoeffs:
            nlte = rttov.getcoefval.rttov_get_coef_val_i0("NLTECOEF")
            if nlte == 1:
                return True
            else:
                return False
        else:
            return False

    def hasInTheName(self, aString, typeCoefffile):
        """ return True if coefficient file as the a string in the name """
        if self.loadCoeffs:
            if self.fileName[typeCoefffile] != "":
                if aString in self.fileName[typeCoefffile]:
                    return True
        return False

    def get_instrument_and_platform_name(self):
        """ return instrument and platform name """
        if self.loadCoeffs:
            id_platform = rttov.getcoefval.rttov_get_coef_val_i0('ID_PLATFORM')
            id_inst = rttov.getcoefval.rttov_get_coef_val_i0('ID_INST')
            id_sat = rttov.getcoefval.rttov_get_coef_val_i0('ID_SAT')
            try:
                platform_name = rttov.getcoefval.platform_name[id_platform - 1]
                inst_name = rttov.getcoefval.instrument_name[id_inst]
                logging.debug("platform_name" +
                              platform_name + str(id_platform))
                logging.debug("inst_name" + inst_name + str(id_inst))
            except Exception:
                logging.warning("cannot recognize instrument and platform")
                inst_name = "unknown"
                platform_name = "unknown"
            return (inst_name, platform_name, id_sat)
        else:
            return (None, None, None)


class pProfile(rttov.profile.Profile):
    """ project Profile Class: add specific functionalities to the
        rttov Profile Class"""

    def __init__(self):
        rttov.profile.Profile.__init__(self)

    def diff(self, profile1):
        """ make difference my.Profile - profile and return the difference """
        diffP = pProfile()
        item = "T"
        diffP[item] = self[item] - profile1[item]
        for item in self.gas_list:
            if (self[item] is not None and profile1[item] is not None):
                diffP[item] = self[item] - profile1[item]
        diffP["P"] = self["P"]

        return diffP

    def ctrlCoherenceClouds(self):
        """ control the cloud coherence """
        cloudCumulDensity = numpy.zeros(self['NLAYERS'], dtype=self['P'].dtype)
        if self.anyCloud():
            print("ctrlCoherenceClouds", self.cloud_list)
            for cloud in self.cloud_list:
                if self[cloud] is not None:
                    cloudCumulDensity = cloudCumulDensity + self[cloud]

            for layer in range(self['NLAYERS']):
                if cloudCumulDensity[layer] == 0:
                    self['CFRAC'][layer] = 0.
                else:
                    if self['CFRAC'][layer] == 0.:
                        self['CFRAC'][layer] = 0.20
                    if self['CFRAC'][layer] >= 1:
                        self['CFRAC'][layer] = 1

    def hasGas(self, gas):
        if gas in self:
            if self[gas] is None:
                return False
            else:
                return True
        else:
            return False

    def hasClouds(self):
        hasclouds = False
        for cloud in self.cloud_list:
            if cloud in self:
                if self[cloud] is not None:
                    hasclouds = True
        return hasclouds

    def hasAerosols(self):
        hasAer = False
        for aer in self.aerosol_list:
            if aer in self:
                if self[aer] is not None:
                    hasAer = True
        return hasAer


class pOption(rttov.option.Option):

    def __init__(self):
        rttov.option.Option.__init__(self)

    def updatePCChoicesList(self, myCoeff):
        """ update the PCChoicesList according to coefficients files"""
        """ FMV_PC_BAND = return Bands availables in the coefficient file """
        """ FMV_PC_SET = return max ipcreg  for each band """
        FMV_PC_BAND = myCoeff.getFMV_PC_BANDS()
        FMV_PC_SET = myCoeff.getFMV_PC_SETS()
        if self["IPCBND"].value > FMV_PC_BAND:
            self["IPCBND"].value = FMV_PC_BAND

        logging.debug("FMV_PC_BAND=" + str(FMV_PC_BAND))
        logging.debug("FMV_PC_SET=" + str(FMV_PC_SET))
        self["IPCBND"].odict = {x: self.fullIPCBND[x] for x in range(
            1, FMV_PC_BAND + 1)}
        maxIPCREG = myCoeff.getMaxIPCREG(self["IPCBND"].value)
        logging.debug("maxIPCREG " + str(maxIPCREG))

        self["IPCREG"].odict = {x: str(x) for x in range(1, maxIPCREG + 1)}
        if self["IPCREG"].value > maxIPCREG:
            self["IPCREG"].value = maxIPCREG


class Project(object):
    '''
    The Project class makes the interface with rttov_hdf_mod and rttov_guy_f2py
    attributes:
    myCoeffs: to store the name of the coefficients files and to know if
    they are loaded or not
    myProfile: a rttov_hdf_mod Profile object
    myOption:  a rttov_hdf_mod Option object
    config: configuration (path ans so on)

    '''
    # take care number of threads should be 1 if the instrument contains
    # both pure solar and pure thermal channels
    # if the add_solar option is false a thread can have only pure solar
    # channels so this thread fails with message
    # "No thermal or solar calculations to carry out"
    nthreads = 1

    def __init__(self, karchivedepth=4):
        '''
        Constructor
        '''
        self.debug = False
        self.myCoeffs = Coefficients()
        self.myProfile = pProfile()
        self.myOption = pOption()
        self.myEmissivity = rttov.emissivity.Emissivity()
        self.myReflectance = rttov.reflectance.Reflectance()
        self.loadProfile = False
        self.useAtlas = False
        self.atlasLon = 0
        self.atlasLat = 0
        self.firstTimeSurface = True
        self.config = config.Config()
        # name files configured by user
        self.savedProfileFileName = self.config.profileDefaultFileName
        self.savedProjectFileName = None
        self.savedSurfaceFileName = None
        # default fileNames
        self.profileFileName = self.config.profileDefaultFileName
        self.radianceFileName = self.config.radianceDefaultFileName
        self.KMatrixFileName = self.config.KMatrixDefaultFileName
        self.transFileName = self.config.transDefaultFileName
        self.surfaceFileName = self.config.surfaceDefaultFileName
        self.pcFileName = self.config.pcDefaultFileName
        self.pcKMatrixFileName = self.config.pcKMatrixFileName

        self.RTTOV_MODE = "DIRECT"
        self.RTTOV_MODE_K = "K"
        # create logger
        self.myLogger = logging.getLogger('rttov gui')
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        self.ListeParamAtlas = ["LATITUDE", "LONGITUDE", "DATE", "AZANGLE",
                                "ZENANGLE", "SUNAZANGLE", "SUNZENANGLE"]
        self.paramAtlas = {}
        self.tmpFileOut = self.config.tmpFileOut
        self.tmpFileErr = self.config.tmpFileErr
        self.npcscores = 200
        self.previousRunIsPc = False
        signal.signal(signal.SIGSEGV, self.sig_handler)
        signal.signal(signal.SIGFPE, self.sig_handler)
        signal.signal(signal.SIGTRAP, self.sig_handler)
        self.aKmats = karchive.Karchive(karchivedepth)

    def sendMessage(self, mess1, mess2=None):
        if wx.VERSION[0] <= 3:
            pub.sendMessage(mess1)
        else:
            pub.sendMessage(mess1, msg=mess2)

    def setFileNameMark(self, mark):
        """ put a mark on filenames for the project """
        """ if the mark is present do nothing """
        """ TODO: create a list for the Names """
        self.profileFileName = self.profileFileName.replace(
            ".h5", mark + ".h5")
        self.profileFileName = self.profileFileName.replace(mark + mark, mark)
        self.config.profileDefaultFileName = \
            self.config.profileDefaultFileName.replace(".h5", mark + ".h5")
        self.config.profileDefaultFileName = \
            self.config.profileDefaultFileName.replace(mark + mark, mark)
        self.radianceFileName = self.radianceFileName.replace(".h5",
                                                              mark + ".h5")
        self.radianceFileName = self.radianceFileName.replace(
            mark + mark, mark)
        self.KMatrixFileName = self.KMatrixFileName.replace(
            ".h5", mark + ".h5")
        self.KMatrixFileName = self.KMatrixFileName.replace(mark + mark, mark)
        self.transFileName = self.transFileName.replace(".h5", mark + ".h5")
        self.transFileName = self.transFileName.replace(mark + mark, mark)
        self.surfaceFileName = self.surfaceFileName.replace(
            ".h5", mark + ".h5")
        self.surfaceFileName = self.surfaceFileName.replace(mark + mark, mark)
        self.pcFileName = self.pcFileName.replace(".h5", mark + ".h5")
        self.pcFileName = self.pcFileName.replace(mark + mark, mark)
        self.pcKMatrixFileName = self.pcKMatrixFileName.replace(".h5",
                                                                mark + ".h5")
        self.pcKMatrixFileName = self.pcKMatrixFileName.replace(mark + mark,
                                                                mark)

    def sig_handler(self, signum, frame):
        # deal with the signal..
        msg = ("ERROR Something wrong occurred with RTTOV:" +
               " see tmpFileErr.log log file and report the bug\n")
        sys.stderr.write(msg)
        logging.error(msg)

    def write(self, msg):
        logging.info(msg)
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), msg)

    def beforeCallFortran(self):
        # redirect stdin and stdout for fortran code
        # open 2 fds
        self.outFile = os.open(
            self.tmpFileOut, os.O_RDWR | os.O_CREAT | os.O_TRUNC)
        self.errFile = os.open(
            self.tmpFileErr, os.O_RDWR | os.O_CREAT | os.O_TRUNC)
        # save the current file stderr and stdout descriptors to a tuple
        self.save_fds = os.dup(1), os.dup(2)
        # put my files fds on 1 and 2O3_DATA      =>    False <type 'bool'>
        os.dup2(self.outFile, 1)
        os.dup2(self.errFile, 2)

    def afterCallFortran(self):
        # restore standard out and err
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # close the temporary fds
        os.close(self.outFile)
        os.close(self.errFile)
        for fileName in (self.tmpFileOut, self.tmpFileErr):
            f = open(fileName, "r")
            lines = f.readlines()
            for line in lines:
                print(line)
                logging.info(str(line))
            f.close()
        os.close(self.save_fds[0])
        os.close(self.save_fds[1])

    def openAsciiProfile(self, filename):
        p = pProfile()
        erreur = p.loadProfileAscii(filename)
        if(erreur != 1):
            self.myProfile = p
            self.myProfile.controlProfile()
            self.loadProfile = True
            self.ctrlCoherence()
            self.sendMessage("Profile CHANGED")
            self.sendMessage("Options CHANGED")
            self.sendMessage("Profile LOADED")
            logging.debug("project send message Profile LOADED")
            self.profileFileName = self.config.profileDefaultFileName
            self.profileInitialFilename = (filename, 1)
            if self.myProfile["GAS_UNITS"] in [-1, 0, 1, 2]:
                self.myProfile.my_gas_units = self.myProfile.gas_units[
                    self.myProfile["GAS_UNITS"]]
            self.write("CFRACTION=" + str(self.myProfile["CFRACTION"]))
            if self.myProfile["CLW_SCHEME"] is not None:
                self.write("CLW_SCHEME=" + str(self.myProfile["CLW_SCHEME"]))
            return 0
        else:
            print((filename, "cannot read profile wrong format"))
            return 1

    def openProfile(self, filename, number=1):
        '''
        Open a Profile from a hdf5 filename
        '''
        self.write("open profile " + str(filename) + " " + str(number))
        # re-initialize an empty profile
        self.myProfile = pProfile()

        # read the profile
        try:
            self.myProfile.loadProfileNumber(filename, number)
        except Exception:
            self.write("cannot load profile")
            return 1
        # control the values of the Profile regarding min and max values
        self.myProfile.controlProfile()
        # read the options
        try:
            f = h5py.File(filename, "r")
            if 'OPTIONS' in list(f.keys()):
                h5 = f['/OPTIONS/']
                self.myOption.loadh5(h5)
            f.close()
        except Exception:
            return 1
        self.ctrlCoherence()
        self.write("CFRACTION=" + str(self.myProfile["CFRACTION"]))
        if self.myProfile["CLW_SCHEME"] is not None:
            self.write("CLW_SCHEME=" + str(self.myProfile["CLW_SCHEME"]))
        self.loadProfile = True
        self.profileFileName = self.config.profileDefaultFileName
        self.profileInitialFilename = (filename, number)
        if self.myProfile["GAS_UNITS"] in [-1, 0, 1, 2]:
            self.myProfile.my_gas_units = self.myProfile.gas_units[
                self.myProfile["GAS_UNITS"]]
        self.sendMessage("Profile CHANGED")
        self.sendMessage("Options CHANGED")
        self.sendMessage("Profile LOADED")
        if self.useAtlas and self.myCoeffs.loadCoeffs:
            self.write("An atlas was used reinit Emissivity and Reflectance")
            self.reinitSurface()
            self.useAtlas = False
        logging.debug("project send message Profile LOADED")
        return 0

    def saveProfile(self, filename=""):
        """ save Profile and Option in filename """
        if (filename == ""):
            filename = self.profileFileName
        self.write("\nSave profile to file: " + os.path.basename(filename))
        self.profileFileName = filename
        myf = h5py.File(filename, 'w')
        self.myProfile.saveh5(myf, "/PROFILES/0001/")
        self.myOption.saveh5(myf, "/OPTIONS/")
        myf.close()
        del myf

    def updateProfile(self, aProfile, onlySurfaceParameters=False,
                      onlyVerticalParameters=False):
        """ method to update self.myProfile sends a message "Profile CHANGED"
            it is preferable to use this method to update the profile """
        logging.debug("project updateProfile")
        if onlySurfaceParameters:
            logging.info("update profile with surface parameters")
            print(self.myProfile.profSurfParamList)
            for sparam in self.myProfile.profSurfParamList:
                self.myProfile[sparam] = aProfile[sparam]
            for sparam in self.myProfile.sskin_list:
                self.myProfile['SKIN'][sparam] = aProfile['SKIN'][sparam]
            for sparam in self.myProfile.s2m_list:
                self.myProfile['S2M'][sparam] = aProfile['S2M'][sparam]
            self.sendMessage("Surface parameters CHANGED")
        elif onlyVerticalParameters:
            logging.info("update profile with vertical contents")
            self.myProfile["T"] = copy.deepcopy(aProfile["T"])
            for item in (["T"] +
                         self.myProfile.gas_list +
                         self.myProfile.cloud_list +
                         self.myProfile.aerosol_list +
                         ["CFRAC", "ICEDE", "CLWDE"]):
                self.myProfile[item] = copy.deepcopy(aProfile[item])
                self.myProfile[item + "_ATTRIBUTE"] = copy.deepcopy(
                    aProfile[item + "_ATTRIBUTE"])
            self.ctrlCoherence()
            self.sendMessage("Profile CHANGED")
            self.sendMessage("Options CHANGED")

        else:
            logging.info("WARNING profile not updated")

    def updateOptions(self, aOption):
        """ method to update self.myOption sends a message "Options CHANGED"
            it is preferable to use this method to update the options """
        self.myOption = aOption
        self.sendMessage("Options CHANGED")

    def updateEmissivity(self, aEmissivity):
        self.myEmissivity = aEmissivity
        logging.debug("sendMessage Surface CHANGED")
        self.sendMessage("Surface CHANGED")

    def updateReflectance(self, aReflectance):
        self.myReflectance = aReflectance
        logging.debug("sendMessage Surface CHANGED")
        self.sendMessage("Surface CHANGED")

    def openSurface(self, filename, ctrl=True):
        """ open an surface file and load it to the myEmissivity and
            myReflectance objects
            must check if the dimension of EMIS and REFL of
            the file are coherent with number
            of channel loaded """
        self.write("openSurface " + filename)
        f = h5py.File(filename, 'r')
        h5 = f['/EMISSIVITY/']
        emissivity = rttov.emissivity.Emissivity()
        emissivity.loadh5(h5)

        if len(emissivity['EMIS_IN'].shape) == 0:
            emisSize = 1
        else:
            emisSize = emissivity['EMIS_IN'].shape[0]

        # check and correct the data shape
        for key, value in list(emissivity.items()):
            if ("ATTRIBUTE" not in key and len(value.shape) == 0):
                logging.debug("convert value " + str(type(value)) +
                              " " + key)
                foo = value
                emissivity[key] = numpy.array([foo])

        logging.debug("nb channels = " + str(self.myCoeffs.nchannels) +
                      " try to load a file with: " + str(emisSize) +
                      " channels")
        if (ctrl):
            if emisSize != self.myCoeffs.nchannels:
                self.write(
                    "incoherent number of channels: can't load surface file" +
                    filename)
                f.close()
                return 1

        self.myEmissivity = emissivity

        if 'REFLECTANCE' in list(f.keys()):
            logging.debug("REFLECTANCE in surface")
            h5 = f['/REFLECTANCE/']
            self.myReflectance.loadh5(h5)

            # check and correct the data shape
            for key, value in list(self.myReflectance.items()):
                if ("ATTRIBUTE" not in key and len(value.shape) == 0):
                    logging.debug("convert value " + str(type(value)) +
                                  " " + key)
                    foo = value
                    self.myReflectance[key] = numpy.array([foo])
        else:
            # create myReflectance with all zero
            self.myReflectance = rttov.reflectance.Reflectance()
            self.myReflectance.setReflectance(emisSize)
        f.close()
        logging.debug("control Surface ")
        self.controleSurface()
        logging.debug("openSurface send message Surface CHANGED")
        self.sendMessage("Surface CHANGED")

    def saveSurface(self, filename):
        """ save an surface file """
        if (filename == ""):
            filename = self.surfaceFileName
        logging.debug("saveSurface filename" + str(filename))
        f = h5py.File(filename, 'w')
        self.myEmissivity.saveh5(f, '/EMISSIVITY/')
        self.myReflectance.saveh5(f, '/REFLECTANCE/')
        self.write("file surface saved" + filename)
        f.close()

    def reinitSurface(self):
        """ reinit the surface surface file with everything
            to 0 and cal to true """
        if not self.myCoeffs.loadCoeffs:
            sys.stderr.write(
                "Cannot reinit surface coefficient file not loaded")
            self.write("Cannot reinit surface coefficient file not loaded")
            return
        logging.debug("nchannels: " + str(self.myCoeffs.nchannels))
        if self.myEmissivity["EMIS_IN"] is not None:
            self.myEmissivity["EMIS_IN"][:] = 0
            self.myEmissivity["EMIS_OUT"][:] = 0
            self.myEmissivity["CALCEMIS"][:] = True
        if self.myReflectance["REFL_IN"] is not None:
            self.myReflectance["REFL_IN"][:] = 0
            self.myReflectance["REFL_OUT"][:] = 0
            self.myReflectance["CALCREFL"][:] = True
            self.saveSurface(self.surfaceFileName)
            logging.debug("project.reinitSurface send message Surface Changed")
            self.sendMessage("Surface CHANGED")

    def createSurface(self):
        """ create a surface file with everything to 0 and cal to true """
        if not self.myCoeffs.loadCoeffs:
            sys.stderr.write(
                "Cannot create surface coefficient file not loaded")
            self.write("Cannot create surface coefficient file not loaded")
            return
        logging.debug("nchannels: " + str(self.myCoeffs.nchannels))
        self.myEmissivity = rttov.emissivity.Emissivity()
        self.myReflectance = rttov.reflectance.Reflectance()
        self.myEmissivity.setEmissivity(self.myCoeffs.nchannels)
        logging.debug('made new emissivity')
        self.myReflectance.setReflectance(self.myCoeffs.nchannels)
        logging.debug("made new reflectance")
        self.saveSurface(self.surfaceFileName)
        self.useAtlas = False
        logging.debug("project.createSurface send message Surface Changed")
        self.sendMessage("Surface CHANGED")

    def controleSurface(self):
        """ control the surface """
        msgtooBig = ("WARNING: there are values > 1: " +
                     "force them to 1")
        msgtooLittle = ("WARNING: there are values < 0: " +
                        "force them to 0")
        return  # TODO: to do or not to do this control ???
        if self.myCoeffs.loadCoeffs and self.loadProfile:
            if self.myEmissivity["EMIS_IN"] is not None:
                self.write("control Emissivity")
                if (numpy.any(self.myEmissivity["EMIS_IN"] > 1)):
                    tooBig = (self.myEmissivity["EMIS_IN"] > 1)
                    self.myEmissivity["EMIS_IN"][tooBig] = 1
                    self.write(msgtooBig)
                if (numpy.any(self.myEmissivity["EMIS_IN"] < 0)):
                    tooLittle = (self.myEmissivity["EMIS_IN"] < 0)
                    self.myEmissivity["EMIS_IN"][tooLittle] = 0
                    self.write(msgtooLittle)
                if (numpy.all(self.myEmissivity["EMIS_IN"] == 0) and
                        numpy.all(self.myEmissivity["CALCEMIS"] is False)):
                    self.write("WARNING: all emisivity values are 0")
                    self.myEmissivity["CALCEMIS"][:] = True

            if (self.myReflectance["REFL_IN"] is not None and
                    self.myCoeffs.hasSolar()):
                self.write("control BRDF")
                if numpy.any(self.myReflectance["REFL_IN"] > 1):
                    tooBig = (self.myReflectance["REFL_IN"] > 1)
                    self.write(msgtooBig)
                    self.myReflectance["REFL_IN"][tooBig] = 1
                if numpy.any(self.myReflectance["REFL_IN"] < 0):
                    self.write(msgtooLittle)
                    tooLittle = (self.myReflectance["REFL_IN"] < 0)
                    self.myReflectance["REFL_IN"][tooLittle] = 0
                if (numpy.all(self.myReflectance["REFL_IN"] == 0) and
                        numpy.all(self.myReflectance["CALCREFL"] is False)):
                    self.write("WARNING: all BRDF values are 0")
                    self.myReflectance["CALCREFL"][:] = True

    def newSurface(self):
        """ delete surface and re-initialise it """
        self.myEmissivity["EMIS_IN"] = None
        self.myReflectance['REFL_IN'] = None
        self.createSurface()

    def loadCoefficients(self, channel_list=None, redirectOutput=True):
        """ on and load a coefficient file coefficient files
         must have been chosen"""
        if not channel_list:
            channel_list = numpy.array([0], dtype=int)
            all_channels = True
            self.write("loadCoefficient with no channel list selection")
        else:
            self.write("loadCoefficient with channel list selection")
            all_channels = False
        self.myCoeffs.WAVENUMBERS = None
        self.myCoeffs.channel_list = channel_list
        previousNbChannels = self.myCoeffs.nchannels
        if redirectOutput:
            self.beforeCallFortran()
        for k in list(self.myCoeffs.fileName.keys()):
            if self.myCoeffs.fileName[k] != "":
                print(self.myCoeffs.fileName[k])
                if not os.path.isfile(self.myCoeffs.fileName[k]):
                    self.write("file:" + self.myCoeffs.fileName[k] +
                               " does not exist ! ")
                    return 1
        self.myCoeffs.nchannels, err = rttov_gui_f2py.rttov_gui_load(
            channel_list,
            self.myCoeffs.fileName["standard"],
            self.myCoeffs.fileName["aerosols"],
            self.myCoeffs.fileName["clouds"],
            self.myCoeffs.fileName["mfasis"],
            self.myCoeffs.fileName["PC"])
        if redirectOutput:
            self.afterCallFortran()
        if err == 0:
            self.myCoeffs.loadCoeffs = True
            self.myCoeffs.channel_list = channel_list
            n = rttov.getcoefval.rttov_get_coef_val_r1('WAVENUMBERS')
            self.myCoeffs.WAVENUMBERS = n
            if len(n) < 20:
                self.write("WAVENUMBERS " + str(n))
            self.write("rttov_gui_load ok nb channels=" +
                       str(self.myCoeffs.nchannels))
            # if we want all the channel in the file we can say how many
            # channels are in the file
            # and we must keep this information for later
            if all_channels:
                self.myCoeffs.Filenchannels = self.myCoeffs.nchannels
            if self.myCoeffs.nchannels != previousNbChannels:
                self.newSurface()
                self.firstTimeSurface = False
                logging.info("create new surface therefore "
                             " firstTimeSurface is False")
            self.ctrlCoherence()
            logging.info("send message options CHANGED")
            self.sendMessage("Options CHANGED")
            n = rttov.getcoefval.rttov_get_coef_val_r1('REF_TEMPERATURE')
            if len(n) < 20:
                self.write("REF_TEMPERATURE", +str(n))
            self.aKmats.reset()
            logging.info("send message Coefficients CHANGED")
            self.sendMessage("Coefficients CHANGED")
        else:
            self.myCoeffs.loadCoeffs = False
            sys.stderr.write("Cannot load coefficients \n")
            self.write("Cannot load coefficients")
            return
        return err

    def dropCoefficients(self):
        """ drop a coefficients files coefficient files must
            have been loaded"""
        err = -1
        if self.myCoeffs.loadCoeffs:
            self.beforeCallFortran()
            err = rttov_gui_f2py.rttov_gui_drop()
            self.afterCallFortran()
            if err == 0:
                self.myCoeffs.loadCoeffs = False
                self.myCoeffs.channel_list = None
                self.myCoeffs.nchannels = 0
                self.write("drop coefficient file OK ")
                self.ctrlCoherence()
                self.sendMessage("Options CHANGED")
                self.sendMessage("Coefficients CHANGED")
                self.myCoeffs.WAVENUMBERS = None
            else:
                self.write("Cannot drop coefficient files ")
                sys.stderr.write("Cannot drop coefficient files \n")
        else:
            self.write("Coefficient files not loaded ")
            sys.stderr.write("Coefficient files not loaded \n")
        return err

    def runDirect(self):
        err = self.runModel(mode=self.RTTOV_MODE)
        return err

    def isMW(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (self.myCoeffs.isMW()):
                return True
        return False

    def isSolar(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (self.hasSolar() and self.myOption["ADDSOLAR"].value):
                return True
        return False

    def hasSolar(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (self.myCoeffs.hasSolar()):
                return True
        return False

    def isPC(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (self.myCoeffs.hasPC() and self.myOption["ADDPC"].value):
                return True
        return False

    def isClouds(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (self.hasClouds() and self.myOption["ADDCLOUDS"].value):
                return True
        return False

    def hasClouds(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (self.myCoeffs.hasClouds() and self.myProfile.hasClouds()):
                return True
        return False

    def hasGas(self, gas):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if self.myCoeffs.hasGas(gas) and self.myProfile.hasGas(gas):
                return True
        return False

    def isAerosols(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (self.hasAerosols() and self.myOption["ADDAEROSL"].value):
                return True
        return False

    def hasAerosols(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (self.myCoeffs.hasAerosols() and self.myProfile.hasAerosols()):
                if (
                        self.myProfile.aerosol_kind == "opac" and
                        not self.myCoeffs.hasInTheName("cams", "aerosols")):
                    return True
                if (
                    self.myProfile.aerosol_kind == "cams" and
                        self.myCoeffs.hasInTheName("cams", "aerosols")):
                    return True
        return False

    def isPCClouds(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (
                    self.myCoeffs.isPCClouds() and
                    self.myOption["ADDPC"].value and
                    self.myOption["ADDCLOUDS"].value):
                return True
        return False

    def isPCtrace(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (
                    self.myCoeffs.hasInTheName("trace", "PC") and
                    self.myOption["ADDPC"].value):
                return True
        return False

    def isCAMS(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (
                    self.myCoeffs.hasInTheName("cams", "aer") and
                    self.myProfile.aerosol_kind == "cams"):
                return True
        return False

    def isMfasis(self):
        if self.loadProfile and self.myCoeffs.loadCoeffs:
            if (
                    self.myCoeffs.hasMfasis() and
                    self.myProfile.hasMfasis()):
                return True
        return False

    def runModel(self, mode):
        """
        save the profile and run the direct RTTOV model (a profile must exist)
        """
        err = -1
        if (not self.loadProfile):
            sys.stderr.write("profile must be loaded \n")
            err = -1
            return err
        self.ctrlCoherence()
        self.sendMessage("Options CHANGED")
        self.saveProfile()
        logging.debug("self.firstTimeSurface:" + str(self.firstTimeSurface))
        if self.previousRunIsPc:
            self.firstTimeSurface = True
            logging.info(">> previously K then now self.firstTimeSurface:" +
                         str(self.firstTimeSurface))
        if self.firstTimeSurface:
            self.newSurface()
            self.fistTimeSurface = False
        else:
            self.saveSurface(self.surfaceFileName)
        logging.info("self.firstTimeSurface= " + str(self.firstTimeSurface))
        self.write("rttov_gui_f2py.rttov_gui_run with: ")
        self.write(str(self.profileFileName))
        self.write(str(self.surfaceFileName))
        self.write(str(self.radianceFileName))
        self.write(self.KMatrixFileName)
        self.write(self.transFileName)
        self.write("mode:" + str(mode))
        self.write("nthreads:" + str(self.nthreads))
        self.beforeCallFortran()
        try:
            err = rttov_gui_f2py.rttov_gui_run(self.profileFileName,
                                               self.surfaceFileName,
                                               self.radianceFileName,
                                               self.KMatrixFileName,
                                               self.transFileName, mode,
                                               self.nthreads)
            self.afterCallFortran()
            self.write("rttov_gui_run completed return code: " + str(err))
        except Exception:
            sys.stderr.write("ERROR with rttov_gui_f2py.rttov_gui_run\n")
        if err == 0:
            self.write("rttov run mode " + mode + " ok \n")
            self.openSurface(self.surfaceFileName)
        else:
            self.write("rttov run mode " + mode + " in error \n")
        self.previousRunIsPc = False
        return err

    def runK(self, numPreviousRunK=None):
        """  save the profile and run the direct RTTOV model
         (a profile must exist) """
        err = self.runModel(mode="K")
        if err == 0:
            if numPreviousRunK is None:
                numberRunK = None
            else:
                numberRunK = numPreviousRunK + 1
            self.aKmats.update(self.KMatrixFileName, numberRunK)
            self.sendMessage("K-Matrix Changed")

        return err

    def runPC(self):
        self.previousRunIsPc = True
        mode = "PCDIRECT"
        err = self.runPCmodel(mode)
        return err

    def runPCK(self, numPreviousRunK=None):
        self.previousRunIsPc = True
        mode = "PCK"
        err = self.runPCmodel(mode)
        if self.myOption["ADDRADREC"].value and err == 0:
            if numPreviousRunK is None:
                numberRunK = None
            else:
                numberRunK = numPreviousRunK + 1
            self.aKmats.update(self.pcKMatrixFileName, numberRunK)
            self.sendMessage("K-Matrix Changed")
        return err

    def runPCmodel(self, mode):
        """ run PC model """
        err = -1
        if (not self.loadProfile):
            self.write("profile must be loaded \n")
            err = -1
            return err
        self.ctrlCoherence()
        self.sendMessage("Options CHANGED")
        if self.debug:
            self.myOption.display()
        self.saveProfile()
        logging.info("self.firstTimeSurface=" + str(self.firstTimeSurface))
        logging.debug("self.firstTimeSurface=" + str(self.firstTimeSurface))
        if self.firstTimeSurface:
            self.createSurface()
            self.firstTimeSurface = False
        else:
            self.saveSurface(self.surfaceFileName)
        logging.debug("self.firstTimeSurface=" + str(self.firstTimeSurface))
        self.npcscores = self.myOption["NPCSCORES"].value
        self.write("rttov_gui_f2py.rttov_gui_pcrun with: ")
        self.write("  profile : " + self.profileFileName)
        self.write("  surface : " + self.surfaceFileName)
        self.write("  pc      : " + self.pcFileName)
        self.write("  pcKMatrix: " + self.pcKMatrixFileName)
        self.write("  mode    : " + mode)
        self.write("  threads : " + str(self.nthreads))
        self.write("  npcscores:" + str(self.npcscores))
        self.beforeCallFortran()
        err = rttov_gui_f2py.rttov_gui_pcrun(self.profileFileName,
                                             self.surfaceFileName,
                                             self.pcFileName,
                                             self.pcKMatrixFileName, mode,
                                             self.nthreads, self.npcscores)
        self.afterCallFortran()
        self.write("rttov_gui_pcrun completed return code: " + str(err))
        if err == 0:
            self.write("rttov run mode " + mode + " ok \n")
            self.openSurface(self.surfaceFileName, ctrl=False)
        else:
            self.write("rttov run pc mode " + mode + " in error")
            sys.stderr.write("rttov run pc mode " + mode + " in error \n")
        return err

    def openAtlas(self, EMIS_ID_ATLAS=1, BRDF_ID_ATLAS=1):
        """ open an Atlas of Emissivity and/or BRDF"""
        """ the IDs of the different type of Atlas (emis, brdf) """
        """ are currently 1 and 1 """
        """ The CAMEL IR and Fatima Karbou Atlas for MW (both ID 2) """
        """ are not used """

        logging.debug("project openAtlas")
        logging.debug("EMIS_ID_ATLAS: %d" % EMIS_ID_ATLAS)
        logging.debug("BRDF_ID_ATLAS: %d" % BRDF_ID_ATLAS)
        if (not self.loadProfile):
            sys.stderr.write("profile must be loaded \n")
            err = -1
            return err
        self.saveProfile()

        for param in self.ListeParamAtlas:
            self.paramAtlas[param] = self.myProfile[param]

        self.write("load Atlas for lat: " + str(self.paramAtlas['LATITUDE']) +
                   " lon: " + str(self.paramAtlas['LONGITUDE']))
        self.write("date: " + str(self.paramAtlas['DATE'][0]) + " " +
                   str(self.paramAtlas['DATE'][1]) + " " +
                   str(self.paramAtlas['DATE'][2]))
        self.write("AZANGLE: " + str(self.paramAtlas['AZANGLE']))
        self.write("ZENANGLE: " + str(self.paramAtlas['ZENANGLE']))
        self.write("SUNAZANGLE: " + str(self.paramAtlas['SUNAZANGLE']))
        self.write("SUNZENANGLE: " + str(self.paramAtlas['SUNZENANGLE']))
        self.createSurface()
        self.write("rttov_gui_f2py.rttov_gui_get_emisbrdf with: ")
        self.write(self.config.ENV['RTTOV_GUI_EMISS_DIR'] + ", " +
                   self.profileFileName +
                   ", " + self.surfaceFileName)
        self.write("EMIS_ID_ATLAS: %d" % EMIS_ID_ATLAS)
        self.write("BRDF_ID_ATLAS: %d" % BRDF_ID_ATLAS)

        self.beforeCallFortran()
        err = rttov_gui_f2py.rttov_gui_get_emisbrdf(self.config.ENV[
            'RTTOV_GUI_EMISS_DIR'],
            self.profileFileName,
            self.surfaceFileName,
            EMIS_ID_ATLAS, BRDF_ID_ATLAS)
        self.afterCallFortran()
        self.write("load atlas completed: " + str(err))
        if err == 0:
            self.write("load atlas ok \n")
            self.openSurface(self.surfaceFileName)
            self.useAtlas = True
            self.firstTimeSurface = False
        else:
            self.write("cannot load atlas")
            sys.stderr.write("cannot load atlas \n")
        return err

    def displayRadiance(self):
        """ display the radiance file """
        if os.path.exists(self.radianceFileName):
            f = h5py.File(self.radianceFileName, 'r')
            h5 = f['/RADIANCE/']
            r1 = rttov.radiance.Radiance()
            r1.loadh5(h5)
            f.close()
            r1.display()

    def ctrlGasInProfAndCoefVsOption(self):
        """ controle gas in profile vs option status and value"""

        for gas in self.myProfile.gas_list[1:]:
            if gas == "O3":
                option = "OZONE_DATA"
            else:
                option = gas + "_DATA"
            if self.hasGas(gas):
                self.enableOption(option)
            else:
                self.disableOption(option)

    def disableGasIfRequired(self):
        """ controle gas in profile vs option status and value"""

        for gas in self.myProfile.gas_list[1:]:
            if gas == "O3":
                option = "OZONE_DATA"
            else:
                option = gas + "_DATA"
            if not self.hasGas(gas):
                self.disableOption(option)

    def initAllStatusAtTrue(self):
        for theme in ('General configuration options',
                      'General RT options',
                      'VIS/IR-only RT options',
                      'Interpolation and vertical grid options'):
            for option in self.myOption.optionsThemesList[theme]:
                self.myOption[option].status = True

    def disableOption(self, option):
        self.myOption[option].status = False
        if option in self.myOption.options_list_logical:
            self.myOption[option].value = False

    def enableOption(self, option):
        self.myOption[option].status = True

    def disableOptionsTheme(self, theme, exclude=[]):
        for option in self.myOption.optionsThemesList[theme]:
            if option not in exclude:
                self.disableOption(option)

    def enableOptionsTheme(self, theme, exclude=[]):
        for option in self.myOption.optionsThemesList[theme]:
            if option not in exclude:
                self.enableOption(option)

    def ctrlOptionMfasis(self):
        """ control if MFASIS can be enabled """
        if self.myCoeffs.loadCoeffs:
            if (
                    self.myOption["ADDCLOUDS"].value and
                    self.myOption["ADDSOLAR"].value and
                    not self.myOption["ADDAEROSL"].value and
                    self.myCoeffs.hasMfasis()):
                self.myOption[
                    "VIS_SCATT_MODEL"].odict = self.myOption.dictOptions[
                    "ALL_VIS_SCATT_MODEL"]
            else:
                self.myOption[
                    "VIS_SCATT_MODEL"].odict = self.myOption.dictOptions[
                    "VIS_SCATT_MODEL"]

    def ctrlOptionsPC(self):
        """ control in PC cases if all other options are OK """
        if self.myOption["IPCBND"].value < 1:
            self.myOption["IPCBND"].value = 1
        if self.myOption["IPCREG"].value < 1:
            self.myOption["IPCREG"].value = 1
        if not self.myCoeffs.loadCoeffs:
            self.disableOptionTheme("PC-RTTOV options")
            return

        if self.myCoeffs.hasPC():
            self.enableOptionsTheme("PC-RTTOV options")
            if self.isPC():
                self.myOption.updatePCChoicesList(self.myCoeffs)
                if self.isPCtrace():
                    exclude_list = ["DO_NLTE_CORRECTION"]
                else:
                    exclude_list = ["DO_NLTE_CORRECTION",
                                    "ADDCLOUDS"]
                    for option in ["OZONE_DATA", "CO2_DATA", "N2O_DATA",
                                   "CO_DATA", "CH4_DATA", "ADDCLOUDS",
                                   "ADDAEROSL"]:
                        self.disableOption(option)
                self.disableOptionsTheme("VIS/IR-only RT options",
                                         exclude=exclude_list)
                self.disableOption("DO_LAMBERTIAN")
                self.disableOption("LAMBERTIAN_FIXED_ANGLE")
                if not self.myCoeffs.hasNLTE():
                    self.disableOption("DO_NLTE_CORRECTION")
            else:
                for option in ("IPCBND", "IPCREG", "ADDRADREC"):
                    self.disableOption(option)
                self.enableOptionsTheme("VIS/IR-only RT options")
                self.ctrlGasInProfAndCoefVsOption()
        else:
            self.disableOptionsTheme("PC-RTTOV options")
            self.myOption["IPCBND"].value = 1
            self.myOption["IPCREG"].value = 1
            return

        if (self.myCoeffs.isPCClouds() and
                self.myProfile.hasClouds()):
            self.enableOption("ADDCLOUDS")
        else:
            self.disableOption("ADDCLOUDS")

        # 'ADDRADREC' and 'SWITCHRAD'  coherence
        if self.isPC():
            if self.myOption['ADDRADREC'].value is False:
                self.disableOption('SWITCHRAD')
            else:
                self.enableOption('SWITCHRAD')

    def isScatt(self):
        return self.isAerosols() or self.isClouds()

    def isLambertian(self):
        return self.myOption["DO_LAMBERTIAN"].value

    def isDomScatt(self):
        if not self.isScatt():
            return False
        if self.isSolar():
            isdom = (self.myOption["IR_SCATT_MODEL"].value == 1 or
                     self.myOption["VIS_SCATT_MODEL"].value == 1)
        else:
            isdom = (self.myOption["IR_SCATT_MODEL"].value == 1)
        return isdom

    def isDomScattSolar(self):
        value = self.isSolar() and self.isDomScatt()
        return value

    def ctrlCoherence(self):
        """ control the coherence between options and profile """
        """ and loaded coefficients"""
        """ options which can't be modified have a False status """
        """ warning: control order is very important !!!"""
        debug = False
        if debug:
            print("DEBUG ctrlCoherence")
            print(self.profileFileName)
            print(self.myCoeffs.fileName)

        for option in self.myOption.options_list:
            self.myOption[option].status = False

        if not self.myCoeffs.loadCoeffs:
            return

        self.initAllStatusAtTrue()
        if debug:
            print("DEBUG ctrlGasInProfAndCoefVsOption ")
        self.ctrlGasInProfAndCoefVsOption()
        if debug:
            print("DEBUG ctrlOptionsPC")
        self.ctrlOptionsPC()

        # Aerosols control
        if debug:
            print("DEBUG ctrl aerosols")
        if self.hasAerosols():
            self.enableOption("ADDAEROSL")
        else:
            self.disableOption("ADDAEROSL")

        # cloud control
        if debug:
            print("DEBUG ctrl clouds")
        if (self.myCoeffs.fileName["clouds"] != "" and
                self.myProfile.anyCloud()):
            self.enableOption("ADDCLOUDS")
            self.enableOption("CLDCOL_THRESHOLD")
            self.enableOption("GRID_BOX_AVG_CLOUD")
        else:
            self.disableOption("ADDCLOUDS")
            self.disableOption("CLDCOL_THRESHOLD")
            self.disableOption("GRID_BOX_AVG_CLOUD")

        # Micro Wave control
        if debug:
            print("DEBUG ctrl MW")
        if self.isMW():
            self.enableOptionsTheme("MW-only RT options")
            self.disableOptionsTheme("VIS/IR-only RT options")
            if self.myOption["FASTEM_VERSION"].value == 0:
                self.disableOption("SUPPLY_FOAM_FRACTION")
            else:
                self.enableOption("SUPPLY_FOAM_FRACTION")
        else:
            self.disableOptionsTheme("MW-only RT options")
        if self.myProfile['CLW'] is None:
            self.disableOption("CLW_DATA")

        # solar
        if debug:
            print("DEBUG ctrl solar")
        if self.hasSolar() and not self.isPC():
            self.enableOption("ADDSOLAR")
        else:
            self.disableOption("ADDSOLAR")

        if self.isSolar():
            self.enableOption("RAYLEIGH_SINGLE_SCATT")
            self.enableOption("RAYLEIGH_MAX_WAVELENGTH")
            self.enableOption("RAYLEIGH_MIN_PRESSURE")
        else:
            self.disableOption("RAYLEIGH_SINGLE_SCATT")
            self.disableOption("RAYLEIGH_MAX_WAVELENGTH")
            self.disableOption("RAYLEIGH_MIN_PRESSURE")

        # ADDSOLAR not compatible with PC
        if self.myOption["ADDSOLAR"].value:
            self.disableOptionsTheme("PC-RTTOV options")

        # SIMPLE CLOUD: not done for the moment
        for option in ("USER_AER_OPT_PARAM", "USER_CLD_OPT_PARAM"):
            self.disableOption(option)

        # NLTE
        if not self.myCoeffs.hasNLTE():
            self.disableOption("DO_NLTE_CORRECTION")
        else:
            self.enableOption("DO_NLTE_CORRECTION")

        if debug:
            print("DEBUG scatt")
        # scatt_options only possible if ADDAERSOSL or ADDCLOUDS checked
        self.ctrlOptionMfasis()
        if self.isScatt():
            for opt in self.myOption.scatt_options:
                self.enableOption(opt)
            if not self.isSolar():
                self.disableOption("VIS_SCATT_MODEL")
            if not self.isDomScatt():
                self.disableOption("DOM_NSTREAMS")
                self.disableOption("DOM_ACCURACY")
                self.disableOption("DOM_OPDEP_THRESHOLD")
                self.disableOption("DOM_RAYLEIGH")
        else:
            for opt in self.myOption.scatt_options:
                self.disableOption(opt)

        if debug:
            print("DEBUG  disable gas if required")
        self.disableGasIfRequired()

        if self.isLambertian():
            self.enableOption("LAMBERTIAN_FIXED_ANGLE")
        else:
            self.disableOption("LAMBERTIAN_FIXED_ANGLE")

        if self.isDomScattSolar() and self.myCoeffs.is_V13():
            self.enableOption("DOM_RAYLEIGH")
        else:
            self.disableOption("DOM_RAYLEIGH")


def OpenAProfile(filename, number):
    # read the profile
    profile = pProfile()
    if number == 1:
        f = h5py.File(filename, "r")
        h5 = f['/PROFILES/0001/']
        profile.loadh5(h5)
        f.close()
    else:
        profile.loadProfileNumber(filename, number)
    profile.controlProfile()
    return profile


def checkCode(err):
    if err != 0:
        print((">>>> ERROR %d" % (err)))


if __name__ == "__main__":
    p = Project()
    profileName = (
        p.config.ENV['RTTOV_GUI_PROFILE_DIR'] + "/cldaer101lev_allgas.H5")
    err = p.openProfile(profileName, 1)
    p.myOption.display()
    checkCode(err)
    coefFile = (p.config.ENV['RTTOV_GUI_COEFF_DIR'] +
                "/rttov7pred54L/rtcoef_noaa_19_hirs.dat")
    cloudCoefFile = (
        p.config.ENV['RTTOV_GUI_COEFF_DIR'] +
        "/cldaer/sccldcoef_noaa_19_hirs.dat")
    p.myCoeffs.fileName["standard"] = coefFile

    p.myOption["IR_SEA_EMIS_MODEL"].value = 1

    err = p.loadCoefficients()
    checkCode(err)
    print("runDirect")
    err = p.runDirect()
    checkCode(err)
    p.displayRadiance()
    err = p.runK()
    checkCode(err)
    err = p.dropCoefficients()
    checkCode(err)
    print("... END OK")
