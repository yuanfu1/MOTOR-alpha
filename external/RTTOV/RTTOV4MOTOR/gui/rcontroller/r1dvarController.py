# -*- coding: utf-8 -*-

from util import GenericController
import rview
import rmodel
import rttov
import copy
import r1Dvar
import time
import argparse
from util import is_valid_file
import sys
import logging
import wx
from pubsub import pub

r1dvarViewTitle = "1DVAR Algorithm"
r1dvarProfViewTitle = "Retrieved profile"
r1dvarRadViewTitle = "Brightness Temperatures differences"


class R1dvarController (GenericController):
    """ controller for the 1Dvar Window """

    def __init__(self, parent, project=None):
        """ Constructor of the R1dvarController
            controls the r1dvarView Windows
            the r1dVarView windows can have 2 children :
              - the r1dvarprofileView window
              - the RBtView
            and the Project1dvar object

        """
        self.theParentFrame = parent
        super().__init__()
        self.myViews[r1dvarViewTitle] = None
        self.myViews[r1dvarProfViewTitle] = None
        self.myViews[r1dvarRadViewTitle] = None
        self.profileFileName = None

        self.controlerName = "R1DvarController"
        self.retrieveProject = None
        self.pTrue = None
        self.stepCounter = 0
        self.pBgOrigin = rmodel.project.Project(0)
        if project:
            self.project = project
        pub.subscribe(self.CoeffsLoadedListener, "Coefficients CHANGED")
        self.hasTrueProfile = False
        self.changeTheBackground = False
        self.oldR1dvarradianceView = []

    def ProfileLoadedControl(self):
        """ the profile have changed """
        """ must prevent to continue with 1dvar if not 54 levels """
        self.prepLog()
        logging.debug("r1dvar Controller ProfileLoadedControl")
        isValid = self.controleProfile(self.project.myProfile)
        mustClose1dVarWindow = False
        logging.debug("is profile valid ? " + str(isValid))
        if self.retrieveProject:
            logging.debug("we had a retrieveProject")
            if isValid:
                isValidwithCoeff = self.controleProject(
                    self.project.myProfile, self.project.myCoeffs)
                if isValidwithCoeff:
                    self.reset()
                else:
                    self.retrieveProject = None
                    self.stepCounter = 0
                    mustClose1dVarWindow = True
            else:
                self.retrieveProject = None
                self.stepCounter = 0
                mustClose1dVarWindow = True
        else:
            if not isValid:
                mustClose1dVarWindow = True
        logging.debug(
            "must close the 1 dvar window " + str(mustClose1dVarWindow))
        if self.myViews[r1dvarViewTitle] and mustClose1dVarWindow:
            self.myViews[r1dvarViewTitle].WarmError(
                "invalid Background profile for 1D-VAR")
            self.myViews[r1dvarViewTitle].Close()

    def CoeffsLoadedListener(self, msg):
        """ the coefficients have changed """
        """ must prevent to continue with the steps """
        self.prepLog()
        logging.debug("r1dvar Controller CoeffseLoadedListener")
        # must control if the project is still valid
        isValid = self.controleProject(
            self.pBgOrigin.myProfile, self.project.myCoeffs)
        if self.retrieveProject:
            isValid = self.controleProject(
                self.pBgOrigin.myProfile, self.project.myCoeffs)
            if isValid:
                (inst, platform, num) = (
                    self.project.myCoeffs.get_instrument_and_platform_name())
                satellite = platform + "-" + str(num)
                self.write("satellite " + satellite + " instrument " + inst)
                self.pBgOrigin.myCoeffs = copy.deepcopy(self.project.myCoeffs)
                self.retrieveProject.reinitCoeff(
                    satellite, inst, self.pBgOrigin.myCoeffs)
            else:
                self.retrieveProject = None
                self.reset()
                self.stepCounter = 0
                if self.myViews[r1dvarViewTitle]:
                    self.myViews[r1dvarViewTitle].Close()

    def InitBgProject(self):
        """ Init a RTTOV BackGround Project from project """
        self.pBgOrigin = rmodel.project.Project(0)
        self.pBgOrigin.myProfile = copy.deepcopy(self.project.myProfile)
        self.pBgOrigin.myOption = self.project.myOption.__deepcopy__()
        self.pBgOrigin.myCoeffs = copy.deepcopy(self.project.myCoeffs)
        self.pBgOrigin.profileInitialFilename = copy.deepcopy(
            self.project.profileInitialFilename)
        self.pBgOrigin.loadProfile = self.project.loadProfile

    def Show(self):
        """ show the r1dvarView  """
        exists = self.ShowWindowIfExists(self.myViews[r1dvarViewTitle])
        if not exists:
            self.myViews[r1dvarViewTitle] = rview.r1dvarView.R1dvarView(
                self.theParentFrame,
                r1dvarViewTitle,
                self.pBgOrigin.config.ENV['RTTOV_GUI_PROFILE_DIR'])
            self.MakeBinding()

    def MakeBinding(self):
        """  Define the different binding of the application  """
        # binding for the MaineView (mv)
        self.myViews[r1dvarViewTitle].Bind(
            wx.EVT_BUTTON, self.OpenATrueProfile, self.myViews[r1dvarViewTitle].btnOpen)
        self.myViews[r1dvarViewTitle].Bind(wx.EVT_MENU, self.OpenATrueProfile,
                                           self.myViews[r1dvarViewTitle].items["openTrueProfile"])
        self.myViews[r1dvarViewTitle].Bind(wx.EVT_BUTTON, self.OnStep,
                                           self.myViews[r1dvarViewTitle].btnStep)
        self.myViews[r1dvarViewTitle].Bind(wx.EVT_BUTTON, self.OnReset,
                                           self.myViews[r1dvarViewTitle].btnReset)

    def controleProject(self, aProfile, aCoeff):
        """ control if it is a valid project for 1D var """
        """ nb level must be 54 """
        """ for profile and coefficients """
        flagProfile = self.controleProfile(aProfile)
        flagCoeff = self.controleCoeff(aCoeff)
        status = flagProfile and flagCoeff

        return (status)

    def controleCoeff(self, aCoeff):
        if aCoeff.loadCoeffs:
            (inst, sat, id_sat) = aCoeff.get_instrument_and_platform_name()
            logging.info("satellite : " + sat +
                         str(id_sat) + " instrument : " + inst)
            if inst in ["airs", "iasi", "amsua", "amsub", "mhs",
                        "hirs", "atms", "ssmis", "cris"]:
                return True
            else:
                return False
        else:
            return False

    def controleProfile(self, aProfile):

        if aProfile:
            nlevels = aProfile["NLEVELS"]

            if nlevels == 54:
                return True
            else:
                return False
        else:
            return False

    def OpenATrueProfile(self, e):
        """ Open the profile , if it contains more than one
            profile ask for the number of the profile """
        self.prepLog()

        self.profileFileName = self.myViews[
            r1dvarViewTitle].OnOpenTrueProfile(e)

        self.write("Open A true profile filename= " +
                   str(self.profileFileName))
        if (self.profileFileName):
            try:
                self.pTrue = rmodel.project.Project(0)
                number = rttov.profile.getNumberOfProfiles(
                    self.profileFileName)
                self.write("number of profile n this file : " + str(number))
                if (number > 1):
                    number = self.myViews[r1dvarViewTitle].ChooseNumber(
                        number, "select a profile number",
                        "profile number selection")
                self.pTrue = rmodel.project.Project(0)
                self.pTrue.openProfile(self.profileFileName, number)
                okProfile = self.controleProfile(self.pTrue.myProfile)
                if okProfile:
                    self.hasTrueProfile = True
                    # if we had a previous retriveProject must create a new one
                    # because we have now a new pTrue
                    if self.retrieveProject:
                        self.retrieveProject = None
                        self.stepCounter = 0
                        if self.myViews[r1dvarProfViewTitle]:
                            self.myViews[r1dvarProfViewTitle].Close()
                            self.myViews[r1dvarProfViewTitle] = None
                        if self.myViews[r1dvarRadViewTitle]:
                            self.myViews[r1dvarRadViewTitle].Close()
                            self.myViews[r1dvarRadViewTitle] = None

                    # we suppose that coeffs are loaded
                    if not self.project.myCoeffs.loadCoeffs:
                        self.myViews[r1dvarViewTitle].WarmError(
                            "Coefficient not loaded ! ")
                    else:
                        okCoeff = self.controleCoeff(self.project.myCoeffs)
                        if okCoeff:
                            if not self.retrieveProject:
                                (inst, platform, num) = \
                                    self.project.myCoeffs.\
                                    get_instrument_and_platform_name()
                                satellite = platform + "-" + str(num)
                                self.write("satellite " +
                                           satellite + " instrument " + inst)
                                self.satellite = satellite
                                self.instrument = inst
                                self.write("Init the Retrieve Project")
                                self.InitBgProject()
                                self.retrieveProject = (
                                    r1Dvar.r1dvar.Project1dvar(
                                        self.pTrue, self.pBgOrigin,
                                        self.satellite, self.instrument))
                                self.write("Retrieve Project initialized")

                            else:
                                self.reset()
                        else:
                            self.myViews[r1dvarViewTitle].WarmError(
                                "the coefficients file is not "
                                "compatible with the 1Dvar algorithm")
                            return
                else:
                    self.myViews[r1dvarViewTitle].WarmError(
                        "the profile must have 54 levels")
                    # profile not OK
                    # warm
                    self.pTrue = None
                    self.hasTrueProfile = False
                    self.retrieveProject = None
                    self.stepCounter = 0
                    if self.myViews[r1dvarProfViewTitle]:
                        self.myViews[r1dvarProfViewTitle].Close()
                        self.myViews[r1dvarProfViewTitle] = None
                    if self.myViews[r1dvarRadViewTitle]:
                        self.myViews[r1dvarRadViewTitle].Close()
                        self.myViews[r1dvarRadViewTitle] = None

                if self.myViews[r1dvarViewTitle]:
                    self.myViews[r1dvarViewTitle].Restore()
            except IOError:
                self.write(
                    "Error while opening and making 1dvar project " +
                    str(self.profileFileName))
                self.pTrue = None

    def OnStep(self, e):
        """ perform a step of retrieval """
        """ max step = 10 """

        if self.stepCounter > 9:
            self.myViews[r1dvarViewTitle].WarmError(
                "maximum step range reached ! ")
            return

        time.sleep(2)
        if not self.hasTrueProfile:
            self.myViews[r1dvarViewTitle].WarmError(
                "you must open a True Profile !")
            return
        if not self.retrieveProject:
            self.reset()
            self.stepCounter = 0

        maxNoise = self.myViews[r1dvarViewTitle].GetMaxNoise()
        obsError = self.myViews[r1dvarViewTitle].GetObsError()
        bgError = self.myViews[r1dvarViewTitle].GetBgError()
        self.retrieveProject.setFactorB(bgError)
        self.retrieveProject.setFactorR(obsError)
        self.retrieveProject.setMaxNoise(maxNoise)

        self.myViews[r1dvarViewTitle].btnStep.Unbind(wx.EVT_BUTTON)

        print(("maxNoise", maxNoise, 'obsError', obsError, 'bgError', bgError))
        try:
            self.myViews[r1dvarViewTitle].BeginBusy()

            if self.stepCounter == 0:
                # initialize the retrieval project
                logging.info("initialize the retrieval project")
                self.retrieveProject.retrieve1d()
            else:
                logging.info("stepCounter= " + str(self.stepCounter))
                self.retrieveProject.stepRetrieve1d()
            self.myViews[r1dvarViewTitle].EndBusy()
            self.stepCounter = self.stepCounter + 1
            exists = self.ShowWindowIfExists(self.myViews[r1dvarProfViewTitle])
            if not exists:
                self.myViews[r1dvarProfViewTitle] = (
                    rview.r1dvarprofileframe.r1dvarProfileView(
                        self.myViews[r1dvarViewTitle],
                        r1dvarProfViewTitle,
                        self.retrieveProject))
                self.myViews[r1dvarProfViewTitle].Show()
            else:
                self.myViews[r1dvarProfViewTitle].RePlot(self.retrieveProject)

            # open the window with BT
            radBgFile = self.retrieveProject.pBg.radianceFileName
            radTrueFile = self.retrieveProject.pTrue.radianceFileName.replace(
                ".h5", "_Noise.h5")

            exists = self.ShowWindowIfExists(self.myViews[r1dvarRadViewTitle])
            if not exists:
                self.myViews[
                    r1dvarRadViewTitle] = rview.rBtView.RBtView(
                    self.myViews[r1dvarViewTitle], r1dvarRadViewTitle,
                    radTrueFile, radBgFile, self.stepCounter)
                self.myViews[r1dvarRadViewTitle].Show()
            else:
                # must control if the r1dvarradianceView can be updated
                # if not make a new one
                # the old one will live its life and will be close by
                # the user or when the gui will be closed
                # the r1dvarradianceView can be updated only if we have the
                # same number of channel
                nbChanFromRadianceView = self.myViews[r1dvarRadViewTitle].nchan
                if nbChanFromRadianceView != self.pBgOrigin.myCoeffs.nchannels:
                    self.oldR1dvarradianceView.append(
                        self.myViews[r1dvarRadViewTitle])
                    self.myViews[r1dvarRadViewTitle] = rview.rBtView.RBtView(
                        self.myViews[r1dvarViewTitle],
                        r1dvarRadViewTitle,
                        radTrueFile,
                        radBgFile,
                        self.stepCounter)
                    self.myViews[r1dvarRadViewTitle].Show()
                else:
                    if self.changeTheBackground:
                        self.myViews[r1dvarRadViewTitle].ReRead(radTrueFile,
                                                                radBgFile)
                        self.myViews[r1dvarRadViewTitle].Show()

        except Exception:
            self.myViews[r1dvarViewTitle].WarmError(
                "Sorry a problem occurred  ")
            self.myViews[r1dvarViewTitle].EndBusy()
            self.myViews[r1dvarViewTitle].Bind(
                wx.EVT_BUTTON,
                self.OnStep,
                self.myViews[r1dvarViewTitle].btnStep)
            return

        if self.stepCounter >= 1:
            if self.myViews[r1dvarViewTitle]:
                self.myViews[r1dvarViewTitle].slider["noise"].Disable()

        self.myViews[r1dvarViewTitle].Bind(
            wx.EVT_BUTTON, self.OnStep, self.myViews[r1dvarViewTitle].btnStep)

    def OnReset(self, e):
        """ perform a reset of the retrieveProject :
            take into account the New Background from the main window"""
        """ we close the r1dvar window we keep the true profile """
        self.reset()

    def reset(self):
        logging.info("reset 1DVAR project")
        if self.myViews[r1dvarProfViewTitle]:
            self.myViews[r1dvarProfViewTitle].Close()
        if self.myViews[r1dvarRadViewTitle]:
            self.myViews[r1dvarRadViewTitle].Close()

        self.myViews[r1dvarProfViewTitle] = None
        self.myViews[r1dvarRadViewTitle] = None
        for r1dvarView in self.oldR1dvarradianceView:
            if r1dvarView:
                r1dvarView.Close()
        isOK = self.controleProject(
            self.project.myProfile, self.project.myCoeffs)
        if isOK:
            self.stepCounter = 0
            if self.myViews[r1dvarViewTitle]:
                self.myViews[r1dvarViewTitle].slider["noise"].Enable()
            self.InitBgProject()
            (inst, platform, num) = (
                self.project.myCoeffs.get_instrument_and_platform_name())
            satellite = platform + "-" + str(num)
            self.write("satellite " + satellite + " instrument " + inst)
            self.satellite = satellite
            self.instrument = inst
            logging.info("initialize a new retrieve project")
            self.InitBgProject()
            self.retrieveProject = r1Dvar.r1dvar.Project1dvar(
                self.pTrue, self.pBgOrigin, self.satellite, self.instrument)
        else:
            logging.warn("cannot initialize a new retrieve project")
            if self.myViews[r1dvarViewTitle]:
                self.myViews[r1dvarViewTitle].WarmError(
                    "cannot initialize a new retrieve project : "
                    "coefficients file or profile is not compatible ! ")


if __name__ == "__main__":

    app = wx.App(False)

    parser = argparse.ArgumentParser(
        description="launch RTTOV GUI 1Dvar", conflict_handler='resolve')
    parser.add_argument('-P', '--profile',
                        dest="profile",
                        help="input file name",
                        required=False,
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-s', '--std-coeff',
                        dest="stdcoeff",
                        help="standard coefficient file",
                        required=False,
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))
    args = parser.parse_args()
    pBg = rmodel.project.Project(0)
    if args.profile is not None:
        filename = args.profile.name
    else:
        filename = pBg.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"

    sys.stderr.write("... load profile 1 from" + filename + "\n")
    try:
        err = pBg.openProfile(filename, 1)
        if err != 0:
            sys.stderr.write("ERROR wrong profile file \n")
            logging.info("wrong profile names ")
            sys.exit(1)
    except Exception:
        sys.stderr.write("ERROR wrong profile names \n")
        logging.info("wrong profile names ")
        sys.exit(1)
    if args.stdcoeff is not None:
        coefFile = args.stdcoeff.name
    else:
        coefFile = pBg.config.ENV["RTTOV_GUI_COEFF_DIR"] + \
            "/rttov7pred54L/rtcoef_noaa_19_hirs.dat"
    sys.stderr.write("... load coefficients from" + coefFile + "\n")
    pBg.myCoeffs.fileName["standard"] = coefFile
    err = pBg.loadCoefficients()
    if err != 0:
        sys.stderr.write("ERROR cannot load coefficients \n")
        sys.exit(1)

    ex = wx.App()

    r1dcController = R1dvarController(None, pBg)
    isValid = r1dcController.controleProject(pBg.myProfile, pBg.myCoeffs)
    if isValid:
        r1dcController.Show()
        ex.MainLoop()
    else:
        sys.stderr.write(
            "ERROR cannot run 1dVar (not the right instrument"
            " or wrong number of levels) \n")
        sys.exit(1)
