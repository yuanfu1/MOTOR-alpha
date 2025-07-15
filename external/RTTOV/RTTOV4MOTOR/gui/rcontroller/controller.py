# -*- coding: utf-8 -*-

try:
    import matplotlib
    matplotlib.use('wxagg')
except Exception:
    import sys
    sys.stderr.write(
        'ERROR : cannot start rttovgui : check your python installation\n')
    sys.stderr.write(
        'ERROR : matplotlib is not installed or is installed'
        ' with an incorrect backend (wxagg is required)\n')
    sys.exit(1)

import rmodel
import wx
import rttov
import sys
import logging
import rview
from util import GenericController
from util import is_valid_file
from optionctrl import OptionController
from surfacectrl import SurfaceController
from profilectrl import ProfileController
from r1dvarController import R1dvarController
import argparse
import warnings
from pubsub import pub

radViewTitle = "Radiance Viewer"
pcViewTitle = "PC Scores"
kpcViewTitle = "KP C Viewer run "
kViewTitle = "K Matrix Viewer run "
kpcProfViewTitle = "KP C Profile View run "


class MainController (GenericController):
    """ Main Controller of the application """

    def __init__(self, app, debug):
        """ Constructor of the MainController
            define a Project object
            define the MainView Object
            define the undo redo Object
            make the binding with the different view

        """
        super().__init__()

        GenericController.project = rmodel.project.Project()
        self.controlerName = "MainController"

        debug = True
        if debug:
            level_logging = logging.DEBUG
        else:
            level_logging = logging.INFO

        logging.basicConfig(
            filename=self.project.config.ENV['GUI_WRK_DIR'] + "/rttovgui.log",
            filemode="w",
            format=("[%(asctime)s] %(levelname)s [%(module)s:%(funcName)s:"
                    "%(lineno)d] %(message)s"),
            level=level_logging,
            datefmt="%Y:%m:%d %H:%M:%S")
        logging.info("start main controller")
        # set up the Mainview
        self.mv = rview.console.MainView(None)

        self.optionController = OptionController(self.mv)
        self.chooseCoeffsController = ChooseCoeffsController(self.mv)
        self.surfaceController = SurfaceController(self.mv)
        self.profileController = ProfileController(self.mv)
        self.r1dvarController = R1dvarController(self.mv)

        # set the binding
        self.MakeBinding()
        self.myViews[radViewTitle] = None
        self.myViews[pcViewTitle] = None
        self.RunkNumber = 0
        pub.subscribe(self.ProfileLoadedListener, "Profile LOADED")
        pub.subscribe(self.CoeffsLoadedListener, "Coefficients CHANGED")
        self.write("Welcome  ! ")

    def ProfileLoadedListener(self, msg):
        self.prepLog()
        logging.debug("mainController ProfileLoadedListener")
        self.mv.EnableMenuItem('profileWindow')
        self.mv.EnableMenuItem('optionsWindow')
        self.mv.EnableMenuItem('surfaceWindow')
        if self.project.myCoeffs.loadCoeffs:
            self.mv.EnableMenuItem('runDirect')
            self.mv.EnableMenuItem('runK')
        self.Ctrl1dvarMenu()
        self.r1dvarController.ProfileLoadedControl()

    def CoeffsLoadedListener(self, msg):
        """ must be called when the coefficient are loaded or dropped
        enable or disable menu rttov                              """
        self.prepLog()
        if self.project.myCoeffs.loadCoeffs:
            if self.project.myCoeffs.fileName["PC"] == "":
                self.mv.EnableMenuItem('selectChannels')
            if self.project.loadProfile:
                self.mv.EnableMenuItem('runDirect')
                self.mv.EnableMenuItem('runK')

            if self.myViews[radViewTitle]:
                # if we have a radiance frame : must close it (if not
                # add further control in radianceframe.py
                # (for channel selection)
                self.mv.write("instrument change : close the radiance window")
                self.myViews[radViewTitle].WarmError(
                    "radiance window will close")
                self.myViews[radViewTitle].Close()
                self.myViews[radViewTitle] = None
        else:
            self.mv.DisableMenuItem('selectChannels')
            self.mv.DisableMenuItem('runDirect')
            self.mv.DisableMenuItem('runK')
        self.Ctrl1dvarMenu()

    def SaveTheProfileAs(self, e):
        """ apply profile change and save profile+option with a new name """
        self.prepLog()
        # re-initialize the name of the saved profile file to None
        self.project.savedProfileFileName = None
        self.SaveTheProfile(e)

    def MakeBinding(self):
        """  Define the different binding of the application  """
        # binding for the MaineView (mv)
        self.mv.Bind(wx.EVT_MENU, self.OpenProfile,
                     self.mv.items["openProfile"])
        self.mv.Bind(wx.EVT_MENU, self.OpenAsciiProfile,
                     self.mv.items["openAsciiProfile"])
        self.mv.Bind(wx.EVT_MENU, self.mv.OnLoadCoefficients,
                     self.mv.items["loadCoefficients"])
        self.mv.Bind(wx.EVT_MENU, self.SaveTheProfile,
                     self.mv.items["saveProfile"])
        self.mv.Bind(wx.EVT_MENU, self.optionController.OnOptions,
                     self.mv.items["optionsWindow"])
        self.mv.Bind(wx.EVT_MENU, self.surfaceController.OnSurface,
                     self.mv.items["surfaceWindow"])
        self.mv.Bind(wx.EVT_MENU, self.profileController.OnProfile,
                     self.mv.items["profileWindow"])
        self.mv.Bind(wx.EVT_MENU, self.SaveTheProfileAs,
                     self.mv.items["saveProfileAs"])
        self.mv.Bind(wx.EVT_MENU, self.LoadCoefficients,
                     self.mv.items["loadCoefficients"])
        self.mv.Bind(wx.EVT_MENU, self.SelectChannels,
                     self.mv.items["selectChannels"])
        self.mv.Bind(wx.EVT_MENU, self.RunDirect, self.mv.items["runDirect"])
        self.mv.Bind(wx.EVT_MENU, self.RunK, self.mv.items["runK"])
        self.mv.Bind(wx.EVT_MENU, self.Config1Dvar,
                     self.mv.items["config1Dvar"])

    def SaveTheProfile(self, e):  self.SaveProfile(e, self.mv)

    def RunDirect(self, e):
        """ Run direct model, if project.isPC() run direct PC"""
        self.prepLog()
        if self.project.isPC():
            self._RunDirectPC(e)
        else:
            self._RunDirect(e)

    def _RunDirectPC(self, e):
        """ run direct model with PC """
        if self.project.myOption["ADDSOLAR"].value:
            addsolar = 1
        else:
            addsolar = 0
        # ask number or pc :
        if self.project.myCoeffs.loadCoeffs and self.project.isPC():
            if self.project.myProfile['SKIN']['SURFTYPE'] != 1:
                if "landsea" not in self.project.myCoeffs.fileName["PC"]:
                    self.mv.WarmError(
                        "pc coefficient file not compatible with"
                        " a land or a sea ice surface")
                    return
            if self.project.myOption["NPCSCORES"].value < 0:
                self.project.npcscores = 200
            else:
                self.project.npcscores = self.project.myOption[
                    "NPCSCORES"].value
            self.mv.BeginBusy()
            err = self.project.runPC()
            self.mv.EndBusy()
            self.write("RunDirectPC return code = " + str(err))
            if err == 0:
                exists = self.ShowWindowIfExists(self.myViews[pcViewTitle])
                if exists:
                    self.myViews[pcViewTitle].ReRead(self.project.pcFileName)
                    self.myViews[pcViewTitle].Show()
                else:
                    self.myViews[pcViewTitle] = rview.pcView.pcView(
                        self.mv, pcViewTitle, self.project.pcFileName)
                if self.project.myOption["ADDRADREC"].value:
                    if self.myViews[radViewTitle]:
                        if self.myViews[radViewTitle].IsIconized():
                            self.myViews[radViewTitle].Restore()
                        self.myViews[radViewTitle].read_rad_h5(
                            self.project.pcFileName, addsolar)
                        self.myViews[radViewTitle].Refresh()
                    else:
                        self.myViews[
                            radViewTitle] = rview.radianceframe.RadianceFrame(
                            self.mv, radViewTitle,
                            self.project.pcFileName, addsolar)
                        self.myViews[radViewTitle].Show()
            else:
                logging.info("pc-rttov runDirect Failed")
                self.mv.WarmError("Run direct Failed ")

    def _RunDirect(self, e):
        """ run the direct model if coefficient """
        """ are loaded and a profile saved"""
        logging.debug("_RunDirect")
        if self.project.myOption["ADDSOLAR"].value:
            addsolar = 1
        else:
            addsolar = 0
        if self.project.myCoeffs.loadCoeffs:
            self.mv.BeginBusy()
            err = self.project.runDirect()
            self.mv.EndBusy()
            if err == 0:
                exists = self.ShowWindowIfExists(self.myViews[radViewTitle])
                if exists:
                    self.myViews[radViewTitle].read_rad_h5(
                        self.project.radianceFileName, addsolar)
                    self.myViews[radViewTitle].Refresh()
                else:
                    self.myViews[
                        radViewTitle] = rview.radianceframe.RadianceFrame(
                        self.mv, radViewTitle,
                        self.project.radianceFileName, addsolar)
                    self.myViews[radViewTitle].Show()
            else:
                logging.info("rttov runDirect Failed")
                self.mv.WarmError("Run direct Failed ")
        else:
            self.mv.WarmError("You must load a Coefficients File ")

    def RunK(self, e):

        self.prepLog()
        if self.project.isPC():
            self._RunKPC(e)
        else:
            self._RunK(e)

    def _RunK(self, e):
        """ run the K model if coefficient are loaded and a profile saved"""

        if self.project.myCoeffs.loadCoeffs:
            self.mv.BeginBusy()
            err = self.project.runK(self.RunkNumber)
            self.mv.EndBusy()
            if err == 0:
                self.RunkNumber = self.RunkNumber + 1
                self.write("RunK successful")
                nkViewTitle = kViewTitle + str(self.RunkNumber)
                print("new K View", nkViewTitle)
                self.myViews[nkViewTitle] = rview.kmatrixframe.KMatrixFrame(
                    self.mv, nkViewTitle,
                    self.project.KMatrixFileName,
                    self.RunkNumber, self.project)
                self.myViews[nkViewTitle].Show()
            else:
                logging.info("rttov runK Failed")
                self.mv.WarmError("RunK Failed ")
        else:
            self.mv.WarmError("You must load a Coefficients File ")

    def _RunKPC(self, e):
        """ run the K PC model if coefficient are loaded and a profile saved"""
        if self.project.myCoeffs.loadCoeffs and self.project.isPC():
            if self.project.myProfile['SKIN']['SURFTYPE'] != 1:
                if "landsea" not in self.project.myCoeffs.fileName["PC"]:
                    dialog = self.mv.WarmError(
                        "pc coefficient file not compatible with"
                        " a land or a sea ice surface")
                    return
                    dialog.destroy()
            if self.project.myOption["NPCSCORES"].value < 0:
                self.project.npcscores = 200
            else:
                self.project.npcscores = self.project.myOption[
                    "NPCSCORES"].value
            self.mv.BeginBusy()
            err = self.project.runPCK(self.RunkNumber)
            self.mv.EndBusy()
            if err == 0:
                self.write("Run PC K successful")
                self.RunkNumber += 1
                fileName = self.project.pcKMatrixFileName
                print("fileName:", fileName, "create kpcmatrixframe ")
                if self.project.myOption["ADDRADREC"].value:
                    nkViewTitle = kViewTitle + str(self.RunkNumber)
                    self.myViews[
                        nkViewTitle] = rview.kmatrixframe.KMatrixFrame(
                            self.mv, nkViewTitle,
                            fileName, self.RunkNumber, self.project)
                else:
                    nkViewTitle = kpcViewTitle + str(self.RunkNumber)
                    self.myViews[
                        nkViewTitle] = rview.kpcmatrixframe.KPCMatrixFrame(
                        self.mv, nkViewTitle,
                        fileName, self.project.npcscores, self.RunkNumber)
                    nkpcProfViewTitle = kpcProfViewTitle + str(self.RunkNumber)
                    self.myViews[
                        nkpcProfViewTitle] = rview.kpcView.kpcView(
                        self.mv, nkpcProfViewTitle, None, None,
                        0, 10, fileName, run=self.RunkNumber)
            else:
                logging.info("pc-rttov runK Failed")
                self.mv.WarmError("RunK Failed ")
        else:
            self.mv.WarmError("You must load a Coefficients File ")

    def OpenAsciiProfile(self, e):
        """ open an ascii profile """
        self.prepLog()
        profile_ascii = self.mv.OnOpenFile(e, self.project.config.ENV[
            'RTTOV_GUI_PROFILE_DIR'] + "/../profile-datasets-py")

        if (profile_ascii):
            self.write("Open ascii profile : " + str(profile_ascii))
            err = self.project.openAsciiProfile(profile_ascii)
            if err != 0:
                self.write("Error while opening " + str(profile_ascii))

    def OpenProfile(self, e):
        """ Open the profile , if it contains more than one profile """
        """ ask for the number of the profile """
        self.prepLog()
        self.profileFileName = self.mv.OnOpenFile(
            e, self.project.config.ENV['RTTOV_GUI_PROFILE_DIR'])
        self.write("Open profile filename : " + str(self.profileFileName))
        if (self.profileFileName):
            try:
                number = rttov.profile.getNumberOfProfiles(
                    self.profileFileName)
                self.write("number of profile n this file : " + str(number))
                if (number > 1):
                    number = self.mv.ChooseNumber(
                        number,
                        "select a profile number", "profile number selection")
                self.project.openProfile(self.profileFileName, number)
            except IOError:
                # ne marche pas pb h5py erreur non remontee
                self.write("Error while opening " + str(self.profileFileName))

    def LoadCoefficients(self, e):
        """ open the chooseCoefsController windows """
        self.prepLog()
        self.chooseCoeffsController.OnChoose(e)

    def SelectChannels(self, e):
        """ select the channels """
        self.prepLog()
        self.chooseCoeffsController.OnSelect(e)

    def Ctrl1dvarMenu(self):
        """ control if the items of the r1dvar menu can be enable """
        """ the condition is profileLoaded and coefficient file loaded """
        """ if the status is true assume that this 2 condition are ok """
        """ so just control the number of levels of the profile """
        if self.project.myCoeffs.loadCoeffs and self.project.loadProfile:
            status = self.r1dvarController.controleProject(
                self.project.myProfile, self.project.myCoeffs)
            if status:
                self.mv.EnableMenuItem("config1Dvar")
            else:
                self.mv.DisableMenuItem("config1Dvar")
        else:
            self.mv.DisableMenuItem("config1Dvar")

    def Config1Dvar(self, e):
        """ said to the r1dvar controller to open the r1Dvar control window """
        self.r1dvarController.Show()


class ChooseCoeffsController (GenericController):
    """ controller for the Option windows """

    def __init__(self, parentFrame):
        """ Constructor of the MainController
            define a Project object
            create the OptionView Object
            make the binding with the options view

        """
        self.theParentFrame = parentFrame
        super().__init__()
        self.chooseCoeffsWindow = None
        self.controlerName = "ChooseCoeffsController"

    def OnSelect(self, e):
        """ select the channel from list """
        if not self.project.myCoeffs.loadCoeffs:
            self.write("error : coefficient file are not loaded")
            return
        nbChannels = self.project.myCoeffs.Filenchannels
        if nbChannels is None:
            self.mv.WarmError("Coefficients not loaded ! ")
            return
        if nbChannels < 30:
            dlg = rview.coeff.PickUpChannelsDialog(
                self.theParentFrame, -1, nbChannels, "Select Channels",
                "Select channels")
        else:
            dlg = rview.coeff.SelectChannelsDialog(
                self.theParentFrame, -1, nbChannels, "Select Channels",
                "Select channels")
        if (dlg.ShowModal() == wx.ID_OK):
            channel_list = dlg.GetSelections()
            self.write("channel selection of {} channels".format(
                len(channel_list)))
            dlg.Close()
            self.theParentFrame.BeginBusy()
            err = self.project.loadCoefficients(channel_list)

            if err == 0:
                self.write(
                    "Coefficients successfully loaded with channel selection ")
                self.theParentFrame.ShowInfoMessage(
                    "Coefficients successfully loaded with channel selection ")
            else:
                self.write("Error Coefficients not loaded ")
                self.theParentFrame.WarmError("Can't Load coefficients : \n ")
            self.theParentFrame.EndBusy()

    def OnChoose(self, e):
        """ create the choose coefficients file window and make the binding
            of the items of this
            window with the ChooseCoeffsControlle methods """
        self.chooseCoeffsWindow = rview.coeff.CoeffFilesView(
            self.theParentFrame, self.project.myCoeffs)
        self.chooseCoeffsWindow.Bind(
            wx.EVT_MENU, self.OnApplyCoefficients,
            self.chooseCoeffsWindow.items["applyCoefficients"])
        self.chooseCoeffsWindow.Bind(
            wx.EVT_MENU, self.OnLoadCoefficients,
            self.chooseCoeffsWindow.items["loadCoefficients"])
        self.chooseCoeffsWindow.Bind(
            wx.EVT_BUTTON, self.OnApplyCoefficients,
            self.chooseCoeffsWindow.applyBtn)
        self.chooseCoeffsWindow.Bind(
            wx.EVT_BUTTON, self.OnLoadCoefficients,
            self.chooseCoeffsWindow.loadBtn)
        self.chooseCoeffsWindow.Bind(
            wx.EVT_BUTTON, self.OnDropCoefficients,
            self.chooseCoeffsWindow.dropBtn)

    def OnLoadCoefficients(self, e):
        """ save selected name in project object
            and try to load coefficients """
        self.project.myCoeffs = self.chooseCoeffsWindow.myCoeffs
        if (self.project.myCoeffs.fileName["standard"] == ""):
            self.chooseCoeffsWindow.WarmError(
                "you must choose the coefficients files")
        else:
            self.write("loadCoefficients " +
                       str(self.project.myCoeffs.fileName))
            try:
                self.chooseCoeffsWindow.BeginBusy()
                returnCode = self.project.loadCoefficients()
                if (returnCode == 0):
                    self.write("Coefficients successfully loaded")
                    self.chooseCoeffsWindow.EndBusy()
                    self.chooseCoeffsWindow.ShowInfoMessage(
                        "Coefficients successfully loaded")
                    self.chooseCoeffsWindow.Close()
                else:
                    self.chooseCoeffsWindow.EndBusy()
                    self.chooseCoeffsWindow.WarmError(
                        "Can't Load coefficients : \n \
                    files are not rttov coefficient files or \n \
                    check file names or compatibility issue \n \
                    between standard  clouds  aerosols and pc  ")

            except IOError:
                self.write("Error while opening " +
                           self.project.myCoeffs.fileName)

    def OnApplyCoefficients(self, e):
        """ save selected name in project object """
        self.project.myCoeffs = self.chooseCoeffsWindow.myCoeffs

    def OnDropCoefficients(self, e):
        """ drop the coeeficients """
        ret = self.project.dropCoefficients()
        if ret == 0:
            self.chooseCoeffsWindow.OnClear(e)
            self.project.myCoeffs = self.chooseCoeffsWindow.myCoeffs
        else:
            pass


def check_prerequisite():

    matplotlib_version = matplotlib.__version__
    if wx.VERSION[0] == 2 and wx.VERSION[1] == 8:
        sys.stderr.write(
            'WARNING: wx python version = 2.8.x  some \
            features of RTTOVGUI might not work properly\n')
    if matplotlib_version < "1.3.1":
        sys.stderr.write(
            'WARNING: matplotlib version < 1.3.1 some features \
             of RTTOVGUI might not work properly\n')
    import h5py.version
    if h5py.version.version_tuple[0] < 2:
        sys.stderr.write(
            'WARNING: h5py version < 2.0.1 some features of \
            RTTOVGUI might not work properly\n')


def fxn():
    warnings.warn("deprecated", DeprecationWarning)
    # todo this does not work
    warnings.warn("wxPyDeprecationWarning", DeprecationWarning)
    warnings.warn("MatplotlibDeprecationWarning",
                  DeprecationWarning, stacklevel=1)


def main():
    check_prerequisite()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fxn()
    app = wx.App(False)

    parser = argparse.ArgumentParser(
        description="launch RTTOV GUI", conflict_handler='resolve')
    parser.add_argument('-P', '--profile', dest="profile",
                        help="input file name", required=False,
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument('-A', '--ascii-profile', dest="ascii_profile",
                        help="input ascii file name", required=False,
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))

    parser.add_argument('-s', '--std-coeff', dest="stdcoeff",
                        help="standard coefficient file",
                        required=False,
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-c', '--cloud-coeff', dest="cldcoeff",
                        help="clouds coefficient file",  required=False,
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-m', '--mfasis-lut', dest="mfasislut",
                        help="mfasis lut file",  required=False,
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-a', '--aer-coeff', dest="aercoeff",
                        help="aerosols coefficient file ", required=False,
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument('-p', '--pc-coeff', dest="pccoeff",
                        help="PC coefficient file ",   required=False,
                        metavar="FILE",
                        type=lambda x: is_valid_file(parser, x))
    parser.add_argument("-d", "--debug", dest="debug",
                        action="store_true",
                        help="debug", required=False)

    args = parser.parse_args()

    controller = MainController(app, args.debug)
    python_version = sys.version
    wx_version = wx.__version__
    matplotlib_version = matplotlib.__version__
    logging.info("*******   RTTOVGUI environment       *********************")
    logging.info("Python version :" + python_version)
    logging.info('wxPython version :' + wx_version)
    logging.info('matplotlib version :' + matplotlib_version)
    for k, v in list(controller.project.config.ENV.items()):
        sys.stderr.write(k + ": " + v + "\n")
    sys.stderr.write
    sys.stderr.write

    logging.info(str(controller.project.config.ENV))
    logging.info("*********************************************************")

    sys.stderr.write("\n")
    sys.stderr.write("Start RTTOVGUI ...\n")
    if args.profile is not None:
        controller.profileFileName = args.profile.name
        sys.stderr.write("... load profile " +
                         controller.profileFileName + "\n")
        logging.info("profile: " + controller.profileFileName)
        print("profile:" + controller.profileFileName)
        err = controller.project.openProfile(
            controller.profileFileName, 1)
        if err != 0:
            sys.stderr.write("ERROR wrong profile file \n")
            logging.info("wrong profile names ")
            sys.exit(1)
    else:
        if args.ascii_profile is not None:
            controller.profileFileName = args.ascii_profile.name
            print("profile:" + controller.profileFileName)
            print("open ascii profile")
            sys.stderr.write("... load profile " +
                             controller.profileFileName + "\n")
            err = controller.project.openAsciiProfile(
                controller.profileFileName)
            if err != 0:
                sys.stderr.write("ERROR wrong profile file \n")
                logging.info("wrong profile names ")
                sys.exit(1)

    if args.cldcoeff is not None:
        controller.project.myCoeffs.fileName["clouds"] = args.cldcoeff.name
    if args.aercoeff is not None:
        controller.project.myCoeffs.fileName["aerosols"] = args.aercoeff.name
    if args.pccoeff is not None:
        controller.project.myCoeffs.fileName["PC"] = args.pccoeff.name
    if args.mfasislut is not None:
        controller.project.myCoeffs.fileName["mfasis"] = args.mfasislut.name
    if args.stdcoeff is not None:
        controller.project.myCoeffs.fileName["standard"] = args.stdcoeff.name
        sys.stderr.write("... load coefficient" +
                         str(controller.project.myCoeffs.fileName) + "\n")
        logging.info(" load coefficients")
        logging.info(str(controller.project.myCoeffs.fileName))

        err = controller.project.loadCoefficients()
        if err != 0:
            print("ERROR wrong coefficient names ")
            logging.info("wrong coefficient names ")
            sys.exit(1)
    logging.info("start MainLoop")

    app.MainLoop()


if __name__ == "__main__":
    main()
