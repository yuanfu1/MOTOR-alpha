# -*- coding: utf-8 -*-

from util import GenericController
import rview
import wx
import logging
from pubsub import pub

surfaceViewTitle = "Surface"
BRDFGridViewTitle = "BRDF Editor"
EmisGridViewTitle = "Emissivity Editor"


class SurfaceController (GenericController):
    """ controller for the surface windows """

    def __init__(self, parent):
        """ Constructor of the SurfaceController
            create the SurfaceView Object
            make the binding with the surface view

        """
        self.theParentFrame = parent
        super().__init__()
        # surface view
        self.myViews[surfaceViewTitle] = None
        self.emisController = None
        self.reflController = None
        # subscribe to project changes
        pub.subscribe(self.ProfileLoadedListener, "Profile CHANGED")
        pub.subscribe(self.CoeffsLoadedListener, "Coefficients CHANGED")
        pub.subscribe(self.SurfaceChangedListener, "Surface CHANGED")
        self.controlerName = "SurfaceController"

    def SurfaceChangedListener(self, msg):
        logging.debug("SurfaceChangedListener")
        self.prepLog()
        if self.myViews[surfaceViewTitle]:
            logging.debug("plotReflectanceEmissivity")
            self.myViews[surfaceViewTitle].SetProfileValues(
                self.project.myProfile)
            self.myViews[surfaceViewTitle].PlotReflectanceEmissivity(
                self.project.myEmissivity, self.project.myReflectance,
                self.project.myCoeffs.WAVENUMBERS)
            logging.debug(str(self.project.myReflectance))

    def ProfileLoadedListener(self, msg):
        """ update window of the surface window if it exists """
        self.prepLog()
        logging.debug("SurfaceController ProfileLoadedListener")
        if self.myViews[surfaceViewTitle]:
            self.myViews[surfaceViewTitle].SetProfileValues(
                self.project.myProfile)

        logging.debug("useAtlas :" + str(self.project.useAtlas))
        needOpenAtlas = False
        if self.project.useAtlas:
            for param in self.project.ListeParamAtlas:
                if param == "DATE":
                    for i in range(3):
                        if self.project.paramAtlas[
                                param][i] != self.project.myProfile[param][i]:
                            needOpenAtlas = True
                elif self.project.paramAtlas[
                        param] != self.project.myProfile[param]:
                    needOpenAtlas = True
            self.write("Need Re-Open Atlas : " + str(needOpenAtlas))
            logging.debug("need to re-open Atlas")
            if needOpenAtlas:
                answer = self.theParentFrame.WarmQuestion(
                    "Profile has changed : do you want to re-open Atlas ?")
                if answer == wx.ID_YES:
                    self.project.openAtlas()

    def CoeffsLoadedListener(self, msg):
        """ update menus of the surface window if it exists """
        self.prepLog()
        if self.myViews[surfaceViewTitle]:
            self.UpdateMenus()

    def UpdateMenus(self):
        """ allow all the menus items """
        self.myViews[surfaceViewTitle].EnableMenuItem(
            'createEmissivityReflectance')
        self.myViews[surfaceViewTitle].EnableMenuItem('loadAtlas')
        self.myViews[surfaceViewTitle].EnableMenuItem('saveSurface')
        self.myViews[surfaceViewTitle].EnableMenuItem('saveSurfaceAs')
        self.myViews[surfaceViewTitle].EnableMenuItem('modifyEmissivity')
        self.myViews[surfaceViewTitle].EnableMenuItem('modifyReflectance')

    def UpdateMenusItemModifyEmissivityReflectance(self):
        pass

    def OnSurface(self, e):
        """ create the surface view """
        logging.debug("surfacectrl OnSurface")
        self.prepLog()
        if self.project.loadProfile:
            self.ShowSurface(self.project)

            if self.project.myCoeffs.loadCoeffs:
                self.UpdateMenus()
            if self.project.myEmissivity["EMIS_IN"] is not None:
                self.myViews[surfaceViewTitle].PlotReflectanceEmissivity(
                    self.project.myEmissivity, self.project.myReflectance,
                    self.project.myCoeffs.WAVENUMBERS)
        else:
            self.theParentFrame.WarmError("You must open a Profile ")

    def ShowSurface(self, project):
        logging.debug("surfacectrl showSurface")
        self.prepLog()
        exists = self.ShowWindowIfExists(self.myViews[surfaceViewTitle])
        if not exists:
            self.myViews[surfaceViewTitle] = rview.surface.SurfaceView(
                self.theParentFrame, surfaceViewTitle, project)
            if project.myEmissivity['EMIS_IN'] is not None:
                self.UpdateMenusItemModifyEmissivityReflectance()
            # create the 2 controllers for emissivity and reflectance edition
            self.emisController = ReflEmisController(
                self.myViews[surfaceViewTitle], EmisGridViewTitle, "EMIS")
            self.reflController = ReflEmisController(
                self.myViews[surfaceViewTitle], BRDFGridViewTitle, "REFL")
            self.MakeSurfaceBinding()

    def applyProfileChange(self, e):
        """ apply the change to the profile """
        self.prepLog()
        aProfile = self.myViews[surfaceViewTitle].OnApplyChange(e)
        self.project.updateProfile(aProfile, onlySurfaceParameters=True)

    def openEmissivityReflectance(self, e):
        """ open an EmissityReflectance file """
        self.prepLog()
        surfaceFileName = self.myViews[surfaceViewTitle].OnOpenFile(
            e, self.project.config.ENV['RTTOV_GUI_PREFIX'])
        self.write("filename= " + surfaceFileName)
        if (surfaceFileName):
            try:
                self.project.openSurface(surfaceFileName)
            except IOError:
                self.write("Error while opening " + surfaceFileName)

    def createEmissivityReflectance(self, e):
        """ create EmissivityReflectance values with defaut values """
        self.prepLog()
        logging.debug("surface controller create EmissivityReflectance ")
        self.project.createSurface()

    def loadAtlasEmissivityReflectance(self, e):
        """ Load EmissivityReflectance values from Atlas"""
        self.prepLog()
        okForOpenAtlas = False
        logging.debug("surface controller loadAtlasEmissivityReflectance ")
        if self.project.myProfile['SKIN']['SURFTYPE'] == 1:
            msg = "No atlas over sea !"
            self.myViews[surfaceViewTitle].WarmError(msg)
        else:
            if self.project.myCoeffs.hasSolar():
                msg = "No BRDF data for sea ice surface"
                if self.project.myProfile['SKIN']['SURFTYPE'] == 2:
                    self.myViews[surfaceViewTitle].WarmError(msg)
                else:
                    okForOpenAtlas = True
            else:
                okForOpenAtlas = True
        if okForOpenAtlas:
            number = 1
            typeAtlas = "IR"
            maxNumber = 3
            if self.project.isMW():
                typeAtlas = "NW"
                maxNumber = 2
            number = self.myViews[surfaceViewTitle].ChooseNumber(
                maxNumber,
                "select a version for %s" % typeAtlas, "atlas version")
            self.myViews[surfaceViewTitle].BeginBusy()
            self.project.openAtlas(number, 1)
            self.myViews[surfaceViewTitle].EndBusy()

    def modifyEmissivity(self, e):
        """ modify Emissivity values """
        self.emisController.OnGridView(e)

    def modifyReflectance(self, e):
        """ modify Emissivity values """
        self.reflController.OnGridView(e)

    def saveProfile(self, e):
        """ apply profile change and save profile+option """
        self.prepLog()
        logging.debug("debug saveProfile surface controller")
        # self.project.myProfile is maybe modified elsewhere (profile editor)
        # (but not necessary with the listener)
        profile = self.myViews[surfaceViewTitle].OnApplyChange(
            e, self.project.myProfile)
        self.project.updateProfile(profile, onlySurfaceParameters=True)
        self.SaveProfile(e, self.myViews[surfaceViewTitle])

    def saveProfileAs(self, e):
        """ apply profile change and save profile+option with a new name """
        logging.debug("debug saveProfileAs surface controller")
        # re-initialize the name of the saved profile file to None
        self.project.savedProfileFileName = None
        self.saveProfile(e)

    def saveSurface(self, e):
        """ apply surface file default name """

        logging.debug("debug saveSurface surface controller")
        # self.project.myProfile is maybe modified elsewhere (profile editor)
        if self.project.savedSurfaceFileName is not None:
            self.project.saveSurface(self.project.savedSurfaceFileName)
        else:
            self.saveSurfaceAs(e)

    def saveSurfaceAs(self, e):
        """ save  surface file with a new name """
        self.prepLog()
        logging.debug("debug saveSurfaceAs surface controller")
        # re-initialize the name of the saved profile file to None
        fileName = self.myViews[surfaceViewTitle].OnSaveFileAs(
            e, title="Save Surface File")
        if fileName is not None:
            self.project.saveSurface(fileName)
            self.project.savedSurfaceFileName = fileName

    def MakeSurfaceBinding(self):
        # binding for the SurfaceView
        self.myViews[surfaceViewTitle].Bind(
            wx.EVT_MENU, self.applyProfileChange,
            self.myViews[surfaceViewTitle].items['applyChange'])
        self.myViews[surfaceViewTitle].Bind(
            wx.EVT_MENU, self.createEmissivityReflectance,
            self.myViews[surfaceViewTitle].items[
                'createEmissivityReflectance'])
        self.myViews[surfaceViewTitle].Bind(
            wx.EVT_MENU, self.modifyEmissivity,
            self.myViews[surfaceViewTitle].items['modifyEmissivity'])
        self.myViews[surfaceViewTitle].Bind(
            wx.EVT_MENU, self.modifyReflectance,
            self.myViews[surfaceViewTitle].items['modifyReflectance'])
        self.myViews[surfaceViewTitle].Bind(
            wx.EVT_MENU,
            self.saveProfileAs,
            self.myViews[surfaceViewTitle].items['saveProfileAs'])
        self.myViews[surfaceViewTitle].Bind(
            wx.EVT_MENU,
            self.loadAtlasEmissivityReflectance,
            self.myViews[surfaceViewTitle].items['loadAtlas'])
        self.myViews[surfaceViewTitle].Bind(
            wx.EVT_MENU,
            self.saveProfile,
            self.myViews[surfaceViewTitle].items['saveProfile'])
        self.myViews[surfaceViewTitle].Bind(
            wx.EVT_MENU,
            self.saveSurface,
            self.myViews[surfaceViewTitle].items['saveSurface'])
        self.myViews[surfaceViewTitle].Bind(
            wx.EVT_MENU,
            self.saveSurfaceAs,
            self.myViews[surfaceViewTitle].items['saveSurfaceAs'])
        self.myViews[surfaceViewTitle].Bind(
            wx.EVT_BUTTON,
            self.applyProfileChange,
            self.myViews[surfaceViewTitle].applyBtn)


class ReflEmisController(GenericController):
    """ controller for the Emissivity or Reflectance editor windows """

    def __init__(self, parent, gridViewTitle, dataName):
        """ Constructor of the Emissivity or Reflectance window
            create the GridView Object
            make the binding with the GridView

        """
        self.gridViewTitle = gridViewTitle
        self.dataName = dataName
        self.theParentFrame = parent
        super().__init__()
        self.controlerName = "ReflEmisController"
        self.myViews[self.gridViewTitle] = None
        logging.debug("Create ReflEmisController")
        pub.subscribe(self.SurfaceChangedListener, "Surface CHANGED")

    def SurfaceChangedListener(self, msg):
        logging.debug("ReflEmisController" + self.dataName +
                      "SurfaceChangedListener")
        if self.myViews[self.gridViewTitle] is not None:
            if self.dataName == 'REFL':
                logging.debug(str(self.project.myReflectance))
                logging.debug(
                    "call self.myViews[self.gridViewTitle].OnUpdateData"
                    "(self.project.myReflectance REFL")
                self.myViews[self.gridViewTitle].OnUpdateData(
                    self.project.myReflectance, "REFL")
            else:
                logging.debug(
                    "call self.myViews[self.gridViewTitle].OnUpdateData("
                    "self.project.myReflectance EMIS")
                self.myViews[self.gridViewTitle].OnUpdateData(
                    self.project.myEmissivity, "EMIS")
        else:
            logging.debug("no gridView")

    def OnGridView(self, e):
        exists = self.ShowWindowIfExists(self.myViews[self.gridViewTitle])
        if not exists:
            if self.dataName == "REFL":
                logging.debug("create GridView REFL")
                self.myViews[self.gridViewTitle] = rview.surfedit.GridView(
                    self.theParentFrame, self.gridViewTitle,
                    self.project.myReflectance,
                    "REFL")
                self.MakeBinding()
            else:
                logging.debug("create GridView EMIS")
                self.myViews[self.gridViewTitle] = rview.surfedit.GridView(
                    self.theParentFrame, self.gridViewTitle,
                    self.project.myEmissivity,
                    "EMIS")
                self.MakeBinding()

    def OnApplyChange(self, e):
        logging.debug("controller Emis/refl " +
                      self.dataName + " OnApplyChange")
        self.prepLog()
        data = self.myViews[self.gridViewTitle].OnApplyChange(e)
        if self.dataName == 'REFL':
            self.project.updateReflectance(data)
        else:
            self.project.updateEmissivity(data)

    def MakeBinding(self):
        # binding for the SurfaceView
        self.myViews[self.gridViewTitle].Bind(
            wx.EVT_MENU,
            self.OnApplyChange,
            self.myViews[self.gridViewTitle].items['applyChange'])
