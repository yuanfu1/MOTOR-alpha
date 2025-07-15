# -*- coding: utf-8 -*-

from util import GenericController

import rview
import wx
import logging
from pubsub import pub

profileViewTitle = "Profile Editor"


class ProfileController (GenericController):
    """ controller for the profile windows """

    def __init__(self, parent):
        """ Constructor of the ProfileController
            create the ProfileView Object
            make the binding with the profile view

        """
        self.theParentFrame = parent
        super().__init__()
        # surface view
        self.myViews[profileViewTitle] = None
        self.Redraw = True

        # subscribe to project changes
        pub.subscribe(self.ProfileChangedListener, "Profile CHANGED")
        pub.subscribe(self.SurfaceParametersChangedListener,
                      "Surface parameters CHANGED")
        self.controlerName = "ProfileController"

    def ProfileChangedListener(self, msg):
        if self.myViews[profileViewTitle]:
            if self.Redraw:
                self.myViews[profileViewTitle].RePlotAll(
                    self.project.myProfile)

    def SurfaceParametersChangedListener(self, msg):
        cfraction = self.project.myProfile["CFRACTION"]
        if self.myViews[profileViewTitle] is not None:
            oldcfraction = self.myViews[profileViewTitle].cfraction
            if cfraction != oldcfraction:
                self.myViews[profileViewTitle].UpdateCfraction(cfraction)

    def OnProfile(self, e):
        """ create the profile view """
        if self.project.loadProfile:
            self.ShowProfile(self.project.myProfile)
        else:
            self.theParentFrame.WarmError("You must open a Profile ")

    def ShowProfile(self, profile):
        exists = self.ShowWindowIfExists(self.myViews[profileViewTitle])
        if not exists:
            self.myViews[profileViewTitle] = rview.profileframe.ProfileView(
                self.theParentFrame,
                profileViewTitle,
                profile)
            self.MakeProfileBinding()

    def applyProfileChange(self, e):
        """ apply the change to the profile """
        self.write("apply profile changes")
        aProfile = self.myViews[profileViewTitle].OnApplyChange(e)
        self.Redraw = False
        # we don't want to redraw everything in this case
        # (project.updateProfile send a message "Profile CHANGED")
        self.project.updateProfile(aProfile, onlyVerticalParameters=True)
        self.project.ctrlCoherence()
        self.Redraw = True

    def saveProfile(self, e):
        """ apply profile change and save profile+option """
        logging.debug("debug saveProfile surface controller")
        self.applyProfileChange(e)
        self.SaveProfile(e, self.myViews[profileViewTitle])

    def saveProfileAs(self, e):
        """ apply profile change and save profile+option with a new name """
        logging.debug("debug saveProfileAs surface controller")
        # re-initialize the name of the saved profile file to None
        self.project.savedProfileFileName = None
        self.saveProfile(e)

    def MakeProfileBinding(self):
        """ binding the binding for the profile view """
        self.myViews[profileViewTitle].Bind(
            wx.EVT_MENU, self.applyProfileChange,
            self.myViews[profileViewTitle].items['applyChange'])
        self.myViews[profileViewTitle].Bind(wx.EVT_MENU, self.saveProfileAs,
                                            self.myViews[profileViewTitle].items['saveProfileAs'])
        self.myViews[profileViewTitle].Bind(wx.EVT_MENU, self.saveProfile,
                                            self.myViews[profileViewTitle].items['saveProfile'])
