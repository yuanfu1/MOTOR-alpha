# -*- coding: utf-8 -*-

import util
import wx
import rview
from pubsub import pub

OptionViewTitle = "Options"


class OptionController (util.GenericController):
    """ controller for the Option windows """

    def __init__(self, parent):
        """ Constructor of the MainController
            define a Project object
            create the OptionView Object
            make the binding with the options view

        """
        self.theParentFrame = parent
        super().__init__()
        self.myViews[OptionViewTitle] = None
        pub.subscribe(self.OptionsChangedListener, "Options CHANGED")

    def OptionsChangedListener(self, msg):
        """ update the Options window if it exist """
        if self.myViews[OptionViewTitle] is not None:
            self.myViews[OptionViewTitle].SetOptions(self.project.myOption)

    def MakeBinding(self):
        """  Define the different binding of the application  """
        # binding for the MainView (mv)
        self.mv.Bind(wx.EVT_MENU, self.OnOptions,
                     self.mv.items["optionsWindow"])

    def _MakeBindingInOptionView(self):

        self.myViews[OptionViewTitle].Bind(wx.EVT_MENU, self.applyOptions,
                                           self.myViews[OptionViewTitle].items[
                                               'applyOptions'])
        self.myViews[OptionViewTitle].Bind(
            wx.EVT_BUTTON, self.applyOptions,
            self.myViews[OptionViewTitle].applyBtn)

    def OnOptions(self, e):

        if self.project.loadProfile:
            self.project.ctrlCoherence()
            self.ShowOptions(self.project)
        else:
            self.theParentFrame.WarmError(
                "You must open a Profile or create a new Profile ")

    def applyOptions(self, e):
        """ apply the options """
        self.myViews[OptionViewTitle].OnApply(e)
        self.project.ctrlCoherence()
        self.project.updateOptions(self.myViews[OptionViewTitle].myOptions)
        self.write("options applied")

    def ShowOptions(self, project):
        exists = self.ShowWindowIfExists(self.myViews[OptionViewTitle])
        if not exists:
            self.myViews[OptionViewTitle] = rview.option.OptionView(
                self.theParentFrame, OptionViewTitle, project)
            self._MakeBindingInOptionView()
