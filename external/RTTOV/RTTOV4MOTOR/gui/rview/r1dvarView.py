# -*- coding: utf-8 -*-

import wx
import rmodel
from rview import util
import os


class R1dvarView(util.GenericView):
    """ Window for launchind and control the Basic 1Dvar Algorithm       """
    helpMessage = """

        """
    helpPage = os.environ["RTTOV_GUI_PREFIX"] + "/doc/helpR1DVAR.html"

    def __init__(self, parent, title, profileDirName="", nlevels=54):

        super().__init__(parent,  title)
        self.profileDirName = profileDirName
        bgSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.nlevels = nlevels

        self.CreateMenuBar()
        self.SetSize((600, 600))
        self.SetMinSize((600, 500))
        self.SetTitle('1DVAR algorithm')
        self.panel1 = wx.Panel(self, -1)

        self._CreateEntry()
        bgSizer.Add((10, 10))
        bgSizer.Add(self.sizer)
        bgSizer.Add((10, 10))
        self.panel1.SetSizer(bgSizer)
        self.panel1.Fit()
        self.Centre()

        self.Show(True)
        self.myConfig = rmodel.config.Config()

    def _CreateEntry(self):

        self.sizer.Add((10, 10))
        self.sizer.Add(wx.StaticText(
            self.panel1, -1,
            "           1DVAR basic Algorithm control"), border=10)
        self.sizer.Add((10, 10))
        self.sizer.Add(wx.StaticText(
            self.panel1, -1,
            "The first opened profile becomes the background profile"),
            border=50)
        self.sizer.Add(wx.StaticText(
            self.panel1, -1,
            ("You must now select a new profile"
             " which will become the true profile")))
        self.sizer.Add(wx.StaticText(
            self.panel1, -1,
            ("warning : restriction will be made on RTTOV options, "
             "geometry and surface values :")))
        self.sizer.Add(wx.StaticText(self.panel1, -1, "see help for details"))

        self.sizer.Add((10, 10))
        btnsizer1 = wx.BoxSizer(wx.HORIZONTAL)
        self.btnOpen = wx.Button(self.panel1, label=" Open a True Profile ")
        self.btnOpen.SetHelpText("Open a True Profile")
        self.btnOpen.SetDefault()
        btnsizer1.Add(self.btnOpen, 0, wx.GROW |
                      wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.TOP, 5)
        self.sizer.Add(btnsizer1, 0, wx.ALIGN_CENTER_VERTICAL |
                       wx.ALL | wx.EXPAND, 5)
        line = wx.StaticLine(self.panel1, -1, size=(20, -1),
                             style=wx.LI_HORIZONTAL)

        self.sizer.Add(line, 0, wx.GROW |
                       wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.TOP, 5)

        slidSizer = wx.BoxSizer(wx.VERTICAL)
        self.slider = {}
        # sliders in sizer2 :
        slidSizer.Add(wx.StaticText(
            self.panel1, -1,
            "Assumed Observation Error Scaling (x10) : "),
            border=10)

        self.slider["ObsError"] = wx.Slider(
            self.panel1,  -1, 10,  1, 100, (30, 60), (250, -1),
            wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS
        )

        self.slider["ObsError"].SetTickFreq(100)
        slidSizer.Add(self.slider["ObsError"], border=10)

        slidSizer.Add((10, 10))
        slidSizer.Add(wx.StaticText(
            self.panel1, -1,
            "Assumed Background Error Scaling (x10) : "),
            border=10)
        self.slider["BgError"] = wx.Slider(
            self.panel1, -1, 10,  1, 100, (30, 60), (250, -1),
            wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS
        )
        self.slider["BgError"].SetTickFreq(100)
        slidSizer.Add(self.slider["BgError"], border=10)
        slidSizer.Add((10, 10))
        slidSizer.Add(wx.StaticText(
            self.panel1, -1,
            "Noise applied to true BT (x10) : "), border=10)
        self.slider["noise"] = wx.Slider(
            self.panel1, -1, 10,  0, 100, (30, 60), (250, -1),
            wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS
        )
        self.slider["noise"].SetTickFreq(100)
        slidSizer.Add(self.slider["noise"], border=10)

        self.sizer.Add(slidSizer, 0, wx.ALIGN_CENTER_VERTICAL |
                       wx.ALL | wx.EXPAND, 5)

        btnsizer2 = wx.BoxSizer(wx.HORIZONTAL)

        self.btnStep = wx.Button(self.panel1, label="RUN 1DVAR ")
        self.btnStep.SetHelpText("Run the 1DVAR Algorithm")
        self.btnStep.SetDefault()
        btnsizer2.Add(self.btnStep, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)

        self.btnReset = wx.Button(self.panel1, label="Reset")
        self.btnReset.SetHelpText("Reset the 1dvar Project")
        self.btnReset.SetDefault()

        btnsizer2.Add(self.btnReset, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)

        self.sizer.Add(btnsizer2, 0, wx.ALIGN_CENTER_VERTICAL |
                       wx.ALL | wx.EXPAND, 5)

    def GetMaxNoise(self):
        return self.slider["noise"].GetValue() / 10.

    def GetObsError(self):
        return self.slider["ObsError"].GetValue() / 10.

    def GetBgError(self):
        return self.slider["BgError"].GetValue() / 10.

    def ShowErrorMessageDialogBox(self, message):
        """ Show a Error Message Dialog Box with message """
        dlg = wx.MessageDialog(
            None, message, caption="Error", style=wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.DeletePendingEvents()
        wx.CallAfter(dlg.Destroy)

    def OnOpenTrueProfile(self, e):
        """ Select an true profile """
        trueProfileFileName = self.OnOpenFile(
            e, self.profileDirName, "Open a true profile")

        return trueProfileFileName

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu
                ("Open a True Profile", "Open a True Profile",
                 self.OnOpenTrueProfile, "openTrueProfile", True),

                ("", "", "", "", True),
                ('&Quit', 'Quit', self.OnQuit, "quit", True)),

               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True),
                ("&Help", "Help", self.OnHelpHTML, "help", True)))


if __name__ == "__main__":

    pBg = rmodel.project.Project()
    profileDirName = pBg.config.ENV['RTTOV_GUI_PROFILE_DIR']
    ex = wx.App()
    r1dvarv = R1dvarView(None, profileDirName)

    ex.MainLoop()
