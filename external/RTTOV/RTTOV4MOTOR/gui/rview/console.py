# -*- coding: utf-8 -*-
import sys
try:
    import wx
except ImportError:
    sys.stderr.write('ERROR: wxPython is not installed\n')
    sys.exit(1)

from . import util
import os


class MainView(util.GenericView):
    """ MainView of the application """
    helpMessage = """    To run the RTTOV model,
    you have to :
    open a Profile,
    load coefficients files,
    modify options, profile and surface parameters
    with the different windows if necessary.
    """

    def __init__(self, parent, logfile=None):
        util.GenericView.__init__(self, parent, "Main View")
        sizer = wx.GridBagSizer(2, 2)

        # make the menu bar
        self.CreateMenuBar()

        # trivial binding
        self._MakeBinding()
        windowSize = (700, 600)
        self.SetSize(windowSize)
        self.SetTitle('RTTOV GUI')
        self.sb = self.CreateStatusBar()
        self.panel = wx.Panel(self, -1)

        self.panel.SetSizer(sizer)
        self.log = wx.TextCtrl(self.panel, -1, style=wx.TE_MULTILINE |
                               wx.TE_READONLY | wx.TE_WORDWRAP | wx.TE_RICH)
        self.log.SetEditable(False)

        self.log.SetBackgroundColour("Light Grey")

        sizer.Add(self.log, pos=(0, 0), flag=wx.EXPAND, border=50)
        sizer.AddGrowableCol(0)
        sizer.AddGrowableRow(0)
        # redirect output in self.log
        redir = RedirectText(self.log)
        if logfile is not None:
            try:
                self.logFile = open(logfile, 'w')
            except Exception:
                print("error open logfile")
        sys.stdout = redir
        # but redirect stderr is not wise ...

        self.Centre()
        self.Show(True)

    def Go2End(self):
        self.log.SetInsertionPointEnd()

    def OnConfig1Dvar(self):
        pass

    def OnHelp1Dvar(self, e):
        helpPage = os.environ["RTTOV_GUI_PREFIX"] + "/doc/helpR1DVAR.html"
        HelpWindow = util.helpframe.HelpWindow(
            self, -1, "Help 1Dvar ", helpPage)
        HelpWindow.Show()

    def changeCursor(self, status):
        """"""
        myCursor = wx.StockCursor(status)
        self.panel.SetCursor(myCursor)

    def OnClose(self, e):
        self.DeletePendingEvents()
        self.Destroy()

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu

                ("Open Profile", "Open a profile",
                 self.OnOpenProfile, "openProfile", True),
                ("Open Ascii Profile", "Open a profile",
                 self.OnOpenAsciiProfile, "openAsciiProfile", True),
                ("Save Profile as", "Save the profile",
                 self.OnSaveProfileAs, "saveProfileAs", True),
                ("Save Profile", "Save the profile",
                 self.OnSaveProfile, "saveProfile", True),
                ("", "", "", "", True),
                ('Quit', 'Quit application', self.OnQuit, "quit", True)),

               ("&Windows",
                ('Profile Editor Window', 'Open the profile editor window',
                 self.OnProfileWindow, "profileWindow", False),
                ('Options Editor Window', 'Open the option editor window',
                 self.OnOptionsWindow, "optionsWindow", False),
                ('Surface Editor Window', 'Open the surface editor window',
                 self.OnSurfaceWindow, "surfaceWindow", False)),
               ("&Rttov",
                ('Load Coefficients', "Load the coefficients files",
                 self.OnLoadCoefficients, "loadCoefficients", True),
                ('Select Channels',
                 "Select channels from the coefficient files",
                 self.OnSelectChannels, "selectChannels", False),
                ('Run RTTOV direct', "run RTTOV direct model",
                 self.OnRunDirect, "runDirect", False),
                ('Run RTTOV K', "run RTTOV jacobian model", self.OnRunK,
                 "runK", False)),
               ("&1Dvar",
                ('Configure 1Dvar', "Configure the basic 1Dvar algorithm",
                 self.OnConfig1Dvar, "config1Dvar", False),
                ('Help 1Dvar', "Help about the basic 1Dvar algorithm",
                 self.OnHelp1Dvar, "help1Dvar", True)),
               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True),
                ("&Help", "Help", self.OnHelp, "help", True)))

    def OnNewProject(self, e): pass

    def OnOpenProject(self, e): pass

    def OnSaveProjectAs(self, e): pass

    def OnSaveProject(self, e): pass

    def OnCreateNewProfile(self, e): pass

    def OnOpenProfile(self, e):
        print("OnOpenProfile MainView")
        print(self.OnOpenFile(e, title="Open a profile"))

    def OnOpenAsciiProfile(self, e):
        print("OnOpenProfile MainView")
        print(self.OnOpenFile(e, title="Open an Ascii profile"))

    def OnSaveProfileAs(self, e):
        print(self.OnSaveFileAs(e))

    def OnSaveProfile(self, e):
        print(self.OnSaveFileAs(e))

    def OnProfileWindow(self, e): pass

    def OnOptionsWindow(self, e): pass

    def OnSurfaceWindow(self, e): pass

    def OnLoadCoefficients(self, e): pass

    def OnSelectChannels(self, e): pass

    def OnRunDirect(self, e): pass

    def OnRunK(self, e): pass

    def _MakeBinding(self):
        """ set the trivial Binding for the View """
        pass


class RedirectText:

    def __init__(self, aWxTextCtrl):
        self.out = aWxTextCtrl

    def write(self, string):
        self.out.WriteText(string)


if __name__ == "__main__":

    ex = wx.App()
    mv = MainView(None, "essai.log")
    # binding must be in the controler
    # here is just some code to test the look of the view

    for i in range(140):
        print("essai essai                                                 "
              "                       essai   essai essai essai essai essai")

    ex.MainLoop()
