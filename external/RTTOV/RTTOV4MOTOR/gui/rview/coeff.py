'''
Created on Jun 20, 2012

@author: pascale
'''
# -*- coding: utf-8 -*-


import wx
import rmodel
from rview import util
import os


class CoeffFilesView(util.GenericView):
    """ Option Dialog window of the application """
    helpMessage = """
Choose your coefficient files :
RTTOV coefficient file is required (54 levels if possible)
optional :
Aerosols coefficient file
Clouds coefficient file
PC coefficient file (only for airs, iasi or iasing)
Click the "Load" button to load the coefficients files.

        """
    dictTypeCoeffFiles = {0: "standard", 1: "aerosols",
                          2: "clouds", 3: "mfasis", 4: "PC"}
    dictTextCoeffFiles = {"standard": " RTTOV coefficient file :",
                          "aerosols": " Aerosols coefficient file : ",
                          "clouds": " Clouds coefficient file :",
                          "mfasis": " Cloud LUT file for mfasis :",
                          "PC": " PC coefficient file :"}

    def __init__(self, parent, coeffs):
        self.myCoeffs = coeffs
        util.GenericView.__init__(self, parent,  "Coefficients Files")

        self.sizer = wx.BoxSizer(wx.VERTICAL)

        self.CreateMenuBar()
        self.SetSize((600, 500))
        self.SetMinSize((600, 500))
        self.SetTitle('Coefficients Files')

        self.panel1 = self
        self.sizer.Add((10, 10))
        self.sizer.Add(wx.StaticText(
            self.panel1, -1,
            "Choose the coefficient files and load them"))
        self.sizer.Add((10, 10))
        self._CreateEntry()
        self._CreateButtons()
        self.SetAutoLayout(True)
        self.SetSizer(self.sizer)
        self.Layout()
        self.Fit()
        self.Centre()
        self.MakeModal(modal=True)
        self.Show(True)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.myConfig = rmodel.config.Config()

    def _MakeBinding(self):
        """ set the trivial Binding for the View """
        pass

    def OnChooseFile(self, e):
        name = self.OnOpenFile(e, directory=self.myConfig.ENV[
                               "RTTOV_GUI_COEFF_DIR"], title="Choose a file")
        if (name is not None):
            for i in range(5):
                if (self.BtnList[i] == e.GetEventObject()):
                    self.myCoeffs.fileName[self.dictTypeCoeffFiles[i]] = name
                    self.TextList[i].SetValue(os.path.basename(name))
        return self.myCoeffs

    def OnClearFile(self, e):
        for i in range(5):
            if (self.BtnClearList[i] == e.GetEventObject()):
                self.TextList[i].SetValue("")
                self.myCoeffs.fileName[self.dictTypeCoeffFiles[i]] = ""
        return self.myCoeffs

    def _CreateEntry(self):

        self.BtnList = []
        self.TextList = []
        self.BtnClearList = []
        for coeff in ("standard", "aerosols", "clouds", "mfasis", "PC"):
            self._CreateAFileChooser(self.panel1, self.dictTextCoeffFiles[
                                     coeff], "Choose...", "help",
                                     self.OnChooseFile,
                                     ctrlText=self.myCoeffs.fileName[coeff])
        self.sizer.Add((10, 10))

    def _CreateAFileChooser(self, parent, theText, theLabel,
                            helpText, theBinding, ctrlText=""):
        """ create a button at the position associated
        at a text  and make the binding
        arguments :
        parent : the parent Frame
        theText : explanation text for the file
        label of the button
        help text
        action
        ctrlText initialisation """
        self.sizer.Add((10, 10))
        self.sizer.Add(wx.StaticText(parent, -1, label=theText))

        mysizer = wx.BoxSizer(wx.HORIZONTAL)
        mysizer.Add((5, 5))

        text = wx.TextCtrl(parent, -1, style=wx.TE_READONLY, size=(400, -1))

        mysizer.Add(text, flag=wx.GROW)
        mysizer.Add((5, 5))
        btn = wx.Button(parent, -1, label=theLabel)
        mysizer.Add(btn)
        self.BtnList.append(btn)
        text.SetValue(os.path.basename(ctrlText))
        self.TextList.append(text)
        btn.Bind(wx.EVT_BUTTON, theBinding)
        self.sizer.Add((10, 10))
        btnClear = wx.Button(parent, -1, label="Clear")
        self.BtnClearList.append(btnClear)
        btnClear.Bind(wx.EVT_BUTTON, self.OnClearFile)
        mysizer.Add((5, 5))
        mysizer.Add(btnClear)

        self.sizer.Add(mysizer, flag=wx.GROW)

    def _CreateButtons(self):

        line = wx.StaticLine(self.panel1, -1, size=(20, -1),
                             style=wx.LI_HORIZONTAL)
        self.sizer.Add(line, 0, wx.GROW |
                       wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.TOP, 5)
        btnsizer = wx.BoxSizer(wx.HORIZONTAL)

        self.loadBtn = wx.Button(self.panel1, label="Load")
        self.loadBtn.SetHelpText("Load rttov with coefficients Files")
        self.loadBtn.SetDefault()
        btnsizer.Add(self.loadBtn, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)

        self.applyBtn = wx.Button(self.panel1, label="Apply")
        self.applyBtn.SetHelpText(
            "Apply the choices (but not load the coefficient files)")
        btnsizer.Add(self.applyBtn, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)

        self.dropBtn = wx.Button(self.panel1, label="Drop")
        self.dropBtn.SetHelpText("Drop coefficient Files")
        btnsizer.Add(self.dropBtn, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)

        self.sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL |
                       wx.ALL | wx.EXPAND, 5)

    def ShowErrorMessageDialogBox(self, message):
        """ Show a Error Message Dialog Box with message """
        dlg = wx.MessageDialog(
            None, message, caption="Error", style=wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.DeletePendingEvents()
        wx.CallAfter(dlg.Destroy)

    def OnApply(self, e):
        """ keep value from the windows and save it in coeffs """
        pass

    def OnLoad(self, e):
        """ take value from the windows and save it in options"""

        self.BusyCursor()

    def OnClear(self, e):
        """ Clear all"""
        for i in range(5):
            self.TextList[i].SetValue("")
            self.myCoeffs.fileName[self.dictTypeCoeffFiles[i]] = ""
        return self.myCoeffs

    def OnClose(self, e):
        """ close the option windows"""
        self.MakeModal(modal=False)
        self.DeletePendingEvents()
        wx.CallAfter(self.Destroy)

    def _initItem(self):
        self.applyItem = 0
        self.saveItem = 0

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu
                ("Apply coefficients", "Apply the coefficients",
                 self.OnApply, "applyCoefficients", True),
                ("Load coefficients", "Load the coefficients",
                 self.OnLoad, "loadCoefficients", True),
                ("clear coefficients names", "Clear the coefficients names",
                 self.OnClear, "clearCoefficients", True),
                ("", "", "", "", True),
                ('&Quit', 'Quit', self.OnQuit, "quit", True)),

               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True),
                ("&Help", "Help", self.OnHelp, "help", True)))


class SelectChannelsDialog(wx.Dialog):

    def __init__(
            self, parent, ID,  nbChannels=10, title="", text="",
            size=wx.DefaultSize, pos=wx.DefaultPosition,
            style=wx.DEFAULT_DIALOG_STYLE,
            useMetal=False, limits_scale=(0.01, 10), limits_offset=(-100, 100)
    ):

        # Instead of calling wx.Dialog.__init__ we precreate the dialog
        # so we can set an extra style that must be set before
        # creation, and then we create the GUI object using the Create
        # method.

        if wx.VERSION[0] <= 3:
            pre = wx.PreDialog()
            pre.SetExtraStyle(wx.DIALOG_EX_CONTEXTHELP)
            pre.Create(parent, ID, title, pos, size, style)

            # This next step is the most important, it turns this Python
            # object into the real wrapper of the dialog (instead of pre)
            # as far as the wxPython extension is concerned.
            self.PostCreate(pre)
        else:
            wx.Dialog.__init__(self)
            self.SetExtraStyle(wx.DIALOG_EX_CONTEXTHELP)
            self.Create(parent, ID, title)

        # This extra style can be set after the UI object has been created.
        if 'wxMac' in wx.PlatformInfo and useMetal:
            self.SetExtraStyle(wx.DIALOG_EX_METAL)

        # Now continue with the normal construction of the dialog
        self.nbChannels = nbChannels

        # contents
        sizer = wx.BoxSizer(wx.VERTICAL)

        label = wx.StaticText(self, -1, text)
        label.SetHelpText("help text")
        sizer.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)

        # box 1 cmin
        box = wx.BoxSizer(wx.HORIZONTAL)

        label = wx.StaticText(self, -1, "first channel:")
        box.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.cmin = wx.Slider(self, -1, 1, 1, self.nbChannels, (30, 60),
                              (400, -1),
                              wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS
                              )
        self.cmin.SetTickFreq(5)

        box.Add(self.cmin, 1, wx.ALIGN_CENTRE | wx.ALL, 5)
        sizer.Add(box, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.cmin.Bind(wx.EVT_SCROLL_CHANGED, self.OnChangedCmin)

        # box 2 cmax
        box = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(self, -1, "last channel :")
        box.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.cmax = wx.Slider(self, -1, self.nbChannels, 1, self.nbChannels,
                              (30, 60), (400, -1),
                              wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS
                              )
        box.Add(self.cmax, 1, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.cmin.SetTickFreq(5)
        sizer.Add(box, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        self.cmax.Bind(wx.EVT_SCROLL_CHANGED, self.OnChangedCmax)

        # box 3 thinning factor
        box = wx.BoxSizer(wx.HORIZONTAL)
        label = wx.StaticText(self, -1, "thinning factor:")
        box.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.thinning_factor = wx.Slider(self, -1, 1, 1, 30, (30, 60),
                                         (400, -1),
                                         wx.SL_HORIZONTAL | wx.SL_AUTOTICKS |
                                         wx.SL_LABELS
                                         )
        box.Add(self.thinning_factor, 1, wx.ALIGN_CENTRE | wx.ALL, 5)
        sizer.Add(box, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20, -1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL |
                  wx.RIGHT | wx.TOP, 5)

        # buttons
        btnsizer = wx.StdDialogButtonSizer()

        btn = wx.Button(self, wx.ID_OK)

        btn.SetDefault()
        btnsizer.AddButton(btn)
        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)
        self.Show()

    def OnChangedCmin(self, e):
        self.cmax.SetMin(self.cmin.GetValue())

    def OnChangedCmax(self, e):
        self.cmin.SetMax(self.cmax.GetValue())

    def computeChannelList(self, cmin, cmax, thinning_factor):
        if cmax == cmin:
            my_liste = [cmin]
        else:
            if thinning_factor > abs(cmax - cmin):
                my_liste = [min(cmin, cmax), max(cmin, cmax)]
            else:
                if cmin <= cmax:
                    my_liste = [x for x in range(
                        cmin, cmax + 1, thinning_factor)]
                else:
                    my_liste = [x for x in range(
                        cmax, cmin + 1, thinning_factor)]
        return my_liste

    def GetSelections(self):
        print(("channel selection cmin=", self.cmin.GetValue()))
        print(("channel selection cmax=", self.cmax.GetValue()))
        print(("channel selection thinning_factor:",
               self.thinning_factor.GetValue()))
        myList = self.computeChannelList(
            self.cmin.GetValue(), self.cmax.GetValue(),
            self.thinning_factor.GetValue())
        return myList


class PickUpChannelsDialog(wx.Dialog):
    """ Dialog window which allows to select the channels"""

    def __init__(
            self, parent, ID,  nbChannels, title="", text="",
            size=wx.DefaultSize, pos=wx.DefaultPosition,
            style=wx.DEFAULT_DIALOG_STYLE,
            useMetal=False, limits_scale=(0.01, 10), limits_offset=(-100, 100)
    ):

        if wx.VERSION[0] <= 3:
            # Instead of calling wx.Dialog.__init__ we precreate the dialog
            # so we can set an extra style that must be set before
            # creation, and then we create the GUI object using the Create
            # method.
            pre = wx.PreDialog()
            pre.SetExtraStyle(wx.DIALOG_EX_CONTEXTHELP)
            pre.Create(parent, ID, title, pos, size, style)

            # This next step is the most important, it turns this Python
            # object into the real wrapper of the dialog (instead of pre)
            # as far as the wxPython extension is concerned.
            self.PostCreate(pre)
        else:
            wx.Dialog.__init__(self)
            self.SetExtraStyle(wx.DIALOG_EX_CONTEXTHELP)
            self.Create(parent, ID, title)

        # This extra style can be set after the UI object has been created.
        if 'wxMac' in wx.PlatformInfo and useMetal:
            self.SetExtraStyle(wx.DIALOG_EX_METAL)

        # Now continue with the normal construction of the dialog

        self.nbChannels = nbChannels

        # contents
        sizer = wx.BoxSizer(wx.VERTICAL)

        label = wx.StaticText(self, -1, text)
        label.SetHelpText("help text")
        sizer.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.cb = {}
        for c in range(1, self.nbChannels + 1):
            self.cb[c] = wx.CheckBox(self, -1, "channel {}".format(c))
            sizer.Add(self.cb[c], 0, wx.GROW |
                      wx.ALIGN_CENTER_VERTICAL | wx.RIGHT | wx.TOP, 5)

        line = wx.StaticLine(self, -1, size=(20, -1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL |
                  wx.RIGHT | wx.TOP, 5)

        # buttons
        btnsizer = wx.StdDialogButtonSizer()
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)
        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)
        self.Show()

    def GetSelections(self):

        return [c for c in range(1, self.nbChannels + 1) if
                self.cb[c].IsChecked()]


if __name__ == "__main__":
    p = rmodel.project.Project()
    p.openProfile(
        p.config.ENV["RTTOV_GUI_PROFILE_DIR"] + '/cldaer101lev_allgas.H5')
    ex = wx.App()
    select_channel = False
    if select_channel:
        sc = PickUpChannelsDialog(
            None, -1, 20, "Select Channels", "Please select the channels")
    else:
        lc = CoeffFilesView(None, p.myCoeffs)
    ex.MainLoop()
