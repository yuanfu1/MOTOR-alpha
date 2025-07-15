# -*- coding: utf-8 -*-
import wx
import rmodel
from rview import util
import wx.lib.agw.floatspin as FS
import os
import logging


class OptionView(util.GenericView):
    """ Option Dialog window of the application """
    helpPage = os.environ["RTTOV_GUI_PREFIX"] + "/doc/helpOptions.html"
    helpTitle = "Options help"

    def __init__(self, parent, title, project):
        self.project = project
        self.myOptions = project.myOption  # linked objects
        self.optionsThemes = self.myOptions.optionsThemes
        self.optionsThemesList = self.myOptions.optionsThemesList
        super().__init__(parent,  title)
        self.CreateMenuBar()

        self.SetMinSize((900, 840))
        self.panel1 = self
        self.masterSizer = wx.BoxSizer(wx.VERTICAL)

        self.cb = {}
        self.combo = {}
        self.boxes = {}
        self.numberParam = {}
        self.boxsizers = {}
        self.gb = wx.GridBagSizer(hgap=5, vgap=5)
        self.numberParamDict = {"RAYLEIGH_MAX_WAVELENGTH": {"step": 0.5,
                                                            "format": "%f",
                                                            "ndigits": 2},
                                "RAYLEIGH_MIN_PRESSURE": {"step": 1,
                                                          "format": "%f",
                                                          "ndigits": 1},
                                "CLDCOL_THRESHOLD": {"step": 0.001,
                                                     "format": "%f",
                                                     "ndigits": 3},
                                "DOM_ACCURACY": {"step": 0.001,
                                                 "format": "%f",
                                                 "ndigits": 3},
                                "DOM_OPDEP_THRESHOLD": {"step": 0.1,
                                                        "format": "%f",
                                                        "ndigits": 2},
                                "DOM_NSTREAMS": {"step": 2,
                                                 "format": "%f",
                                                 "ndigits": 0},
                                "NPCSCORES": {"step": 1,
                                              "format": "%f",
                                              "ndigits": 0}}

        mysizer = {}
        # position of option theme in the GridBagSizer
        posInGridBagSizer = {0:(0,0),
                             4:(0,1),
                             2:(1,0),
                             3:(0,2),
                             1:(1,1),
                             5:(1,2)}
        for i in range(6):
            self.boxes[self.optionsThemes[i]] = wx.StaticBox(
                self.panel1, label=self.optionsThemes[i])
            self.boxsizers[self.optionsThemes[i]] = wx.StaticBoxSizer(
                self.boxes[self.optionsThemes[i]], wx.VERTICAL)

            for parametre in self.optionsThemesList[self.optionsThemes[i]]:
                if not self.myOptions[parametre].hidden:
                    if self.myOptions[parametre].otype == dict:
                        odict = self.myOptions[parametre].odict
                        mysizer[parametre] = wx.BoxSizer(wx.HORIZONTAL)
                        self.combo[parametre] = wx.ComboBox(
                            self.panel1,
                            choices=list(odict.values()),
                            size=(190, 26))
                        mysizer[parametre].Add(self.combo[parametre], border=5)
                        mysizer[parametre].Add(wx.StaticText(
                            self.panel1, -1, parametre.swapcase()), border=5)

                        self.boxsizers[self.optionsThemes[i]].Add(
                            mysizer[parametre],
                            flag=wx.LEFT | wx.TOP | wx.EXPAND,
                            border=5)
                    else:
                        if (parametre in list(self.numberParamDict.keys())):
                            self.numberParam[parametre] = FS.FloatSpin(
                                self.panel1, -1,
                                min_val=self.myOptions[
                                    parametre].min_max_values[0],
                                max_val=self.myOptions[
                                    parametre].min_max_values[1],
                                increment=self.numberParamDict[parametre][
                                    "step"],
                                agwStyle=FS.FS_LEFT)
                            self.numberParam[parametre].SetFormat(
                                self.numberParamDict[parametre]["format"])
                            self.numberParam[parametre].SetDigits(
                                self.numberParamDict[parametre]["ndigits"])
                            mysizer[parametre] = wx.BoxSizer(wx.HORIZONTAL)
                            mysizer[parametre].Add(
                                self.numberParam[parametre], border=5)
                            mysizer[parametre].Add(wx.StaticText(
                                self.panel1, -1,
                                parametre.swapcase()),
                                border=5)
                            self.boxsizers[self.optionsThemes[i]].Add(
                                mysizer[parametre], flag=wx.LEFT | wx.TOP,
                                border=5)
                        else:
                            self.cb[parametre] = wx.CheckBox(
                                self.panel1, label=parametre.swapcase())
                            self.boxsizers[self.optionsThemes[i]].Add(
                                self.cb[parametre], flag=wx.LEFT | wx.TOP,
                                border=5)


            self.gb.Add(self.boxsizers[self.optionsThemes[i]],
                            pos=posInGridBagSizer[i], border=20,
                            flag=wx.ALIGN_CENTER|wx.EXPAND )


        self.masterSizer.Add(self.gb, flag=wx.ALIGN_CENTER, border=10)
        self.masterSizer.Add((10, 10), flag=wx.EXPAND)

        self._CreateButtons()

        self.masterSizer.Add(self.btnSizer, flag=wx.ALIGN_CENTER, border=10)
        self.masterSizer.Add((10, 10), flag=wx.EXPAND)
        self.SetValueItems()
        self.sb = self.CreateStatusBar()

        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.combo['IPCBND'].Bind(wx.EVT_COMBOBOX, self.UpdateOptions)
        self.cb['ADDPC'].Bind(wx.EVT_CHECKBOX, self.UpdateOptions)
        self.cb['DO_LAMBERTIAN'].Bind(wx.EVT_CHECKBOX, self.UpdateOptions)
        self.cb['ADDRADREC'].Bind(wx.EVT_CHECKBOX, self.UpdateOptions)
        self.cb['ADDCLOUDS'].Bind(wx.EVT_CHECKBOX, self.UpdateOptions)
        self.cb['ADDAEROSL'].Bind(wx.EVT_CHECKBOX, self.UpdateOptions)
        self.combo["FASTEM_VERSION"].Bind(wx.EVT_COMBOBOX, self.UpdateOptions)
        self.combo["IR_SCATT_MODEL"].Bind(wx.EVT_COMBOBOX, self.UpdateOptions)
        self.combo["VIS_SCATT_MODEL"].Bind(wx.EVT_COMBOBOX, self.UpdateOptions)

        self.cb['ADDSOLAR'].Bind(wx.EVT_CHECKBOX, self.UpdateOptions)
        self.UpdateOptions(None)

        self.SetSizer(self.masterSizer)
        self.Fit()
        self.Layout()
        self.Centre()
        self.Show(True)

    def UpdateOptions(self, e):
        """ check options values versus status and update the vue """
        self.GetValuesItems()
        self.project.ctrlCoherence()
        self.SetValueItems()

    def _MakeBinding(self):
        """ set the trivial Binding for the View """
        self.Bind(wx.EVT_BUTTON, self.OnCancel, self.cancelBtn)

    def OnItemFocus(self, e):  # TODO A REVOIR
        for (i, value) in list(self.cb.items()):
            if (value == e.GetEventObject()):
                self.sb.PushStatusText(
                    self.myOptions[i].comment, 1)

    def _CreateButtons(self):

        self.btnSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.cancelBtn = wx.Button(self.panel1, wx.ID_CANCEL, label="Revert")
        self.cancelBtn.SetHelpText("Revert to previous options")
        self.btnSizer.Add(
            self.cancelBtn, flag=wx.ALIGN_RIGHT | wx.RIGHT, border=10)

        self.applyBtn = wx.Button(self.panel1, wx.ID_OK, label="Apply")
        self.applyBtn.SetHelpText("Apply options")
        self.applyBtn.SetDefault()
        self.btnSizer.Add(self.applyBtn, flag=wx.ALIGN_RIGHT |
                          wx.RIGHT, border=10)

        # binding cancel button
        self.cancelBtn.Bind(wx.EVT_BUTTON, self.OnCancel)

    def SetValueItems(self):
        for i in self.myOptions.options_list_logical:
            if not self.myOptions[i].hidden:
                self.cb[i].Enable(self.myOptions[i].status)
                self.cb[i].SetValue(self.myOptions[i].value)

        for param in list(self.combo.keys()):
            if not self.myOptions[param].hidden:
                odict = self.myOptions[param].odict
                logging.debug("SetValueItems " + param +
                              str(self.myOptions[param].value) +
                              str(list(odict.values())))

                myval = self.myOptions[param].getValue()
                logging.debug("SetValueItems " + param +
                              str(odict) + " myval " + str(myval))
                self.combo[param].Enable(self.myOptions[param].status)
                self.combo[param].SetItems(list(odict.values()))
                self.combo[param].SetValue(myval)

        for param in list(self.numberParamDict.keys()):
            if not self.myOptions[param].hidden:
                self.numberParam[param].SetValue(self.myOptions[param].value)
                self.numberParam[param].Enable(self.myOptions[param].status)

    def SetOptions(self, options):
        self.myOptions = options
        self.SetValueItems()

    def GetOptions(self):
        return self.myOptions

    def GetValuesItems(self):
        for i in self.myOptions.options_list_logical:
            if not self.myOptions[i].hidden:
                self.myOptions[i].value = self.cb[i].GetValue()
        for i in list(self.combo.keys()):
            if not self.myOptions[i].hidden:
                value = self.combo[i].GetValue()
                odict = self.myOptions[i].odict
                self.myOptions[i].value = [k for k, v in odict.items()
                                           if v == value][0]
                logging.debug("value" + str(value) +
                              " result " + str(self.myOptions[i].value))

        try:
            for parameter in ['CLDCOL_THRESHOLD',
                              "DOM_OPDEP_THRESHOLD",
                              "DOM_ACCURACY",
                              "RAYLEIGH_MAX_WAVELENGTH",
                              "RAYLEIGH_MIN_PRESSURE"]:
                self.myOptions[parameter].value = float(
                    self.numberParam[parameter].GetValue())
            for parameter in ['DOM_NSTREAMS',
                              'NPCSCORES']:
                self.myOptions[parameter].value = int(
                    self.numberParam[parameter].GetValue())

        except ValueError:
            self.ShowErrorMessageDialogBox(
                var=parameter, type='number')

    def ShowErrorMessageDialogBox(self, varName, varType):
        message = "variable " + varName + " must be of " + varType + " type"
        dlg = wx.MessageDialog(
            None, message, caption="Error", style=wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.DeletePendingEvents()
        wx.CallAfter(dlg.Destroy)

    def OnCancel(self, e):
        """ Cancel modifications made in the frame
            set all widget with initial option values"""
        self.SetValueItems()

    def OnApply(self, e):
        """ take value from the windows and save it in options"""
        self.GetValuesItems()

    def OnSave(self, e):
        """ take value from the windows and save it in options"""
        self.GetValuesItems()

    def _initItem(self):
        self.applyItem = 0
        self.saveItem = 0

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu
                ("Apply options", "Apply the options",
                 self.OnApply, "applyOptions", True),
                ("Save options",
                 "Apply and Save the options and the profile in a file",
                 self.OnSave, "saveOptions", True),
                ("", "", "", "", True),
                ('&Quit', 'Quit', self.OnQuit, "quit", True)),
               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True),
                ("&Help", "Help", self.OnHelpHTML, "help", True)))


if __name__ == "__main__":
    p = rmodel.project.Project()
    p.openProfile(p.config.ENV["RTTOV_GUI_PROFILE_DIR"] +
                  '/cldaer101lev_allgas.H5')

    for option in p.myOption.options_list:
        p.myOption[option].status = True

    coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
        "/rttov9pred101L/rtcoef_metop_2_iasi.H5"
    pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
        "/pc/pccoef_metop_2_iasi_landsea_trace_aer.H5"
    p.myCoeffs.fileName["standard"] = coefFile
    p.loadCoefficients()
    ex = wx.App()
    p.myOption.display()
    mv = OptionView(None, "Options", p)

    ex.MainLoop()
