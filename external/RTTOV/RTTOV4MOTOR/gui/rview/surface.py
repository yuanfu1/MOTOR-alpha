# -*- coding: utf-8 -*-

try:
    import wx.lib.masked
    import wx.lib.agw.floatspin as FS
except ImportError:
    import sys
    sys.stderr.write('ERROR: wxPython is not installed\n')
    sys.exit(1)

import rmodel
from rview import util, colors
import collections
import locale
locale.setlocale(locale.LC_ALL, '')

try:
    from matplotlib.backends.backend_wxagg import \
        FigureCanvasWxAgg as FigureCanvas
    from matplotlib.backends.backend_wxagg import \
        NavigationToolbar2WxAgg as ToolBar
    from matplotlib.figure import Figure
except ImportError:
    import sys
    sys.stderr.write('ERROR: matplotlib is not installed\n')
    sys.exit(1)
try:
    import numpy
except ImportError:
    import sys
    sys.stderr.write('ERROR: numpy is not installed\n')
    sys.exit(1)

if wx.VERSION[0] > 3:
    import wx.adv


class PlotPage(wx.Panel):

    def __init__(self, parent, theValues, theName):
        wx.Panel.__init__(self, parent, style=wx.BORDER_SIMPLE)
        self.theName = theName
        self.theValues = theValues
        self.fig = Figure()

        self.canvas = FigureCanvas(self, -1, self.fig)
        self.canvas.mpl_connect('motion_notify_event', self.onMouseMotion)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        tlb = ToolBar(self.canvas)
        tlb.Realize()
        self.sizer.Add(tlb, 0, wx.GROW)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW, 1)
        self.SetSizer(self.sizer)

        self.text = wx.StaticText(self, -1, label="")
        self.sizer.Add(self.text)
        self.Fit()

    def onMouseMotion(self, event):
        # if wn < 400     : lambda = xxx micrometre
        #        else          : nu = yyy GHz

        txt = ""
        spl = 29979245800.0  # en cm
        if event.inaxes:
            x = int(round(event.xdata))
            if x < 1:
                x = 1
            if x > self.nbChannels:
                x = self.nbChannels
            y_in = self.theValuesIn[x - 1]
            y_out = self.theValuesOut[x - 1]
            if self.WAVENUMBERS is not None:

                chan = self.WAVENUMBERS[x - 1]
                if chan >= 400.0:
                    txt0 = '%s=%.3f %sm]' % (
                        u'\u03BB', (10000 / chan), u'\u03BC')
                else:
                    txt0 = '%s=%.4f GHz]' % (
                        u'\u03BD', (spl * chan / 1000000000))

                labelchan = txt0

            else:
                labelchan = ""

            txt = 'Channel #%02g (%s ) %s: in=%.3f out=%.3f' % (
                x, labelchan, self.theName, y_in, y_out)
            self.text.SetLabel(txt)

    def OnPlot(self, theValuesIn, theValueOut, WAVENUMBERS=None):
        self.WAVENUMBERS = WAVENUMBERS
        self.theValuesIn = theValuesIn
        self.theValuesOut = theValueOut
        from matplotlib.ticker import MultipleLocator, FormatStrFormatter
        if len(theValuesIn.shape) == 0:
            nbChannels = 1
        else:
            nbChannels = theValuesIn.shape[0]
        theMarker = None
        self.nbChannels = nbChannels
        if nbChannels < 10:
            majorLocator = MultipleLocator(1)
            minorLocator = MultipleLocator(1)
        majorFormatter = FormatStrFormatter('%d')
        self.fig.clear()
        self.subP = self.fig.add_subplot(1, 1, 1)
        channels = numpy.arange(1, nbChannels + 1)
        if nbChannels < 50:
            theMarker = '+'

        self.subP.plot(channels, theValuesIn,
                       color=colors.surfaceIn, label="in", marker=theMarker)
        self.subP.set_ylabel("Y Axis (" + self.theName + ")")
        self.subP.set_xlabel("X Axis (channels)")
        self.subP.set_ylim((-0.05, 1.05))
        self.subP.set_xlim((0, nbChannels + 1))
        if nbChannels < 10:
            majorLocator = MultipleLocator(2)
            minorLocator = MultipleLocator(1)
            self.subP.xaxis.set_major_locator(majorLocator)
            self.subP.xaxis.set_minor_locator(minorLocator)
        self.subP.xaxis.set_major_formatter(majorFormatter)
        self.subP.xaxis.set_minor_formatter(majorFormatter)
        if theValueOut is not None:
            self.subP.plot(channels, theValueOut,
                           color=colors.surfaceOut, label="out",
                           marker=theMarker)
        self.subP.legend(prop={'size': 10}, shadow=True,
                         fancybox=True, loc='best')

        self.subP.grid(True)
        self.fig.canvas.draw()


class SurfaceView(util.GenericView):
    """ Surface window of the application """
    helpMessage = """

        This window allows you to modify the surface parameters
        of the profile.
        You can also load an Atlas or modify the values of
        reflectance or emissivity of the surface
        if a coefficients file is loaded.
        Do not forget to apply your changes.
        """

    valFastem = collections.OrderedDict()

    profileParamroList = ['NLAYERS', 'NLEVELS']

    def initValFastem(self):
        # Surface type 1 2 3 4 5
        # default rttov values
        self.valFastem['Default'] = [3., 5., 15., 0.1, 0.30000001]
        # Summer land surface
        self.valFastem['Forest'] = [1.7, 1.0, 163.0, 0.0, 0.5]
        self.valFastem['Open Grass'] = [2.2, 1.3, 138.0, 0.0, 0.42]
        self.valFastem['Bare soil'] = [2.3, 1.9, 21.8, 0.0, 0.5]
        # Winter land surface
        self.valFastem['Forest and snow'] = [2.9, 3.4, 27.0, 0.0, 0.0]
        self.valFastem['Deep dry snow'] = [3.0, 24.0, 60.0, 0.1, 0.15]
        self.valFastem['Frozen soil'] = [117.8, 2.0, 0.19, 0.2, 0.35]
        # Sea ice
        self.valFastem['Grease ice'] = [23.7, 7.7, 17.3, 0.0, 0.15]
        self.valFastem['Baltic nilas'] = [1.6, 3.3, 2.2, 0.0, 0.0]
        self.valFastem['New ice (no snow)'] = [2.9, 3.4, 27.0, 0.0, 0.0]
        self.valFastem['New ice (snow)'] = [2.2, 3.7, 122.0, 0.0, 0.15]
        self.valFastem['Brash ice'] = [3.0, 5.5, 183.0, 0.0, 0.0]
        self.valFastem['Compact pack ice'] = [2.0, 1700000, 49000000, 0.0, 0.0]
        self.valFastem['Fast ice'] = [1.5, 77.8, 703, 0.1, 0.35]
        self.valFastem['Lake ice + snow'] = [1.8, 67.1, 534, 0.1, 0.15]
        self.valFastem['Multi-year ice'] = [1.5, 85000, 4700000, 0.0, 0.0]
        # Table 3. Coefficients for FASTEM-2 for different surface types
        # (adapted from English and Hewison, 1998).
        # English S.J. and T.J. Hewison 1998 A fast generic millimetre
        # wave emissivity model. Microwave Remote Sensing of the
        # Atmosphere and Environment Proc. SPIE 3503 22-30 '

    def __init__(self, parent, title, project):
        super().__init__(parent, title)
        self.myProfile = project.myProfile
        self.debug = False
        self.initValFastem()

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        self.CreateMenuBar()
        self.SetSize((1020, 750))
        self.SetMinSize((1020, 750))
        self.panel1 = wx.Panel(
            self, -1, style=wx.BORDER_SIMPLE, size=(200, 500))
        sizer.Add(self.panel1, 1, wx.EXPAND)
        self.gbs = wx.GridBagSizer(5, 4)
        self.panel1.SetSizer(self.gbs)
        self.sskin_list = []
        for label in self.myProfile.sskin_list:
            if "SKIN_" + label not in self.myProfile.notUsedParams:
                self.sskin_list.append(label)
        self.s2m_list = []
        for label in self.myProfile.s2m_list:
            if "S2M_" + label not in self.myProfile.notUsedParams:
                self.s2m_list.append(label)
        self.ShowSurfaceParameters()
        self._MakeBinding()
        self.panel2 = wx.Panel(self, -1, style=wx.BORDER_SIMPLE)
        sizer.Add(self.panel2, 1, wx.EXPAND)
        # creation of the notebook for the panel 2
        nb = wx.Notebook(self.panel2, -1)
        # creation of the 2 pages
        self.emissivityPageGraphic = PlotPage(
            nb, project.myEmissivity['EMIS_IN'], theName="Emissivity")
        self.reflectancePageGraphic = PlotPage(
            nb, project.myReflectance['REFL_IN'], theName="Surface BRDF")
        # add the page to the notebook
        nb.AddPage(self.emissivityPageGraphic, "Emissivity")
        nb.AddPage(self.reflectancePageGraphic, "Surface BRDF")
        # create a second sizer for the notebook
        sizer2 = wx.BoxSizer()
        sizer2.Add(nb, 1, wx.EXPAND)
        self.panel2.SetSizer(sizer2)
        self.sb = self.CreateStatusBar()
        self.Centre()
        self.panel1.Fit()
        self.panel2.Fit()
        self.Fit()
        self.Show(True)

    def _MakeBinding(self):
        """ set the trivial Binding for the View """
        # binding revert button
        self.revertBtn.Bind(wx.EVT_BUTTON, self.OnRevertChange)
        # binding on ApplyButton will be made by the surface controler

    def ShowSurfaceParameters(self):
        """ Show surface parameters of the profile in the panel1
            control limits have been extracted of rttov_const.F90 """

        self.ct = {}
        self.gbs.Add(wx.StaticText(self.panel1, -1, 'ID'),
                     pos=(0, 0), flag=wx.ALL | wx.ALIGN_RIGHT, border=4)
        self.ct['ID'] = wx.TextCtrl(
            self.panel1, -1, self.myProfile['ID'], size=(350, 25))
        self.gbs.Add(self.ct['ID'], pos=(0, 1), span=(
            1, 3), flag=wx.ALIGN_LEFT | wx.TOP, border=8)

        self.gbs.Add(wx.StaticText(self.panel1, -1, 'DATE'),
                     pos=(1, 0), flag=wx.LEFT | wx.ALIGN_RIGHT, border=4)
        self.gbs.Add(wx.StaticText(self.panel1, -1, 'TIME'),
                     pos=(2, 0), flag=wx.LEFT | wx.ALIGN_RIGHT, border=4)

        self.ct["DATE"] = wx.adv.DatePickerCtrl(self.panel1)

        myDate = wx.DateTime.Today()
        myDate.SetDay(int(self.myProfile['DATE'][2]))
        myDate.SetMonth(int(self.myProfile['DATE'][1]))
        myDate.SetYear(int(self.myProfile['DATE'][0]))
        myDate.FormatISODate()
        self.ct["TIME"] = wx.adv.TimePickerCtrl(self.panel1)

        self.gbs.Add(self.ct['DATE'], pos=(1, 1), flag=wx.ALIGN_LEFT)
        self.gbs.Add(self.ct['TIME'], pos=(2, 1), flag=wx.ALIGN_LEFT)
        line = 3
        for label in self.myProfile.profSurfParamList[3:]:
            self.gbs.Add(wx.StaticText(self.panel1, -1, label),
                         pos=(line, 0),
                         flag=wx.LEFT | wx.ALIGN_RIGHT, border=4)

            if (label == 'ICEDE_PARAM' or label == 'ICE_SCHEME' or label ==
                    'CLW_SCHEME'):
                self.ct[label] = wx.ComboBox(
                    self.panel1,
                    choices=list(self.myProfile.typeListe[label].values()))
            else:
                self.ct[label] = FS.FloatSpin(
                    self.panel1, -1,
                    min_val=self.myProfile.minValue[label],
                    max_val=self.myProfile.maxValue[label],
                    increment=(self.myProfile.maxValue[label] -
                               self.myProfile.minValue[label]) / 100,
                    agwStyle=FS.FS_LEFT)
                self.ct[label].SetFormat("%f")
                self.ct[label].SetDigits(3)
            self.gbs.Add(self.ct[label], pos=(line, 1))
            line += 1
        line = 1
        self.gbs.Add(wx.StaticText(self.panel1, -1, 'S2M parameters'),
                     pos=(line, 2),
                     span=(1, 2),
                     flag=wx.ALIGN_CENTER | wx.ALIGN_BOTTOM, border=8)
        for label in self.s2m_list:
            line += 1
            self.gbs.Add(wx.StaticText(self.panel1, -1, label),
                         pos=(line, 2),
                         flag=wx.LEFT | wx.ALIGN_RIGHT, border=4)
            self.ct['S2M_' + label] = FS.FloatSpin(
                self.panel1, -1,
                min_val=self.myProfile.minValue['S2M_' + label],
                max_val=self.myProfile.maxValue['S2M_' + label],
                increment=(self.myProfile.maxValue['S2M_' + label] -
                           self.myProfile.minValue['S2M_' + label]) / 100,
                agwStyle=FS.FS_LEFT)
            self.ct['S2M_' + label].SetFormat("%f")
            self.ct['S2M_' + label].SetDigits(3)

            self.gbs.Add(self.ct['S2M_' + label], pos=(line, 3))

        line += 1
        self.gbs.Add(wx.StaticText(self.panel1,
                                   -1, 'SKIN parameters'), pos=(line, 2),
                     span=(1, 2),
                     flag=wx.ALIGN_CENTER | wx.ALIGN_BOTTOM, border=8)
        for label in self.sskin_list[:-1]:
            line += 1
            self.gbs.Add(wx.StaticText(self.panel1, -1, label),
                         pos=(line, 2),
                         flag=wx.LEFT | wx.ALIGN_RIGHT, border=4)
            if (label == 'SURFTYPE' or label == 'WATERTYPE'):
                self.ct['SKIN_' + label] = wx.ComboBox(
                    self.panel1, choices=list(self.myProfile.typeListe[
                        label].values()))
            else:
                self.ct['SKIN_' + label] = FS.FloatSpin(
                    self.panel1, -1,
                    min_val=self.myProfile.minValue['SKIN_' + label],
                    max_val=self.myProfile.maxValue['SKIN_' + label],
                    increment=(self.myProfile.maxValue['SKIN_' + label] -
                               self.myProfile.minValue['SKIN_' + label]) / 100,
                    agwStyle=FS.FS_LEFT)
                self.ct['SKIN_' + label].SetFormat("%f")
                self.ct['SKIN_' + label].SetDigits(3)
            self.gbs.Add(self.ct['SKIN_' + label], pos=(line, 3))
        line += 1
        self.gbs.Add(wx.StaticText(self.panel1, -1, 'Fastem parameters'),
                     pos=(line, 2),
                     span=(1, 2),
                     flag=wx.ALIGN_CENTER | wx.ALIGN_BOTTOM, border=8)
        self.cb = wx.ComboBox(self.panel1, choices=list(self.valFastem.keys()))
        line += 1
        self.gbs.Add(self.cb, pos=(line, 2), span=(
            1, 2), flag=wx.LEFT | wx.ALIGN_RIGHT)
        line += 1

        self.applyBtn = wx.Button(self.panel1, label="Apply")
        self.applyBtn.SetHelpText("Apply changes")
        self.applyBtn.SetDefault()
        line += 1
        self.gbs.Add(self.applyBtn, pos=(line, 3))
        self.revertBtn = wx.Button(self.panel1, label="Revert")
        self.revertBtn.SetHelpText("Revert to previous values")
        self.gbs.Add(self.revertBtn, pos=(line, 2),
                     flag=wx.RIGHT | wx.BOTTOM, border=4)
        self.gbs.AddGrowableRow(line)
        self.panel1.SetSizerAndFit(self.gbs)
        self.SetProfileValues()

    def OnApplyChange(self, e, aProfile=None):
        """ apply the choice to the profile
            may re-initialize the profile with
            aProfile (if modified elsewhere) """
        self.GetProfileValues(aProfile)
        return self.myProfile

    def OnRevertChange(self, e):
        """ revert to the precedent value of the profile """

        self.SetProfileValues()

    def GetProfileValues(self, aProfile=None):
        """ get the value from the frame and return a profile
            may re-initialize the profile with aProfile
            (if modified elsewhere) """
        if aProfile is not None:
            self.myProfile = aProfile
        try:
            self.myProfile['DATE'][0] = self.ct['DATE'].GetValue().GetYear()
            self.myProfile['DATE'][1] = self.ct[
                'DATE'].GetValue().GetMonth() + 1
            self.myProfile['DATE'][2] = self.ct['DATE'].GetValue().GetDay()
            self.myProfile['TIME'][0] = self.ct['TIME'].GetTime()[0]
            self.myProfile['TIME'][1] = self.ct['TIME'].GetTime()[1]
            self.myProfile['TIME'][2] = self.ct['TIME'].GetTime()[2]
            self.myProfile['ID'] = str(self.ct['ID'].GetValue())
            for param in self.myProfile.profSurfParamList[3:]:
                if param in ('ICEDE_PARAM', 'ICE_SCHEME', 'CLW_SCHEME'):
                    val = self.ct[param].GetValue()
                    self.myProfile[param] = [
                        k for (k, v) in list(self.myProfile.typeListe[
                            param].items()) if v == val][0]
                else:
                    self.myProfile[param] = self.ct[param].GetValue()
            for param in self.s2m_list:
                self.myProfile['S2M'][param] = self.ct[
                    'S2M_' + param].GetValue()
            if self.cb.GetValue() != "":
                self.myProfile['SKIN']['FASTEM'] = numpy.array(
                    self.valFastem[self.cb.GetValue()])
            for param in self.sskin_list[:-1]:
                if (param == 'SURFTYPE' or param == 'WATERTYPE'):
                    val = self.ct['SKIN_' + param].GetValue()
                    self.myProfile['SKIN'][param] = [
                        k for (k, v) in list(self.myProfile.typeListe[
                            param].items()) if v == val][0]
                else:
                    self.myProfile['SKIN'][param] = self.ct[
                        'SKIN_' + param].GetValue()
        except ValueError:
            print("error getProfileValues")
        return self.myProfile

    def SetProfileValues(self, aProfile=None):
        """ set all values of the profile to the windows """
        if self.debug:
            print("surface SetProfileValues")
        if aProfile is not None:
            self.myProfile = aProfile
        myDate = wx.DateTime.Today()
        myDate.SetDay(1)
        myDate.SetDay(int(self.myProfile['DATE'][2]))
        myDate.SetMonth(int(self.myProfile['DATE'][1] - 1))
        myDate.SetYear(int(self.myProfile['DATE'][0]))
        self.ct['DATE'].SetValue(myDate)
        self.ct['TIME'].SetTime(self.myProfile['TIME'][0],
                                self.myProfile['TIME'][1],
                                self.myProfile['TIME'][2])
        for param in self.myProfile.profSurfParamList[2:]:
            if (param == 'ID'):
                self.ct[param].SetValue(self.myProfile[param])
            else:
                if param in ('ICEDE_PARAM', 'ICE_SCHEME', 'CLW_SCHEME'):
                    self.ct[param].SetValue(self.myProfile.typeListe[
                                            param][self.myProfile[param]])
                else:
                    self.ct[param].SetValue(float(self.myProfile[param]))
        for param in self.s2m_list:
            self.ct[
                'S2M_' + param].SetValue(float(self.myProfile['S2M'][param]))
        for param in self.sskin_list[:-1]:
            # all parameters except Fastem
            if (param == 'SURFTYPE' or param == 'WATERTYPE'):
                self.ct['SKIN_' + param].SetValue(
                    self.myProfile.typeListe[param][
                        self.myProfile['SKIN'][param]])
            else:
                self.ct[
                    'SKIN_' + param].SetValue(float(
                        self.myProfile['SKIN'][param]))
        for key in list(self.valFastem.keys()):
            value = 'Default'
            myA = numpy.array(self.valFastem[key])
            myB = abs(myA - self.myProfile['SKIN']['FASTEM']) < 0.000001
            if (myB.all):
                value = key
                break
        self.cb.SetHelpText(self.myProfile['SKIN'][
                            'FASTEM_ATTRIBUTE']['COMMENT'])
        self.cb.SetValue(value)

    def ShowErrorMessageDialogBox(self, varName, varType):
        message = "variable " + varName + "must be of " + varType + " type"
        dlg = wx.MessageDialog(
            None, message, caption="Error", style=wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.DeletePendingEvents()
        wx.CallAfter(dlg.Destroy)

    def PlotReflectanceEmissivity(self, emis, refl, WAVENUMBERS=None):
        if self.debug:
            print("plot Emissivity ")
        self.emissivityPageGraphic.OnPlot(
            emis['EMIS_IN'], emis['EMIS_OUT'], WAVENUMBERS)
        if self.debug:
            print("plot reflectance")
        if refl['REFL_IN'] is not None:
            if self.debug:
                print(refl['REFL_IN'], refl['REFL_OUT'])
            self.reflectancePageGraphic.OnPlot(
                refl['REFL_IN'], refl['REFL_OUT'], WAVENUMBERS)

    def OnCreateEmissivityReflectance(self, e):
        """ create Emissivity/Reflectance values for each
            channels from scratch"""
        pass

    def OnLoadAtlas(self, e):
        """ load Emissivity/Reflectance values for each chanels from Atlas"""
        pass

    def OnSaveProfile(self, e):
        """ return the name for a profile File """
        fileSaved = self.OnSaveFile(e, "Save a profile")
        return fileSaved

    def OnSaveProfileAs(self, e):
        """ return the name for a profile File """
        fileSaved = self.OnSaveFile(e, "Save a profile")
        return fileSaved

    def OnSaveSurface(self, e):
        """ return the name of a Surface File """
        fileSaved = self.OnSaveFile(e, "Save a Surface File")
        return fileSaved

    def OnSaveSurfaceAs(self, e):
        """ return the name of a Surface File """
        fileSaved = self.OnSaveFile(e, "Save a Surface File")
        return fileSaved

    def _initItem(self):
        self.applyItem = 0
        self.saveItem = 0

    def OnEditEmissivity(self, e): pass

    def OnEditReflectance(self, e): pass

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu

                ("Apply change", "Apply change for the profile ",
                 self.OnApplyChange, "applyChange", True),
                ("Initialize Emissivity/Reflectance",
                 "Create Emissivity/Reflectances values for each channel",
                 self.OnCreateEmissivityReflectance,
                 "createEmissivityReflectance", False),
                ("Modify values of Emissivity",
                 "Modify Emissivity values for each channel",
                 self.OnEditEmissivity, "modifyEmissivity", False),
                ("Modify values of Surface BRDF",
                 "Modify Surface BRDF values for each channel",
                 self.OnEditReflectance, "modifyReflectance", False),
                ("Load Atlas",
                 "Load value of Emissivity/Surface BRDF from atlas",
                 self.OnLoadAtlas, "loadAtlas", False),
                ("", "", "", "", True),

                ("Save profile", "Save the profile file",
                 self.OnSaveProfile, "saveProfile", True),
                ("Save profile as", "Save the profile file",
                 self.OnSaveProfileAs, "saveProfileAs", True),
                ("Save surface", "Save the surface file",
                 self.OnSaveSurface, "saveSurface", False),
                ("Save surface as", "Save the surface file",
                 self.OnSaveSurfaceAs, "saveSurfaceAs", False),
                ('&Quit', 'Quit', self.OnQuit, "quit", True)),
               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True),
                ("&Help", "Help", self.OnHelp, "help", True)))


if __name__ == "__main__":
    import sys
    import os
    p = rmodel.project.Project()
    prefix = p.config.ENV['RTTOV_GUI_PREFIX']
    print("Configuration : ", dir)

    err = p.openProfile(p.config.ENV["RTTOV_GUI_PROFILE_DIR"] +
                        "/cldaer101lev_allgas.H5", 1)
    if err != 0:
        print("error open profile")
        sys.exit(1)
    p.myCoeffs.fileName["standard"] = (p.config.ENV["RTTOV_GUI_COEFF_DIR"] +
                                       "/rttov9pred54L/rtcoef_eos_2_modis.dat")
    p.loadCoefficients()
    print(p.myCoeffs.nchannels)
    p.myOption["ADDSOLAR"].value = True
    p.runDirect()
    surfaceFile = os.path.join(p.config.ENV["GUI_WRK_DIR"],
                               "surface.h5")
    p.openSurface(surfaceFile)
    ex = wx.App()
    sv = SurfaceView(None, "Surface", p)
    sv.PlotReflectanceEmissivity(p.myEmissivity, p.myReflectance)
    # these binding must be in the controller
    sv.applyBtn.Bind(wx.EVT_BUTTON, sv.OnApplyChange)
    sv.Bind(wx.EVT_MENU, sv.OnSaveSurfaceAs, sv.items["saveSurfaceAs"])
    sv.Bind(wx.EVT_MENU, sv.OnSaveSurface, sv.items["saveSurface"])
    sv.Bind(wx.EVT_MENU, sv.OnSaveProfileAs, sv.items["saveProfileAs"])
    sv.Bind(wx.EVT_MENU, sv.OnSaveProfile, sv.items["saveProfile"])
    sv.items["modifyEmissivity"].Enable(True)
    ex.MainLoop()
