# -*- coding: utf-8 -*-

import wx
import rmodel
from rview import util
import locale
import numpy
import matplotlib
from rview.profileframeutils import GenericPlotItemPanel, MyNotebook
from rview.profileframeutils import PlotItemPanelAll
from rview import profileframeutils as pfu
import logging
from rview import colors
from rview import myunits
itemColor = colors.profileItemColors
oldKColors = colors.oldKColors
itemMarker = colors.profileItemMarkers

locale.setlocale(locale.LC_ALL, '')

if wx.VERSION[1] == 8:
    from wx.lib.pubsub import Publisher as pub
else:
    if wx.VERSION[0] <= 3:
        # wxPython 2.9 or 3.0:
        # wx.lib.pubsub: Pusub now defaults to the new "kwarg" version of
        # the API.
        # In order to continue using the original "arg1" API you
        # will need to import wx.lib.pubsub.setuparg1 before importing any
        # other
        # pubsub modules.
        import wx.lib.pubsub.setuparg1
        from wx.lib.pubsub import pub
    else:
        # wxPython 4
        from pubsub import pub


class PlotItemPanel(GenericPlotItemPanel):
    """ plot on a PlotPanel multi curves """

    def __init__(self, parent, value, pression, theName,
                 kind="GASES", xlegend="(mw/cm-1/ster/sq.m)/ppmv",
                 layerstyle=False, layer=None, yInPressions=True, tskin=None,
                 previousValues=None, runNumber=None, xlabel=None):
        edit = False

        '''
        Constructor
        '''
        layerstyle = False
        layer = None
        self.StyleForPreviousValues = colors.StyleForPreviousValues
        self.previousValues = previousValues
        GenericPlotItemPanel. __init__(self, parent, value,
                                       pression, theName,
                                       kind,
                                       xlegend,
                                       edit,
                                       layerstyle,
                                       layer,
                                       yInPressions,
                                       tskin,
                                       xlabel=xlabel)

        self.SetTickSize(8)
        self.runNumber = runNumber
        self.data.label = "last run (" + str(self.runNumber) + ")"
        if self.previousValues is not None:
            self.addPrevious2myChannelList()
        self.PlotCurves()

    def addPrevious2myChannelList(self):
        if self.previousValues is not None:
            nkrun = self.runNumber
            ncurve = 0
            for previousV in self.previousValues:
                if ncurve < len(self.StyleForPreviousValues):
                    myStyle = self.StyleForPreviousValues[ncurve]
                    myColor = oldKColors[ncurve]
                else:
                    myStyle = self.StyleForPreviousValues[ncurve][-1]
                    myColor = oldKColors[ncurve]
                nkrun = nkrun - 1
                ncurve = ncurve + 1
                logging.debug("add data nkrun" + str(nkrun) + " ncurve " + str(
                    ncurve) + myStyle)
                data = pfu.Data(previousV, self.y,
                                theName=self.theName,
                                theColor=myColor,
                                theMarker=itemMarker[self.theName],
                                theStyle=myStyle,
                                theLabel="run " + str(nkrun))
                self.myChannelList.append(data)

    def PlotCurves(self):
        logging.debug("PlotCurves for " + self.theName)
        logging.debug(
            "PlotCurves len(self.myCurves)" + str(len(self.myCurves)))
        logging.debug("len self.axes.lines" + str(len(self.axes.lines)))
        while(len(self.axes.lines) > 0):
            self.axes.lines.remove(self.axes.lines[0])
        while(len(self.myCurves) > 0):
            self.myCurves.pop()
        for data in self.myChannelList:
            c = self.axes.plot(
                data.x, data.y, color=data.color, marker=data.marker,
                label=data.label, linestyle=data.style)
            self.myCurves.append(c)
        self.fig.canvas.draw_idle()

        self.axes.legend(loc=0, fontsize=10)


class KProfileView (util.GenericViewRadio):
    """ Profile window of the application """
    helpTitle = "Help Profile"
    helpMessage = """
    Select and visualize a component profile on the right panel
    Click left button to modify the profile.
    Click left and drag a zone to zoom in.
    Click right button to zoom out.
    Apply your changes or save the profile
    for the next run of RTTOV.
        """

    def __init__(self, parent, title, channel=1,
                 edit=False, yInPressions=True, runNumber=1,
                 baseProfile=None, kProfile=None, project=None):

        super().__init__(parent,  title)
        self.channel = int(channel)
        self.descChannel = ""
        self.edit = edit
        self.myProject = project
        self.oldValues = None
        self.runNumber = runNumber
        self.yInPressions = yInPressions
        self.baseProfile = baseProfile
        isOK = False
        self.initializedByProject = False

        if self.canUpdateFromProject():
            nchan = project.aKmats.aKmat[-1].getnchannels()
            print("kProfileView: initialization with Project")
            print(("nchannels:", nchan, "channel number:", self.channel))
            isOK = self.InitValuesFromProject()
            if isOK:
                self.initializedByProject = True
        else:
            # initialization with baseProfile and kProfile
            print("kProfileView: initialization with kProfile")
            if baseProfile is not None and kProfile is not None:
                self.myKProfile = kProfile
                self.myKProfile['P'] = baseProfile['P']
                isOK = True

        if not isOK:
            print("ERROR cannot create kProfileView")
            return
        self.InitUnitsForThePlots()
        self.setOldValues()

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        self.CreateMenuBar()

        self.items["ypressions"].Enable(True)
        self.items["ypressions"].Check(True)
        self.items["ylevels"].Check(False)

        self.SetSize((1000, 700))
        self.SetMinSize((1200, 700))
        # panel 1 notebook with all curves (GASES, AEROSOLS, CLOUDS)
        self.panel1 = wx.Panel(self, -1, style=wx.BORDER_SIMPLE)
        self.panel1.SetSize((200, 500))
        sizer.Add(self.panel1, 1, wx.EXPAND)
        # creation of notebook for the panel 1
        self.nb_all = MyNotebook(self.panel1, isRightPage=False)

        sizer1 = wx.BoxSizer()
        sizer1.Add(self.nb_all, 1, wx.EXPAND)
        self.panel1.SetSizer(sizer1)

        # panel 2 notebook with one curve
        self.panel2 = wx.Panel(self, -1, style=wx.BORDER_SIMPLE)
        sizer.Add(self.panel2, 1, wx.EXPAND)
        # creation of the notebook for the panel 2
        self.nb = MyNotebook(self.panel2)
        self.axesDef = []

        # creation des graphiques
        self.Plot(self.myKProfile)

        # create a second sizer for the notebook
        sizer2 = wx.BoxSizer()
        sizer2.Add(self.nb, 1, wx.EXPAND)
        self.panel2.SetSizer(sizer2)
        self.sb = self.CreateStatusBar()
        self.sb.SetBackgroundColour('WHITE')
        txt = self.descChannel
        self.sb.SetStatusText(txt)
        self.Centre()
        self.Show(True)
        pub.subscribe(self.OnKMatChanged, "K-Matrix Changed")

        if self.initializedByProject:
            self.EnableMenuItem("resetk")
            self.EnableMenuItem("scaleTen")

    def InitUnitsForThePlots(self):
        self.Tunit = myunits.makePretty(
            "T",
            self.myKProfile["T_ATTRIBUTE"]['UNITS'],
            self.myKProfile['GAS_UNITS']
        )
        self.GasUnit = myunits.makePretty(
            "Q",
            self.myKProfile["Q_ATTRIBUTE"]['UNITS'],
            self.myKProfile['GAS_UNITS']
        )
        self.TunitforAxe = myunits.makePrettyForMatplotlib(
            "T",
            self.myKProfile["T_ATTRIBUTE"]['UNITS'],
            self.myKProfile['GAS_UNITS']
        )
        self.GasUnitForAxe = myunits.makePrettyForMatplotlib(
            "Q",
            self.myKProfile["Q_ATTRIBUTE"]['UNITS'],
            self.myKProfile['GAS_UNITS']
        )

    def InitValuesFromProject(self):
        """ init the kprofile from the Project k-matrix archive """
        logging.info("InitValuesFromProject")
        self.oldValues = None
        if self.myProject is not None:
            aKmat = self.myProject.aKmats.aKmat
            logging.debug("get last run from aKmat " + str(len(aKmat)))
            if len(aKmat) > 0:
                nchan = aKmat[-1].getnchannels()
                if (self.channel <= nchan and self.channel > 0):
                    self.myKProfile = aKmat[-1].getchanprof(self.channel - 1)
                    self.baseProfile = aKmat[-1].profile
                    wavenumber = aKmat[-1].getwavenumber(self.channel)
                    instrument = aKmat[-1].getinstrument()
                    sat = aKmat[-1].getsatname()
                    nchannels = aKmat[-1].getnchannels()
                    self.descChannel = (sat.upper() + " / " +
                                        instrument.upper() +
                                        " C=#" + str(self.channel) +
                                        " / " + str(nchannels) + " (" +
                                        "wavenumber: " +
                                        str(wavenumber) + "(cm" +
                                        myunits.exponentMinus1 +
                                        "))")
                    logging.debug("profile" + str(self.myKProfile['T'].shape))
                    logging.debug("baseProfile" +
                                  str(self.baseProfile['T'].shape))
                    self.myKProfile['P'] = self.baseProfile['P']
                    self.runNumber = self.myProject.aKmats.runNumber
                    self.InitUnitsForThePlots()
                    return True
                else:
                    logging.warning("Cannot init Values From Project")
                    return False
            else:
                logging.warning("Cannot init Values From Project")
                return False
        else:
            return False

    def setOldValues(self):
        """ init old values from the Project k-matrix archive """
        logging.debug(">>>setOldValues")
        self.oldValues = {}
        self.oldValues['T'] = None
        for gas in self.myKProfile.gas_list:
            self.oldValues[gas] = None

        archiveLen = 0
        hasOldValues = False

        if self.canUpdateFromProject():
            logging.debug("yes canUpdateFromProject")
            archiveLen = self.myProject.aKmats.getLengh()
            logging.debug("DEBUG: archilLen" + str(archiveLen))
            if archiveLen >= 2:
                hasOldValues = True

        if hasOldValues:
            self.oldValues['T'] = []
            for gas in self.myKProfile.gas_list:
                if self.myKProfile[gas] is not None:
                    self.oldValues[gas] = []

            for ikmat in range(archiveLen - 2, -1, -1):
                kmat = self.myProject.aKmats.aKmat[ikmat]
                ov = kmat.getchanprof(self.channel - 1)['T']
                self.oldValues['T'].append(ov)
                for gas in self.myKProfile.gas_list:
                    if self.myKProfile[gas] is not None:
                        ov = kmat.getchanprof(self.channel - 1)[gas]
                        self.oldValues[gas].append(ov)

    def canUpdateFromProject(self):
        can = False
        if self.myProject is not None:
            project = self.myProject
            if len(project.aKmats.aKmat) > 0:
                nchan = project.aKmats.aKmat[-1].getnchannels()
                if self.channel > 0 and self.channel <= nchan:
                    can = True
        return can

    def OnKMatChanged(self, msg):
        if not self.canUpdateFromProject():
            # cannot update the kprofileframe : not Kmat matrice to read
            # in Project :
            # it happens when using the 1DVAR
            return
        for kmat in self.myProject.aKmats.aKmat:
            kmat.kscale(1.0)
        self.InitValuesFromProject()
        self.setOldValues()
        self.RePlotAll()
        self.sb.SetBackgroundColour('WHITE')
        txt = self.descChannel
        self.sb.SetStatusText(txt)

    def PlotLeft(self, profile=None):
        # plot panel 1 with all gas
        self.allGraphicsPages = {}
        self.panel1.ClearBackground()
        self.allGraphicsPages["last"] = PlotItemPanelAll(
            self.nb_all,
            self.myKProfile,
            kind='GASES',
            xlegendT=self.TunitforAxe,
            xlegend=self.GasUnitForAxe,
            yInPressions=self.yInPressions,
            addTskin=True,
            runNumber=self.runNumber)
        self.allGraphicsPages['last'].SetTickSize(8)
        self.nb_all.AddPage(self.allGraphicsPages['last'],
                            'RUN ' + str(self.runNumber))

        if self.myProject is not None:
            for kmat in self.myProject.aKmats.aKmat[-2::-1]:
                runNumber = kmat.runNumber
                kprofile = kmat.getchanprof(self.channel - 1)
                kprofile['P'] = kmat.profile['P']
                self.allGraphicsPages[runNumber] = PlotItemPanelAll(
                    self.nb_all,
                    kprofile,
                    kind='GASES',
                    XinLog=False,
                    xlegendT=self.TunitforAxe,
                    xlegend=self.GasUnitForAxe,
                    yInPressions=self.yInPressions,
                    addTskin=True,
                    runNumber=runNumber)
                self.allGraphicsPages[runNumber].SetTickSize(8)
                self.nb_all.AddPage(self.allGraphicsPages[runNumber],
                                    'RUN ' + str(runNumber))
        self.panel1.Refresh()
        self.panel1.Update()

    def Plot(self, profile=None):
        if profile is not None:
            self.myKProfile = profile
            self._ComputeLayers(self.myKProfile['P'])
        self.graphicPages = {}
        self.graphicPages['T'] = PlotItemPanel(
            self.nb,
            self.myKProfile['T'],
            self.myKProfile['P'],
            theName='T',
            xlegend=self.TunitforAxe,
            xlabel=self.Tunit,
            yInPressions=self.yInPressions,
            tskin=self.myKProfile[
                'SKIN']['T'],
            previousValues=self.oldValues['T'],
            runNumber=self.runNumber)
        self.nb.AddPage(self.graphicPages['T'], 'T')

        for gas in self.myKProfile.gas_list:
            if self.myKProfile[gas] is not None:
                self.graphicPages[gas] = PlotItemPanel(
                    self.nb,
                    self.myKProfile[gas],
                    self.myKProfile['P'],
                    theName=gas,
                    xlegend=self.GasUnitForAxe,
                    xlabel=self.GasUnit,
                    yInPressions=(self.yInPressions),
                    previousValues=self.oldValues[gas],
                    runNumber=self.runNumber)
                self.nb.AddPage(self.graphicPages[gas], gas)

        # delete empty graphicPages
        for key in list(self.graphicPages.keys()):
            if self.myKProfile[key] is None:
                del self.graphicPages[key]

        # plot panel 1 with all gas
        self.PlotLeft()

    def _ComputeLayers(self, pression):
        """ Compute the mean value of pression in a layer """
        foo = numpy.empty(pression.shape[0] - 1)
        for i in range(foo.shape[0]):
            foo[i] = (pression[i + 1] + pression[i]) / 2
        self.layer = foo

    def _MakeBinding(self):
        """ set the trivial Binding for the View """
        # binding cancel button

    def RePlotAll(self, profile=None):
        """ Plot the 2 panels with (new) profile
            (delete everything before redraw) """

        if profile is not None:
            self.myKProfile = profile
            self._ComputeLayers(self.myKProfile['P'])
        # remove all pages of the notebook
        self.nb.DeleteAllPages()
        self.nb_all.DeleteAllPages()
        self.Plot()

    def RePlotAllLeftPanel(self, profile=None):
        """ Plot the 2 panels with (new) profile
            (delete everything before redraw) """

        if profile is not None:
            self.myKProfile = profile
            self._ComputeLayers(self.myKProfile['P'])
        # remove all pages of the notebook

        self.nb_all.DeleteAllPages()
        self.PlotLeft()

    def OnMouseMove(self, e):
        """ print x y value of the left plot in the status bar  """
        pass

    def OnYPressions(self, e):
        self.yInPressions = True
        self.RePlotAll()

    def OnYLevels(self, e):
        self.yInPressions = False
        self.RePlotAll()

    def Onscale5k(self, e):
        self.Scalexk(0.05)

    def Onscale10k(self, e):
        self.Scalexk(0.1)

    def Onscale20k(self, e):
        self.Scalexk(0.2)

    def Scalexk(self, scalingFactor):
        if self.myProject is not None:
            for kmat in self.myProject.aKmats.aKmat:
                kmat.kscale(scalingFactor)
            self.InitValuesFromProject()
            self.setOldValues()
            self.RePlotAll()

    def OnResetK(self, e):
        self.Scalexk(1.0)

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu
                ('&Quit', 'Quit', self.OnQuit, "quit", True, False)),
               ("&Edit",  # Edit Menu
                ("Yaxis in pressure units", "put y in pressure units",
                 self.OnYPressions, "ypressions", False, True),
                ("Yaxis in level units", "put y in level unit",
                 self.OnYLevels, "ylevels", True, True),
                ("", "", "", "", True, False),
                ("Scale x 10% profile", "scale by 10% of profile",
                 self.Onscale10k, "scaleTen", False, False),
                ("Reset kprofile", "reset the kprofile",
                 self.OnResetK, "resetk", False, False)),

               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True, False)))


def checkCode(err):
    if err != 0:
        print((">>>> ERROR %d" % (err)))


if __name__ == "__main__":

    print("version matplotlib :", matplotlib.__version__)
    p = rmodel.project.Project()
    ex = wx.App()

    profileName = (p.config.ENV['RTTOV_GUI_PROFILE_DIR'] +
                   "/cldaer101lev_allgas.H5")
    err = p.openProfile(profileName, 1)
    checkCode(err)

    coefFile = (p.config.ENV['RTTOV_GUI_COEFF_DIR'] +
                "/rttov7pred54L/rtcoef_noaa_19_hirs.dat")

    p.myCoeffs.fileName["standard"] = coefFile
    p.myOption["INTERP_MODE"].value = 5
    err = p.loadCoefficients()
    p.runK()

    p.myOption["INTERP_MODE"].value = 2

    p.runK()
    p.myOption["INTERP_MODE"].value = 4

    p.runK()
    p.myOption["INTERP_MODE"].value = 3

    p.runK()

    checkCode(err)
    frame = KProfileView(None, "Kprofile", channel=18,
                         yInPressions=True, runNumber=3, project=p)
    frame.Show()

    ex.MainLoop()

    print("loop")
    ex.MainLoop()
