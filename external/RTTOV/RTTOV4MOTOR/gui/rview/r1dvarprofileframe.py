# -*- coding: utf-8 -*-

import locale
import rmodel
from rview import util
from rview import wxmpl
import numpy
import matplotlib
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
import r1Dvar
import collections
import wx.lib.scrolledpanel
from rview import colors
from rview.util import kindOfItem
from matplotlib.backends.backend_wxagg import\
    NavigationToolbar2WxAgg as ToolBar

stepColor = colors.RttovGui1DvarStepColor
item1DvarColors = colors.RttovGui1DvarItemColors
itemColor = colors.profileItemColors
itemMarker = colors.profileItemMarkers

axesDef = {"GASES":     {"xlimits": None, "xscale": "linear",
                         "ylimits": None, "yscale": "log"},

           'Q':   {"xlimits": None, "xscale": "linear",
                   "ylimits": None, "yscale": "log"},

           "T":   {"xlimits": None, "xscale": "linear",
                   "ylimits": None, "yscale": "log"},


           }

locale.setlocale(locale.LC_ALL, '')


class Data(wxmpl.Channel):
    """ Class to keep data for the PlotItemPanel data """

    def __init__(self, x, y, theName="", theColor=None, theStyle=None,
                 theMarker=None):
        wxmpl.Channel.__init__(
            self, name=theName, color=theColor,
            style=theStyle, marker=theMarker)
        self.x = x
        self.y = y
        self.label = theName

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def myUpdate(self, x, y):
        self.x = numpy.zeros(x.shape[0]) + x
        self.y = numpy.zeros(y.shape[0]) + y
        self.setChanged('True')


class PlotItemPanel(wx.Panel):
    """ plot on a PlotPanel one curve """

    def __init__(self, parent, profile, background, retrieved, xlegend,
                 theName, yInpressions=True):
        """ profile : is the True profile """
        """ background : is the initial profile """
        """ retrieved is a list of retrieved profiles """
        self.theName = theName
        self.theParent = parent
        self.xlegend = xlegend
        self.profile = profile
        self.background = background
        self.retrieved = retrieved
        self.yInPression = yInpressions

        self.pression = profile["P"]
        self.value = profile[theName]

        wx.Panel.__init__(self, parent, style=wx.BORDER_SIMPLE)
        self.myCurves = []
        self.fig = Figure()
        self.canvas = FigureCanvas(self, -1, self.fig)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.GROW, 1)
        self.OnPlot()
        self.tlb = ToolBar(self.canvas)
        self.sizer.Add(self.tlb, 0, wx.GROW)
        self.tlb.Realize()
        self.text = wx.StaticText(self, -1, label="")
        self.sizer.Add(self.text)
        self.Fit()
        self.SetSizer(self.sizer)
        self.valueHistory = []
        self.valueHistoryRedo = []
        self.canvas.mpl_connect('motion_notify_event', self.onMouseMotion)

    def onMouseMotion(self, event):
        """ set text when moving mousse """

        if event.inaxes:

            xdata = event.xdata
            ydata = event.ydata
            xstr = "%0.4g" % xdata
            ystr = "%0.4g" % ydata

            value = str(self.axes.get_ylabel()) + "=" + ystr + \
                "  " + str(self.axes.get_xlabel()) + "=" + xstr

            self.text.SetLabel(value)

    def SetXlimits(self, theName=None, xmin=None, xmax=None):
        """ set x limits """,

        if axesDef[self.theName]["xlimits"] is not None:
            self.axes.set_xlim(axesDef[self.theName]["xlimits"])

        self.axes.set_xscale(axesDef[self.theName]["xscale"])

    def DrawCurves(self):

        while(len(self.myCurves) > 0):
            self.axes.lines.remove(self.axes.lines[0])
            self.myCurves.pop()

        for data in self.myChannelList:
            c = self.axes.plot(
                data.x, data.y, color=data.color, marker=data.marker,
                label=data.label)
            self.myCurves.append(c)

        self.axes.legend(loc=0, fontsize=10)
        self.fig.canvas.draw_idle()

    def UpdateData(self, dataX):
        self.x = dataX
        self.data.setChanged(True)
        self.data.myUpdate(self.x, self.y)
        self.stripCharter.update()

    def OnPlot(self):
        matplotlib.rc('xtick', labelsize=8)
        matplotlib.rc('ytick', labelsize=8)
        fig = self.fig
        self.axes = fig.gca()
        self.x = self.value[::1]
        if self.yInPression:
            self.axes.set_yscale("log")
            self.axes.set_yticks((0.00005, 0.0001, 0.0002, 0.0005, 0.001,
                                  0.002, 0.005, 0.01,
                                  0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10,
                                  25, 50, 100, 200, 300, 500, 1000))
            label = ('5e-5', '1e-4', '2e-4', '5e-4', '1e-3',
                     '2e-3', '5e-3', '0.01', '0.02', '0.05',
                     '0.1', '0.2', '0.5', '1', '2', '5', '10',
                     '25', '50', '100', '200', '300', '500', '1000')
            for tick in self.axes.yaxis.get_major_ticks():
                tick.label.set_fontsize(10)
            self.axes.set_yticklabels(label)
            self.axes.set_ylabel('pressure (hPa)')
            self.axes.set_ylim((self.pression[-1] + 150, self.pression[0]))
        else:
            self.axes.set_ylim(self.value.shape[0] + 2, 1)
            self.axes.set_ylabel('level')

        if self.yInPression:
            self.y = self.pression[::1]
        else:
            self.y = numpy.arange(1, self.myProfile["T"].shape[0] + 1, 1)

        self.data = collections.OrderedDict()
        self.data["TrueProfile"] = Data(self.profile[self.theName], self.y,
                                        theName=self.theName + " True",
                                        theColor=item1DvarColors["true"],
                                        theMarker=itemMarker[self.theName])
        self.data["BgProfile"] = Data(self.background[self.theName], self.y,
                                      theName=self.theName + " Bg",
                                      theColor=item1DvarColors["background"],
                                      theMarker=itemMarker[self.theName])

        for i in range(1, len(self.retrieved) + 1):
            self.data["retrieved" + str(i)] = Data(
                self.retrieved[i - 1][self.theName],
                self.y, theName=self.theName + " run " + str(i),
                theColor=stepColor[i - 1],
                theMarker=itemMarker[self.theName], theStyle="-")

        self.axes.set_xlabel(self.xlegend)
        formatter = self.axes.xaxis.get_major_formatter()
        formatter.set_powerlimits((-3, 4))
        self.axes.xaxis.set_major_formatter(formatter)

        self.axes.grid(True, axis='both')

        self.myChannelList = []
        for k, v in list(self.data.items()):
            self.myChannelList.append(v)

        self.DrawCurves()


class MyNotebook(wx.Notebook):
    """ MyNotebook Class inherits from Notebook
        Bind the event EVT_NOTEBOOK_PAGE_CHANGED
        with Update method of the page """

    def __init__(self, parent, isRightPage=True):
        wx.Notebook.__init__(self, parent, id=wx.ID_ANY, style=wx.BK_DEFAULT
                             )
        self.parent = parent
        self.isRightPage = isRightPage
        self.BindEvent()

    def BindEvent(self):
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)

    def OnPageChanged(self, event):
        sel = self.GetSelection()
        thePage = self.GetPage(sel)
        thePage.Update()
        if self.isRightPage:
            item = thePage.theName
        event.Skip()


class r1dvarProfileView (util.GenericViewRadio):
    """ r1dvarProfile View : display True, Background and retrieved
        profile from different tries """
    helpTitle = "Help Profile"
    helpMessage = """
    This windows display True, Background and retrieved profiles
        """

    def __init__(self, parent, title, the1dvarProject):

        self.yInPression = True
        self.myProfile = the1dvarProject.pTrue.myProfile
        self.myBgInitialProfile = the1dvarProject.BgInitialProfile
        self.myRetrievedProfiles = the1dvarProject.retrievedProfiles
        super().__init__(parent,  title)
        self.infos = collections.OrderedDict()
        self.infos["True Profile File Name"] = (
            the1dvarProject.pTrue.profileInitialFilename[0] +
            " number {}".format(
                the1dvarProject.pTrue.profileInitialFilename[1]))
        self.infos["Bg Profile File Name"] = (
            the1dvarProject.pBg.profileInitialFilename[0] +
            " number {}".format(
                the1dvarProject.pBg.profileInitialFilename[1]))
        self.step_counter = 1
        self.infos["Instrument run 1"] = (
            the1dvarProject.pBg.myCoeffs.get_instrument_and_platform_name()[0])
        self.infos["Satellite run 1"] = (
            the1dvarProject.pBg.myCoeffs.get_instrument_and_platform_name()[
                1] +
            str(
                the1dvarProject.pBg.myCoeffs.get_instrument_and_platform_name()[2]))

        cl = the1dvarProject.pBg.myCoeffs.getFF_ORI_CHN()
        lon = len(cl)
        if lon < 20:
            self.infos['Channel selection run 1'] = "{}".format(cl)
        else:
            self.infos['Channel selection run 1'] = (
                "[ {}".format(cl[0]) + " {}".format(cl[1]) +
                " ... {}".format(cl[lon - 2]) + " {} ]".format(cl[lon - 1]))
        self.infos["Number of channels run 1"] = lon

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        self.CreateMenuBar()
        self.SetSize((1300, 700))
        self.SetMinSize((1100, 700))
        self.SetTitle("Retrieved profile")
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
        self.infoPages = {}
        self.infoPages['INFO'] = None
        self.Plot(self.myProfile, self.myBgInitialProfile,
                  self.myRetrievedProfiles)

        # create a second sizer for the notebook
        sizer2 = wx.BoxSizer()
        sizer2.Add(self.nb, 1, wx.EXPAND)
        self.panel2.SetSizer(sizer2)
        self.sb = self.CreateStatusBar()
        self.sb.SetBackgroundColour('WHITE')
        txt = ''
        self.sb.SetStatusText(txt)
        self.Centre()
        self.Show(True)

    def RePlot(self, the1dvarProject):
        self.myProfile = the1dvarProject.pTrue.myProfile
        self.myBgInitialProfile = the1dvarProject.BgInitialProfile
        self.myRetrievedProfiles = the1dvarProject.retrievedProfiles
        self.nb.DeleteAllPages()
        self.nb_all.DeleteAllPages()
        self.step_counter = self.step_counter + 1
        self.infos["Instrument run {}".format(self.step_counter)] = (
            the1dvarProject.pBg.myCoeffs.get_instrument_and_platform_name()[0])
        self.infos["Satellite run {}".format(self.step_counter)] = (
            the1dvarProject.pBg.myCoeffs.get_instrument_and_platform_name()[
                1] +
            str(
                the1dvarProject.pBg.myCoeffs.get_instrument_and_platform_name()[
                    2]))

        cl = the1dvarProject.pBg.myCoeffs.getFF_ORI_CHN()
        lon = len(cl)
        if lon < 20:
            self.infos['Channel selection run {}'.format(
                self.step_counter)] = "{}".format(cl)
        else:
            self.infos['Channel selection run {}'.format(
                self.step_counter)] = "[ {}".format(
                cl[0]) + " {}".format(cl[1]) + " ... {}".format(
                cl[lon - 2]) + " {} ]".format(cl[lon - 1])
        self.infos["Number of channels run {}".format(self.step_counter)] = lon
        self.Plot(self.myProfile, self.myBgInitialProfile,
                  self.myRetrievedProfiles)

    def Plot(self, profile, background, retrieved):

        self._ComputeLayers(self.myProfile['P'])

        self.graphicPages = {}
        self.graphicPages['T'] = PlotItemPanel(
            self.nb, profile, background, retrieved, theName='T',
            xlegend=self.myProfile['T_ATTRIBUTE']['UNITS'])
        self.nb.AddPage(self.graphicPages['T'], 'T')

        for gas in "Q":
            if self.myProfile[gas] is not None:
                self.graphicPages[gas] = PlotItemPanel(
                    self.nb, profile, background, retrieved,
                    theName=gas, xlegend=self.myProfile[
                        gas + '_ATTRIBUTE']['UNITS'],
                    yInpressions=self.yInPression)
                self.nb.AddPage(self.graphicPages[gas], gas)

        # delete empty graphicPages
        for key in list(self.graphicPages.keys()):
            if self.myProfile[key] is None:
                del self.graphicPages[key]

        # plot panel 1 with all gas
        if not self.infoPages['INFO']:
            self.infoPages['INFO'] = wx.lib.scrolledpanel.ScrolledPanel(
                self.nb_all, -1)
            self.infoPages['INFO'].SetupScrolling()
        self.FillInfoPanel(
            self.infoPages['INFO'], profile, background, retrieved)

        self.nb_all.AddPage(self.infoPages['INFO'], 'INFO')

    def FillInfoPanel(self, panel, profile, background, retrieved):
        sizer = wx.BoxSizer(wx.VERTICAL)
        panel.SetSizer(sizer)
        panel.Fit()
        # print info from Background :
        for k, v in list(self.infos.items()):
            sizer.Add(wx.StaticText(panel, -1, k + " : " + str(v)))
        sizer.Add((10, 20))
        sizer.Add(wx.StaticText(
            panel, -1, "Surface Information on Background Profile : "))
        sizer.Add(wx.StaticText(panel, -1, " Tskin =" + "%.2f" %
                                background["SKIN"]["T"] + " " + background[
                                    "T_ATTRIBUTE"]["UNITS"]))
        sizer.Add(wx.StaticText(panel, -1, " Q 2m =" + "%.2f" %
                                background["S2M"]["Q"] + " " + background[
                                    "Q_ATTRIBUTE"]["UNITS"]))
        sizer.Add(wx.StaticText(panel, -1, " T 2m =" + "%.2f" %
                                background["S2M"]["T"] + " " + background[
                                    "T_ATTRIBUTE"]["UNITS"]))
        sizer.Add((10, 20))
        sizer.Add(wx.StaticText(
            panel, -1, "Surface Information on True Profile : "))
        sizer.Add(wx.StaticText(panel, -1, " Tskin =" + "%.2f" %
                                profile["SKIN"]["T"] + " " + profile[
                                    "T_ATTRIBUTE"]["UNITS"]))
        sizer.Add(wx.StaticText(panel, -1, " Q 2m =" + "%.2f" %
                                profile["S2M"]["Q"] + " " + profile[
                                    "Q_ATTRIBUTE"]["UNITS"]))
        sizer.Add(wx.StaticText(panel, -1, " T 2m =" + "%.2f" %
                                profile["S2M"]["T"] + " " + profile[
                                    "T_ATTRIBUTE"]["UNITS"]))
        sizer.Add((10, 20))
        sizer.Add(wx.StaticText(
            panel, -1, "Surface Information on Retrieved Profiles : "))
        for i in range(1, len(retrieved) + 1):
            sizer.Add(wx.StaticText(panel, -1, " run " + "%d" % i))
            sizer.Add(wx.StaticText(
                panel, -1, "   Tskin =" + "%.2f" % retrieved[i - 1][
                      "SKIN"]["T"] + " " + retrieved[i - 1]["T_ATTRIBUTE"][
                    "UNITS"]))
            sizer.Add(wx.StaticText(
                panel, -1,
                "   Q 2m =" + "%.2f" % retrieved[
                    i - 1]["S2M"]["Q"] + " " + retrieved[i - 1][
                    "Q_ATTRIBUTE"]["UNITS"]))
            sizer.Add(wx.StaticText(
                panel, -1,
                "   T 2m =""%.2f" % retrieved[i - 1][
                      "S2M"]["T"] + " " + retrieved[i - 1][
                    "T_ATTRIBUTE"]["UNITS"]))
            sizer.Add((10, 20))

    def _ComputeLayers(self, pression):
        """ Compute the mean value of pression in a layer """
        foo = numpy.empty(pression.shape[0] - 1)
        for i in range(foo.shape[0]):
            foo[i] = (pression[i + 1] + pression[i]) / 2
        self.layer = foo

    def _MakeBinding(self):
        """ set the trivial Binding for the View """
        # binding cancel button

    def RePlotAllLeftPanel(self, profile=None):
        """ Plot the 2 panels with (new) profile (delete everything
            before redraw) """

        if profile is not None:
            self.myProfile = profile
            self._ComputeLayers(self.myProfile['P'])
        # remove all pages of the notebook
        self.nb_all.DeleteAllPages()
        self.PlotLeft()

    def addRightPage(self, item):
        """ add an new item page  """
        kind = kindOfItem[item]
        if kind == "GASES":
            myY = self.myProfile['P']
        else:
            myY = self.layer
        self.graphicPages[item] = PlotItemPanel(
            self.nb, self.myProfile[item], self.myProfile['P'],
            theName=item, kind=kind,
            xlegend=self.myProfile[item + '_ATTRIBUTE']['UNITS'],
            layer=myY, yInPression=self.yInPression)
        self.nb.AddPage(self.graphicPages[item], item)

    def OnMouseMove(self, e):
        """ print x y value of the left plot in the status bar  """
        pass

    def OnYPressions(self, e):
        self.yInPression = True
        self.RePlotAll()

    def OnYLevels(self, e):
        self.yInPression = False
        self.RePlotAll()

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu
                ('&Quit', 'Quit', self.OnQuit, "quit", True, False)),
               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True, False)))


if __name__ == "__main__":

    print("version matplotlib :", matplotlib.__version__)

    ex = wx.App()
    print(">>>>>>>>>> testrunRetrieve_HIRS ")
    pTrue = rmodel.project.Project()
    pBg = rmodel.project.Project()
    pTrue.setFileNameMark("TrueTestHIRS")
    pBg.setFileNameMark("BgTestHIRS")

    print(("pBg.profileFileName", pBg.profileFileName))
    filename = pTrue.config.ENV[
        "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
    pTrue.openProfile(filename, 1)
    pBg.openProfile(filename, 6)
    coefFile = pTrue.config.ENV["RTTOV_GUI_COEFF_DIR"] + \
        "/rttov7pred54L/rtcoef_noaa_19_hirs.dat"
    pBg.myCoeffs.fileName["standard"] = coefFile

    Bmatrix = (pTrue.config.ENV["RTTOV_GUI_PREFIX"] +
               "/r1Dvar/data/Sample_Bmatrices/Bmatrix_54L")
    retrieveProject = r1Dvar.r1dvar.Project1dvar(
        pTrue, pBg, matrixBfile=Bmatrix,
        satellite="noaa-19", instrument="hirs")
    retrieveProject.setFactorB(1)
    retrieveProject.setFactorR(1)

    Xr1 = retrieveProject.retrieve1d()

    Xr2 = retrieveProject.stepRetrieve1d()

    frame = r1dvarProfileView(None, "Retrieved profile", retrieveProject)

    frame.Show()

    print("loop")
    ex.MainLoop()
    # ex.MainLoop()
