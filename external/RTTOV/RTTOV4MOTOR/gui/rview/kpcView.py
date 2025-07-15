# -*- coding: utf-8 -*-
import wx.lib.scrolledpanel
from rview import util
import rmodel
import os
import numpy
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import \
    NavigationToolbar2WxAgg as ToolBar
from matplotlib.figure import Figure
import rttov
import h5py
from rview import colors


class PlotPage(wx.Panel):

    def __init__(self, parent, theValues, level, pression, theName,
                 nbpc, axeXLegend="", inPressure=True, startPC=0):
        wx.Panel.__init__(self, parent, style=wx.BORDER_SIMPLE)
        self.colors = colors.kpcviewColors
        self.linestyles = colors.kpcviewLinestyles
        self.linewidths = colors.kpcviewLinewidths
        self.yInPressure = inPressure
        self.theName = theName
        self.theValues = theValues
        self.pression = pression
        self.level = level
        self.fig = Figure()
        self.nbpc = nbpc  # number of pc to draw

        self.sizer = wx.BoxSizer(wx.VERTICAL)

        self.SetSizer(self.sizer)

        self.axeXLegend = axeXLegend

        self.canvas = FigureCanvas(self, -1, self.fig)
        self.canvas.mpl_connect('motion_notify_event', self.onMouseMotion)
        tlb = ToolBar(self.canvas)
        tlb.Realize()

        self.sizer.Add(tlb, 0, wx.GROW)
        self.sizer.Add(self.canvas, 1, wx.LEFT |
                       wx.TOP | wx.GROW, 1, wx.EXPAND)

        self.text = wx.StaticText(self, -1, label="")
        self.sizer.Add(self.text)

        self.Fit()
        self.OnPlot(theValues, startPC)

    def onMouseMotion(self, event):

        txt = ""
        if event.inaxes:
            x = event.xdata
            y = event.ydata
            txt = 'X=%.1f Y=%.2f' % (x, y)
            self.text.SetLabel(txt)

    def OnPlot(self, theValues, startPC=0):

        self.fig.clear()
        self.axes = self.fig.add_subplot(1, 1, 1)
        theMarker = ''
        if self.yInPressure:
            self.axes.set_ylim((self.pression[-1], self.pression[0]))
            self.axes.set_yscale("log")
            self.axes.set_yticks((0.00005, 0.0001, 0.0002, 0.0005,
                                  0.001, 0.002, 0.005, 0.01,
                                  0.02, 0.05, 0.1, 0.2, 0.5,
                                  1, 2, 5, 10, 25, 50, 100,
                                  200, 300, 500, 900, 1000))
            label = ('5e-5', '1e-4', '2e-4', '5e-4',
                     '1e-3', '2e-3', '5e-3', '0.01',
                     '0.02', '0.05', '0.1', '0.2', '0.5',
                     '1', '2', '5', '10', '25', '50', '100',
                     '200', '300', '500', '900', '1000')
            for tick in self.axes.yaxis.get_major_ticks():
                tick.label.set_fontsize(10)
            self.axes.set_yticklabels(label)
            self.axes.set_ylabel('pressure (hPa)')
            self.axes.set_ylim((self.pression[-1], self.pression[0]))
            y = self.pression

        else:
            self.axes.set_ylim(self.level.shape[0] + 1, 1)
            self.  axes.set_ylabel("level")
            y = self.level

        for i in range(0, len(theValues)):
            self.axes.plot(theValues[i], y, color=self.colors[i],
                           label="PC" + str(startPC + i + 1),
                           marker=theMarker, linestyle=self.linestyles[i])

        self.axes.set_xlabel(self.axeXLegend)
        self.axes.legend(prop={'size': 10})

        self.axes.grid(True)
        self.fig.canvas.draw()


class kpcView(util.GenericViewRadio):
    """ Surface window of the application """
    helpMessage = """

       On this window you can visualize profiles from the KPC matrix
       you can choose the fist profile to be displayed
       and how many profiles you can display.
       Be patient while waiting for your drawings...
      """
    helpPage = os.environ["RTTOV_GUI_PREFIX"] + "/doc/helpKPC.html"

    def __init__(self, parent, title="", profileList=None, baseProfile=None,
                 startProfile=0, nbProfile=10, kpcmatrixFileName=None,
                 inPressure=True, run=1):

        super().__init__(parent, title)
        self.YinPressure = inPressure
        self.profileList = profileList
        self.baseProfile = baseProfile
        if kpcmatrixFileName is not None:
            self.OpenKPCMATRIX(kpcmatrixFileName)
        self.pression = self.baseProfile['P']

        self.items = {}
        self.items["T"] = []
        self.items["Q"] = []
        self.items["CO2"] = []
        self.items['O3'] = []
        self.level = numpy.arange(1, self.profileList[0]['T'].shape[0] + 1, 1)
        self.nbpc = len(self.profileList)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        self.CreateMenuBar()
        self.SetSize((1020, 770))
        self.SetMinSize((1020, 770))
        self.SetTitle('K PC profile RUN ' + str(run))
        self.panel1 = wx.Panel(self, -1, style=wx.BORDER_SIMPLE)
        sizer.Add(self.panel1, 1, wx.EXPAND)

        for item in ("T", "Q", "O3", "CO2"):
            if item in self.profileList[0]:
                for pc in range(0, self.nbpc):
                    self.items[item].append(self.profileList[pc][item])

        self.startProfile = 0
        self.endProfile = startProfile + 10
        self.nb = None
        self.plot()
# create a second sizer for the notebook

        sizer1 = wx.BoxSizer()
        sizer1.Add(self.nb, 1, wx.EXPAND)

        self.panel1.SetSizer(sizer1)
        sizerRight = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(sizerRight)
        # panel 2 : sliders
        self.panel2 = wx.Panel(self, -1, style=wx.BORDER_SIMPLE)
        self.sizer2 = wx.BoxSizer(wx.VERTICAL)
        self.panel2.SetSizer(self.sizer2)

        sizerRight.Add(self.panel2)
        self.slider1 = wx.Slider(
            self.panel2, 10, 10, 1, 10, (30, 60), (250, -1),
            wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS
        )
        self.slider1.SetTickFreq(5)
        self.sizer2.Add(wx.StaticText(self.panel2, -1,
                                      "Number of K PC Profiles to show : "),
                        border=10)
        self.sizer2.Add(self.slider1, border=10)
        self.slider1.Bind(wx.EVT_SCROLL_THUMBRELEASE, self.RePlot)
        self.slider2 = wx.Slider(
            self.panel2, 10, 1, 1, self.nbpc, (30, 60), (250, -1),
            wx.SL_HORIZONTAL | wx.SL_AUTOTICKS | wx.SL_LABELS
        )
        self.slider2.SetTickFreq(5)
        self.slider2.Bind(wx.EVT_SCROLL_THUMBRELEASE, self.RePlot)
        # self.slider2.Bind(wx.EVT_SCROLL_CHANGED,self.RePlot)
        self.sizer2.Add(wx.StaticText(
            self.panel2, -1, "Start from :"), border=10)
        self.sizer2.Add(self.slider2, border=10)

        # panel3 + sizer3 : TSKIN
        self.panel3 = wx.lib.scrolledpanel.ScrolledPanel(
            self, -1, size=wx.Size(200, 600),
            style=wx.BORDER_SIMPLE | wx.EXPAND)
        self.panel3.SetupScrolling()
        sizerRight.Add(self.panel3, flag=wx.EXPAND)
        self.panel3.SetAutoLayout(1)
        self.sizer3 = wx.BoxSizer(wx.VERTICAL)
        self.panel3.SetSizer(self.sizer3)
        self.sizer3.Add(wx.StaticText(self.panel3, -1, "TSKIN"), border=5)
        self.printTskin()
        self.panel3.SetupScrolling()

        self.sb = self.CreateStatusBar()
        self.sb.SetBackgroundColour('WHITE')
        self.Centre()
        self.Show(True)
        self.StartShow = 0

    def plot(self):
        """ create the notebook and the graphics pages"""
        self.PCPageGraphic = {}
        if self.nb is None:
            self.nb = wx.Notebook(self.panel1, -1)
        else:
            self.nb.DeleteAllPages()
        for item in ("T", "Q", "O3", "CO2", "CO", "N2O", "CH4"):
            if item in self.profileList[0]:
                if self.profileList[0][item] is not None:
                    self.PCPageGraphic[item] = PlotPage(
                        self.nb,
                        self.items[item][self.startProfile:self.endProfile],
                        self.level, self.pression, item, self.nbpc,
                        axeXLegend=self.profileList[0][
                            item + '_ATTRIBUTE']['UNITS'],
                        inPressure=self.YinPressure, startPC=self.startProfile)
                    self.nb.AddPage(self.PCPageGraphic[item], item)

    def RePlot(self, e):
        self.BeginBusy()
        nb = self.slider1.GetValue()
        self.startProfile = self.slider2.GetValue() - 1
        self.endProfile = min(self.startProfile + nb, self.nbpc)
        self.write("Show  kp profile from " +
                   str(self.startProfile) + " to " + str(self.endProfile))
        self.ShowStartEnd(self.startProfile, self.endProfile)
        self.EndBusy()

    def ShowStartEnd(self, start, end):
        for item in list(self.PCPageGraphic.keys()):
            self.PCPageGraphic[item].OnPlot(
                self.items[item][start:end], startPC=start)

    def printTskin(self):
        for i in range(0, len(self.profileList)):
            tskin = self.profileList[i]["SKIN"]['T']
            self.sizer3.Add(wx.StaticText(
                self.panel3, -1,
                "tskin PC" + str(i + 1) + ": " + str(tskin)))

    def OnYPressions(self, e):
        self.YinPressure = True
        self.plot()

    def OnYLevels(self, e):
        self.YinPressure = False
        self.plot()

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu
                ('&Quit', 'Quit', self.OnQuit, "quit", True, False)),
               ("&Edit",  # Edit Menu
                ("Yaxis in pressure units", "put y in pressure units",
                 self.OnYPressions, "ypressions", True, True),
                ("Yaxis in level units", "put y in level unit",
                 self.OnYLevels, "ylevels", True, True)),
               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True, False),
                ("&Help", "Help", self.OnHelpHTML, "help", True, False)))

    def OpenKPCMATRIX(self, fileName):
        """ Open a KPC Matrix File and prepare
            the profile list and the baseProfile to be displayed"""

        f = h5py.File(fileName, 'r')
        h5 = f['/']
        kpc = rttov.kpcmatrix.Kpcmatrix()
        # Load kmatrix
        kpc.loadh5(h5)
        nbpc = kpc.kpcmatrix['T'].shape[1]
        self.write("display profiles " + fileName +
                   ' loaded. nbpc =' + str(nbpc))
        profile_list = []
        profile_kpc = rttov.profile.Profile()

        for pc in range(0, nbpc):
            profile_kpc = rttov.profile.Profile()
            profile_kpc = kpc.getkpcprof(pc)
            profile_list.append(profile_kpc)
        f.close()
        self.profileList = profile_list

        self.baseProfile = rmodel.project.OpenAProfile(fileName, 1)


if __name__ == "__main__":
    from util import rttov_gui_data_test_dir
    data_test_dir = rttov_gui_data_test_dir()
    fileName = os.path.join(data_test_dir, "pckmat.h5")
    p = rmodel.project.Project()
    ex = wx.App()
    sv = kpcView(None, "KPC View", None, None, 0, 10, fileName)
    ex.MainLoop()
