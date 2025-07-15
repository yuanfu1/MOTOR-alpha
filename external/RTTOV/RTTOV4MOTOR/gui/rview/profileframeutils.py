# -*- coding: utf-8 -*-

from . import wxmpl
import numpy
import matplotlib
from . import layeritem
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import\
    NavigationToolbar2WxAgg as ToolBar
from matplotlib.figure import Figure
import wx.lib.agw.flatnotebook
from . import colors
from .util import axesDef
import copy

itemColor = colors.profileItemColors

itemMarker = colors.profileItemMarkers


class Data(wxmpl.Channel):
    """ Class to keep data and style for the PlotItemPanel data """

    def __init__(self, x, y, theName="", theColor=None, theStyle="-",
                 theMarker=None, theLabel=None):
        wxmpl.Channel.__init__(
            self, name=theName, color=theColor, style=theStyle,
            marker=theMarker)
        self.x = x
        self.y = y
        self.name = theName
        self.color = theColor
        self.style = theStyle
        self.marker = theMarker
        if theLabel is None:
            self.label = theName
        else:
            self.label = theLabel
        self.changed = False

    def getX(self):
        return self.x

    def getY(self):
        return self.y

    def myUpdate(self, x, y=None):
        self.x = numpy.zeros(x.shape[0]) + x
        if y is not None:
            self.y = numpy.zeros(y.shape[0]) + y
        self.changed = True


class GenericPlotItemPanel(wx.Panel):
    """ plot on a PlotPanel one curve """

    def __init__(self, parent, value, pression, theName,
                 kind="GASES", xlegend="ppmv", edit=False,
                 layerstyle=False, layer=None, yInPressions=True, tskin=None,
                 tickSize=10, xlabel=""):

        self.theName = theName
        self.theParent = parent
        self.xlegend = xlegend
        self.xlabel = xlabel
        self.edit = edit
        self.kind = kind
        self.yInPressions = yInPressions
        self.layer = layer
        self.layerstyle = layerstyle
        self.tickSize = tickSize
        self.pression = pression
        self.value = value
        self.myLayeritem = None
        self.tskin = []
        self.ytskin = []
        if tskin:
            self.tskin.append(tskin)

        wx.Panel.__init__(self, parent, style=wx.BORDER_SIMPLE)

        # define object for matplotlib
        self.fig = Figure()
        self.canvas = FigureCanvas(self, -1, self.fig)
        self.canvas.mpl_connect('motion_notify_event', self.onMouseMotion)

        self.text = wx.StaticText(self, -1, label="")
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.GROW, 1)
        self.tlb = ToolBar(self.canvas)
        self.sizer.Add(self.tlb, 0, wx.GROW)
        self.tlb.Realize()
        self.SetSizer(self.sizer)

        self.text = wx.StaticText(self, -1, label="")
        self.sizer.Add(self.text)
        self.Fit()

        self.onInsert = True
        self.myCurves = []
        self.PrepareAxeAndData()
        self.valueHistory = []
        self.valueHistoryRedo = []
        self.cid = None

    def onMouseMotion(self, event):
        """ set text when moving mouse """

        if event.inaxes:

            xdata = event.xdata
            ydata = event.ydata
            xstr = "%0.4g" % xdata
            ystr = "%0.4g" % ydata

            value = self.axes.get_ylabel() + "=" + ystr + \
                "  " + self.xlabel + "=" + xstr

            self.text.SetLabel(value)

    def MyReleaseMouse(self):
        try:
            if self.canvas.HasCapture():
                self.canvas.ReleaseMouse()
            if self.HasCapture():
                self.ReleaseMouse()
        except:
            pass

    def PrepareAxeAndData(self):
        """ effectively perform the graphics """
        self.SetTickSize(self.tickSize)

        self.fig.clear()
        self.axes = self.fig.add_subplot(1, 1, 1)

        self.x = self.value[::1]

        if self.yInPressions:
            self.axes.set_yscale("log")
            self.axes.set_yticks((0.00005, 0.0001, 0.0002, 0.0005, 0.001,
                                  0.002, 0.005, 0.01, 0.02, 0.05,
                                  0.1, 0.2, 0.5, 1, 2, 5, 10, 25, 50,
                                  100, 200, 300, 500, 1000))
            label = ('5e-5', '1e-4', '2e-4', '5e-4', '1e-3',
                     '2e-3', '5e-3', '0.01', '0.02', '0.05',
                     '0.1', '0.2', '0.5', '1', '2', '5', '10', '25', '50',
                     '100', '200', '300', '500', '1000')
            self.axes.set_yticklabels(label)
            self.axes.set_ylabel('Pressure (hPa)')
            self.axes.set_ylim((self.pression[-1] + 150, self.pression[0]))
        else:
            self.axes.set_ylim(self.value.shape[0] + 2, 1)
            self.axes.set_ylabel('Level')

        if self.kind == "GASES":
            if self.yInPressions:
                self.y = self.pression[::1]
            else:
                self.y = numpy.arange(1, self.value.shape[0] + 1, 1)
        else:
            if self.yInPressions:
                self.y = self.layer[::1]
            else:
                self.y = numpy.arange(1.5, self.value.shape[0], 1)

        if not self.layerstyle:
            self.data = Data(self.x, self.y, theName=self.theName,
                             theColor=itemColor[self.theName],
                             theMarker=itemMarker[self.theName])
        else:
            if self.yInPressions:
                self.myLayeritem = layeritem.Layeritem(
                    layeritem=self.x, pression=self.pression[::1])
            else:
                self.myLayeritem = layeritem.Layeritem(
                    layeritem=self.x, pression=numpy.arange(
                        1, self.value.shape[0] + 1, 1))
            (self.xlayeritem, self.ylayeritem) = (
                self.myLayeritem.computeLayerLine(layers=self.y))
            self.data = Data(self.xlayeritem, self.ylayeritem,
                             theName=self.theName,
                             theColor=itemColor[self.theName],
                             theMarker=itemMarker[self.theName],
                             theLabel=self.theName)

        self.axes.set_xlabel(self.xlegend)
        formatter = self.axes.xaxis.get_major_formatter()
        formatter.set_powerlimits((-3, 4))
        self.axes.xaxis.set_major_formatter(formatter)
        self.axes.grid(True, axis='both')

        self.myChannelList = []
        self.myChannelList.append(self.data)

        if self.theName == "T":
            if len(self.tskin) > 0:
                if self.yInPressions:
                    self.ytskin.append(self.pression[-1] + 50)
                else:
                    self.ytskin.append(self.value.shape[0] + 1)
                datatskin = Data(self.tskin, self.ytskin,
                                 theName='TSKIN', theColor="red",
                                 theMarker="*", theLabel='Tskin')
                self.myChannelList.append(datatskin)

        if wx.Platform == '__WXMAC__':
            self.PlotCurves()

    def SetTickSize(self, size):
        matplotlib.rc('xtick', labelsize=size)
        matplotlib.rc('ytick', labelsize=size)

    def ConnectCanvasEVT_POINT(self, methode):
        self.cid = self.fig.canvas.mpl_connect("button_press_event", methode)
        self.MyReleaseMouse()

    def DisconnectCanvasEVT_POINT(self):
        if self.cid is not None:
            self.fig.canvas.mpl_disconnect(self.cid)
            self.MyReleaseMouse()

    def SetXlimits(self,  xmin=None, xmax=None, data=None):
        """ set x limits """
        if xmin is not None and xmax is not None:
            self.axes.set_xlim((xmin, xmax))
        elif data is not None:
            xmin = data.min()
            xmax = data.max()
            gap = xmax - xmin
            self.axes.set_xlim((xmin, xmax + gap / 10.))

        self.MyReleaseMouse()

    def PlotCurves(self):
        """ erase the curve if necessary and redraw """
        if len(self.myCurves) == 1:
            if len(self.axes.lines) == 1:
                self.axes.lines.remove(self.axes.lines[0])
            self.myCurves.pop()

        for data in self.myChannelList:
            c = self.axes.plot(
                data.x, data.y, color=data.color, marker=data.marker)
            self.myCurves.append(c)
        self.fig.canvas.draw_idle()
        self.MyReleaseMouse()

    def UpdateDataAndRePlot(self, dataX):
        self.x = dataX
        self.data.setChanged(True)
        if not self.layerstyle:
            self.data.myUpdate(self.x, self.y)
        else:
            (self.xlayeritem, self.ylayeritem) = (
                self.myLayeritem.computeLayerLine(layeritem=dataX))
            self.data.myUpdate(self.xlayeritem, self.ylayeritem)
        self.PlotCurves()

    def OnRedo(self):
        if self.valueHistoryRedo != []:
            if not self.layerstyle:
                X = numpy.zeros(self.x.shape[0]) + self.x
                self.valueHistory.append(X)
                X = self.valueHistoryRedo.pop()
                self.x = numpy.zeros(X.shape[0]) + X
                self.data.myUpdate(self.x, self.y)
            else:
                X = numpy.zeros(self.xlayeritem.shape[0]) + self.xlayeritem
                self.valueHistory.append(X)
                X = self.valueHistoryRedo.pop()
                self.xlayeritem = numpy.zeros(X.shape[0]) + X
                self.x = self.myLayeritem.getLayeritem(self.xlayeritem)
                self.myLayeritem.update(self.xlayeritem, self.ylayeritem)
                self.data.myUpdate(self.xlayeritem, self.ylayeritem)
            self.PlotCurves()
        self.MyReleaseMouse()

    def OnUndo(self):

        if self.valueHistory != []:
            if not self.layerstyle:
                X = numpy.zeros(self.x.shape[0]) + self.x
                self.valueHistoryRedo.append(X)
                X = self.valueHistory.pop()
                self.x = numpy.zeros(X.shape[0]) + X
                self.data.myUpdate(self.x, self.y)
            else:
                X = numpy.zeros(self.xlayeritem.shape[0]) + self.xlayeritem
                self.valueHistoryRedo.append(X)
                X = self.valueHistory.pop()
                self.xlayeritem = numpy.zeros(X.shape[0]) + X
                self.x = self.myLayeritem.getLayeritem(self.xlayeritem)
                self.data.myUpdate(self.xlayeritem, self.ylayeritem)

            self.PlotCurves()
        self.MyReleaseMouse()

    def OnPoint(self, e):
        """ OnPoint Methods """
        self.MyReleaseMouse()
        if (e.button == 1) or (e.dblclick):
            self.MyReleaseMouse()
            return(False)
        if e.xdata is None or e.ydata is None:
            self.MyReleaseMouse()
            return False
        if (e.ydata < self.y.min() or e.ydata > self.y.max()):
            self.MyReleaseMouse()
            return(False)

        if not self.layerstyle:
            y = numpy.zeros(self.x.shape[0]) + self.x
            self.valueHistory.append(y)
            mini = 1000
            for index in range(self.y.shape[0]):
                dist = abs(self.y[index] - e.ydata)
                if dist < mini:
                    imin = index
                    mini = dist
            if self.kind != "GASES" and not self.onInsert:
                self.x[imin] = 0
            else:
                self.x[imin] = e.xdata
            self.data.setChanged(True)
            self.data.myUpdate(self.x, self.y)
            self.PlotCurves()

        else:

            y = numpy.zeros(self.xlayeritem.shape[0]) + self.xlayeritem
            self.valueHistory.append(y)
            mini = 1000
            for index in range(self.ylayeritem.shape[0]):
                dist = self.ylayeritem[index] - e.ydata
                if dist < mini and dist > 0:
                    imin = index
                    mini = dist

            if not self.onInsert:
                self.xlayeritem[imin] = 0
                # we have 2 points to move and its depends if imin is odd
                if imin % 2 != 0:
                    if imin != self.xlayeritem.shape[0]:
                        self.xlayeritem[imin + 1] = 0
                else:
                    if imin != 0:
                        self.xlayeritem[imin - 1] = 0
            else:
                self.xlayeritem[imin] = e.xdata
                # we have 2 points to move and its depends if imini is odd
                if imin % 2 != 0:
                    if imin != self.xlayeritem.shape[0]:
                        self.xlayeritem[imin + 1] = e.xdata
                else:
                    if imin != 0:
                        self.xlayeritem[imin - 1] = e.xdata

            self.data.setChanged(True)
            self.data.myUpdate(self.xlayeritem, self.ylayeritem)
            self.x = self.myLayeritem.getLayeritem(self.xlayeritem)
            self.PlotCurves()

        self.MyReleaseMouse()
        return True

    def GetItem(self):
        """ get value from curve (=data) and return a profile for the item """
        self.MyReleaseMouse()
        myX = self.data.getX()
        if self.layerstyle:
            layerX = self.myLayeritem.getLayeritem(myX)
            return layerX
        else:
            return myX


class PlotItemPanelAll(wx.Panel):
    """ Plot all gas in the same graphic """

    def __init__(self, parent, theProfile, kind="GASES", layer=None,
                 xlegendT=None,
                 xlegend=None, yInPressions=True, addTskin=False,
                 XinLog=False, tickSize=8, runNumber=None):
        wx.Panel.__init__(self, parent, style=wx.BORDER_SIMPLE)

        # define object for matplotlib
        self.fig = Figure()
        self.canvas = FigureCanvas(self, -1, self.fig)
        self.text = wx.StaticText(self, -1, label="")
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.GROW, 1)
        self.tlb = ToolBar(self.canvas)
        self.sizer.Add(self.tlb, 0, wx.GROW)
        self.tlb.Realize()
        self.text = wx.StaticText(self, -1, label="")
        self.sizer.Add(self.text)
        self.SetSizer(self.sizer)
        self.Fit()
        self.runNumber = runNumber
        self.myProfileRef = theProfile
        self.myProfile = copy.deepcopy(theProfile)
        self.pression = self.myProfile['P']
        self.layer = layer
        self.layerlevel = numpy.arange(1.5, self.myProfile['T'].shape[0], 1)
        self.tickSize = tickSize

        self.kind = kind
        self.theName = "all " + kind
        self.yInPressions = yInPressions
        self.itemsList = {}
        self.itemsList['GASES'] = self.myProfile.gas_list
        self.itemsList['AEROSOLS'] = self.myProfile.aerosol_list
        self.itemsList['CLOUDS'] = self.myProfile.cloud_list
        self.labelx = {}
        self.labelx[kind] = xlegend
        self.labelxT = xlegendT

        self.tskin = []
        self.ytskin = []
        if addTskin:
            self.tskin.append(self.myProfile['SKIN']['T'])
        self.XinLog = XinLog
        self.myCurves = []
        self.PrepareAxeAndData()

    def SetStatusText(self, text, color=None):
        self.text.SetLabel(text)
        if color is not None:
            self.text.SetBackgroundColour(color)

    def SetTickSize(self, size):
        matplotlib.rc('xtick', labelsize=size)
        matplotlib.rc('ytick', labelsize=size)

    def SetXlimits(self,  xmin, xmax):
        """ set x limits """
        if axesDef[self.kind]["xlimits"] is not None:
            self.axes.set_xlim(axesDef[self.theKind]["xlimits"])
        else:
            gap = xmax - xmin
            newxmax = xmax + gap / 10.
            self.axes.set_xlim((xmin, newxmax))

    def PrepareAxeAndData(self, theProfile=None):

        if theProfile is not None:
            self.myProfileRef = theProfile
            self.myProfile = copy.deepcopy(theProfile)

        self.SetTickSize(self.tickSize)
        fig = self.fig
        fig.clear()
        self.axes = fig.gca()
        if self.XinLog:
            self.axes.set_xscale("log")

        self.data = {}

        if self.yInPressions:
            if len(self.tskin) >= 1:
                self.ytskin.append(self.pression[-1] + 50)
            self.axes.set_yscale("log")
            self.axes.set_yticks((0.00005, 0.0001, 0.0002, 0.0005, 0.001,
                                  0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2,
                                  0.5, 1, 2, 5, 10, 25, 50, 100,
                                  200, 300, 500, 1000))
            label = ('5e-5', '1e-4', '2e-4', '5e-4', '1e-3',
                     '2e-3', '5e-3', '0.01', '0.02', '0.05', '0.1', '0.2',
                     '0.5', '1', '2', '5', '10', '25', '50', '100',
                     '200', '300', '500', '1000')

            self.axes.set_yticklabels(label)
            self.axes.set_ylabel('Pressure (hPa)')
            self.axes.set_ylim((self.pression[-1], self.pression[0]))
            if self.kind == 'GASES':
                y = self.myProfile['P'][::1]
            else:
                y = self.layer[::1]
            ytempe = self.myProfile['P'][::1]
        else:
            self.axes.set_ylabel('Level')
            y = numpy.arange(1, self.myProfile['T'].shape[0] + 1, 1)
            if len(self.tskin) >= 1:
                self.ytskin.append(self.myProfile['T'].shape[0] + 1)
            if self.kind != 'GASES':
                y = self.layerlevel
            ytempe = numpy.arange(1, self.myProfile['T'].shape[0] + 1, 1)
            self.axes.set_ylim(self.myProfile['T'].shape[0] + 2, 1)

        self.axes.set_xlabel(self.labelx[self.kind])

        if self.kind == "GASES":
            if self.labelx["GASES"] is None:
                self.axes.set_xlabel(self.myProfile['Q_ATTRIBUTE']['UNITS'])
            markerFlag = False

        else:
            markerFlag = True

        if self.kind == "AEROSOLS":
            if self.labelx["AEROSOLS"] is None:
                for aero in self.myProfile.aerosol_list:
                    if self.myProfile[aero] is not None and self.myProfile[
                            aero + '_ATTRIBUTE']['UNITS'] != 'n/a':
                        self.axes.set_xlabel(
                            self.myProfile[aero + '_ATTRIBUTE']['UNITS'])
            else:
                self.axes.set_xlabel(self.labelx["AEROSOLS"])

        if self.kind == "CLOUDS":
            if self.labelx["CLOUDS"] is None:
                for cloud in self.myProfile.cloud_list:
                    if self.myProfile[cloud] is not None and self.myProfile[
                            cloud + '_ATTRIBUTE']['UNITS'] != 'n/a':
                        self.axes.set_xlabel(
                            self.myProfile[cloud + '_ATTRIBUTE']['UNITS'])
            else:
                self.axes.set_xlabel(self.labelx["CLOUDS"])

        self.axesT = self.axes.twiny()
        if self.labelxT is None:
            self.axesT.set_xlabel(self.myProfile['T_ATTRIBUTE']['UNITS'])
        else:
            self.axesT.set_xlabel(self.labelxT)

        # must redefine ylim for the levels case
        if self.yInPressions:
            self.axes.set_ylim((self.pression[-1] + 150, self.pression[0]))
        else:
            self.axes.set_ylim(self.myProfile['T'].shape[0] + 2, 1)

        self.stripCharterT = wxmpl.StripCharter(
            self.axesT, loc_legend="upper right")

        self.myChannelList = []
        self.myChannelListT = []
        xmax = 0.
        xmin = 9.e20
        self.axes.grid(True, axis='both')

        for item in self.itemsList[self.kind]:
            if self.myProfile[item] is not None:
                x = self.myProfile[item][::1]
                xmax = max(xmax, x.max())
                xmin = min(xmin, x.min())
                if markerFlag:
                    self.data[item] = Data(x, y, theName=item,
                                           theColor=itemColor[item],
                                           theMarker=itemMarker[item])
                else:
                    self.data[item] = Data(
                        x, y, theName=item, theColor=itemColor[item])
                self.myChannelList.append(self.data[item])

        x = self.myProfile['T'][::1]
        self.data['T'] = Data(x, ytempe, theName='T', theColor=itemColor['T'])
        self.myChannelListT.append(self.data['T'])

        if len(self.tskin) >= 1:
            self.data['TSKIN'] = Data(
                self.tskin, self.ytskin, theName='TSKIN', theColor="red",
                theMarker="*")
            self.myChannelListT.append(self.data['TSKIN'])

        self.SetXlimits(xmin, xmax)

        if not self.XinLog:
            formatter = self.axes.xaxis.get_major_formatter()
            formatter.set_powerlimits((-3, 4))
            self.axes.xaxis.set_major_formatter(formatter)

            formatter = self.axesT.xaxis.get_major_formatter()
            formatter.set_powerlimits((-3, 4))
            self.axesT.xaxis.set_major_formatter(formatter)

        self.PlotCurves()

    def PlotCurves(self):
        """ erase the curve if necessary and redraw """

        lenAll = len(self.myCurves)
        if lenAll > 0:
            # update existing curves
            for l in range(0, lenAll):
                c = self.myCurves[l]["line"]
                data = self.myCurves[l]["data"]
                c.set_xdata(data.x)
                self.fig.canvas.draw_idle()
        else:
            # first time plot
            for data in self.myChannelList:
                c, = self.axes.plot(
                    data.x, data.y, color=data.color, marker=data.marker,
                    label=data.label)
                self.myCurves.append({"line": c, "data": data})
            self.axes.legend(loc="upper right", fontsize=10)

            for data in self.myChannelListT:
                c, = self.axesT.plot(
                    data.x, data.y, color=data.color, marker=data.marker,
                    label=data.label)
                self.myCurves.append({"line": c, "data": data})
            self.axesT.legend(loc="upper left", fontsize=10)
            self.fig.canvas.draw_idle()

    def UpdateDataAndRePlot(self, profile, itemName):
        """ PlotCurves the data with the profile and update the graphic
            this method can't be used if there is a modification in
            number of gas or number of layers
            use PrepareAxeAndData in this case """
        self.myProfileRef = profile
        self.myProfile = copy.deepcopy(profile)
        if self.myProfile[itemName] is not None:
            if (
                    self.kind != 'GASES' and itemName != 'T' and
                    itemName not in ("CFRAC", "ICEDE", "CLW", "CLWDE")):
                self.data[itemName].myUpdate(
                    self.myProfile[itemName], self.layer)
            else:
                if itemName not in ("CFRAC", "ICEDE", "CLW", "CLWDE"):
                    self.data[itemName].myUpdate(
                        self.myProfile[itemName], self.myProfile['P'])
        self.PlotCurves()


class MyNotebook(wx.lib.agw.flatnotebook.FlatNotebook):
    """ MyNotebook Class inherits from Notebook
        Bind the event EVT_NOTEBOOK_PAGE_CHANGED
        with PlotCurves method of the page """

    def __init__(self, parent, isRightPage=True):
        #        wx.Notebook.__init__(self, parent, id=wx.ID_ANY, style=
        #                             wx.BK_DEFAULT
        #                             )
        wx.lib.agw.flatnotebook.FlatNotebook.__init__(
            self, parent, id=wx.ID_ANY,
            agwStyle=(wx.lib.agw.flatnotebook.FNB_NO_X_BUTTON |
                      wx.lib.agw.flatnotebook.FNB_FANCY_TABS))
        self.parent = parent
        self.isRightPage = isRightPage
        self.BindEvent()

    def BindEvent(self):
        self.Bind(wx.EVT_NOTEBOOK_PAGE_CHANGED, self.OnPageChanged)

    def OnPageChanged(self, event):
        sel = self.GetSelection()
        thePage = self.GetPage(sel)
        thePage.PlotCurves()
        event.Skip()
