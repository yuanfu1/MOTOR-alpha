# -*- coding: utf-8 -*-

import wx
from rview import util
import rmodel
import os
import h5py
import rttov
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import \
    NavigationToolbar2WxAgg as ToolBar
from matplotlib.figure import Figure
from rview import colors


class PlotPage(wx.Panel):

    def __init__(self, parent, theValues, theName, theColor="Blue"):
        wx.Panel.__init__(self, parent, style=wx.BORDER_SIMPLE)
        self.colors = colors.btViewColors
        self.listChoiceRun = ["   try  1   "]

        self.listRefRun = ["ref try 1"]
        self.listDifRun = ["diff try"]
        self.dicDifRun = {}
        self.theName = theName
        self.theColor = theColor
        self.theValues = []
        self.theValues.append(theValues)
        self.fig = Figure()
        self.nbRun = 1

        self.canvas = FigureCanvas(self, -1, self.fig)
        self.sizer = wx.BoxSizer(wx.VERTICAL)

        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.GROW, 1)

        tlb = ToolBar(self.canvas)
        self.sizer.Add(tlb, 0, wx.GROW)
        tlb.Realize()
        self.SetSizer(self.sizer)
        self.Fit()
        self.OnPlot(self.theValues, 0)

    def onMouseMotion(self, event):

        txt = ""
        if event.inaxes:
            ax = event.inaxes.title.get_text()
            x = event.xdata
            y = event.ydata
            txt = 'X=%.1f Y=%.2f' % (x, y)
            self.text.SetLabel(txt)

    def OnPlot(self, theValues, i, theLabel="try 1"):
        if len(theValues[i].shape) == 0:
            n = 1
        else:
            n = theValues[i].shape[0]
        theMarker = None
        if n < 10:
            theMarker = "+"
        self.fig.clear()
        self.subP = self.fig.add_subplot(1, 1, 1)
        x = list(range(1, n + 1))

        theColor = self.theColor

        self.subP.plot(x, theValues[i], color=theColor,
                       label=theLabel, marker=theMarker)
        self.subP.set_ylabel("Y Axis (" + self.theName + ")")
        self.subP.set_xlabel("X Axis (Channels)")
        self.subP.grid(True)
        self.fig.canvas.draw()

    def PlotDif(self, diffRun):
        theMarker = None
        self.fig.clear()
        self.subP = self.fig.add_subplot(1, 1, 1)
        if len(self.theValues[0].shape) == 0:
            n = 1
        else:
            n = self.theValues[0].shape[0]
        if n < 10:
            theMarker = "+"
        x = list(range(1, n + 1))
        diffValues = self.theValues[
            diffRun[0] - 1] - self.theValues[diffRun[1] - 1]
        self.subP.plot(x, diffValues, color="red", label="", marker=theMarker)

        self.subP.set_ylabel("Y Axis (Differences)")
        self.subP.set_xlabel("X Axis (Principal Components)")
        self.subP.set_title(
            "try " + str(diffRun[0]) + " - " + "try " + str(diffRun[1]))
        self.subP.grid(True)
        self.fig.canvas.draw()

    def updateCbDifRun(self, refRun):

        self.listDifRun = []
        self.dicDifRun = {}
        self.listDifRun.append("diff try")
        nbRun = len(self.listChoiceRun)

        if nbRun > 1:
            for r in range(1, nbRun + 1):
                if r != refRun:
                    txt = "try " + str(r) + " - try " + str(refRun)
                    self.listDifRun.append(txt)
                    self.dicDifRun[txt] = (r, refRun)
        foo = self.cbDifRun.GetValue()
        self.cbDifRun.SetItems(self.listDifRun)
        self.cbDifRun.SetValue(foo)

    def RePlot(self, theValues, theRun=1):

        self.listChoiceRun.append("   try " + str(theRun) + "  ")
        self.listRefRun.append("ref try " + str(theRun) + " ")
        self.cbChoiceRun.SetItems(self.listChoiceRun[::-1])
        self.cbChoiceRun.SetValue(self.listChoiceRun[-1])

        foo = self.cbRefRun.GetValue()
        self.cbRefRun.SetItems(self.listRefRun)
        self.cbRefRun.SetValue(foo)

        refRun = self.listRefRun.index(str(self.cbRefRun.GetValue())) + 1
        self.updateCbDifRun(refRun)

        self.theValues.append(theValues)
        self.OnPlot(self.theValues, theRun - 1, "Run " + str(theRun))

    def OnChoiceRun(self, e):
        choiceRun = self.listChoiceRun.index(str(self.cbChoiceRun.GetValue()))

        self.OnPlot(self.theValues, choiceRun)

    def OnRefRun(self, e):
        refRun = self.listRefRun.index(str(self.cbRefRun.GetValue())) + 1
        self.updateCbDifRun(refRun)

    def OnDifRun(self, e):
        difRun = self.cbDifRun.GetValue()
        if difRun == "diff try":
            # redraw normal try
            self.OnChoiceRun(e)
        else:

            self.PlotDif(self.dicDifRun[difRun])


class RBtView(util.GenericView):
    """ Display difference in Brightness Temperature for
        the 1dvar controller """
    helpMessage = """

       TODO : write help message
      """
    helpPage = os.environ["RTTOV_GUI_PREFIX"] + "/doc/helpDiffRad.html"

    def __init__(self, parent, title, radFileNameT, radFileNameBg, counter=1):

        super().__init__(parent,  title)
        self.nbRun = 1

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        self.CreateMenuBar()
        self.SetSize((800, 600))
        self.SetMinSize((800, 600))
        self.SetTitle('True BT- Background BT')
        self.btdiff = None
        self.misc = None
        self.instrument = None
        self.radT = self.OpenRad(radFileNameT)
        self.radBg = self.OpenRad(radFileNameBg)
        if len(self.radT["BT"].shape) == 0:
            self.nchan = 1
        else:
            self.nchan = self.radT["BT"].shape[0]
        self.btdiff = self.radT["BT"] - self.radBg["BT"]

        self.instrument = self.radT.misc['INSTRUMENT']

        self.write("instrument " + str(self.instrument))
        self.panel1 = wx.Panel(
            self, -1, style=wx.BORDER_SIMPLE, size=(200, 500))
        sizer.Add(self.panel1, 1, wx.EXPAND)
        nb = wx.Notebook(self.panel1, -1)

        self.PCPageGraphic = PlotPage(nb, self.btdiff, "True -Bg",
                                      theColor=colors.btViewColors["T-Bg"])
        self.PCPageGraphicT = PlotPage(nb, self.radT['BT'], "True",
                                       theColor=colors.btViewColors["True"])
        self.PCPageGraphicBg = PlotPage(nb, self.radBg["BT"], "Bg",
                                        theColor=colors.btViewColors["Bg"])
        # add the page to the notebook
        nb.AddPage(self.PCPageGraphic, "True-Bg")
        nb.AddPage(self.PCPageGraphicT, "True")
        nb.AddPage(self.PCPageGraphicBg, "Bg")

        # create a second sizer for the notebook

        sizer2 = wx.BoxSizer()
        sizer2.Add(nb, 1, wx.EXPAND)
        self.panel1.SetSizer(sizer2)
        self.sb = self.CreateStatusBar()
        self.Centre()
        self.Show(True)

    def getNbChan(self):
        return self.nchan

    def ReRead(self, radFileName, radFileNameYXb):

        self.radTrue = self.OpenRad(radFileName)
        self.radBg = self.OpenRad(radFileNameYXb)
        self.btdiff = self.radTrue["BT"] - self.radBg["BT"]

        if len(self.btdiff.shape) == 0:
            nchan = 1
        else:
            nchan = self.btdiff.shape[0]

        if nchan != self.nchan:
            txt = 'WARNING : nb channels : %s --> %s' % (self.nchan, nchan)
            self.writeSB(txt, "red", 10, 1)
            self.nbRun = 1
            return
        if str(self.radTrue.misc['INSTRUMENT']) != self.instrument:
            txt = 'WARNING :isntrument : %s --> %s' % (
                self.instrument, str(self.misc['INSTRUMENT']))
            self.writeSB(txt, "red", 10, 1)
            self.nbRun = 1
            return
        self.nbRun += 1
        self.PCPageGraphic.RePlot(self.btdiff, self.nbRun)
        self.PCPageGraphicT.RePlot(self.radTrue["BT"], self.nbRun)
        self.PCPageGraphicBg.RePlot(self.radBg["BT"], self.nbRun)

    def OpenRad(self, radFileName):
        self.write("Open radiance file : " + radFileName)
        try:
            rad = rttov.radiance.Radiance()
            rad.read(radFileName)
            return rad

        except:
            txt = 'error access file : %s' % radFileName
            self.writeSB(txt, 'RED', 10, 1)
            self.nbRun = -1
            return None

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu

                ('&Quit', 'Quit', self.OnQuit, "quit", True)),
               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True),
                ("&Help", "Help", self.OnHelpHTML, "help", True)))


if __name__ == "__main__":
    p = rmodel.project.Project()
    print("Configuration : ", p.config.ENV['RTTOV_GUI_PREFIX'])
    ex = wx.App()

    fileNameT = '/home/pascale/.rttov/radrTrue.h5'
    fileNameBg = '/home/pascale/.rttov/radrBg.h5'

    sv = RBtView(None,  fileNameT, fileNameBg, 1)

    ex.MainLoop()
