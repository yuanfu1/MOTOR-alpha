# -*- coding: utf-8 -*-
import wx
from rview import util
import rmodel
import os
import numpy
import h5py
import rttov
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wxagg import \
    NavigationToolbar2WxAgg as ToolBar
from matplotlib.figure import Figure
from rview import colors


class PlotPage(wx.Panel):

    def __init__(self, parent, theValues, theName):
        wx.Panel.__init__(self, parent, style=wx.BORDER_SIMPLE)
        self.colors = colors.pcViewColors
        self.listChoiceRun = ["   run  1   "]

        self.listRefRun = ["ref run 1"]
        self.listDifRun = ["diff run"]
        self.dicDifRun = {}
        self.theName = theName
        self.theValues = []
        self.theValues.append(theValues)
        self.fig = Figure()
        self.nbRun = 1

        self.canvas = FigureCanvas(self, -1, self.fig)
        self.canvas.mpl_connect('motion_notify_event', self.onMouseMotion)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        horSizer = wx.BoxSizer(wx.HORIZONTAL)
        self.cbChoiceRun = wx.ComboBox(
            self, choices=self.listChoiceRun, size=(180, 26))
        self.cbChoiceRun.SetValue(self.listChoiceRun[0])
        self.cbRefRun = wx.ComboBox(
            self, choices=self.listRefRun, size=(180, 26))
        self.cbRefRun.SetValue(self.listRefRun[0])
        self.cbDifRun = wx.ComboBox(
            self, choices=self.listDifRun, size=(180, 26))
        self.cbDifRun.SetValue(self.listDifRun[0])
        self.cbChoiceRun.Bind(wx.EVT_COMBOBOX, self.OnChoiceRun)
        self.cbRefRun.Bind(wx.EVT_COMBOBOX, self.OnRefRun)
        self.cbDifRun.Bind(wx.EVT_COMBOBOX, self.OnDifRun)
        horSizer.Add(self.cbChoiceRun, border=1)
        horSizer.Add(self.cbRefRun, border=1)
        horSizer.Add(self.cbDifRun, border=1)
        self.sizer.Add(horSizer, border=1)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.GROW, border=2)
        self.sizer.Add
        tlb = ToolBar(self.canvas)
        self.sizer.Add(tlb, 0, wx.GROW)
        tlb.Realize()
        self.SetSizer(self.sizer)

        self.text = wx.StaticText(self, -1, label="")
        self.sizer.Add(self.text)
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

    def OnPlot(self, theValues, i, theLabel="run 1"):
        n = theValues[i].shape[0]
        theMarker = None
        self.fig.clear()
        self.subP = self.fig.add_subplot(1, 1, 1)
        x = numpy.arange(1, n + 1)
        if i >= len(self.colors):
            theColor = "black"
        else:
            theColor = self.colors[i]
        self.subP.plot(x, theValues[i], color=theColor,
                       label=theLabel, marker=theMarker)
        self.subP.set_ylabel("Y Axis (" + self.theName + ")")
        self.subP.set_xlabel("X Axis (Principal Components)")
        self.subP.set_title("run " + str(i + 1))
        self.subP.set_yscale('symlog', linthreshy=100.0)
        self.subP.plot([0, n], [100, 100], 'k--')
        self.subP.plot([0, n], [-100, -100], 'k--')
        self.subP.grid(True)
        self.fig.canvas.draw()

    def PlotDif(self, diffRun):
        theMarker = None
        self.fig.clear()
        self.subP = self.fig.add_subplot(1, 1, 1)
        n = self.theValues[0].shape[0]
        x = numpy.arange(1, n + 1)
        diffValues = self.theValues[
            diffRun[0] - 1] - self.theValues[diffRun[1] - 1]
        self.subP.plot(x, diffValues, color="red", label="", marker=theMarker)

        self.subP.set_ylabel("Y Axis (Differences)")
        self.subP.set_xlabel("X Axis (Principal Components)")
        self.subP.set_title(
            "run " + str(diffRun[0]) + " - " + "run " + str(diffRun[1]))
        self.subP.grid(True)
        self.fig.canvas.draw()

    def updateCbDifRun(self, refRun):

        self.listDifRun = []
        self.dicDifRun = {}
        self.listDifRun.append("diff run")
        nbRun = len(self.listChoiceRun)

        if nbRun > 1:
            for r in range(1, nbRun + 1):
                if r != refRun:
                    txt = "run " + str(r) + " - run " + str(refRun)
                    self.listDifRun.append(txt)
                    self.dicDifRun[txt] = (r, refRun)
        foo = self.cbDifRun.GetValue()
        self.cbDifRun.SetItems(self.listDifRun)
        self.cbDifRun.SetValue(foo)

    def RePlot(self, theValues, theRun=1):
        self.listChoiceRun.append("   run " + str(theRun) + "  ")
        self.listRefRun.append("ref run " + str(theRun) + " ")
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
        # print "plot",self.cbChoiceRun.GetValue(), choiceRun
        self.OnPlot(self.theValues, choiceRun)

    def OnRefRun(self, e):
        refRun = self.listRefRun.index(str(self.cbRefRun.GetValue())) + 1
        self.updateCbDifRun(refRun)

    def OnDifRun(self, e):
        difRun = self.cbDifRun.GetValue()
        if difRun == "diff run":
            # redraw normal run
            self.OnChoiceRun(e)
        else:
            # print "diff ",difRun
            self.PlotDif(self.dicDifRun[difRun])


class pcView(util.GenericView):
    """ Surface window of the application """
    helpMessage = """

       TODO : write help message
      """
    helpPage = os.environ["RTTOV_GUI_PREFIX"] + "/doc/helpPC.html"

    def __init__(self, parent, title, pcFileName):

        super().__init__(parent,  title)
        self.nbRun = 1

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.SetSizer(sizer)
        self.CreateMenuBar()
        self.SetSize((800, 600))
        self.SetMinSize((800, 600))
        self.pc = None
        self.misc = None
        self.instrument = None
        self.OpenPC(pcFileName)
        self.instrument = self.misc['INSTRUMENT']
        self.write("instrument " + str(self.instrument))
        self.panel1 = wx.Panel(
            self, -1, style=wx.BORDER_SIMPLE, size=(200, 500))
        sizer.Add(self.panel1, 1, wx.EXPAND)
        nb = wx.Notebook(self.panel1, -1)
        self.nbpc = self.pc['TOTAL_PCSCORES'].shape[0]
        self.PCPageGraphic = PlotPage(
            nb, self.pc['TOTAL_PCSCORES'], "TOTAL_PCSCORES")

        # add the page to the notebook
        nb.AddPage(self.PCPageGraphic, "TOTAL_PCSCORES")
        # nb.AddPage(self.BTPageGraphic,"BT_PCCOMP")
        # create a second sizer for the notebook

        sizer2 = wx.BoxSizer()
        sizer2.Add(nb, 1, wx.EXPAND)
        self.panel1.SetSizer(sizer2)
        self.sb = self.CreateStatusBar()
        self.Centre()
        self.Show(True)

    def ReRead(self, pcfileName):
        self.OpenPC(pcfileName)
        nbpc = self.pc['TOTAL_PCSCORES'].shape[0]
        if nbpc != self.nbpc:
            txt = 'WARNING : nb pc scores : %s --> %s' % (self.nbpc, nbpc)
            self.writeSB(txt, "red", 10, 1)
            self.nbRun = 1
            return
        if str(self.misc['INSTRUMENT']) != str(self.instrument):
            txt = 'WARNING :isntrument : %s --> %s' % (
                self.instrument, str(self.misc['INSTRUMENT']))
            self.writeSB(txt, "red", 10, 1)
            self.nbRun = 1
            return
        self.nbRun += 1
        print("number of PC run=", self.nbRun)
        self.PCPageGraphic.RePlot(self.pc['TOTAL_PCSCORES'], self.nbRun)

    def Refresh(self):
        pass

    def OpenPC(self, pcFileName):
        self.write("Open pc file : " + pcFileName)
        try:
            f = h5py.File(pcFileName, 'r')
            # get the Dataset
            h5 = f['/PCCOMP/']
            self.pc = rttov.pccomp.PCCOMP()
            # Load
            self.pc.loadh5(h5)

            # misc
            h5 = f['/MISC/']
            self.misc = rttov.pccomp.PCMISC()
            self.misc.loadh5(h5)

            # Close HDF file
            f.close()

        except:
            txt = 'error access file : %s' % pcFileName
            self.writeSB(txt, 'RED', 10, 1)
            self.nbRun = -1
            return

    def MenuData(self):
        """ define the data for the menu
        """
        return(("&File",  # File Menu

                ('&Quit', 'Quit', self.OnQuit, "quit", True)),
               ("&Help",  # Help Menu
                ("About", "About screen", self.OnAbout, "about", True),
                ("&Help", "Help", self.OnHelpHTML, "help", True)))


if __name__ == "__main__":
    from util import rttov_gui_data_test_dir
    data_test_dir = rttov_gui_data_test_dir()
    pcfile1 = os.path.join(data_test_dir, "pc1.h5")
    pcfile2 = os.path.join(data_test_dir, "pc2.h5")
    pcfile3 = os.path.join(data_test_dir, "pc3.h5")
    p = rmodel.project.Project()
    ex = wx.App()

    f = h5py.File(pcfile1, 'r')

    # get the Dataset

    # Display option
    print("PC-----------------------")

    sv = pcView(None, "PC SCORES", pcfile1)
    sv.ReRead(pcfile2)
    sv.ReRead(pcfile3)
    ex.MainLoop()
