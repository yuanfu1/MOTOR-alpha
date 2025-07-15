# -*- coding: utf-8 -*-
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_wx import _load_bitmap
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigureCanvas
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import h5py
import numpy as np
import os
import sys
import rttov
import rview
import time
from rview import util
import wx
import rmodel
from rview.util import MyCustomToolbar

# 2014-01-29 : cas avec T et Q (pas de gaz)
# 2014-02-17 : messagebox de saisie du no canal pour profil
# 2014-02-22 : lecture KMATRIX_K_PC si existe
# 2014-02-24 : equivalent channels = size kmat(P)
# 2014-05-19 : meme taille graphes pour 44/101 niveaux

# constants
chnx = 500
# all images will be interpolated on levx
levx = 101
# tick positions for y axes (change it if you change levx)
ticky_positions = [100, 75, 50, 25, 0]
wix = 800
wnx = 400.0
spl = 29979245800.0  # en cm
gaslist = ['CH4', 'CO', 'CO2', 'N2O', 'O3']


class figPlot(plt.Figure):

    def __init__(self):
        plt.Figure.__init__(self)


class KPCMatrixFrame(util.GenericView):

    helpTitle = "Help"
    helpMessage = """Run K"""
    helpPage = os.environ["RTTOV_GUI_PREFIX"] + "/doc/helpKpcMatrixFrame.html"

    def __init__(self, parent, title, fname, npcscores=200, runNumber=1):
        util.GenericView.__init__(self, parent,  title)

        self.CreateMenuBar()
        self.runNumber = runNumber
        self.mykProfileView = None
        self.SetSize((900, 650))
        self.SetPosition((10, 10))
        self.sbgen = self.CreateStatusBar()
        self.sbgen.SetBackgroundColour('WHITE')
        txt = 'read KPCMatrix file %s' % fname
        self.sbgen.SetStatusText(txt)
        self.npcscores = npcscores

# timer affichage StatusBar
        self.txtsb = ''
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.refreshSB, self.timer)

        """ read kpcmat.h5 file """
        self.fh5 = fname
        """ read profile in kpcmat.h5 """
        self.baseProfile = rttov.profile.Profile()
        self.baseProfile = rmodel.project.OpenAProfile(self.fh5, 1)

        t0 = time.time()
        frad = h5py.File(self.fh5, 'r')
        h5 = frad['/']
        self.kmat = rttov.kpcmatrix.Kpcmatrix()
        self.kmat.loadh5(h5)
        frad.close()
        t1 = time.time()
        txt = 'read KPCMatrix file %s : %f sec.' % (self.fh5, (t1 - t0))
        self.sbgen.SetStatusText(txt)

        # self.kmat.kpcmatrix.display()

        self.sat = self.kmat.misc['SATELLITE']
        self.ins = self.kmat.misc['INSTRUMENT']

        self.chn = self.kmat.kpcmatrix['T'].shape[1]

        self.lev = self.kmat.profile['NLEVELS']
#
#  growth the levels  < 101
#
# table d'equivalence des indices
        self.tabilevels = np.zeros(levx, int)
        for l in range(0, levx):
            self.tabilevels[l] = round(l * self.lev / levx)

        self.tabP = self.kmat.profile['P']
        # print 'CHANNELS/NLEVELS :',self.chn, self.lev
        if self.lev < levx:
            tabNP = np.zeros(levx)
            for l in range(0, levx):
                tabNP[l] = self.tabP[self.tabilevels[l]]
            self.tabP = tabNP

        self.tabW = self.kmat.misc['WAVENUMBERS']
        self.tabY = {}

        """
        gasplot : 0.0 not shown ?
        """
        self.tcomm = {}
        self.tunit = {}
        self.tcbar = {}
        self.gasplot = []
        for gas in ['T', 'Q'] + gaslist:
            if self.kmat.kpcmatrix[gas] is None:
                continue
            kmin = np.amin(self.kmat.kpcmatrix[gas])
            kmax = np.amax(self.kmat.kpcmatrix[gas])
            self.tcbar[gas] = 1
            if kmin == 0.0 and kmax == 0.0:
                self.tcbar[gas] = 0
                continue

            attr = '%s_ATTRIBUTE' % gas
            self.tcomm[gas] = u"{}".format(
                self.kmat.kpcmatrix[attr]['COMMENT'])
            self.tunit[gas] = self.kmat.kpcmatrix[attr]['UNITS']
            # print gas, attr, self.tcomm[gas],' :',kmin,'<==>',kmax

            if not gas == 'T' and not gas == 'Q':
                self.gasplot.append(gas)
#
# init figure
#
        self.fig = figPlot()
        self.cnv = FigureCanvas(self, -1, self.fig)
        self.cnv.mpl_connect('motion_notify_event', self.onMMotion)
        self.cnv.mpl_connect('key_press_event', self.onkeypress)
#
# subplot width
#
        self.barwidth = 20

        self.yMax = self.lev

        self.xMax = self.barwidth * self.chn

        self.txtsb = '%s / %s' % (
            self.sat.replace(' ', '').upper(), self.ins.upper())
#
# matplotlib ToolBar selon nombre canaux
#

#
# choice gas in subplot 313
#
        self.glabel = {}
        self.gasID = {}

        self.tlb = MyCustomToolbar(self.cnv)
        iconsdir = os.environ["RTTOV_GUI_PREFIX"] + '/icons/'

        if len(self.gasplot) > 1:
            for ig in self.gasplot:
                self.gasID[ig] = wx.Window.NewControlId()
                self.glabel[self.gasID[ig]] = ig
                ico = iconsdir + '%s.png' % ig
                self.tlb.AddSimpleTool(
                    self.gasID[ig], _load_bitmap(ico), ig, '')
                wx.EvtHandler.Bind(self,
                                   event=wx.EVT_TOOL,
                                   handler=self.gaschoice,
                                   id=self.gasID[ig])
        self.LEFTLEFT_ID = wx.Window.NewControlId()
        self.tlb.AddSimpleTool(self.LEFTLEFT_ID, _load_bitmap(
            iconsdir + 'hand.xpm'), 'One screen left', '')
        wx.EvtHandler.Bind(self,
                           event=wx.EVT_TOOL,
                           handler=self.onleftleft,
                           id=self.LEFTLEFT_ID)

        self.LEFT_ID = wx.Window.NewControlId()
        self.tlb.AddSimpleTool(self.LEFT_ID, _load_bitmap(
            iconsdir + 'stock_left.xpm'), 'left', '')
        wx.EvtHandler.Bind(self,
                           event=wx.EVT_TOOL,
                           handler=self.onleft,
                           id=self.LEFT_ID)

        self.RIGHT_ID = wx.Window.NewControlId()
        self.tlb.AddSimpleTool(self.RIGHT_ID, _load_bitmap(
            iconsdir + 'stock_right.xpm'), 'right', '')
        wx.EvtHandler.Bind(self,
                           event=wx.EVT_TOOL,
                           handler=self.onright,
                           id=self.RIGHT_ID)

        self.RIGHTRIGHT_ID = wx.Window.NewControlId()
        self.tlb.AddSimpleTool(self.RIGHTRIGHT_ID, _load_bitmap(
            iconsdir + 'hand.xpm'), 'One screen right', '')
        wx.EvtHandler.Bind(self,
                           event=wx.EVT_TOOL,
                           handler=self.onrightright,
                           id=self.RIGHTRIGHT_ID)

        self.UP_ID = wx.Window.NewControlId()
        self.tlb.AddSimpleTool(self.UP_ID, _load_bitmap(
            iconsdir + 'stock_up.xpm'), 'scroll up', '')
        wx.EvtHandler.Bind(self,
                           event=wx.EVT_TOOL,
                           handler=self.onup,
                           id=self.UP_ID)

        self.DOWN_ID = wx.Window.NewControlId()
        self.tlb.AddSimpleTool(self.DOWN_ID, _load_bitmap(
            iconsdir + 'stock_down.xpm'), 'scroll down', '')
        wx.EvtHandler.Bind(self,
                           event=wx.EVT_TOOL,
                           handler=self.ondown,
                           id=self.DOWN_ID)
#
# add kscale icon : K / K%
#
        self.scalek = 0

        self.tlb.Realize()

# canvas and toolbar in sizer
        sizer = wx.BoxSizer(wx.VERTICAL)

        sizer.Add(self.tlb, 0, wx.GROW)
        sizer.Add(self.cnv, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(sizer)
        self.Fit()
#
        self.kprof = None
#

        self.subP = {}
        self.cax = {}
#
# subplots T and Q
#
        self.subP['T'] = self.fig.add_subplot(3, 1, 1)
        self.tabY['T'] = self.extract_kgas('T')
        img = self.subP['T'].imshow(self.tabY['T'], origin='upper')
        self.subP['T'].set_ylabel('P levels')
        divider = make_axes_locatable(self.subP['T'])
        self.cax['T'] = divider.append_axes("right", size="1%", pad=0.05)
        self.fig.colorbar(img, cax=self.cax['T'])
        self.subP['T'].set_xlim([0, min(self.xMax, wix)])
        self.draw_plot('T', 'T', self.kmat.kpcmatrix[
                       '%s_ATTRIBUTE' % 'T']['COMMENT'])

        self.subP['Q'] = self.fig.add_subplot(
            3, 1, 2, sharex=self.subP['T'], sharey=self.subP['T'])
        self.tabY['Q'] = self.extract_kgas('Q')
        img = self.subP['Q'].imshow(self.tabY['Q'], origin='upper')
        self.subP['Q'].set_ylabel('P levels')
        divider = make_axes_locatable(self.subP['Q'])
        self.cax['Q'] = divider.append_axes("right", size="1%", pad=0.05)
        self.fig.colorbar(img, cax=self.cax['Q'])
        self.subP['Q'].set_xlim([0, min(self.xMax, wix)])
        self.draw_plot('Q', 'Q', self.kmat.kpcmatrix[
                       '%s_ATTRIBUTE' % 'Q']['COMMENT'])

#
# subplot other gases if exist
#
        if self.gasplot:
            for gas in self.gasplot:
                self.tabY[gas] = self.extract_kgas(gas)

            gas = self.gasplot[0]
            self.subP['G'] = self.fig.add_subplot(
                3, 1, 3, sharex=self.subP['T'], sharey=self.subP['T'])
            img = self.subP['G'].imshow(self.tabY[gas], origin='upper')
            self.subP['G'].set_ylabel('P levels')

            divider = make_axes_locatable(self.subP['G'])
            self.cax['G'] = divider.append_axes("right", size="1%", pad=0.05)
            self.fig.colorbar(img, cax=self.cax['G'])

            self.subP['G'].set_xlim([0, min(self.xMax, wix)])
            self.draw_plot('G', gas, self.kmat.kpcmatrix[
                           '%s_ATTRIBUTE' % gas]['COMMENT'])

        self.fig.tight_layout()
        self.fig.canvas.draw()
        self.Show()

    def gaschoice(self, evt):
        gas = self.glabel[evt.GetId()]
        # print "choix gas = ", gas, self.tcomm[gas]
        # print self.subP['G']
        self.draw_plot('G', gas, self.tcomm[gas])
        return

    def onkeypress(self, evt):
        if (evt.key == 'p' or evt.key == 'P') and evt.inaxes:

            if self.chn > chnx:
                c = int(round(evt.xdata)) + 1
            else:
                c = int(round(evt.xdata) / self.barwidth) + 1

            # print 'x=', evt.xdata ,' c=',c
            if c < 0 or c > self.chn:
                return

            txt = 'display profile channel : %d from %s' % (c, self.fh5)
            self.writeSB(txt, 'green', 5, 1)
            pchan = self.kmat.getchanprof(c - 1)
            self.mykProfileView = rview.kprofileframe.KProfileView(
                self, pchan, channel=c, baseProfile=self.baseProfile,
                runNumber=self.runNumber)

            # pchan.display()

    def extract_kgas(self, gas):

        # growth levels ?
        if self.lev < levx:
            tab101 = np.zeros((levx, self.chn))
            for l in range(0, levx):
                tab101[l][:] = self.kmat.kpcmatrix[gas][self.tabilevels[l]][:]

            self.kmat.kpcmatrix[gas] = tab101

        tY = np.zeros(shape=(self.barwidth * self.chn, levx))
#
# min/max pour couleur bande separatrice
#
        kmin = np.amin(self.kmat.kpcmatrix[gas])
        kmax = np.amax(self.kmat.kpcmatrix[gas])
        kmoy = (kmax + kmin) / 2
        pc = 0
        for ic in range(0, self.chn):

            for ir in range(0, self.barwidth - 1):
                tY[pc] = np.transpose(self.kmat.kpcmatrix[gas])[ic]
                pc += 1
#
# separation bandes
#
            tY[pc] = kmoy
            pc += 1

        return np.transpose(tY)

    # pan the graph to the left
    def onleftleft(self, evt):
        if evt.GetId() == self.LEFTLEFT_ID:
            axes = self.cnv.figure.axes[0]
            x1, x2 = axes.get_xlim()
            shift = x2 - x1
            if x1 - shift <= 0:
                axes.set_xlim(0, shift)
            else:
                axes.set_xlim(x1 - shift, x2 - shift)
            self.redraw_fig()

    def onleft(self, evt):
        if evt.GetId() == self.LEFT_ID:
            axes = self.cnv.figure.axes[0]
            x1, x2 = axes.get_xlim()
            shift = (x2 - x1) / 10
            axes.set_xlim(x1 - shift, x2 - shift)
            self.redraw_fig()

    def onrightright(self, evt):
        if evt.GetId() == self.RIGHTRIGHT_ID:
            axes = self.cnv.figure.axes[0]
            x1, x2 = axes.get_xlim()
            shift = x2 - x1
            if x2 + shift >= self.xMax:
                axes.set_xlim(self.xMax - shift, self.xMax)
            else:
                axes.set_xlim(x1 + shift, x2 + shift)
            self.redraw_fig()

    def onright(self, evt):
        if evt.GetId() == self.RIGHT_ID:
            axes = self.cnv.figure.axes[0]
            x1, x2 = axes.get_xlim()
            shift = (x2 - x1) / 10
            axes.set_xlim(x1 + shift, x2 + shift)
            self.redraw_fig()

    def onup(self, evt):
        if evt.GetId() == self.UP_ID:
            axes = self.cnv.figure.axes[0]
            y1, y2 = axes.get_ylim()
            axes.set_ylim(y1 - 1, y2 - 1)
            self.redraw_fig()

    def ondown(self, evt):
        if evt.GetId() == self.DOWN_ID:
            axes = self.cnv.figure.axes[0]
            y1, y2 = axes.get_ylim()
            axes.set_ylim(y1 + 1, y2 + 1)
            self.redraw_fig()

    # scale K-K%
    def onscale5k(self, evt):
        if evt.GetId() == self.SCALE5K_ID:
            self.scalexk(0.05)

    def onscale10k(self, evt):
        if evt.GetId() == self.SCALE10K_ID:
            self.scalexk(0.1)

    def onscale20k(self, evt):
        if evt.GetId() == self.SCALE20K_ID:
            self.scalexk(0.2)

    def scalexk(self, coeff):
        if self.scalek == 1:
            return
            # self.onresetk(None)
        self.scalek = 1
        text = 'init kmatrix scale'
        self.writeSB(text, 'yellow', 5, 1)
        t1 = time.clock()
        # self.kmat.kpcmatrix.kscale(coeff)
        t2 = time.clock()
        text1 = 'end k scale : %s seconds' % (t2 - t1)
        self.writeSB(text, 'green', 5, 1)

#
# maj pour tous gaz presents
#
        for ig in ['Q'] + self.gasplot:
            self.tabY[ig] = self.extract_kgas(ig)
            self.tcomm[ig] = '%s scaled by %3.2f input profile' % (self.tcomm[
                                                                   ig], coeff)

#
# redraw subplots
#
        for a in self.cnv.figure.axes:
            gas = a.get_title().split(' :')[0]
            if gas == 'T':
                continue
            elif gas:
                if gas == 'Q':
                    plot = 'Q'
                else:
                    plot = 'G'

                img = self.subP[plot].imshow(self.tabY[gas], origin='upper')
                self.fig.colorbar(img, cax=self.cax[plot])
                self.draw_plot(plot, gas, self.tcomm[gas])

        self.redraw_fig()
        t3 = time.clock()
        text2 = ' / update graph % sec.' % (t3 - t2)
        self.writeSB(text1 + text2, 'green', 5, 1)

    def onresetk(self, evt):
        if self.scalek == 0:
            return
        self.scalek = 0
        text = 'init kmatix reset'
        self.writeSB(text, 'yellow', 5, 1)
        t1 = time.clock()
        t1 = time.clock()
        frad = h5py.File(self.fh5, 'r')
        h5 = frad['/']
        self.kmat.loadh5(h5)
        frad.close()
        t2 = time.clock()
        text1 = 'end k reset : %s seconds' % (t2 - t1)


#
# maj pour tous gaz presents
#
        for ig in ['Q'] + self.gasplot:
            self.tabY[ig] = self.extract_kgas(ig)
            attr = '%s_ATTRIBUTE' % ig
            self.tcomm[ig] = self.kmat.kpcmatrix[attr]['COMMENT']

        for a in self.cnv.figure.axes:
            gas = a.get_title().split(' :')[0]
            if gas == 'T':
                continue
            elif gas:
                if gas == 'Q':
                    self.draw_plot('Q', 'Q', self.tcomm[gas])
                else:
                    self.draw_plot('G', gas, self.tcomm[gas])

        self.redraw_fig()
        t3 = time.clock()
        text2 = ' / update graph % sec.' % (t3 - t2)
        self.writeSB(text1 + text2, 'green', 5, 1)

    def onselchan(self, evt):
        if evt.GetId() != self.SELCHAN_ID:
            return

        txtbox = 'Select a channel number [%01d..%d] : ' % (1, self.chn)
        chnbox = util.SimpleFSDialogBox(
            self, -1, "Select channel", text=txtbox, text2="Channel", minval=1,
            maxval=self.chn, increment=1)

        if (chnbox.ShowModal() == wx.ID_OK):
            c = chnbox.GetSelections()
            txt = 'display profile channel : %d from %s' % (c, self.fh5)
            self.writeSB(txt, 'green', 5, 1)
            pchan = self.kmat.getchanprof(c - 1)
            self.mykProfileView = rview.kprofileframe.KProfileView(
                self, pchan, channel=c, baseProfile=self.baseProfile,
                runNumber=self.runNumber)
        chnbox.DeletePendingEvents()
        wx.CallAfter(chnbox.Destroy)

    def draw_plot(self, plot, gas, comment):
        tt = u'{} : {}'.format(gas, comment.decode("utf-8"))
        self.subP[plot].set_title(tt)
        img = self.subP[plot].imshow(self.tabY[gas], origin='upper')
        if self.tcbar[gas]:
            self.fig.colorbar(img, cax=self.cax[plot])

# x ticks = channel number
        if self.chn <= chnx:
            xtl = []
            xtp = []
            if gas == 'T':
                for i in range(0, self.chn + 1):
                    xtl.append(i + 1)
                    xtp.append((xtl[i] - 1) * self.barwidth +
                               int(self.barwidth) / 2)
                    # print i,' ==> lab :', xtl[i],'pos :',xtp[i]

                self.subP['T'].set_xticks(xtp)
                self.subP['T'].set_xticklabels(xtl)

        else:
            xa = self.subP[plot].get_xaxis()
            xa.set_major_locator(MaxNLocator(nbins=10, integer=True))

# astuce pour cacher le 0
        my_ticks = self.subP[plot].get_xticks()
        xlabel = np.arange(len(my_ticks) + 1)

        for i in range(0, len(my_ticks)):
            xlabel[i] = i + 1
        self.subP[plot].set_xticklabels(xlabel, rotation=-45)

        self.redraw_fig()

    def redraw_fig(self):
        #
        # for each axes compute and set yticklables
        #
        taxes = self.fig.get_axes()
        ytickslabels = []
        for pos in ticky_positions:
            ytickslabels.append(str(self.tabilevels[pos] + 1))
        for ax in range(0, len(taxes)):
            if taxes[ax].get_title():
                taxes[ax].set_yticks(ticky_positions)
                taxes[ax].set_yticklabels(ytickslabels)
        self.fig.canvas.draw_idle()
        return

    def MenuData(self):
        return (
            ("&File",  # File Menu
             ('&Quit', 'Quit', self.OnQuit, "quit", True)),
            ("&Help",  # Help Menu
             ("About", "About screen", self.OnAbout, "about", True),
             ("&Help", "Help", self.OnHelpHTML, "help", True)))

    def onMMotion(self, evt):
        if evt.inaxes:
            gas = evt.inaxes.title.get_text().split(' :')[0]

            if self.chn > chnx:
                c = int(round(evt.xdata)) + 1
            else:
                c = int(round(evt.xdata) / self.barwidth) + 1

            # print 'x=', evt.xdata ,' c=',c
            if c < 0 or c > self.xMax:
                return

            y = int(evt.ydata)
# separation canaux
            if (
                    c > 1 and self.barwidth > 1 and
                    round(evt.xdata) % self.barwidth == 0.0):
                self.writeSB('channel separator', 'grey', 5, 0)
                return
#
# pression, canal, k
#
            try:
                k = self.tabY[gas][y][int(round(evt.xdata))]
                text = "%s : C=#%d P=%.3f hPa k=%g %s" % (
                    gas, c, self.tabP[y], k, self.tunit[gas])
#
# ecriture message dans status bar
#
                self.writeSB(text, 'grey', 15, 0)
            except Exception:
                return


if __name__ == "__main__":
    print("******************************************************************")
    print("***** KPCMatrixFrame environment *********************************")
    print("******************************************************************")
    print('Executing on', os.uname())
    print('Python version', sys.version)
    print('wxPython version', wx.__version__)
    print('matplotlib version', matplotlib.__version__)
    print("******************************************************************")
    fh5 = []

    for i in range(len(sys.argv)):
        fh5.append(sys.argv[i])

    print("f=", fh5)
    if len(fh5) < 2:
        from util import rttov_gui_data_test_dir
        data_test_dir = rttov_gui_data_test_dir()
        kpcmatrixfile = os.path.join(data_test_dir, "pckmat.h5")
    else:
        kpcmatrixfile = fh5[1]
    print("kpcmatrixfile:", kpcmatrixfile)
    app = wx.App(0)
    frame = KPCMatrixFrame(None, "K Matrix Viewer", kpcmatrixfile)
    frame.Show()

    app.MainLoop()
