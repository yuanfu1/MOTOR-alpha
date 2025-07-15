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
from rview import util
from rview.util import MyCustomToolbar
import wx
import rmodel
from rview import myunits
from pubsub import pub

# constants
chnx = 50
# all images will be interpolated on levx
levx = 101
# tick positions for y axes (change it if you change levx)
ticky_positions = [100, 75, 50, 25, 0]
wix = 783
wnx = 400.0
spl = 29979245800.0  # en cm
gaslist = ['CH4', 'CO', 'CO2', 'N2O', 'O3', 'SO2']

obsoleteFrameMessage = (
    "kP invalid function:" +
    "\n- use the last k Matrix Viewer to display the K profile" +
    "\n- or run RTTOV K if necessary to display a new k Matrix Viewer")
kProfViewTitle = "run K profile channel "


class figPlot(plt.Figure):

    def __init__(self):
        plt.Figure.__init__(self)


class KMatrixFrame(util.GenericView):

    helpTitle = "Help"
    helpMessage = """Run K"""
    helpPage = os.environ["RTTOV_GUI_PREFIX"] + "/doc/helpKMatrixFrame.html"

    def __init__(self, parent, title, fname, runNumber=1, project=None):
        super().__init__(parent, title)
        self.title = title
        plt.rcParams['figure.autolayout'] = True
        self.myParent = parent
        self.myProject = project
        self.CreateMenuBar()
        self.runNumber = runNumber
        self.mykProfileView = None
        self.SetSize((900, 650))
        self.SetPosition((10, 10))
        self.sbgen = self.CreateStatusBar()
        self.sbgen.SetBackgroundColour('WHITE')
        txt = 'KMatrix file = %s' % fname
        self.sbgen.SetStatusText(txt)

        """  timer affichage StatusBar """
        self.txtsb = ''
        self.timer = wx.Timer(self)

        """ read kmat.h5 file """
        self.fh5 = fname
        self.kmat = rttov.kmatrix.Kmatrix()
        self.kmat101 = rttov.kmatrix.Kmatrix()

        """ read profile in kmat.h5 """
        self.baseProfile = rmodel.project.OpenAProfile(self.fh5, 1)

        frad = h5py.File(self.fh5, 'r')
        h5 = frad['/']
        self.kmat.loadh5(h5)
        frad.close()

        txt = 'lecture fichier %s' % (self.fh5)
        self.sat = self.kmat.misc['SATELLITE']
        self.ins = self.kmat.misc['INSTRUMENT']

        self.lev = self.kmat.profile['NLEVELS']
        self.tabP = self.kmat.profile['P']
        self.txtsb = '%s / %s' % (self.sat.replace(
            ' ', '').upper(), self.ins.upper())

        """ growth the levels  < 101 and map them """

        self.tabilevels = np.zeros(levx, int)
        for l in range(0, levx):
            self.tabilevels[l] = round(l * self.lev // levx)

        if self.lev < levx:
            tabNP = np.zeros(levx)
            for l in range(0, levx):
                tabNP[l] = self.tabP[self.tabilevels[l]]
            self.tabP = tabNP
        self.tabW = self.kmat.misc['WAVENUMBERS']
        self.tabY = {}

        self.chn = self.kmat.kmatrix['T'].shape[1]
        self.barwidth = max(1, int(wix / self.chn))
        self.yMax = self.lev
        if self.chn >= chnx:
            self.xMax = self.chn
        else:
            self.xMax = self.barwidth * self.chn
        self.tcomm = {}
        self.tunit = {}
        self.tcbar = {}
        self.gasplot = []
        for gas in ['T', 'Q'] + gaslist:
            if self.kmat.kmatrix[gas] is None:
                continue
            kmin = np.amin(self.kmat.kmatrix[gas])
            kmax = np.amax(self.kmat.kmatrix[gas])
            self.tcbar[gas] = 1
            if kmin == 0.0 and kmax == 0.0:
                self.tcbar[gas] = 0

            attr = '%s_ATTRIBUTE' % gas
            self.tcomm[gas] = u"{}".format((
                self.kmat.kmatrix[attr]['COMMENT']).decode("utf-8"))
            self.tunit[gas] = myunits.makePretty(
                gas,
                self.kmat.kmatrix[attr]['UNITS'],
                self.kmat.kmatrix['GAS_UNITS'])

            if not gas == 'T' and not gas == 'Q':
                self.gasplot.append(gas)

        """ extract gas """
        self.tabY['T'] = self.extract_kgas('T')
        self.tabY['Q'] = self.extract_kgas('Q')
        if self.gasplot:
            for gas in self.gasplot:
                self.tabY[gas] = self.extract_kgas(gas)
        self.init_canvas()
        self.chosen_gas = None
        if self.gasplot:
            self.chosen_gas = self.gasplot[0]
        """ plot references """
        self.sub_orig_P = None
        self.sub_scaled_P = None
        self.create_plots(self.sub_orig_P)
        self.create_toolbar()
        """ canvas and toolbar in sizer """
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.tlb, 0, wx.GROW)
        sizer.Add(self.cnv, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(sizer)
        self.Fit()
        pub.subscribe(self.ClosekProfViewListener, "CLOSE")
        self.myViews = {}
        self.MakeBinds()
        self.Show()

    def MakeBinds(self):
        self.Bind(wx.EVT_TIMER, self.refreshSB, self.timer)

    def ClosekProfViewListener(self, msg):
        viewTitle = msg
        if viewTitle == self.title:
            pub.unsubscribe(self.ClosekProfViewListener, "CLOSE")
        if viewTitle in self.myViews:
            print("remove reference to Window " + msg)
            self.myViews[viewTitle] = None

    def init_canvas(self):
        """ init figure """
        self.fig = figPlot()
        self.cnv = FigureCanvas(self, -1, self.fig)
        self.cnv.mpl_connect('motion_notify_event', self.onMMotion)
        self.cnv.mpl_connect('key_press_event', self.onkeypress)

    def create_plots(self, plot_refs):
        """ create the sublpots """

        #
        self.kprof = None
        # plot references
        plot_refs = {}
        self.cax = {}
        self.cb = {}
        #
        # subplots T and Q
        #
        plot_refs['T'] = self.fig.add_subplot(3, 1, 1)
        tt = u'%s : %s' % ('T', self.tcomm['T'])
        plot_refs['T'].set_title(tt)
        img = plot_refs['T'].imshow(self.tabY['T'], origin='upper')
        plot_refs['T'].set_ylabel('P levels')
        divider = make_axes_locatable(plot_refs['T'])
        self.cax['T'] = divider.append_axes(
            "right", size="1%", pad=0.05)
        self.cb['T'] = self.fig.colorbar(img, cax=self.cax['T'])
        plot_refs['T'].set_xlim([0, min(self.xMax, wix)])
        self.draw_plot('T', 'T', self.tcomm['T'], plot_refs)

        plot_refs['Q'] = self.fig.add_subplot(
            3, 1, 2, sharex=plot_refs['T'], sharey=plot_refs['T'])
        tt = u'%s : %s' % ('Q', self.tcomm['Q'])
        plot_refs['Q'].set_title(tt)
        img = plot_refs['Q'].imshow(self.tabY['Q'], origin='upper')
        plot_refs['Q'].set_ylabel('P levels')
        divider = make_axes_locatable(plot_refs['Q'])
        self.cax['Q'] = divider.append_axes(
            "right", size="1%", pad=0.05)
        self.cb['Q'] = self.fig.colorbar(img, cax=self.cax['Q'])
        plot_refs['Q'].set_xlim([0, min(self.xMax, wix)])
        self.draw_plot('Q', 'Q', self.tcomm['Q'], plot_refs)

#
# subplot other gases if exist
#
        if self.gasplot:
            gas = self.chosen_gas
            plot_refs['G'] = self.fig.add_subplot(
                3, 1, 3, sharex=plot_refs['T'],
                sharey=plot_refs['T'])
            tt = u'%s : %s' % (gas, self.tcomm[gas])
            plot_refs['G'].set_title(tt)
            img = plot_refs['G'].imshow(
                self.tabY[gas], origin='upper')
            plot_refs['G'].set_ylabel('P levels')

            divider = make_axes_locatable(plot_refs['G'])
            self.cax['G'] = divider.append_axes(
                "right", size="1%", pad=0.05)
            self.cb['G'] = self.fig.colorbar(img, cax=self.cax['G'])

            plot_refs['G'].set_xlim([0, min(self.xMax, wix)])
            self.draw_plot(
                'G', gas, self.tcomm[gas], plot_refs)

        self.make_yticks()
        # self.fig.tight_layout()
        self.fig.canvas.draw()
        self.plot_refs = plot_refs

    def create_toolbar(self):
        """ toolbar : custom buttons """
        self.tlb = MyCustomToolbar(self.cnv)
        iconsdir = os.environ["RTTOV_GUI_PREFIX"] + '/icons/'
        self.glabel = {}
        self.gasID = {}

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
        self.SCALE10K_ID = wx.Window.NewControlId()
        self.tlb.AddSimpleTool(self.SCALE10K_ID, _load_bitmap(
            iconsdir + 'k10.png'),
            'scale Kx10%', '')
        wx.EvtHandler.Bind(self,
                           event=wx.EVT_TOOL,
                           handler=self.onscale10k,
                           id=self.SCALE10K_ID)

        self.RESETK_ID = wx.Window.NewControlId()
        self.tlb.AddSimpleTool(self.RESETK_ID, _load_bitmap(
            iconsdir + 'reset.png'),
            'reset K', '')
        wx.EvtHandler.Bind(self,
                           event=wx.EVT_TOOL,
                           handler=self.onresetk,
                           id=self.RESETK_ID)

#
# add box : select channel
#
        self.SELCHAN_ID = wx.Window.NewControlId()
        self.tlb.AddSimpleTool(self.SELCHAN_ID, _load_bitmap(
            iconsdir + 'kP.png'),
            'select channnel', '')
        wx.EvtHandler.Bind(self,
                           event=wx.EVT_TOOL,
                           handler=self.onselchan,
                           id=self.SELCHAN_ID)

        self.tlb.Realize()

    def zoomp(self, evt, plot_refs):
        zp = int(self.szp.val)
        print("z :", zp)
        axes = plot_refs['T']
        y1, y2 = axes.get_ylim()
        ym = (y1 - y2) / 2
        print('avant :', y1, y2, ym)
        yy1 = ym + 50 / zp
        yy2 = ym - 50 / zp
        axes.set_ylim(yy1, yy2)
        print('apres :', yy1, yy2)
        self.make_yticks()

    def gaschoice(self, evt):
        gas = self.glabel[evt.GetId()]
        self.chosen_gas = gas
        plot_refs = self.plot_refs
        plot_refs['G'].clear()
        self.cb['G'].remove()
        plot_refs['G'] = self.fig.add_subplot(
            3, 1, 3, sharex=plot_refs['T'],
            sharey=plot_refs['T'])
        tt = u'%s : %s' % (gas, self.tcomm[gas])
        print(tt)
        plot_refs['G'].set_title(tt)
        img = plot_refs['G'].imshow(
            self.tabY[gas], origin='upper')
        plot_refs['G'].set_ylabel('P levels')

        divider = make_axes_locatable(plot_refs['G'])
        self.cax['G'] = divider.append_axes(
            "right", size="1%", pad=0.05)
        self.cb['G'] = self.fig.colorbar(img, cax=self.cax['G'])

        plot_refs['G'].set_xlim([0, min(self.xMax, wix)])
        self.draw_plot('G', gas, self.tcomm[gas], plot_refs)
        self.fig.canvas.draw_idle()
        self.make_yticks()

    def canDisplayKprofile(self):
        if self.myProject is None:
            return True
        else:
            if len(self.myProject.aKmats.aKmat) > 0:
                if self.runNumber != self.myProject.aKmats.aKmat[-1].runNumber:
                    return False
                else:
                    return True
            else:
                return False

    def onkeypress(self, evt):
        if (evt.key == 'p' or evt.key == 'P') and evt.inaxes:
            if not self.canDisplayKprofile():
                self.WarmError(obsoleteFrameMessage)
                return
            if self.chn > chnx:
                c = int(round(evt.xdata)) + 1
            else:
                c = int(round(evt.xdata) / self.barwidth) + 1

            # print 'x=', evt.xdata ,' c=',c
            if c < 0 or c > self.chn:
                return
            self.showkprof(c)

    def showkprof(self, c):
        """ show the kprofile frame """
        txt = 'display profile channel : %d from %s' % (c, self.fh5)
        self.writeSB(txt, 'green', 5, 1)
        pchan = self.kmat.getchanprof(c - 1)
        if self.myParent is not None:
            parent = self.myParent
        else:
            parent = self
        mykProfViewTitle = kProfViewTitle + str(c)
        self.myViews[mykProfViewTitle] = rview.kprofileframe.KProfileView(
            parent, mykProfViewTitle, channel=c, runNumber=self.runNumber,
            kProfile=pchan,
            baseProfile=self.baseProfile,
            project=self.myProject)

    def extract_kgas(self, gas):
        # levels : 101 niveaux qqsoit le profil

        tab101 = np.zeros((levx, self.chn))
        for l in range(0, levx):
            tab101[l][:] = self.kmat.kmatrix[gas][self.tabilevels[l]][:]

        self.kmat101.kmatrix[gas] = tab101

        if self.chn >= chnx:
            return self.kmat101.kmatrix[gas]
        else:
            tY = np.zeros(shape=((self.barwidth) * self.chn, levx))

            kmin = np.amin(self.kmat101.kmatrix[gas])
            kmax = np.amax(self.kmat101.kmatrix[gas])
            kmoy = (kmax + kmin) / 2
            pc = 0
            for ic in range(0, self.chn):

                for ir in range(0, self.barwidth - 1):
                    tY[pc] = np.transpose(self.kmat101.kmatrix[gas])[ic]
                    pc += 1
            #
            # bands separation
            #
                tY[pc] = kmoy
                pc += 1

            return np.transpose(tY)

    def onleftleft(self, evt):
        """ pan the graph to the left """
        if evt.GetId() == self.LEFTLEFT_ID:
            axes = self.cnv.figure.axes[0]
            x1, x2 = axes.get_xlim()
            shift = x2 - x1
            if x1 - shift <= 0:
                axes.set_xlim(0, shift)
            else:
                axes.set_xlim(x1 - shift, x2 - shift)
            self.make_yticks()

    def onleft(self, evt):
        if evt.GetId() == self.LEFT_ID:
            axes = self.cnv.figure.axes[0]
            x1, x2 = axes.get_xlim()
            shift = (x2 - x1) / 10
            axes.set_xlim(x1 - shift, x2 - shift)
            self.make_yticks()

    def onrightright(self, evt):
        if evt.GetId() == self.RIGHTRIGHT_ID:
            axes = self.cnv.figure.axes[0]
            x1, x2 = axes.get_xlim()
            shift = x2 - x1
            if x2 + shift >= self.xMax:
                axes.set_xlim(self.xMax - shift, self.xMax)
            else:
                axes.set_xlim(x1 + shift, x2 + shift)
            self.make_yticks()

    def onright(self, evt):
        if evt.GetId() == self.RIGHT_ID:
            axes = self.cnv.figure.axes[0]
            x1, x2 = axes.get_xlim()
            shift = (x2 - x1) / 10
            axes.set_xlim(x1 + shift, x2 + shift)
            self.make_yticks()

    def onup(self, evt):
        if evt.GetId() == self.UP_ID:
            axes = self.cnv.figure.axes[0]
            y1, y2 = axes.get_ylim()
            axes.set_ylim(y1 - 1, y2 - 1)
            self.make_yticks()

    def ondown(self, evt):
        if evt.GetId() == self.DOWN_ID:
            axes = self.cnv.figure.axes[0]
            y1, y2 = axes.get_ylim()
            axes.set_ylim(y1 + 1, y2 + 1)
            self.make_yticks()

    def redraw_y(self):
        axes = self.cnv.figure.axes[0]
        y1, y2 = axes.get_ylim()
        axes.set_ylim(y1, y2)

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
        """ scale gas and plot """
        if self.scalek == 1:
            return

        self.scalek = 1
        text = 'init kmatrix scale'
        self.writeSB(text, 'yellow', 5, 1)
        self.kmat.kscale(coeff)
        text1 = u'end k scale'
        self.writeSB(text, 'green', 5, 1)
#
# update for all present gases
#
        for ig in ['Q'] + self.gasplot:
            self.tabY[ig] = self.extract_kgas(ig)
            self.tcomm[ig] = u'%s scaled by %3.2f input profile' % (
                self.tcomm[ig], coeff)

        for item in self.cax.keys():
            for axis in ['top', 'bottom', 'left', 'right']:
                self.cax[item].spines[axis].set_linewidth(0)
                self.cax[item].set_xticks([])
                self.cax[item].set_yticks([])
        for a in self.cnv.figure.axes:
            a.clear()
        for item in self.cax.keys():
            self.cax[item].clear()
        for item in self.cb.keys():
            self.cb[item].remove()
        self.create_plots(self.sub_scaled_P)
        text2 = u' / update graph '
        self.writeSB(text1 + text2, 'green', 5, 1)

    def onresetk(self, evt):
        """ reread kmat,  clear all axes and recreate plots """
        if self.scalek == 0:
            return
        self.scalek = 0
        text = 'init kmatix reset'
        self.writeSB(text, 'yellow', 5, 1)
        frad = h5py.File(self.fh5, 'r')
        h5 = frad['/']
        self.kmat.loadh5(h5)
        frad.close()
        text1 = 'end k reset'
        for ig in ['Q'] + self.gasplot:
            self.tabY[ig] = self.extract_kgas(ig)
            attr = '%s_ATTRIBUTE' % ig
            self.tcomm[ig] = str(self.kmat.kmatrix[attr]['COMMENT'])

        for a in self.cnv.figure.axes:
            a.clear()
        for item in self.cax.keys():
            self.cax[item].clear()
        for item in self.cb.keys():
            self.cb[item].remove()
        self.create_plots(self.sub_orig_P)

        text2 = ' / update graph '
        self.writeSB(text1 + text2, 'green', 5, 1)

    def onselchan(self, evt):
        if evt.GetId() != self.SELCHAN_ID:
            return
        if not self.canDisplayKprofile():
            self.WarmError(obsoleteFrameMessage)
            return
        txtbox = 'Select a channel number [%01d..%d] : ' % (1, self.chn)
        chnbox = util.SimpleFSDialogBox(
            self, -1, "Select channel", text=txtbox, text2="Channel",
            minval=1, maxval=self.chn, increment=1)
        chnbox.CenterOnParent()
        if (chnbox.ShowModal() == wx.ID_OK):
            c = chnbox.GetSelections()
            self.showkprof(c)
            chnbox.DeletePendingEvents()
            wx.CallAfter(chnbox.Destroy)

    def draw_plot(self, plot, gas, comment, plot_refs):
        """ draw the plot """
        tt = u'%s : %s' % (gas, comment)
        plot_refs[plot].set_title(tt)
        plot_refs[plot].imshow(self.tabY[gas], origin='upper')

        # x ticks = channel number
        if self.chn <= chnx:
            if self.chn <= 10:
                xtl = []
                xtp = []
                if gas == "T":
                    for i in range(0, self.chn + 1):
                        xtl.append(i + 1)
                        xtp.append((xtl[i] - 1) *
                                   self.barwidth + int(self.barwidth) // 2)

                    plot_refs[gas].set_xticks(xtp)
                    plot_refs[gas].set_xticklabels(xtl)
            else:
                xtl = []
                xtp = []
                if gas == "T":
                    for i in range(0, (self.chn + 1) // 2):
                        xtl.append(2 * i + 1)
                        xtp.append(int(self.barwidth // 2) +
                                   (xtl[i] - 1) * self.barwidth)

                    plot_refs[gas].set_xticks(xtp)
                    plot_refs[gas].set_xticklabels(
                        xtl, horizontalalignment='center')

        else:

            xa = plot_refs[plot].get_xaxis()
            xa.set_major_locator(MaxNLocator(nbins=10, integer=True))

    def make_yticks(self):
        """ for each axes compute and set yticklables """
        taxes = self.fig.get_axes()
        ytickslabels = []
        for pos in ticky_positions:
            ytickslabels.append(str(self.tabilevels[pos] + 1))
        for ax in range(0, len(taxes)):
            if taxes[ax].get_title():
                taxes[ax].set_yticks(ticky_positions)
                taxes[ax].set_yticklabels(ytickslabels)
        # self.fig.tight_layout()
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

            # print 'x=', evt.xdata ,' c=',c, 'xmax=', self.xMax
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
                tabwc = self.tabW[c - 1]
            except Exception:
                tabwc = self.tabW

            try:
                k = self.tabY[gas][y][int(round(evt.xdata))]

                if tabwc >= wnx:
                    txt0 = ' %s=%.3f %sm]' % (
                        u'\u03BB', (10000 / tabwc), u'\u03BC')
                else:
                    txt0 = ' %s=%.4f GHz]' % (
                        u'\u03BD', (spl * tabwc / 1000000000))
                mylevel = self.tabilevels[y] + 1
                text = u"%s : C=#%d P(%d)=%.3f hPa k=%g %s [wn=%.3f cm%s" % (
                    gas, c, mylevel, self.tabP[y], k,
                    self.tunit[gas], tabwc, myunits.exponentMinus1) + txt0
#
# ecriture message dans status bar
#
                self.writeSB(text, 'grey', 15, 0)
            except Exception:
                return


if __name__ == "__main__":
    print("******************************************************************")
    print("*******   KMatrixFrame environment   *****************************")
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

    app = wx.App(0)
    if len(fh5) < 2:
        from util import rttov_gui_data_test_dir
        data_test_dir = rttov_gui_data_test_dir()
        kmatrixfile = os.path.join(data_test_dir, "kmat.h5")
    else:
        kmatrixfile = fh5[1]
    print("kmatrixfile:", kmatrixfile)
    frame = KMatrixFrame(None, "K Matrix Viewer", kmatrixfile)

    app.MainLoop()
