# -*- coding: utf-8 -*-

try:
    import wx
except ImportError:
    import sys
    sys.stderr.write('ERROR: wxPython is not installed\n')
    sys.exit(1)

import os
from rview import helpframe
import wx.lib.agw.floatspin as FS
from pubsub import pub
try:
    import datetime
except ImportError:
    import sys
    sys.stderr.write('ERROR: d is not installed\n')
    sys.exit(1)

import logging
from matplotlib.backends.backend_wxagg import \
    NavigationToolbar2WxAgg as ToolBar

kindOfItem = {'Q': "GASES", 'O3': "GASES", 'CO2': "GASES", 'CH4': "GASES",
              'CO': "GASES", 'N2O': "GASES", "T": "GASES", 'SO2': "GASES",
              'INSO': "AEROSOLS", 'WASO': "AEROSOLS", 'SOOT': "AEROSOLS",
              'SSAM': "AEROSOLS", 'SSCM': "AEROSOLS", 'MINM': "AEROSOLS",
              'MIAM': "AEROSOLS", 'MICM': "AEROSOLS", 'MITR': "AEROSOLS",
              'SUSO': "AEROSOLS", 'VOLA': "AEROSOLS", 'VAPO': "AEROSOLS",
              'ASDU': "AEROSOLS", 'BCAR': "AEROSOLS", 'DUS1': "AEROSOLS",
              'DUS2': "AEROSOLS", 'DUS3': "AEROSOLS", 'SULP': "AEROSOLS",
              'SSA1': "AEROSOLS", 'SSA2': "AEROSOLS", 'SSA3': "AEROSOLS",
              'OMAT': "AEROSOLS",
              'STCO': "CLOUDS", 'STMA': "CLOUDS", 'CUCC': "CLOUDS",
              'CUCP': "CLOUDS", 'CUMA': "CLOUDS", 'CIRR': "CLOUDS",
              'CFRAC': "CLOUDS", "CLW": "CLOUDS", 'ICEDE': "CLOUDS",
              'CLWDE': "CLOUDS"
              }

axesDef = {"GASES":     {"xlimits": None},
           "AEROSOLS":  {"xlimits": None},
           "CLOUDS":    {"xlimits": None}
           }


def rttov_gui_data_test_dir():
    import sys
    import os
    print(os.environ["RTTOV_GUI_PREFIX"])
    try:
        prefix = os.environ["RTTOV_GUI_PREFIX"]
    except Exception:
        print("RTTOV_GUI_PREFIX is not defined")
        sys.exit(1)
    return os.path.join(prefix, "rttov_gui_data_test")


class MyCustomToolbar(ToolBar):
    toolitems = [t for t in ToolBar.toolitems if
                 t[0] in ("Home", "Back", "Forward", 'Subplots', "Save")]

    def __init__(self, plotCanvas):
        # create the default toolbar
        ToolBar.__init__(self, plotCanvas)

    def AddSimpleTool(self, id, bitmap, mess1, mess2):
        if wx.VERSION[0] <= 3:
            ToolBar.AddSimpleTool(self, id, bitmap, mess1, mess2)
        else:
            ToolBar.AddTool(self, id, mess2, bitmap, mess1)


class SimpleFSDialogBox(wx.Dialog):

    def __init__(
            self, parent, ID, title, text="", text2="", helptext="",
            minval=1, maxval=100, increment=None, defval=1,
            size=wx.DefaultSize, pos=wx.DefaultPosition,
            style=wx.DEFAULT_DIALOG_STYLE,
            useMetal=False,
    ):

        # Instead of calling wx.Dialog.__init__ we precreate the dialog
        # so we can set an extra style that must be set before
        # creation, and then we create the GUI object using the Create
        # method.
        if increment is None:
            increment = 10
        if wx.VERSION[0] <= 3:
            pre = wx.PreDialog()
            pre.SetExtraStyle(wx.DIALOG_EX_CONTEXTHELP)
            pre.Create(parent, ID, title, pos, size, style)

            # This next step is the most important, it turns this Python
            # object into the real wrapper of the dialog (instead of pre)
            # as far as the wxPython extension is concerned.
            self.PostCreate(pre)
        else:
            wx.Dialog.__init__(self)
            self.SetExtraStyle(wx.DIALOG_EX_CONTEXTHELP)
            self.Create(parent, ID, title)

        # This extra style can be set after the UI object has been created.
        if 'wxMac' in wx.PlatformInfo and useMetal:
            self.SetExtraStyle(wx.DIALOG_EX_METAL)

        # Now continue with the normal construction of the dialog
        # contents
        sizer = wx.BoxSizer(wx.VERTICAL)

        label = wx.StaticText(self, -1, text)
        label.SetHelpText(helptext)
        sizer.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)

        # box 1 xmin
        box = wx.BoxSizer(wx.HORIZONTAL)

        label = wx.StaticText(self, -1, text2)
        box.Add(label, 0, wx.ALIGN_CENTRE | wx.ALL, 5)
        self.channelFS = FS.FloatSpin(self, -1, min_val=minval, max_val=maxval,
                                      digits=0, value=defval,
                                      increment=increment, agwStyle=FS.FS_LEFT)
        box.Add(self.channelFS, 1, wx.ALIGN_CENTRE | wx.ALL, 5)
        sizer.Add(box, 0, wx.GROW | wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        # buttons
        btnsizer = wx.StdDialogButtonSizer()

        btn = wx.Button(self, wx.ID_OK)

        btn.SetDefault()
        btnsizer.AddButton(btn)
        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)
        self.Show()

    def GetSelections(self):
        return (self.channelFS.GetValue())


class GenericView(wx.Frame):
    '''
    Generic Frame for RTTOV GUI with basic methods like OnAbout function
    '''
    aboutMessage = "Welcome to RTTOV GUI for RTTOV 13\n(C)opyright NWPSAF"
    aboutTitle = "About RTTOV GUI"
    helpMessage = ""
    helpTitle = "Help"

    def __init__(self, parent, title, style=wx.DEFAULT_FRAME_STYLE):
        '''
        Constructor
        '''

        wx.Frame.__init__(self, parent, -1, title, style=style)
        self.items = {}

        # timer affichage StatusBar
        self.txtsb = ''
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self.refreshSB, self.timer)
        self.log = None
        self.title = title
        self.Bind(wx.EVT_CLOSE, self.OnClose)

    def MakeModal(self, modal=True):
        if wx.VERSION[0] <= 3:
            wx.Frame.MakeModal(self, modal)

    def write(self, msg):
        win = wx.GetTopLevelWindows()[0]
        # if the topLevelParent is the console then log is defined (log is the
        # wx.TextCtrl of the mian window)
        logging.info("<" + self.GetTitle() + " : >" + msg)
        if win.log is not None:
            win.log.SetInsertionPointEnd()
        print((datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
               self.title, msg))

    def BeginBusy(self):
        """ set a busy cursor """
        wx.Yield()
        wx.BeginBusyCursor()

    def EndBusy(self):
        """ end the busy cursor """
        wx.EndBusyCursor()

    def refreshSB(self, evt):
        sbgen = self.GetStatusBar()
        if sbgen is not None:
            sbgen.SetBackgroundColour('WHITE')
            sbgen.SetStatusText(self.txtsb)

    def writeSB(self, msg, col, nbs, pr):
        msg = self.txtsb + ' : ' + msg
        mainW = self.GetTopLevelParent()
        if mainW.StatusBar:
            mainW.StatusBar.SetBackgroundColour(col)
            mainW.SetStatusText(msg)
        if pr > 0:
            self.write(msg)
        if nbs > 0:
            mainW.Refresh()
            self.timer.Stop()
            self.timer.Start(nbs * 1000)
        else:
            self.timer.Stop()

    def ChooseNumber(self, maxNumber, question, title):
        """ Open a dialog box for choosing a number between 1 and maxNumber """
        number = 1
        choice = []
        for i in range(1, maxNumber + 1):
            ch = "%i" % i
            choice.append(ch)
        dlg = wx.SingleChoiceDialog(None, question, title, choice)
        if dlg.ShowModal() == wx.ID_OK:
            number = int(dlg.GetStringSelection())
        dlg.DeletePendingEvents()
        wx.CallAfter(dlg.Destroy)
        return number

    def CreateMenuBar(self):
        """
        Create a Menu Bar use MenuData function which return the list of
        the menu this function must be declared for each specific view
        """
        menuBar = wx.MenuBar()
        for eachMenuData in self.MenuData():
            menuLabel = eachMenuData[0]
            menuItems = eachMenuData[1:]
            menuBar.Append(self.CreateMenu(menuItems), menuLabel)
        self.SetMenuBar(menuBar)

    def CreateMenu(self, menuItems):
        """ create a menu from menuData
            menuData is e tuple which must contains :
            - the label of the menu
            - the status of the menu
            - the handler of the menu
        """
        menu = wx.Menu()
        for label, status, handler, item, enable in menuItems:
            if not label:
                menu.AppendSeparator()
                continue
            self.items[item] = menu.Append(-1, label, status)
            self.Bind(wx.EVT_MENU, handler, self.items[item])
            self.items[item].Enable(enable)
        return menu

    def EnableMenuItem(self, aMenuItemName):
        """ enable a menu item """
        self.items[aMenuItemName].Enable(True)

    def DisableMenuItem(self, aMenuItemName):
        """ disable a menu item """
        self.items[aMenuItemName].Enable(False)

    # handler must/may be redefined in specific View
    def OnSave(self, e): print("OnSave GenericView")

    def OnRedo(self, e): pass

    def OnUndo(self, e): pass

    def OnDestroy(self, e):
        """ Close the window and send a message for the controller """
        self.Unbind(wx.EVT_TIMER)
        self.Unbind(wx.EVT_TOOL)
        self.Unbind(wx.EVT_MENU)
        self.Unbind(wx.EVT_CLOSE)
        self.Unbind(wx.EVT_COMBOBOX)
        self.Unbind(wx.EVT_NOTEBOOK_PAGE_CHANGED)
        self.Unbind(wx.EVT_CHECKBOX)
        self.Unbind(wx.EVT_SCROLL_THUMBRELEASE)
        self.Unbind(wx.EVT_BUTTON)
        self.Unbind(wx.EVT_KEY_DOWN)
        self.Unbind(wx.EVT_SCROLL_CHANGED)
        pub.sendMessage("CLOSE", msg=self.GetTitle())
        self.DeletePendingEvents()
        # self.Destroy may cause a segmentation fault even if Called after
        wx.CallAfter(self.Destroy)

    def OnClose(self, e):
        """ just hide the window """
        self.DeletePendingEvents()
        self.Show(False)

    def OnQuit(self, event):
        print("Close window", self.Title)
        # send event EVT_CLOSE
        self.Close()

    def OnAbout(self, event):
        dialog = wx.MessageDialog(self, self.aboutMessage,
                                  self.aboutTitle, wx.OK | wx.ICON_INFORMATION)
        dialog.CentreOnParent()
        dialog.ShowModal()
        dialog.DeletePendingEvents()
        wx.CallAfter(dialog.Destroy)

    def OnHelp(self, e):
        """ show a Help Message dialog box
            self.helpMessage must be defined in specific View
        """
        helpDialog = wx.MessageDialog(
            None, self.helpMessage, caption=self.helpTitle, style=wx.OK)
        helpDialog.ShowModal()
        helpDialog.DeletePendingEvents()
        wx.CallAfter(helpDialog.Destroy)

    def OnHelpHTML(self, e):
        """ show a Help Message in a HTML window
            self.helpPage and self.helpTitle must be defined in specific View
         """
        self.HelpWindow = helpframe.HelpWindow(
            self, -1, self.helpTitle, self.helpPage)

    def ShowDialog(self, parent, message, title, flags):
        dialog = wx.MessageDialog(parent, message, title, flags)
        dialog.CentreOnParent()
        dialog.ShowModal()
        dialog.DeletePendingEvents()
        wx.CallAfter(dialog.Destroy)

    def WarmError(self, message):
        """ Open a wx.MessageDialog with an error message """
        self.ShowDialog(self, message, "Error", wx.ICON_ERROR)

    def WarmQuestion(self, message):
        """ Open a wx.MessageDialog with an error message """
        self.ShowDialog(
            self, message, "Warning", wx.YES_NO | wx.ICON_QUESTION)

    def ShowInfoMessage(self, message):
        """ Open a wx.MessageDialog with an information message """
        self.ShowDialog(self, message, "Info", wx.OK)

    def OnOpenFile(self, event, directory="", title="Choose a file"):
        """ Open a wx.FileDialog window  and return the filename chosen
            else return None  """
        if directory == "":
            directory = os.getcwd()
        if wx.VERSION[0] >= 4:
            myflag = wx.FD_OPEN
        else:
            myflag = wx.OPEN
        dlg = wx.FileDialog(self, title, directory, "", "*.*", myflag)
        self.profileFile = None
        mypath = None
        if dlg.ShowModal() == wx.ID_OK:
            mypath = dlg.GetPath()
        dlg.DeletePendingEvents()
        wx.CallAfter(dlg.Destroy)
        return mypath

    def OnSaveFileAs(self, event, title="Save a file"):
        """ Open a Save file wx.FileDialog window return the chosen file
            name else return None"""
        dlg = wx.FileDialog(self, title, os.getcwd(), "", "*.*", wx.FD_SAVE)
        path = None
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
        dlg.DeletePendingEvents()
        wx.CallAfter(dlg.Destroy)
        return path

    def OnSaveFile(self, event, title="Save a file"):
        """Open a Save file wx.FileDialog window return the chosen file name
           else return None"""
        dlg = wx.FileDialog(self, title, os.getcwd(), "", "*.*", wx.FD_SAVE)
        path = None
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
        dlg.DeletePendingEvents()
        wx.CallAfter(dlg.Destroy)
        return path


class GenericViewRadio(GenericView):
    '''
    Generic Frame for RTTOV GUI with basic methods like OnAbout function
    This frame is inherited for GenericView
    It allows menu radio item
    '''

    def __init__(self, parent, title):
        '''
        Constructor
        '''
        super().__init__(parent, title)
        self.menu = None

    def CreateMenu(self, menuItems):
        """ create a menu from menuData
        menuData is e tuple which must contains :
        - the label of the menu
        - the status of the menu
        - the handler of the menu
        """
        menu = wx.Menu()
        self.menu = menu

        for label, status, handler, item, enable, menuRadio in menuItems:
            if not label:
                menu.AppendSeparator()
                continue
            if menuRadio:
                self.items[item] = menu.AppendRadioItem(-1, label)
            else:
                self.items[item] = menu.Append(-1, label, status)
            self.Bind(wx.EVT_MENU, handler, self.items[item])
            self.items[item].Enable(enable)
        return menu
