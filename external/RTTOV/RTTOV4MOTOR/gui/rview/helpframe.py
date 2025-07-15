'''
Created on Apr 25, 2013

@author: pascale
'''
# -*- coding: utf-8 -*-

# helpwindow.py

import wx
import wx.html as html
import os


class HelpWindow(wx.Frame):

    def __init__(self, parent, myid, title, page,
                 icondir=os.environ["RTTOV_GUI_PREFIX"] + "/icons"):
        """ create an help window with page page
            icondir is the directory where to find the icons """
        wx.Frame.__init__(self, parent, myid, title, size=(570, 400))
        self.page = page
        toolbar = self.CreateToolBar()
        if wx.VERSION[0] <= 3:
            toolbar.AddLabelTool(1, 'Exit', wx.Bitmap(icondir + '/exit.png'))
            toolbar.AddLabelTool(2, 'Help', wx.Bitmap(icondir + '/help.png'))
        else:
            toolbar.AddTool(1, 'Exit', wx.Bitmap(icondir + '/exit.png'))
            toolbar.AddTool(2, 'Help', wx.Bitmap(icondir + '/help.png'))

        toolbar.Realize()

        self.splitter = wx.SplitterWindow(self, -1)
        self.panelLeft = wx.Panel(self.splitter, -1, style=wx.BORDER_SUNKEN)

        self.panelRight = wx.Panel(self.splitter, -1)
        vbox2 = wx.BoxSizer(wx.VERTICAL)
        header = wx.Panel(self.panelRight, -1, size=(-1, 20))
        header.SetBackgroundColour('#6f6a59')
        header.SetForegroundColour('WHITE')
        hbox = wx.BoxSizer(wx.HORIZONTAL)

        st = wx.StaticText(header, -1, 'Help', (5, 5))
        font = st.GetFont()
        font.SetPointSize(9)
        st.SetFont(font)
        hbox.Add(st, 1, wx.TOP | wx.BOTTOM | wx.LEFT, 5)

        close = wx.BitmapButton(header, -1,
                                wx.Bitmap(icondir + '/fileclose.png',
                                          wx.BITMAP_TYPE_PNG),
                                style=wx.NO_BORDER)
        close.SetBackgroundColour('#6f6a59')
        hbox.Add(close, 0)
        header.SetSizer(hbox)

        vbox2.Add(header, 0, wx.EXPAND)

        myHelp = html.HtmlWindow(self.panelRight, -1, style=wx.NO_BORDER)

        myHelp.LoadPage(page)

        vbox2.Add(myHelp, 1, wx.EXPAND)

        self.panelRight.SetSizer(vbox2)

        self.panelLeft.SetFocus()

        self.splitter.SplitVertically(self.panelLeft, self.panelRight)
        self.splitter.Unsplit()

        self.Bind(wx.EVT_BUTTON, self.CloseHelp, id=close.GetId())
        self.Bind(wx.EVT_TOOL, self.OnClose, id=1)
        self.Bind(wx.EVT_TOOL, self.OnHelp, id=2)

        self.Bind(wx.EVT_KEY_DOWN, self.OnKeyPressed)
        self.splitter.SplitVertically(self.panelLeft, self.panelRight)
        self.panelLeft.SetFocus()
        self.CreateStatusBar()

        self.Centre()
        self.Show(True)

    def OnClose(self, event):
        self.Close()

    def OnHelp(self, event):
        self.splitter.SplitVertically(self.panelLeft, self.panelRight)
        self.panelLeft.SetFocus()

    def CloseHelp(self, event):
        self.splitter.Unsplit()
        self.panelLeft.SetFocus()

    def OnKeyPressed(self, event):
        keycode = event.GetKeyCode()
        if keycode == wx.WXK_F1:
            self.splitter.SplitVertically(self.panelLeft, self.panelRight)
            self.panelLeft.SetFocus()


if __name__ == "__main__":
    app = wx.App()
    import os
    HelpWindow(None, -1, 'HelpWindow',
               os.environ["RTTOV_GUI_PREFIX"] + "/doc/helpOptions.html")
    app.MainLoop()
