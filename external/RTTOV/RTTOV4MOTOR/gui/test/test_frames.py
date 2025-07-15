# -*- coding: utf-8 -*-
import wx
from pubsub import pub
import sys
from rview.kmatrixframe import KMatrixFrame
from rview.radianceframe import RadianceFrame
from rview.profileframe import ProfileView
import rmodel


class MySubFrame(wx.Frame):
    def __init__(self, parent, title):
        super().__init__(parent, title=title)
        self.title = title
        self.panel = wx.Panel(self)
        self.SetSize(200, 200)
        self.Bind(wx.EVT_CLOSE, self.OnClose)
        self.Show()

    def OnClose(self, e):
        pub.sendMessage("CLOSE", msg=self.title)
        print("Destroy ", self.title)
        self.Destroy()


class MyFrame(wx.Frame):
    def __init__(self, listefile=None):
        super().__init__(None, title="MAIN")
        self.fh5 = listefile
        self.panel = wx.Panel(self)
        self.panel.SetSize(600, 600)
        button1 = wx.Button(self.panel, label='new radiance')
        main_sizer = wx.BoxSizer(wx.HORIZONTAL)
        main_sizer.Add(button1)
        button2 = wx.Button(self.panel, label='new kmatrix')
        main_sizer.Add(button2)
        button3 = wx.Button(self.panel, label='new simple frame')
        main_sizer.Add(button3)
        button4 = wx.Button(self.panel, label='new profileframe')
        main_sizer.Add(button4)
        self.panel.SetSizer(main_sizer)

        button1.Bind(wx.EVT_BUTTON, self.NewWindow1)
        button2.Bind(wx.EVT_BUTTON, self.NewWindow2)
        button3.Bind(wx.EVT_BUTTON, self.NewWindow3)
        button4.Bind(wx.EVT_BUTTON, self.NewWindow4)
        self.myViews = {}
        self.counter = 0
        self.log = None

        pub.subscribe(self.listener, "CLOSE")
        self.Show()

    def listener(self, msg):
        print("in listener")
        if msg in self.myViews:
            print("  remove ref to", msg)
            # wx.CallAfter(self.myViews[msg].Destroy) ## with that we have a segmentation fault
            self.myViews[msg] = None

    def NewWindow1(self, e):

        title = "sub window 1 number " + str(self.counter)
        #self.myViews[title] = MySubFrame(self,title)
        myfile = self.fh5[0]

        self.myViews[title] = RadianceFrame(self, title, myfile, False)
        self.myViews[title].Show()
        self.counter = self.counter + 1
        print(self.counter)

    def NewWindow2(self, e):

        title = "sub window 2 number " + str(self.counter)
        #self.myViews[title] = MySubFrame(self,title)
        myfile = self.fh5[1]
        self.myViews[title] = KMatrixFrame(self, title, myfile)
        self.myViews[title].Show()
        self.counter = self.counter + 1
        print(self.counter)

    def NewWindow3(self, e):

        title = "sub window 2 number " + str(self.counter)
        self.myViews[title] = MySubFrame(self, title)
        self.myViews[title].Show()
        self.counter = self.counter + 1
        print(self.counter)

    def NewWindow4(self, e):
        p = rmodel.project.Project()
        p.openProfile(p.config.ENV["RTTOV_GUI_PROFILE_DIR"] +
                      "/cld50lev_co2o3.H5", 1)
        title = "sub window 4 number " + str(self.counter)
        print("here")
        self.myViews[title] = ProfileView(self, title, p.myProfile)
        self.myViews[title].Show()
        self.counter = self.counter + 1
        print(self.counter)


if __name__ == "__main__":
    print(sys.argv)
    fh5 = []
    for i in range(1, len(sys.argv)):
        fh5.append(sys.argv[i])

    print("f=", fh5)
    ex = wx.App(redirect=False)
    myView = MyFrame(listefile=fh5)
    ex.MainLoop()
