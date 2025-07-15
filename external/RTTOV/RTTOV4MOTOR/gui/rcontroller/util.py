# -*- coding: utf-8 -*-

try:
    import wx

except ImportError:
    import sys
    sys.stderr.write('ERROR: check your python installation\n')
    sys.stderr.write('ERROR: wxPython is not installed\n')
    sys.exit(1)
try:
    import datetime
except ImportError:
    import sys
    sys.stderr.write('ERROR: datetime is not installed\n')
    sys.exit(1)


import logging
import os

from pubsub import pub


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file %s does not exist!" % arg)
    else:
        return open(arg, 'r')  # return an open file handle


class UndoRedo:
    """ object keep to list for options and profile
        an undo list every modification of an object profile or options
        is saved in the list
        when the list are too long they must be cleaned """

    forUndo = {"profiles": [], "options": []}
    forRedo = {"profiles": [], "options": []}

    def __init__(self):
        pass

    def clean(self, typeItem):
        pass  # TODO

    def saveforUndo(self, item, typeItem):
        self.forUndo[typeItem].append(item)
        self.clean()

    def undo(self, typeItem):
        item = self.forUndo[typeItem].pop()
        del self.forUndo[typeItem][-1]
        self.forRedo[typeItem].append(item)
        return item

    def redo(self, typeItem):
        item = self.forRedo[typeItem].pop()
        del self.forRedo[typeItem][-1]
        self.forUndo[typeItem].append(item)
        return item


class GenericController (object):
    """ define Generic actions for all controllers
        SaveProfile etc
        project is a static variable : known for all instances
        undoredo is also a static variable """
    project = None
    undoRedo = None

    def __init__(self):
        if (self.undoRedo is None):
            self.undoRedo = UndoRedo()
        pub.subscribe(self.CloseWindowListener, "CLOSE")
        self.controlerName = "generic"
        self.myViews = {}

    def prepLog(self):
        """ before any call to project methods put the cursor
            of the TextEdit of the main window at the end
            in order to prevent messed output """

        mainW = wx.GetTopLevelWindows()[0]
        if mainW is not None:
            if mainW.log is not None:
                mainW.log.SetInsertionPointEnd()

    def ShowWindowIfExists(self, view, needRefresh=False):
        if view is not None:
            logging.debug("show Window " + view.GetTitle())
            if view.IsIconized():
                view.Restore()
            else:
                if not view.IsShown():
                    view.Show()
            if needRefresh:
                view.Refresh()
            return True
        else:
            return False

    def CloseWindowListener(self, msg):
        self.prepLog()
        viewTitle = msg
        if viewTitle in self.myViews:
            logging.debug("remove reference to Window " + msg)
            self.myViews[viewTitle] = None

    def OnChanged(self, msg):
        print("message", msg)

    def write(self, msg):
        self.prepLog()
        logging.info(
            self.controlerName + " " + msg + " in " + self.controlerName)
        print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), msg)

    def SaveProfile(self, e, parentFrame):
        """ Save the Profile (Profile+Option)  """
        try:
            fileName = None
            if (self.project.savedProfileFileName) is not None:
                fileName = self.project.savedProfileFileName
            else:
                fileName = parentFrame.OnSaveFile(e)
                self.project.savedProfileFileName = fileName
            if(fileName):
                self.project.saveProfile(fileName)
                print(fileName, "saved")
        except IOError:
            print("Error while saving the profile file ", fileName)

    def SaveSurface(self, e, parentFrame):
        """ Save the Surface File (Emissivity+Reflectance)  """
        print("SaveSurface GenericController")
        try:
            fileName = None
            if (self.project.savedSurfaceFileName):
                fileName = self.project.savedSurfaceFileName
            else:
                fileName = self.project.surfaceFileName
            if(fileName):
                self.project.saveSurface(fileName)
                print(fileName, "saved")
        except IOError:
            print("Error while saving the surface file ")
