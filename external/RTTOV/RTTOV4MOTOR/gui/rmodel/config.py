# -*- coding: utf-8 -*-

import os
import sys


class Config (object):
    '''
    Configuration class for RTTOV GUI
    (ie profile directory and so on which must have been defined as environment
     variables)
    '''

    def __init__(self):
        '''
        Constructor
        '''
        self.ENV = {}

        try:
            self.ENV['RTTOV_GUI_PREFIX'] = os.environ["RTTOV_GUI_PREFIX"]
        except KeyError:
            self.ENV['RTTOV_GUI_PREFIX'] = "."

        try:
            self.ENV['GUI_WRK_DIR'] = os.environ["RTTOV_GUI_WRK_DIR"]
        except KeyError:
            self.ENV['GUI_WRK_DIR'] = os.environ["HOME"] + "/.rttov/"

        try:
            self.ENV['RTTOV_GUI_COEFF_DIR'] = os.environ["RTTOV_GUI_COEFF_DIR"]
        except KeyError:
            self.ENV['RTTOV_GUI_COEFF_DIR'] = self.ENV['RTTOV_GUI_PREFIX']

        try:
            self.ENV['RTTOV_GUI_PROFILE_DIR'] = os.environ[
                "RTTOV_GUI_PROFILE_DIR"]
        except KeyError:
            self.ENV['RTTOV_GUI_PROFILE_DIR'] = self.ENV['RTTOV_GUI_PREFIX']

        try:
            self.ENV['RTTOV_GUI_EMISS_DIR'] = os.environ["RTTOV_GUI_EMISS_DIR"]
        except KeyError:
            print("ERROR : RTTOV_GUI_EMIS_DIR environment variable is not set")
            sys.exit(1)
        self.workDir = self.ENV['GUI_WRK_DIR']
        self.profileDefaultFileName = self.workDir + "/profile.h5"
        self.radianceDefaultFileName = self.workDir + "/radr.h5"
        self.KMatrixDefaultFileName = self.workDir + "/kmat.h5"
        self.transDefaultFileName = self.workDir + "/trns.h5"
        self.surfaceDefaultFileName = self.workDir + "/surface.h5"
        self.tmpFileOut = self.workDir + "/tmpFileOut.log"
        self.tmpFileErr = self.workDir + "/tmpFileErr.log"
        self.pcDefaultFileName = self.workDir + "/pc.h5"
        self.pcKMatrixFileName = self.workDir + "/pckmat.h5"
        self.logFile = self.workDir + "/rttovgui.log"
        if not os.access(self.workDir, os.W_OK):
            sys.stdout.write("cannot write in " + self.workDir)
            sys.exit(1)


if __name__ == "__main__":
    config = Config()
    print('RTTOV_GUI_PREFIX = ', config.ENV['RTTOV_GUI_PREFIX'])
    print(config.workDir)
