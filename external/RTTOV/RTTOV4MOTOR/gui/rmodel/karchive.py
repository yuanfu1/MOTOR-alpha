# -*- coding: utf-8 -*-
import rttov
import h5py
import numpy as np
import logging
'''
Created on Feb 16, 2016

@author: pascale
'''


class Kmatrixn(rttov.kmatrix.Kmatrix):
    '''
    classdocs :  K-matrix + runKnumber
    '''

    def __init__(self, runNumber=0):
        rttov.kmatrix.Kmatrix.__init__(self)
        self.runNumber = runNumber


class Karchive(object):
    '''
    classdocs : Keep an archive of K-matrixes
    '''

    def __init__(self, depth=4, tag="default"):
        '''
        Constructor : set de depth of the archive
        '''
        # create logger
        self.myLogger = logging.getLogger('rttov gui')
        logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        self.archiveDepth = depth
        self.aKmat = []
        self.runNumber = 0
        self.tag = tag

    def update(self, fileK, numberRun=None):
        """ add a K-matrix to the archive """
        """ read kmat.h5 file """
        """ in some context (GUI for example we want to force the """
        """ runKnumber in the archive """
        """ runKnumber could be a number or an ID or a hash """
        if self.archiveDepth == 0:
            return
        if numberRun is not None:
            self.runNumber = numberRun
        else:
            self.runNumber = self.runNumber
        fh5 = fileK
        kmat = Kmatrixn(self.runNumber)
        fk = h5py.File(fh5, 'r')
        h5 = fk['/']
        kmat.loadh5(h5)
        """ read profile in kmat.h5 """

        baseProfile = kmat.profile
        if (len(self.aKmat) > 0):
            switchrad_was = self.aKmat[-1].option['SWITCHRAD'].value
            switchrad = kmat.option['SWITCHRAD'].value
            if (switchrad != switchrad_was):
                self.reset()

        if len(self.aKmat) != 0:
            pvKmat = self.aKmat[-1]
            pvProfile = pvKmat.profile
            if baseProfile['P'].shape[0] != pvProfile['P'].shape[0]:
                logging.info(
                    "reset kmat archive not same pression shape" + self.tag)
                self.reset()
            else:
                if not np.allclose(baseProfile['P'],
                                   pvProfile['P'],
                                   rtol=0,
                                   atol=1e-4):
                    logging.info(
                        "reset kmat : not same pression values " + self.tag)
                    self.reset()
                else:
                    # check if same gas are present
                    for item in kmat.kmatrix.gas_list:
                        if kmat.kmatrix[item] is not None:
                            if pvKmat.kmatrix[item] is None:
                                logging.info("reset kmat new gas" + self.tag)
                                self.reset()
                        if pvKmat.kmatrix[item] is not None:
                            if kmat.kmatrix[item] is None:
                                logging.info(
                                    "reset kmat missing gas " + self.tag)
                                self.reset()

        if len(self.aKmat) == self.archiveDepth:
            self.aKmat.pop(0)

        self.aKmat.append(kmat)
        fk.close()

    def getLengh(self):
        return len(self.aKmat)

    def reset(self):
        """ purge the archive """
        logging.info(" >> reset kMat archive for run " + str(
            self.runNumber) + " " + self.tag)
        self.aKmat = []

    def setDepth(self, depth):
        if (depth >= 0):
            if (depth < self.archiveDepth):
                while (len(self.aKmat) > depth):
                    self.aKmat.pop(0)
        self.archiveDepth = depth
