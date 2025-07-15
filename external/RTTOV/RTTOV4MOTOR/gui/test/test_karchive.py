'''
Created on Feb 16, 2016

@author: pascale
'''

import unittest
import rmodel
import logging
from test import rttovgui_unittest_class


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):
        level_logging = logging.DEBUG
        self.p = rmodel.project.Project()
        logging.basicConfig(filename=(self.p.config.ENV['GUI_WRK_DIR'] +
                                      "/rttovgui_unittest_test_karchive.log"),
                            format=("[%(asctime)s] %(levelname)s [%(module)s:"
                                    "%(funcName)s:%(lineno)d] %(message)s"),
                            level=level_logging,
                            datefmt="%Y:%m:%d %H:%M:%S",
                            filemode="w")

    def tearDown(self):
        pass

    def testName(self):
        pass

    def test1(self):

        self.p.debug = True
        print("Configuration : ", self.p.config.ENV['RTTOV_GUI_PREFIX'])
        print("open profile")

        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        print(profileName)
        err = self.p.openProfile(profileName)
        self.assertEqual(err, 0)

        coefFile = self.p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metop_2_iasi.H5"

        print(coefFile)

        self.p.myCoeffs.fileName["standard"] = coefFile
        self.p.loadCoefficients()
        self.p.myOption["IR_SEA_EMIS_MODEL"].value = 1
        self.p.ctrlCoherence()
        self.check_option(self.p)
        err = self.p.runDirect()
        self.assertEqual(err, 0)
        err = self.p.runK()
        self.assertEqual(err, 0)
        karchive = rmodel.karchive.Karchive(2, tag="test")
        karchive.update(self.p.KMatrixFileName)
        lenArchive = karchive.getLengh()
        self.assertEqual(lenArchive, 1)
        self.assertEqual(err, 0)
        err = self.p.runK()
        self.assertEqual(err, 0)
        karchive.update(self.p.KMatrixFileName)
        lenArchive = karchive.getLengh()
        self.assertEqual(lenArchive, 2)
        err = self.p.runK()
        self.assertEqual(err, 0)
        karchive.update(self.p.KMatrixFileName)
        lenArchive = karchive.getLengh()
        self.assertEqual(lenArchive, 2)

        # try another profile
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/us76_50lev_allgas.H5"
        print(profileName)
        err = self.p.openProfile(profileName)
        self.assertEqual(err, 0)
        err = self.p.runK()
        self.assertEqual(err, 0)
        karchive.update(self.p.KMatrixFileName)
        lenArchive = karchive.getLengh()
        self.assertEqual(lenArchive, 1)

        karchive.setDepth(0)
        lenArchive = karchive.getLengh()
        self.assertEqual(lenArchive, 0)
        karchive.setDepth(3)
        lenArchive = karchive.getLengh()
        self.assertEqual(lenArchive, 0)

        print(">>>> CO_DATA was ", self.p.myOption["CO_DATA"])
        self.p.myOption["CO_DATA"].value = True
        err = self.p.runK()
        self.assertEqual(err, 0)
        karchive.update(self.p.KMatrixFileName)
        lenArchive = karchive.getLengh()
        self.assertEqual(lenArchive, 1)

        self.p.myOption["CO_DATA"].value = False
        err = self.p.runK()
        self.assertEqual(err, 0)
        karchive.update(self.p.KMatrixFileName)
        lenArchive = karchive.getLengh()
        self.assertEqual(lenArchive, 1)

        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_o3ref.H5"
        err = self.p.openProfile(profileName)
        self.assertEqual(err, 0)

        err = self.p.runK()
        self.assertEqual(err, 0)
        karchive.update(self.p.KMatrixFileName)
        lenArchive = karchive.getLengh()
        self.assertEqual(lenArchive, 1)

        karchive.setDepth(0)
        err = self.p.runK()
        self.assertEqual(err, 0)
        karchive.update(self.p.KMatrixFileName)
        lenArchive = karchive.getLengh()
        self.assertEqual(lenArchive, 0)


if __name__ == "__main__":

    unittest.main()
