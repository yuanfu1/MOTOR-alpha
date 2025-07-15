import unittest
import rmodel
import logging
import glob
from test import rttovgui_unittest_class
import re
import os


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):
        level_logging = logging.DEBUG
        self.p = rmodel.project.Project()

        logging.basicConfig(
            filename=(self.p.config.ENV['GUI_WRK_DIR'] +
                      "/rttovgui_unittest_test_aerosol.log"),
            format=("[%(asctime)s] %(levelname)s [%(module)s:"
                    "%(funcName)s:%(lineno)d] %(message)s"),
            level=level_logging,
            datefmt="%Y:%m:%d %H:%M:%S",
            filemode="w")

    def test_aerosols(self):
        print(">>>>>>>>>>>>start test_aerosols")
        list_aer = []
        profile_with = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"
        profile_without = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard101lev_allgasref.H5"

        for coefFile in glob.glob(self.p.config.ENV["RTTOV_GUI_COEFF_DIR"] +
                                  "/cldaer_ir/scaer*.dat"):
            print("add ", coefFile)
            list_aer.append(coefFile)

        for coefFile in list_aer:
            self.p.myCoeffs.fileName["aerosols"] = coefFile
            std = re.sub('scaer', 'rt', re.sub(
                "cldaer_ir", "rttov7pred54L", coefFile))
            if os.path.exists(std):
                self.p.myCoeffs.fileName["standard"] = std

                print(self.p.myCoeffs.fileName)
                err = self.p.loadCoefficients()
                self.assertEqual(err, 0)
                print(">>>>coefFile;", self.p.myCoeffs.fileName, " loaded")
                err = self.p.openProfile(profile_with, 1)
                self.p.ctrlCoherence()
                self.check_option(self.p)
                err = self.p.runDirect()
                self.assertEqual(err, 0)
                err = self.p.runK()
                self.assertEqual(err, 0)
                self.check_option(self.p)
                self.assertTrue(self.p.myOption["ADDAEROSL"].status)
                self.p.myOption["ADDAEROSL"].value = True
                err = self.p.runDirect()
                self.assertEqual(err, 0)
                err = self.p.runK()
                self.assertEqual(err, 0)
                self.check_option(self.p)
                err = self.p.openProfile(profile_without, 1)
                self.assertEqual(err, 0)
                self.p.ctrlCoherence()
                self.check_option(self.p)

                self.assertFalse(self.p.myOption["ADDAEROSL"].status)
                self.assertFalse(self.p.myOption["ADDAEROSL"].value)

                print(">>>>>>>>>>>>>>>>>>OK for ", coefFile, profile_with, 1)

        print(">>>>>>>>>>>>>>>>> end test_aerosols")


if __name__ == "__main__":
    unittest.main()
