import unittest
import rmodel
import logging
from test import rttovgui_unittest_class
import os


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):
        level_logging = logging.DEBUG
        self.p = rmodel.project.Project()

        logging.basicConfig(
            filename=(self.p.config.ENV['GUI_WRK_DIR'] +
                      "/rttovgui_unittest_test_update12_2.log"),
            format=("[%(asctime)s] %(levelname)s [%(module)s:%(funcName)s"
                    ":%(lineno)d] %(message)s"),
            level=level_logging,
            datefmt="%Y:%m:%d %H:%M:%S",
            filemode="w")

    def testCamsAerosols(self):
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/aercams50lev_co2o3.H5"
        err = self.p.openProfile(profileName, 2)
        self.assertEqual(err, 0)

        stdcoefFile = self.p.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + "/rttov9pred54L/rtcoef_noaa_19_avhrr.dat"
        self.p.myCoeffs.fileName["standard"] = stdcoefFile
        aerfile = "/cldaer_visir/scaercoef_noaa_19_avhrr_cams.dat"
        aercoefFile = self.p.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + aerfile
        self.p.myCoeffs.fileName["aerosols"] = aercoefFile

        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)

        self.p.ctrlCoherence()
        self.check_option(self.p)
        self.assertEqual(self.p.myOption["ADDAEROSL"].status, True)
        self.p.myOption["ADDAEROSL"].value = True
        self.p.ctrlCoherence()
        self.check_option(self.p)
        self.assertEqual(self.p.myOption["ADDAEROSL"].value, True)

        err = self.p.runDirect()
        self.assertEqual(err, 0)
        err = self.p.runK()
        self.assertEqual(err, 0)
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/aer50lev_co2o3.H5"
        err = self.p.openProfile(profileName, 2)
        self.p.ctrlCoherence()
        self.assertEqual(self.p.myOption["ADDAEROSL"].value, False)
        self.assertEqual(self.p.myOption["ADDAEROSL"].status, False)

        print("OK")

    def testMfasis(self):
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cld50lev_co2o3.H5"
        err = self.p.openProfile(profileName, 4)
        self.assertEqual(err, 0)
        self.assertTrue(self.p.isMfasis)
        coefdir = os.environ["RTTOV_GUI_COEFF_DIR"] + "/rttov9pred54L/"

        stdcoef = coefdir + "rtcoef_msg_3_seviri.dat"
        print("stdcoef:", stdcoef)
        cldcoefdir = os.environ["RTTOV_GUI_COEFF_DIR"] + "/cldaer_visir/"
        cldcoef = cldcoefdir + "sccldcoef_msg_3_seviri.dat"
        lutdir = os.environ["RTTOV_GUI_COEFF_DIR"] + "/mfasis_lut/"
        lut = lutdir + "rttov_mfasis_cld_msg_3_seviri_deff.H5"

        self.p.myCoeffs.fileName["standard"] = stdcoef
        self.p.myCoeffs.fileName["mfasis"] = lut
        self.p.myCoeffs.fileName["clouds"] = cldcoef
        print(self.p.myCoeffs.fileName)
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        self.assertTrue(self.p.isMfasis())
        self.assertTrue(self.p.myOption["ADDSOLAR"].status)
        self.assertTrue(self.p.myOption["ADDCLOUDS"].status)
        optVisScatt = self.p.myOption["VIS_SCATT_MODEL"]
        print("optVisScatt", optVisScatt.odict)
        self.assertEqual(len(optVisScatt.odict), 2)
        self.p.myOption["ADDSOLAR"].value = True
        self.p.myOption["ADDCLOUDS"].value = True
        print("optVisScatt", optVisScatt.odict)
        self.p.ctrlCoherence()
        optVisScatt = self.p.myOption["VIS_SCATT_MODEL"]
        self.assertTrue(optVisScatt.status)
        print("optVisScatt", optVisScatt.odict)
        self.assertEqual(len(optVisScatt.odict), 3)
        self.p.myOption["VIS_SCATT_MODEL"].value = 3
        err = self.p.runDirect()
        print("runDirect returns:", err)
        self.assertEqual(err, 0)
        self.p.dropCoefficients()
        self.p.myCoeffs.fileName["mfasis"] = ""
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        self.p.ctrlCoherence()
        self.assertEqual(len(optVisScatt.odict), 2)


if __name__ == "__main__":
    unittest.main()
