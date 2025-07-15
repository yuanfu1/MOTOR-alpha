'''
Created on May 7, 2014

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
                                      "/rttovgui_inittest.log"),
                            format=("[%(asctime)s] %(levelname)s "
                                    "[%(module)s:%(funcName)s:%"
                                    "(lineno)d] %(message)s"),
                            level=level_logging,
                            datefmt="%Y:%m:%d %H:%M:%S",
                            filemode="w")
        logging.info("start main controller")

        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"

        self.p.openProfile(profileName)
        self.p.myOption["IR_SEA_EMIS_MODEL"].value = 1
        self.stdCoeffListToTest = self.mk_coeff_list()

    def xtest_solar1(self):
        self.p.myCoeffs.fileName["standard"] = self.p.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + "/rttov9pred54L/rtcoef_eos_2_modis.dat"
        err = self.p.loadCoefficients()

        self.p.ctrlCoherence()
        self.p.myOption.display()
        self.p.myOption["IR_SEA_EMIS_MODEL"].value = 1
        self.assertEqual(err, 0)
        self.assertTrue(self.p.myOption["ADDSOLAR"].status)
        self.assertFalse(self.p.myOption["ADDSOLAR"].value)
        self.assertFalse(self.p.myOption["RAYLEIGH_SINGLE_SCATT"].status)
        self.assertFalse(self.p.myOption["ADDCLOUDS"].value)
        self.assertFalse(self.p.myOption["ADDAEROSL"].value)
        print("OZONE?", self.p.myOption["OZONE_DATA"].value)
        self.controleOptionsGas()

        print("OZONE?", self.p.myOption["OZONE_DATA"].value)

        self.assertFalse(self.p.myOption["CO2_DATA"].value)
        self.assertFalse(self.p.myOption["N2O_DATA"].value)
        self.assertFalse(self.p.myOption["CO_DATA"].value)
        self.assertFalse(self.p.myOption["ADDPC"].status)
        self.assertFalse(self.p.myOption["ADDPC"].value)
        print("OZONE?", self.p.myOption["OZONE_DATA"].value)
        self.p.myProfile.removeGas("O3")
        print("OZONE?", self.p.myOption["OZONE_DATA"].value)
        self.p.ctrlCoherence()
        print("OZONE?", self.p.myOption["OZONE_DATA"].value)
        self.assertFalse(self.p.myOption["OZONE_DATA"].value)
        self.assertFalse(self.p.myOption["OZONE_DATA"].status)
        self.assertFalse(self.p.myProfile.hasGas("O3"))
        err = self.p.runDirect()
        self.assertEqual(err, 0)
        err = self.p.runK()
        self.assertEqual(err, 0)
        print("OZONE?", self.p.myOption["OZONE_DATA"].value)
        err = self.p.runDirect()
        self.assertEqual(err, 0)
        self.assertFalse(self.p.myOption["OZONE_DATA"].value)
        self.assertFalse(self.p.myOption["OZONE_DATA"].status)
        err = self.p.dropCoefficients()
        self.assertEqual(err, 0)
        print("***************************************************test iasi *")
        pcfile = "/rttov9pred101L/rtcoef_metop_2_iasi_pcrttov_compat.H5"
        self.p.myCoeffs.fileName["standard"] = self.p.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + pcfile
        err = self.p.loadCoefficients()
        print("OZONE?", self.p.myOption["OZONE_DATA"].value)
        self.assertEqual(err, 0)
        self.p.ctrlCoherence()
        self.p.myOption.display()
        self.controleOptionsGas()
        self.assertFalse(self.p.myOption["OZONE_DATA"].value)
        self.assertFalse(self.p.myOption["OZONE_DATA"].status)

        self.assertFalse(self.p.myOption["ADDSOLAR"].value)
        self.assertTrue(self.p.myOption["ADDSOLAR"].status)

        self.assertFalse(self.p.myOption["ADDPC"].value)
        self.controlePC()

        print(("CH4?=", self.p.myOption["CH4_DATA"].value, "status?",
               self.p.myOption["CH4_DATA"].status))
        self.p.myProfile.removeGas("CH4")

        self.p.ctrlCoherence()
        print(("CH4?=", self.p.myOption["CH4_DATA"].value, "status?",
               self.p.myOption["CH4_DATA"].status))
        self.controleOptionsGas()

        self.assertFalse(self.p.myOption["CH4_DATA"].value)
        self.assertFalse(self.p.myOption["CH4_DATA"].status)

        self.assertFalse(self.p.myOption["ADDPC"].value)
        self.assertTrue(self.p.myOption["ADDSOLAR"].status)
        err = self.p.runDirect()
        self.assertEqual(err, 0)
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ok ")
        pccoef = "/pc/pccoef_metop_2_iasi_landsea_trace_nlte.H5"
        self.p.myCoeffs.fileName["PC"] = self.p.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + pccoef
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        self.p.ctrlCoherence()

        self.assertTrue(self.p.myOption["ADDPC"].status)
        self.controlePC()

        self.assertTrue(self.p.myOption["ADDSOLAR"].status)

        self.assertFalse(self.p.myOption["ADDCLOUDS"].status)
        self.assertFalse(self.p.myOption["ADDCLOUDS"].value)

        self.p.myOption["ADDPC"].value = True
        self.p.ctrlCoherence()
        self.controlePC()

        self.assertFalse(self.p.myOption["ADDSOLAR"].status)
        self.assertFalse(self.p.myOption["ADDSOLAR"].value)
        self.assertFalse(self.p.myOption["ADDCLOUDS"].status)
        self.assertFalse(self.p.myOption["ADDCLOUDS"].value)
        err = self.p.runPC()
        self.assertEqual(err, 0)
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ok 1")
        self.p.myOption["ADDPC"].value = False
        err = self.p.runDirect()
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ok 2")
        self.assertEqual(err, 0)
        self.p.ctrlCoherence()
        self.assertTrue(self.p.myOption["ADDSOLAR"].status)
        self.controlePC()
        err = self.p.runDirect()
        self.assertEqual(err, 0)
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ok 3")
        self.p.myOption["ADDSOLAR"].value = True
        self.p.ctrlCoherence()
        self.assertFalse(self.p.myOption["ADDPC"].status)
        self.controlePC()

        err = self.p.runDirect()
        self.assertEqual(err, 0)
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>ok 4")
        self.p.myOption["N2O_DATA"].value = True
        self.p.ctrlCoherence()

        err = self.p.runDirect()
        self.assertEqual(err, 0)
        print(">>>>>>>>>>>>>>>>>>>>>>>>ok 5")
        self.p.myOption["ADDPC"].value = True

        self.p.ctrlCoherence()
        self.p.myOption["ADDSOLAR"].value = True
        err = self.p.runPC()
        self.assertEqual(err, 0)

        print(">>>>>>>>>>>>>>>>>>>>>>>>>End test solar1")

    def test_removeO3(self):
        print(">>>>>>>>>>>>>>>>>>test_removeO3")

        self.p.myProfile.removeGas("O3")
        list_std = self.stdCoeffListToTest

        for coefFile in list_std:
            print(">>>>>>>>>>>>>>>>>>>>>coefFile;", coefFile)
            self.p.myCoeffs.fileName["standard"] = coefFile
            err = self.p.loadCoefficients()
            self.assertEqual(err, 0)
            print(">>>>>>>>>>>>>>>>>>>>>coefFile;", coefFile)

            self.p.ctrlCoherence()
            self.check_option(self.p)
            err = self.p.runDirect()
            self.assertEqual(err, 0)
            err = self.p.runK()
            self.assertEqual(err, 0)
            self.check_option(self.p)
            print(">>>>>>>>>>>>>>>>>>>>>>OK", coefFile)

        print(">>>>>>>>>>>>>>>>> end test remove O3")

    def xtest_N2O(self):
        print(">>>>>>>>>>>>>>>>>test_N2O")
        coefFile = self.p.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + "/rttov8pred101L/rtcoef_metop_2_iasi.H5"
        self.p.myCoeffs.fileName["standard"] = coefFile
        print("coefFile", coefFile)
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        print("has N2O?")
        has_N2O = self.p.hasGas("N2O")
        self.assertEqual(has_N2O, False)
        print("############has_N2O", has_N2O)
        coefFile = self.p.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + "/rttov9pred54L/rtcoef_msg_2_seviri.dat"
        self.p.myCoeffs.fileName["standard"] = coefFile
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        print("has N2O?")
        has_N2O = self.p.hasGas("N2O")
        print("############has_N2O", has_N2O)
        self.assertEqual(has_N2O, False)
        print(">>>>>>>>>>>>>>>>>end test_N2O")


if __name__ == "__main__":
    unittest.main()
