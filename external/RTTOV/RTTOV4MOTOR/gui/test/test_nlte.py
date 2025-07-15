import unittest
import rttov
import rmodel
import logging
import os

from test import rttovgui_unittest_class


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):
        level_logging = logging.DEBUG
        self.p = rmodel.project.Project()

        logging.basicConfig(
            filename=(self.p.config.ENV['GUI_WRK_DIR'] +
                      "/rttovgui_unittest_test_nlte.log"),
            format=("[%(asctime)s] %(levelname)s [%(module)s:%(funcName)s"
                    ":%(lineno)d] %(message)s"),
            level=level_logging,
            datefmt="%Y:%m:%d %H:%M:%S",
            filemode="w")
        self.stdCoeffListToTest = self.mk_coeff_list()

    def controleOptionsGas(self):

        for gas in ["CO2", "N2O", "CO", "CH4", "SO2"]:
            print("controle gas", gas)
            if self.p.myCoeffs.hasGas(gas) and self.p.myProfile.hasGas(gas):
                self.assertTrue(self.p.myOption[gas + "_DATA"].status)
            else:
                self.assertFalse(self.p.myOption[gas + "_DATA"].value)
                self.assertFalse(self.p.myOption[gas + "_DATA"].status)
        gas = "OZONE"
        if self.p.myCoeffs.hasGas("O3") and self.p.myProfile.hasGas("O3"):
            print("gas=OZONE")
            self.assertTrue(self.p.myOption[gas + "_DATA"].status)
        else:
            print(("self.p.myProfile.hasGas('O3')",
                   self.p.myProfile.hasGas("O3")))
            print(self.p.myProfile["O3"])
            print("self.p.myCoeff.hasGas('O3')", self.p.myCoeffs.hasGas("O3"))
            self.assertFalse(self.p.myOption[gas + "_DATA"].value)
            self.assertFalse(self.p.myOption[gas + "_DATA"].status)

    def controlePC(self):
        if (self.p.myOption["ADDPC"].value):
            self.assertTrue(self.p.myOption["ADDRADREC"].status)
            self.assertTrue(self.p.myOption["IPCBND"].status)
            self.assertTrue(self.p.myOption["IPCREG"].status)
        else:
            self.assertFalse(self.p.myOption["ADDRADREC"].value)
            self.assertFalse(self.p.myOption["IPCBND"].status)
            self.assertFalse(self.p.myOption["IPCREG"].status)
            self.assertFalse(self.p.myOption["ADDRADREC"].status)

    def test_full(self):

        list_std = self.stdCoeffListToTest
        prof_dir = self.p.config.ENV["RTTOV_GUI_PROFILE_DIR"]
        list_prof = [os.path.join(prof_dir, "standard101lev_allgas.H5"),
                     os.path.join(prof_dir, "cldaer50lev_co2o3.H5")]

        for coefFile in list_std:

            self.p.myCoeffs.fileName["standard"] = coefFile
            err = self.p.loadCoefficients()
            self.assertEqual(err, 0)
            print(">>>>>>>>>>>>>>>>>>>>>coefFile;", coefFile, " loaded")
            for prof in list_prof:
                print(">>>>>>>>>>>>>>>>>>>open profile ", prof)
                nb = rttov.profile.getNumberOfProfiles(prof)
                for n in range(1, nb + 1):
                    print(">>>>>>>>>>>>>>>profile ", prof, "number ", n)
                    self.p.openProfile(prof, n)
                    self.p.myOption["DO_NLTE_CORRECTION"].value = True
                    self.p.ctrlCoherence()
                    print(("################has NLTE ?",
                           self.p.myCoeffs.hasNLTE()))
                    self.check_option(self.p)
                    if (
                            self.p.myCoeffs.hasNLTE() and
                            (not self.p.isPC() or not self.p.myCoeffs.isMW())):
                        print("##################cas NLTE")
                        self.assertTrue(self.p.myOption[
                                        "DO_NLTE_CORRECTION"].value)
                        self.assertTrue(self.p.myOption[
                                        "DO_NLTE_CORRECTION"].status)
                    else:
                        self.assertFalse(self.p.myOption[
                                         "DO_NLTE_CORRECTION"].value)
                        self.assertFalse(self.p.myOption[
                                         "DO_NLTE_CORRECTION"].status)
                    print("run direct with ")
                    print(("NLTE", self.p.myOption[
                        "DO_NLTE_CORRECTION"].value,
                        self.p.myOption[
                        "DO_NLTE_CORRECTION"].status))
                    err = self.p.runDirect()
                    self.assertEqual(err, 0)
                    err = self.p.runK()
                    self.assertEqual(err, 0)
                    self.check_option(self.p)
                    print(">>>>>>>>>>>>>>>>>>>>>>OK for ", coefFile, prof, n)

        print(">>>>>>>>>>>>>>>>> end test_nlte")


if __name__ == "__main__":
    unittest.main()
