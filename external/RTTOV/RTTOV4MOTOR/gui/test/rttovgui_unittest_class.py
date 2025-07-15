'''
Created on May 14, 2014

@author: pascale
'''
import unittest
import rmodel
import os


class RttovGuiUnitTest(unittest.TestCase):

    def mk_coeff_list(self):
        coeff_dir = os.environ["RTTOV_GUI_COEFF_DIR"]
        coeff_dir_v7_54 = os.path.join(coeff_dir, "rttov7pred54L")
        coeff_dir_v7_101 = os.path.join(coeff_dir, "rttov7pred101L")
        coeff_dir_v8_54 = os.path.join(coeff_dir, "rttov8pred54L")
        coeff_dir_v8_101 = os.path.join(coeff_dir, "rttov8pred101L")
        coeff_dir_v9_54 = os.path.join(coeff_dir, "rttov9pred54L")
        coeff_dir_v9_101 = os.path.join(coeff_dir, "rttov9pred101L")
        coeff_dir_v13_54 = os.path.join(coeff_dir, "rttov13pred54L")
        coeff_dir_v13_101 = os.path.join(coeff_dir, "rttov13pred101L")
        stdCoeffListToTest = [os.path.join(coeff_dir_v7_54,
                                           "rtcoef_metop_2_amsua.dat"),
                              os.path.join(coeff_dir_v7_101,
                                           "rtcoef_metop_2_iasi.H5"),
                              os.path.join(coeff_dir_v7_54,
                                           "rtcoef_dmsp_19_ssmis.dat"),
                              os.path.join(coeff_dir_v8_54,
                                           "rtcoef_metop_2_hirs.dat"),
                              os.path.join(coeff_dir_v8_101,
                                           "rtcoef_metop_2_iasi.H5"),
                              os.path.join(coeff_dir_v9_54,
                                           "rtcoef_metop_2_avhrr.dat"),
                              os.path.join(coeff_dir_v9_101,
                                           "rtcoef_jpss_0_cris.H5")]
        return stdCoeffListToTest

    def test_dummy(self):
        pass

    def controleOptionsGas(self):

        for gas in ["CO2", "N2O", "CO", "CH4"]:
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
            self.assertFalse(self.p.myOption["ADDRADREC"].status)
            self.assertFalse(self.p.myOption["IPCBND"].status)
            self.assertFalse(self.p.myOption["IPCREG"].status)
            self.assertFalse(self.p.myOption["ADDRADREC"].value)
            self.assertFalse(self.p.myOption["ADDPC"].value)

    def check_option(self, p):

        if p is None:
            return
        if p.isPC():
            self.assertTrue(p.myOption["ADDPC"].value)
            self.assertTrue(p.myOption["ADDPC"].status)
        else:
            self.assertFalse(p.myOption["ADDPC"].value)
            self.assertTrue(p.myOption["DO_LAMBERTIAN"].status)

        if p.myCoeffs.hasPC():
            self.assertTrue(p.myOption["ADDPC"].status)
        else:
            self.assertFalse(p.myOption["ADDPC"].status)
            self.assertFalse(p.myOption["ADDPC"].value)

        if p.myCoeffs.isPCClouds() and p.myProfile.hasClouds():
            self.assertTrue(p.myOption["ADDPC"].status)

        if p.myCoeffs.hasSolar():
            if not p.isPC():
                self.assertTrue(p.myOption["ADDSOLAR"].status)
        else:
            self.assertFalse(p.myOption["ADDSOLAR"].status)
            self.assertFalse(p.myOption["ADDSOLAR"].value)

        if not p.myCoeffs.isMW():
            self.assertFalse(p.myOption["FASTEM_VERSION"].status)
            self.assertFalse(p.myOption["CLW_DATA"].status)
            self.assertFalse(p.myOption["CLW_DATA"].value)

        if p.myCoeffs.isMW():
            self.assertTrue(p.myOption["FASTEM_VERSION"].status)
            self.assertTrue(p.myOption["SUPPLY_FOAM_FRACTION"].status)

            if p.myProfile["CLW"] is not None:
                self.assertTrue(p.myOption["CLW_DATA"].status)
            self.assertTrue(p.myOption["DO_LAMBERTIAN"].status)
            self.assertFalse(p.myOption["ADDSOLAR"].status)
            self.assertFalse(p.myOption["DO_NLTE_CORRECTION"].status)
            self.assertFalse(p.myOption["ADDAEROSL"].status)
            self.assertFalse(p.myOption["ADDCLOUDS"].status)
            self.assertFalse(p.myOption["CLDCOL_THRESHOLD"].status)
            self.assertFalse(p.myOption["OZONE_DATA"].status)
            self.assertFalse(p.myOption["CO2_DATA"].status)
            self.assertFalse(p.myOption["N2O_DATA"].status)
            self.assertFalse(p.myOption["CO_DATA"].status)
            self.assertFalse(p.myOption["CH4_DATA"].status)
            self.assertFalse(p.myOption["ADDPC"].status)
            self.assertFalse(p.myOption["IPCBND"].status)
            self.assertFalse(p.myOption["IPCREG"].status)
            self.assertFalse(p.myOption["ADDRADREC"].status)

            self.assertFalse(p.myOption["ADDSOLAR"].status)
            self.assertFalse(p.myOption["DO_NLTE_CORRECTION"].status)
            self.assertFalse(p.myOption["ADDAEROSL"].status)
            self.assertFalse(p.myOption["ADDCLOUDS"].status)

            self.assertFalse(p.myOption["OZONE_DATA"].value)
            self.assertFalse(p.myOption["CO2_DATA"].value)
            self.assertFalse(p.myOption["N2O_DATA"].value)
            self.assertFalse(p.myOption["CO_DATA"].value)
            self.assertFalse(p.myOption["CH4_DATA"].value)
            self.assertFalse(p.myOption["ADDPC"].value)
            self.assertFalse(p.myOption["ADDRADREC"].value)

        else:  # not MW
            self.assertFalse(p.myOption["FASTEM_VERSION"].status)
            self.assertFalse(p.myOption["CLW_DATA"].status)
            self.assertFalse(p.myOption["CLW_DATA"].value)
            if p.isPC() and p.isPCtrace():
                for gas in ("CO", "CO2", "N2O", "CH4"):
                    if p.myCoeffs.hasGas(gas) and p.myProfile.hasGas(gas):
                        self.assertFalse(p.myOption[gas + "_DATA"].value)
                        self.assertTrue(
                            p.myOption[gas + "_DATA"].status)
                self.assertFalse(p.myOption["ADDAEROSL"].value)
                self.assertFalse(p.myOption["ADDAEROSL"].status)
                self.assertFalse(p.myOption["SO2_DATA"].value)
                self.assertFalse(p.myOption["SO2_DATA"].status)

            if p.isPC() and not p.isPCtrace():
                for gas in ("CO", "CO2", "N2O", "CH4", "SO2"):
                    if p.myCoeffs.hasGas(gas) and p.myProfile.hasGas(gas):
                        self.assertFalse(p.myOption[gas + "_DATA"].value)
                        self.assertFalse(
                            p.myOption[gas + "_DATA"].status)
                self.assertFalse(p.myOption["OZONE_DATA"].value)
                self.assertFalse(p.myOption["OZONE_DATA"].status)
                self.assertFalse(p.myOption["ADDAEROSL"].value)
                if p.myCoeffs.hasNLTE():
                    self.assertTrue(p.myOption["DO_NLTE_CORRECTION"].status)
                else:
                    self.assertFalse(p.myOption["DO_NLTE_CORRECTION"].status)
            else:
                for gas in ("CO", "CO2", "N2O", "CH4", "SO2"):
                    if p.myCoeffs.hasGas(gas) and p.myProfile.hasGas(gas):
                        print(">>> projet has gas ?", gas, p.hasGas(gas))
                        self.assertTrue(p.myOption[gas + "_DATA"].status)
                    else:
                        print(">>>", gas, p.myOption[gas + "_DATA"].value)
                        print(p.myOption[gas + "_DATA"].status)
                        print(">>> profile hasGas?", p.myProfile.hasGas(gas))
                        print(">>> coeff hasGas?", p.myCoeffs.hasGas(gas))
                        print(">>> project hasGas?", p.hasGas(gas))
                        self.assertFalse(p.myOption[gas + "_DATA"].value)
                        self.assertFalse(
                            p.myOption[gas + "_DATA"].status)
                if p.myCoeffs.hasAerosols() and p.myProfile.hasAerosols():
                    self.assertTrue(p.myOption["ADDAEROSL"].status)
                else:
                    self.assertFalse(p.myOption["ADDAEROSL"].value)
                    self.assertFalse(p.myOption["ADDAEROSL"].status)

                if p.myCoeffs.hasGas("O3") and p.myProfile.hasGas("O3"):
                    print("OZONE_DATA status:")
                    print(p.myOption["OZONE_DATA"].status)
                    print("%%%%%has gas ?", p.hasGas("O3"))
                    self.assertTrue(p.myOption["OZONE_DATA"].status)
                else:
                    self.assertFalse(p.myOption["OZONE_DATA"].value)
                    self.assertFalse(p.myOption["OZONE_DATA"].status)
        if not p.myCoeffs.hasNLTE():
            self.assertFalse(p.myOption["DO_NLTE_CORRECTION"].status)
            self.assertFalse(p.myOption["DO_NLTE_CORRECTION"].value)
        if p.myOption["ADDCLOUDS"].value:
            self.assertTrue(p.myOption["CLDCOL_THRESHOLD"].status)
        else:
            self.assertFalse(p.myOption["CLDCOL_THRESHOLD"].status)
        if not p.myOption["ADDCLOUDS"].status:
            self.assertFalse(p.myOption["CLDCOL_THRESHOLD"].status)


class test_test(RttovGuiUnitTest):

    def test_check_option(self):
        print("test rttovgui_unittest_class")
        p = rmodel.project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/us76_50lev_allgas.H5"
        print("openProfile", profileName)
        p.openProfile(profileName)

        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metop_2_iasi_pcrttov_compat.H5"
        pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/pc/pccoef_metop_2_iasi_landsea_trace_nlte.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        p.myCoeffs.fileName["PC"] = pcFile

        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        p.ctrlCoherence()
        self.check_option(p)


if __name__ == "__main__":
    unittest.main()
