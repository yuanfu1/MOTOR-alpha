# -*- coding: utf-8 -*-
import unittest
import rmodel
from test import rttovgui_unittest_class


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):
        self.p = rmodel.project.Project()
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        self.p.openProfile(profileName)

    def tearDown(self):
        pass

    def test_read_profile(self):
        self.assertGreater(self.p.myProfile['BE'],  0.20)

    def test_run_pc_1(self):
        self.p.dropCoefficients()
        pccompatfile = "/rttov9pred101L/rtcoef_metop_2_iasi_pcrttov_compat.H5"
        coefFile = self.p.config.ENV['RTTOV_GUI_COEFF_DIR'] + pccompatfile
        pccoef = "/pc/pccoef_metop_2_iasi_landsea_trace_nlte.H5"
        pcFile = self.p.config.ENV['RTTOV_GUI_COEFF_DIR'] + pccoef
        self.p.myCoeffs.fileName["standard"] = coefFile
        self.p.myCoeffs.fileName["PC"] = pcFile
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        self.assertFalse(self.p.isPC())
        self.p.myOption["ADDPC"].value = True
        self.assertTrue(self.p.isPC())

        self.p.myProfile["SKIN"]["SURFTYPE"] = 1
        err = self.p.runPC()
        self.assertEqual(err, 0)

    def test_run_pc_2(self):
        print("-------------------AIRS-------------------------")
        self.p.myProfile["SKIN"]["SURFTYPE"] = 1
        pccompatfile = "/rttov9pred101L/rtcoef_eos_2_airs_pcrttov_compat.H5"
        coefFile = self.p.config.ENV['RTTOV_GUI_COEFF_DIR'] + pccompatfile
        pcFile = self.p.config.ENV[
            'RTTOV_GUI_COEFF_DIR'] + "/pc/pccoef_eos_2_airs_landsea.H5"
        self.p.myCoeffs.fileName["standard"] = coefFile
        self.p.myCoeffs.fileName["PC"] = pcFile
        print(self.p.myCoeffs.fileName)
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        err = self.p.runPC()
        self.assertEqual(err, 0)
        self.check_option(self.p)
        print("-------------------IASI NG-------------------------")
        pcfile = "/rttov9pred101L/rtcoef_metopsg_1_iasing_pcrttov_compat.H5"
        coefFile = self.p.config.ENV['RTTOV_GUI_COEFF_DIR'] + pcfile
        pcFile = self.p.config.ENV[
            'RTTOV_GUI_COEFF_DIR'] + "/pc/pccoef_metopsg_1_iasing_sea.H5"
        self.p.myCoeffs.fileName["standard"] = coefFile
        self.p.myCoeffs.fileName["PC"] = pcFile
        print(self.p.myCoeffs.fileName)
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        err = self.p.runPC()
        self.assertEqual(err, 0)
        self.check_option(self.p)

    def test_read_profile_ascii_1(self):
        print(">>>>>>>>test_read_profile_ascii_1")
        profileName = (self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] +
            "/../profile-datasets-py/us76_50lev_allgas/001.py")
        err = self.p.openAsciiProfile(profileName)
        self.assertEqual(err, 0)

    def test_read_profile_ascii_2(self):
        print(">>>>>>>>test_read_profile_ascii_2")
        profileName = (self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] +
            "/../profile-datasets-py/varying101lev_clw/005.py")
        err = self.p.openAsciiProfile(profileName)
        self.assertEqual(err, 0)

    def test_gas_for_iasi(self):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>test_gas_for_pc")
        p = rmodel.project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metop_2_iasi.H5"

        p.myCoeffs.fileName["standard"] = coefFile
        print(p.myCoeffs.fileName)
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        print("getFMV_GAS_POS:")
        print(self.p.myCoeffs.getFMV_GAS_POS())
        if self.p.myCoeffs.getFMV_GAS_POS() is None:
            self.assertEqual(self.p.myCoeffs.hasGas("O3"), False)
        print("has O3 ?",  self.p.myCoeffs.hasGas("O3"))
        print("end >>>> test_gas_for_pc")

    def test_pc(self):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>test pc")
        p = rmodel.project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("---------------------TEST PC ---------------------------------")
        p.dropCoefficients()
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metop_2_iasi_pcrttov_compat.H5"
        pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/pc/pccoef_metop_2_iasi_landsea_trace_nlte.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        p.myCoeffs.fileName["PC"] = pcFile
        print(p.myCoeffs.fileName)
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        print("isPC? ", p.isPC())
        self.assertFalse(p.isPC())
        print("isPCtrace? ", p.isPCtrace())
        self.assertFalse(p.isPCtrace())
        p.myOption["ADDPC"].value = True
        p.ctrlCoherence()
        self.assertTrue(p.isPC())
        self.assertTrue(p.isPCtrace())
        self.check_option(p)
        print("isPC ?", p.isPC())
        self.assertTrue(p.isPC())
        p.myProfile['SKIN']['SURFTYPE'] = 1
        print("-------------------runPC------------------------")
        p.myOption["ADDRADREC"].value = True
        p.myOption["ADDSOLAR"].value = True
        p.myOption.display()
        err = p.runPC()
        self.assertEqual(err, 0)
        self.assertFalse(p.myOption["ADDSOLAR"].value)
        self.check_option(p)
        print("-------------------AIRS-------------------------")
        p.dropCoefficients()
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_eos_2_airs_pcrttov_compat.H5"
        pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/pc/pccoef_eos_2_airs_landsea.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        p.myCoeffs.fileName["PC"] = pcFile
        print(p.myCoeffs.fileName)
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        err = p.runPC()
        self.assertEqual(err, 0)
        err = p.runPCK()
        self.assertEqual(err, 0)
        self.check_option(p)
        print("-------------------IASI NG-------------------------")
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metopsg_1_iasing_pcrttov_compat.H5"
        pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/pc/pccoef_metopsg_1_iasing_sea.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        p.myCoeffs.fileName["PC"] = pcFile
        print(p.myCoeffs.fileName)
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        err = p.runPC()
        self.assertEqual(err, 0)
        err = p.runPCK()
        self.assertEqual(err, 0)
        self.check_option(p)


if __name__ == "__main__":
    unittest.main()
