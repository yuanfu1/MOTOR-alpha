#! /usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
import rttov
import h5py
import rmodel
import logging
from test import rttovgui_unittest_class


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):

        self.p = rmodel.project.Project()

        level_logging = logging.DEBUG
        logging.basicConfig(
            filename=(self.p.config.ENV['GUI_WRK_DIR'] +
                      "/rttovgui_test_clouds.log"),
            format=("[%(asctime)s] %(levelname)s [%(module)s:%(funcName)s:"
                    "%(lineno)d] %(message)s"),
            level=level_logging,
            datefmt="%Y:%m:%d %H:%M:%S",
            filemode="w")
        logging.info("start test cloud")

        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"
        err = self.p.openProfile(profileName)
        self.assertEqual(err, 0)

    def tearDown(self):
        pass

    def test_read_profile(self):
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        err = self.p.openProfile(profileName)
        self.assertEqual(err, 0)
        self.assertGreater(self.p.myProfile['BE'],  0.20)

    def test_wv_clw(self):
        print("------------------test CLW-------------------------------")
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard101lev_clw.H5"
        self.p.openProfile(profileName)
        coefFile = self.p.config.ENV[
            'RTTOV_GUI_COEFF_DIR'] + "/rttov7pred54L/rtcoef_dmsp_18_ssmis.dat"
        self.p.myCoeffs.fileName["standard"] = coefFile
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        self.p.ctrlCoherence()
        # not yet implemented
        # self.assertTrue(self.p.myOption["CLW_DATA"].status)
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"
        self.p.openProfile(profileName)
        self.p.ctrlCoherence()
        self.assertFalse(self.p.myOption["CLW_DATA"].status)
        self.assertFalse(self.p.myOption["CLW_DATA"].value)

    def test_pc_cloud1(self):
        print("-------------------PC CLOUD-1----------------------")

        self.p.myOption["ADDPC"].value = True
        pcfile = "/rttov9pred101L/rtcoef_metop_2_iasi_pcrttov_compat.H5"
        coefFile = self.p.config.ENV['RTTOV_GUI_COEFF_DIR'] + pcfile
        pccoef = "/pc/pccoef_metop_2_iasi_landsea_trace_nlte.H5"
        pcFile = self.p.config.ENV['RTTOV_GUI_COEFF_DIR'] + pccoef

        self.p.myCoeffs.fileName["standard"] = coefFile
        self.p.myCoeffs.fileName["PC"] = pcFile
        self.p.myProfile['SKIN']['SURFTYPE'] = 1

        print(self.p.myCoeffs.fileName)
        self.p.myOption.display()
        err = self.p.loadCoefficients()
        self.assertEqual(err, 0)
        self.p.myProfile['SKIN']['SURFTYPE'] = 1
        print("is cloud ?", self.p.myCoeffs.isPCClouds())

        err = self.p.runPC()
        self.assertEqual(err, 0)

        print("ADDCLOUD ?:", self.p.myOption["ADDCLOUDS"].value)
        self.p.myOption["ADDCLOUDS"].value = True
        err = self.p.runPC()
        self.assertEqual(err, 0)
        self.check_option(self.p)
        if self.p.myCoeffs.isPCClouds():
            self.assertTrue(self.p.myOption["ADDCLOUDS"].value)
            self.assertTrue(self.p.myOption["ADDCLOUDS"].status)
        else:
            self.assertFalse(self.p.myOption["ADDCLOUDS"].value)
            self.assertFalse(self.p.myOption["ADDCLOUDS"].status)
        print(">>>>>>>>>>>>> ADDCLOUD ?:", self.p.myOption["ADDCLOUDS"].value)
        err = self.p.runPCK()
        self.assertEqual(err, 0)

    def test_pc_cloud2(self):
        print(">>>>>>>>>>>test_pc_cloud2")
        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        self.p.openProfile(profileName)

        print(("profile", profileName, " has cloud ?",
               self.p.myProfile.hasClouds()))
        self.assertFalse(self.p.myProfile.hasClouds())

        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"
        self.p.openProfile(profileName)

        print(("profile", profileName, " has cloud ?",
               self.p.myProfile.hasClouds()))
        self.assertTrue(self.p.myProfile.hasClouds())


if __name__ == "__main__":
    unittest.main()
