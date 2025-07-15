'''
Created on Dec 13, 2013

@author: pascale
'''
import unittest
import rttov
import h5py
import rmodel
import rmodel
from test import rttovgui_unittest_class


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def test_controleCoherencePC(self):
        self.p = rmodel.project.Project()
        p = self.p
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("p.loadProfile", p.loadProfile)
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metop_2_iasi_pcrttov_compat.H5"
        pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/pc/pccoef_metop_2_iasi_landsea_trace_nlte.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        p.myCoeffs.fileName["PC"] = pcFile
        print(p.myCoeffs.fileName)
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        p.myOption["ADDPC"].value = True
        self.assertTrue(p.isPC())
        p.myOption["IPCBND"].value = 6
        p.myOption["IPCREG"].value = 8
        p.ctrlCoherence()
        self.assertEqual(p.myOption["IPCBND"].value, 1)
        self.assertEqual(p.myOption["IPCREG"].value, 4)

    def test_runpc(self):
        self.p = rmodel.project.Project()
        p = self.p
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("p.loadProfile", p.loadProfile)
        # p.myProfile.display()

        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metop_2_iasi_pcrttov_compat.H5"
        pcFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/pc/pccoef_metop_2_iasi_landsea_trace_nlte.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        p.myCoeffs.fileName["PC"] = pcFile
        print(p.myCoeffs.fileName)
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        self.assertTrue(p.myCoeffs.hasPC())
        self.assertFalse(p.myOption["ADDPC"].value)
        self.assertFalse(p.isPC())
        self.assertTrue(p.myCoeffs.loadCoeffs)
        self.assertTrue(p.myCoeffs.hasSolar())
        print(p.myCoeffs.fileName)
        p.myOption["ADDPC"].value = True
        p.myOption["OZONE_DATA"].value = False
        p.myOption["CO2_DATA"].value = False
        p.myOption["DO_CHECKINPUT"].value = False
        p.myOption["IPCBND"].value = 1
        p.myOption["IPCREG"].value = 1
        p.myOption.print_value()
        # p.myProfile.display()
        p.myProfile['SKIN']['SURFTYPE'] = 1
        p.ctrlCoherence()
        print(("surftype :", type(p.myProfile['SKIN']['SURFTYPE']),
               p.myProfile['SKIN']['SURFTYPE']))
        print("p.myCoeffs.hasPC()", p.myCoeffs.hasPC())
        print("p.myOption['ADDPC']", p.myOption["ADDPC"].value)
        print("p.loadProfile", p.loadProfile)
        print("p.myCoeffs.loadCoeffs", p.myCoeffs.loadCoeffs)
        self.assertTrue(p.isPC())
        err = p.runPC()
        print("p.myCoeffs.hasNLTE", p.myCoeffs.hasNLTE())
        self.check_option(p)
        self.assertEqual(err, 0)
        err = p.runPCK()
        self.assertEqual(err, 0)
        print("End test_runpc")


if __name__ == "__main__":
    unittest.main()
