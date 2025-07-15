'''
Created on Mar 3, 2014

@author: pascale
'''
from rmodel import project
import unittest
from test import rttovgui_unittest_class


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def ctrl(self, p):
        import numpy
        print((p.myProfile['LATITUDE'], p.myProfile['LONGITUDE'],
               p.myProfile["SKIN"]["SURFTYPE"]))
        print("firstTimeSurface", p.firstTimeSurface)
        print("use atlas", p.useAtlas)
        print("###########call controleSurface !!!###########")
        p.controleSurface()
        print("######################")
        print(p.myEmissivity["EMIS_IN"])
        print(p.myCoeffs.fileName["standard"])
        print(type(p.myEmissivity["EMIS_IN"]))
        self.assertTrue(numpy.all(p.myEmissivity["EMIS_IN"] <= 1))
        self.assertTrue(numpy.all(p.myEmissivity["EMIS_IN"] >= 0))

    def test_atlas_landsat(self):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>test atlas landsat")
        p = project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("---------------------TEST ATLAS ---------------------------")
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov7pred54L/rtcoef_landsat_4_tm.dat"

        p.myCoeffs.fileName["standard"] = coefFile

        print(p.myCoeffs.fileName)

        err = p.loadCoefficients()
        self.assertEqual(err, 0)

        p.myProfile['LATITUDE'] = 45
        p.myProfile['LONGITUDE'] = 10
        p.myProfile["SKIN"]["SURFTYPE"] = 0

        err = p.openAtlas()
        self.assertEqual(err, 0)
        print("****************")
        print(type(p.myEmissivity["EMIS_IN"]))
        print("*************")
        import numpy
        self.assertTrue(numpy.any(p.myEmissivity["EMIS_IN"] > 0))
        self.ctrl(p)

    def test_atlas_1(self):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>test atlas 1")
        p = project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("---------------------TEST ATLAS 1-------------------------")
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred54L/rtcoef_eos_2_modis.dat"

        p.myCoeffs.fileName["standard"] = coefFile

        print(p.myCoeffs.fileName)

        err = p.loadCoefficients()
        self.assertEqual(err, 0)

        p.myProfile['LATITUDE'] = 45
        p.myProfile['LONGITUDE'] = 10
        p.myProfile["SKIN"]["SURFTYPE"] = 0

        err = p.openAtlas()
        self.assertEqual(err, 0)
        self.ctrl(p)
        import numpy
        self.assertTrue(numpy.any(p.myEmissivity["EMIS_IN"] > 0))

    def test_atlas_2(self):
        print(">>>>>>>>>>>>>test atlas 2")
        print("test_atlas_2 surftype=1")

        p = project.Project()
        profileName = p.config.ENV["RTTOV_GUI_PROFILE_DIR"] + \
            "/cldaer101lev_allgas.H5"

        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("---------------------TEST ATLAS 2 --------------------------")
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred54L/rtcoef_eos_2_modis.dat"

        p.myCoeffs.fileName["standard"] = coefFile

        print(p.myCoeffs.fileName)

        err = p.loadCoefficients()
        self.assertEqual(err, 0)

        p.myProfile['LATITUDE'] = 45
        p.myProfile['LONGITUDE'] = 10
        p.myProfile["SKIN"]["SURFTYPE"] = 1
        print((p.myProfile['LATITUDE'], p.myProfile['LONGITUDE'],
               p.myProfile["SKIN"]["SURFTYPE"]))

        err = p.openAtlas()
        self.assertEqual(err, 0)
        self.assertFalse(p.firstTimeSurface)
        self.ctrl(p)

    def test_atlas_3(self):
        print(">>>>>>>>>test_atlas_3")
        p = project.Project()
        profileName = p.config.ENV["RTTOV_GUI_PROFILE_DIR"] + \
            "/cldaer101lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("---------------------TEST ATLAS ---------------------------")
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred54L/rtcoef_eos_2_modis.dat"

        p.myCoeffs.fileName["standard"] = coefFile

        print(p.myCoeffs.fileName)

        err = p.loadCoefficients()
        self.assertEqual(err, 0)

        p.myProfile['LATITUDE'] = 45
        p.myProfile['LONGITUDE'] = 10
        p.myProfile["SKIN"]["SURFTYPE"] = 1

        err = p.openAtlas()
        self.assertEqual(err, 0)
        self.assertFalse(p.firstTimeSurface)
        self.ctrl(p)

    def test_atlas_4(self):
        print(">>>>>>>>>test_atlas_4")
        p = project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("---------------------TEST ATLAS 4 ----------------------------")
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred54L/rtcoef_eos_2_modis.dat"

        p.myCoeffs.fileName["standard"] = coefFile

        print(p.myCoeffs.fileName)

        err = p.loadCoefficients()
        self.assertEqual(err, 0)

        p.myProfile['LATITUDE'] = 45
        p.myProfile['LONGITUDE'] = 10
        p.myProfile["SKIN"]["SURFTYPE"] = 0

        err = p.openAtlas()
        self.assertEqual(err, 0)
        self.assertFalse(p.firstTimeSurface)
        self.ctrl(p)

    def test_atlas_5(self):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>test_atlas_5")
        p = project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("---------------------TEST ATLAS 5 ----------------------------")
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred54L/rtcoef_mtg_1_fci.dat"

        p.myCoeffs.fileName["standard"] = coefFile

        print(p.myCoeffs.fileName)

        err = p.loadCoefficients()
        self.assertEqual(err, 0)

        p.myProfile['LATITUDE'] = 45
        p.myProfile['LONGITUDE'] = 10
        p.myProfile["SKIN"]["SURFTYPE"] = 1
        err = p.openAtlas()
        self.assertEqual(err, 0)

        err = p.dropCoefficients()
        self.assertEqual(err, 0)

        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred54L/rtcoef_eos_1_modis.dat"
        p.myCoeffs.fileName["standard"] = coefFile
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        err = p.openAtlas()
        self.assertEqual(err, 0)
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred101L/rtcoef_metop_2_iasi.H5"
        p.myCoeffs.fileName["standard"] = coefFile
        err = p.loadCoefficients()
        self.assertEqual(err, 0)
        print("has MW ?")
        self.assertFalse(p.myCoeffs.isMW())
        import rttov
        wavenumbers = rttov.getcoefval.rttov_get_coef_val_r1('WAVENUMBERS')
        print(min(wavenumbers), max(wavenumbers))
        self.assertEqual(min(wavenumbers), 645.0)
        self.assertEqual(max(wavenumbers), 2760.0)
        id_sensor = rttov.getcoefval.rttov_get_coef_val_i0('ID_SENSOR')
        print("id_sensor :", rttov.getcoefval.rttov_get_coef_val_i0(
            'ID_SENSOR'))
        self.assertEqual(id_sensor, 3)
        p.ctrlCoherence()
        self.check_option(p)
        err = p.runDirect()
        self.assertEqual(err, 0)

    def test_atlas_6(self):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>test_atlas_6")
        p = project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("---------------------TEST ATLAS 6----------------------------")
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov9pred54L/rtcoef_mtg_1_fci.dat"

        p.myCoeffs.fileName["standard"] = coefFile

        print(p.myCoeffs.fileName)

        err = p.loadCoefficients()
        self.assertEqual(err, 0)

        p.myProfile['LATITUDE'] = 45
        p.myProfile['LONGITUDE'] = 10
        p.myProfile["SKIN"]["SURFTYPE"] = 1
        print(p.myProfile['DATE'])
        err = p.openAtlas(2, 1)
        self.assertEqual(err, 0)
        err = p.openAtlas(3, 1)
        self.assertEqual(err, 0)

    def test_atlas_7(self):
        print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>test_atlas_7")
        p = project.Project()
        profileName = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/cldaer101lev_allgas.H5"
        print(profileName)
        err = p.openProfile(profileName)
        self.assertEqual(err, 0)
        print("---------------------TEST ATLAS 7----------------------------")
        coefFile = p.config.ENV['RTTOV_GUI_COEFF_DIR'] + \
            "/rttov7pred54L/rtcoef_metop_1_mhs.dat"

        p.myCoeffs.fileName["standard"] = coefFile

        print(p.myCoeffs.fileName)

        err = p.loadCoefficients()
        self.assertEqual(err, 0)

        p.myProfile['LATITUDE'] = 45
        p.myProfile['LONGITUDE'] = 10
        p.myProfile["SKIN"]["SURFTYPE"] = 1
        print(p.myProfile['DATE'])
        err = p.openAtlas(1, 1)
        self.assertEqual(err, 0)
        err = p.openAtlas(2, 1)
        self.assertEqual(err, 0)


if __name__ == "__main__":
    unittest.main()
