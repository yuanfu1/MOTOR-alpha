import unittest
import rmodel
import logging
from test import rttovgui_unittest_class
import r1Dvar
import os


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):
        level_logging = logging.DEBUG
        self.p = rmodel.project.Project()

        logging.basicConfig(filename=(self.p.config.ENV['GUI_WRK_DIR'] +
                                      "/rttovgui_unittest_test_retrieve.log"),
                            format=("[%(asctime)s] %(levelname)s [%(module)s:"
                                    "%(funcName)s:%(lineno)d] %(message)s"),
                            level=level_logging,
                            datefmt="%Y:%m:%d %H:%M:%S",
                            filemode="w")

    def test_runRetrieve_HIRS(self):
        print(">>>>>>>>>> testrunRetrieve_HIRS ")
        pTrue = rmodel.project.Project()
        pBg = rmodel.project.Project()
        pTrue.setFileNameMark("TrueTestHIRS")
        pBg.setFileNameMark("BgTestHIRS")

        print(("pBg.profileFileName", pBg.profileFileName))
        filename = pTrue.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        pTrue.openProfile(filename, 1)
        pTrue.myProfile["LATITUDE"] = -59.796
        pBg.openProfile(filename, 6)
        coefFile = pTrue.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + "/rttov7pred54L/rtcoef_noaa_19_hirs.dat"
        pBg.myCoeffs.fileName["standard"] = coefFile
        pBg.loadCoefficients()
        (inst, sat, id_sat) = pBg.myCoeffs.get_instrument_and_platform_name()
        print(">>>>>>", (sat, inst, id_sat))
        self.assertEqual(sat, "noaa")
        self.assertEqual(inst, "hirs")
        self.assertEqual(id_sat, 19)
        print("id_sat=", id_sat)
        Bmatrix = os.environ[
            "RTTOV_GUI_PREFIX"] + "/r1Dvar/data/Sample_Bmatrices/Bmatrix_54L"
        pTrue.myOption["IR_SEA_EMIS_MODEL"].value = 1
        pBg.myOption["IR_SEA_EMIS_MODEL"].value = 1
        retrieveProject = r1Dvar.r1dvar.Project1dvar(
            pTrue, pBg, matrixBfile=Bmatrix,
            satellite="noaa-19", instrument="hirs")
        retrieveProject.setFactorB(1)
        retrieveProject.setFactorR(1)
        Xr1 = retrieveProject.retrieve1d()
        self.assertEqual(
            pTrue.myProfile["LATITUDE"], pBg.myProfile["LATITUDE"])
        self.assertEqual(pTrue.myProfile["CFRACTION"], 0)
        self.assertEqual(pTrue.myOption["OZONE_DATA"].value, 0)
        # retrieveProject.plot()
        Xr2 = retrieveProject.stepRetrieve1d()
        # retrieveProject.plot()
        Xr3 = retrieveProject.stepRetrieve1d()
        # retrieveProject.plot()
        Xr4 = retrieveProject.stepRetrieve1d()
        retrieveProject.plot()

    def test_runRetrieve_AMSUA(self):
        print(">>>>>>>>>> testrunRetrieve_AMSUA ")

        pTrue = rmodel.project.Project()
        pBg = rmodel.project.Project()
        pTrue.setFileNameMark("TrueTestAMSUA")
        pBg.setFileNameMark("BgTestAMSUA")
        print(("pBg.profileFileName", pBg.profileFileName))
        filename = pTrue.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        pTrue.openProfile(filename, 1)
        pBg.openProfile(filename, 6)
        coefFile = pTrue.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + "/rttov7pred54L/rtcoef_noaa_19_amsua.dat"
        pBg.myCoeffs.fileName["standard"] = coefFile
        pTrue.myOption["IR_SEA_EMIS_MODEL"].value = 1
        pBg.myOption["IR_SEA_EMIS_MODEL"].value = 1
        print(">>>>>>>>>>>> create the Project1dvar")
        retrieveProject = r1Dvar.r1dvar.Project1dvar(
            pTrue, pBg, satellite="noaa-19", instrument="amsua")
        retrieveProject.setFactorB(1)
        retrieveProject.setFactorR(1)
        Xr1 = retrieveProject.retrieve1d()

        # retrieveProject.plot()
        Xr2 = retrieveProject.stepRetrieve1d()
        # retrieveProject.plot()
        Xr3 = retrieveProject.stepRetrieve1d()
        # retrieveProject.plot()
        Xr4 = retrieveProject.stepRetrieve1d()
        print(">>>>> test pBg.radianceFileName change", pBg.radianceFileName)

        retrieveProject.plot()


if __name__ == "__main__":
    unittest.main()
