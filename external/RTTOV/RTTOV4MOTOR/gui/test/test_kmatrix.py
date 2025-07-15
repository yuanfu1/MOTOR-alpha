'''
Created on Feb 23, 2016

@author: pascale
'''
import unittest

import rmodel
import logging
import rttov
from test import rttovgui_unittest_class


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):
        level_logging = logging.DEBUG
        self.p = rmodel.project.Project()
        logging.basicConfig(filename=(self.p.config.ENV['GUI_WRK_DIR'] +
                                      "/rttovgui_unittest_test_karchive.log"),
                            format=("[%(asctime)s] %(levelname)s [%(module)s:"
                                    "%(funcName)s:%(lineno)d] %(message)s"),
                            level=level_logging,
                            datefmt="%Y:%m:%d %H:%M:%S",
                            filemode="w")

    def testKmatrix1(self):
        self.p.debug = True
        print("Configuration : ", self.p.config.ENV['RTTOV_GUI_PREFIX'])
        print("open profile")

        profileName = self.p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        print(profileName)
        err = self.p.openProfile(profileName)
        self.assertEqual(err, 0)

        coefFile = (self.p.config.ENV['RTTOV_GUI_COEFF_DIR'] +
                    "/rttov7pred54L/rtcoef_noaa_19_hirs.dat")
        print(coefFile)

        self.p.myCoeffs.fileName["standard"] = coefFile
        self.p.loadCoefficients()
        self.p.myOption["IR_SEA_EMIS_MODEL"].value = 1

        err = self.p.runK()
        self.assertEqual(err, 0)

        kmat = rttov.kmatrix.Kmatrix()
        kmat.read(self.p.KMatrixFileName)
        print("profile Q")
        print(kmat.profile["Q"])
        print("profile Q * 0.1")
        print(kmat.profile["Q"] * 0.1)

        print("kprofile channel 1 original")
        prof0 = kmat.getchanprof(0)
        print("Q:")
        print(prof0["Q"])

        print("calcul manuel kprofile chanel 1 * 0.1 profile")
        prof_calcule = prof0["Q"] * kmat.profile["Q"] * 0.1
        print(prof_calcule)

        print("kprofile avec methode kscale 0.1")
        kmat.kscale(0.1)
        prof1 = kmat.getchanprof(0)
        print("Q:")
        print(prof1["Q"])

        print("After reset kprofile channel 1")
        kmat.kscale(1.0)
        prof2 = kmat.getchanprof(0)
        print("Q:")
        print(prof2["Q"])

        print("Calcul Manuel reset")
        prof_reset = prof_calcule / (0.1 * kmat.profile["Q"])
        print("Q:")
        print(prof_reset)

        print("hypothese erreur kmat *2 fois par profile ?")
        prof_erreur = (prof_calcule / 0.1) * kmat.profile["Q"]
        print("Q:")
        print(prof_erreur)


if __name__ == "__main__":

    unittest.main()
