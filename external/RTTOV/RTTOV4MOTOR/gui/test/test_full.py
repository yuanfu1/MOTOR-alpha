'''
Created on May 7, 2014

@author: pascale
'''

import unittest
import rttov
import h5py
import rmodel
import logging
from test import rttovgui_unittest_class
import sys
import os


class Test(rttovgui_unittest_class.RttovGuiUnitTest):

    def setUp(self):
        level_logging = logging.DEBUG
        self.p = rmodel.project.Project()

        logging.basicConfig(filename=(self.p.config.ENV['GUI_WRK_DIR'] +
                                      "/rttovgui_unittest_test_full.log"),
                            format=("[%(asctime)s] %(levelname)s [%(module)s:"
                                    "%(funcName)s:%(lineno)d] %(message)s"),
                            level=level_logging,
                            datefmt="%Y:%m:%d %H:%M:%S",
                            filemode="w")
        self.stdCoeffListToTest = self.mk_coeff_list()

    def test_full(self):

        list_std = self.stdCoeffListToTest
        prof_dir = self.p.config.ENV["RTTOV_GUI_PROFILE_DIR"]
        list_profiles = [os.path.join(prof_dir, "standard101lev_allgas.H5"),
                         os.path.join(prof_dir, "cldaer50lev_co2o3.H5")]

        print(list_profiles)
        for coefFile in list_std:

            self.p.myCoeffs.fileName["standard"] = coefFile
            err = self.p.loadCoefficients()
            self.assertEqual(err, 0)
            print(">>>>>>>>>>>>>>>>>>>>>coefFile;", coefFile, " loaded")
            for prof in list_profiles:
                print(">>>>>>>>>>>>>>>>>>>open profile ", prof)
                nb = rttov.profile.getNumberOfProfiles(prof)
                print("nb profiles", nb, "testing the first")
                for n in range(1, 2):
                    print(">>>>>>>>>>>>>>>profile ", prof, "number ", n)
                    self.p.openProfile(prof, n)

                    self.p.ctrlCoherence()
                    self.check_option(self.p)
                    err = self.p.runDirect()
                    self.assertEqual(err, 0)
                    err = self.p.runK()
                    self.assertEqual(err, 0)
                    self.check_option(self.p)
                    print(">>>>>>>>>>>>>>>>>>>>>>OK for ", coefFile, prof, n)

        print(">>>>>>>>>>>>>>>>> end test remove O3")


if __name__ == "__main__":
    unittest.main()
