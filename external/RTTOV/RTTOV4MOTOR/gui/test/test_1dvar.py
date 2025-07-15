'''
Created on May 7, 2014

@author: pascale
'''
import unittest


import rmodel
import logging

from test import rttovgui_unittest_class

import r1Dvar
import os

prefix = os.environ["RTTOV_GUI_PREFIX"]


class Test(rttovgui_unittest_class.RttovGuiUnitTest):
    def setUp(self):
        level_logging = logging.DEBUG
        self.p = rmodel.project.Project()

        logging.basicConfig(filename=(self.p.config.ENV['GUI_WRK_DIR'] +
                                      "/rttovgui_unittest_test_1dvar.log"),
                            format="[%(asctime)s] %(levelname)s "
                            "[%(module)s:%(funcName)s:%(lineno)d] %(message)s",
                            level=level_logging,
                            datefmt="%Y:%m:%d %H:%M:%S",
                            filemode="w")

    def test_1dvar_backround(self):
        fname = prefix + "/r1Dvar/data/Sample_Background/BACKGROUND_43L.dat"
        PBg = r1Dvar.r1dvar.Background()
        PBg.read(fname)
        PBg.print_data()
        prof = PBg.toProfile()
        nlev = prof["NLEVELS"]
        self.assertGreater(prof["P"][10], prof["P"][0], "comparaison P")
        self.assertGreater(prof["Q"][nlev - 1], prof["Q"][0], "comparaison Q")

        PBg.print_data()
        print("-> profile : ")
        prof.display()

    def test_2profile(self):
        print(">>> test_2profile")
        p = rmodel.project.Project()
        filename = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        p.openProfile(filename, 1)
        nlevels = p.myProfile["NLEVELS"]
        pbis = rmodel.project.Project()
        pbis.openProfile(filename, 2)
        vector = p.myProfile.to1DvarVect()
        newprofile = r1Dvar.r1dvar.vector2profile(vector, pbis.myProfile)
        vector2 = newprofile.to1DvarVect()
        print("vector")
        print(vector)
        print("vector2")
        print(vector2)
        print("vector-vector2")
        print(vector - vector2)
        for i in range(len(vector)):
            self.assertEqual(vector[i], vector2[i])
        # verify vector2 and vector are not the same memory
        vector[:] = 1
        for i in range(nlevels):
            self.assertNotEqual(newprofile["T"][i], vector[i])
        for i in range(nlevels - 29, nlevels):
            self.assertNotEqual(newprofile["Q"][i], vector[
                                nlevels + i - nlevels + 29])
        self.assertNotEqual(newprofile["S2M"]["T"], vector[nlevels + 29])
        self.assertNotEqual(newprofile["S2M"]["Q"], vector[nlevels + 29 + 1])
        self.assertNotEqual(newprofile["SKIN"]["T"], vector[nlevels + 29 + 2])
        print("newprofile[Q]")
        print(newprofile["Q"])

        print((newprofile["S2M"]["T"], vector[nlevels + 29]))
        # verify vector2 and pbis are not the same memory
        pbis.myProfile["T"][:] = 2
        pbis.myProfile["Q"][:] = 0
        pbis.myProfile["S2M"]["T"] = 0
        pbis.myProfile["S2M"]["Q"] = 0
        pbis.myProfile["SKIN"]["T"] = 0
        for i in range(nlevels):
            self.assertNotEqual(newprofile["T"][i], pbis.myProfile["T"][i])
        for i in range(nlevels):
            self.assertNotEqual(newprofile["Q"][i], pbis.myProfile["Q"][i])
        self.assertNotEqual(newprofile["S2M"]["T"], pbis.myProfile["S2M"]["T"])
        self.assertNotEqual(newprofile["S2M"]["Q"], pbis.myProfile["S2M"]["Q"])
        self.assertNotEqual(newprofile["SKIN"][
                            "T"], pbis.myProfile["SKIN"]["T"])

    def test_2vector(self):
        print(">>> test_2vector")
        print(" a vector contains only T values, lnq bottom lnq vales ,"
              " Tsurf, lnq surf and Tskin")
        print(" temperature, and ln q are stored"
              " from top of atmosphere to ground ")
        p = rmodel.project.Project()
        filename = p.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        p.openProfile(filename, 1)
        vector = p.myProfile.to1DvarVect()
        print("T:")
        print((p.myProfile["T"]))
        print("Q:")
        print((p.myProfile["Q"]))
        print((p.myProfile["S2M"]["T"]))
        print((p.myProfile["S2M"]["Q"]))
        print((p.myProfile["SKIN"]["T"]))
        print('vector')
        print(vector)
        nlevels = p.myProfile["NLEVELS"]
        self.assertEqual(vector[0], p.myProfile["T"][0])
        self.assertEqual(vector[nlevels - 1], p.myProfile["T"][nlevels - 1])
        self.assertEqual(vector[nlevels + 29 - 1],
                         p.myProfile["Q"][nlevels - 1])
        self.assertEqual(vector[nlevels], p.myProfile["Q"][nlevels - 29])
        # verify p.myProfile and vector are not linked
        p.myProfile["T"][0] = -1
        p.myProfile["Q"][nlevels - 1] = -1
        self.assertNotEqual(vector[0], p.myProfile["T"][0])
        self.assertNotEqual(vector[nlevels + 29 - 1],
                            p.myProfile["Q"][nlevels - 1])

    def test_1dvar_Rmatrices(self):
        print(">>> test_1dvar_Rmatrices")
        # test Rmatrices
        rmatrix = prefix + "/r1Dvar/data/IASI_COEFFS_DIR/Rmatrix_orig"

        print("rmatrix")
        R = r1Dvar.r1dvar.Rmatrix()
        R.read(rmatrix)

    def test_1dvar_Bmatrices(self):
        print(">>> test_1dvar_Bmatrices")
        filename = prefix + "/r1Dvar/data/Sample_Bmatrices/Bmatrix_43L"
        B = r1Dvar.r1dvar.Bmatrix()
        B.read_matrices(filename)


if __name__ == "__main__":
    unittest.main()
