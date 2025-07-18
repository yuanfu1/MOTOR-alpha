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

        logging.basicConfig(filename=self.p.config.ENV['GUI_WRK_DIR'] +
                            "/rttovgui_unittest_test_retrieve_iasi.log",
                            format=("[%(asctime)s] %(levelname)s [%(module)s:"
                                    "%(funcName)s:%(lineno)d] %(message)s"),
                            level=level_logging,
                            datefmt="%Y:%m:%d %H:%M:%S",
                            filemode="w")


#    def tearDown(self):
#        pass

    def test_runRetrieve_IASI(self):
        channel_list = [16, 38, 49, 51, 55, 57, 59, 61, 63, 66, 70, 72, 74, 79,
                        81, 83, 85, 87, 89, 92, 95, 97, 99, 101, 104, 106, 109,
                        111, 113, 116, 119, 122, 125, 128, 131, 133, 135, 138,
                        141, 144, 146, 148, 151, 154, 157, 159, 161, 163, 167,
                        170,
                        173, 176, 179, 180, 185, 187, 193, 199, 205, 207, 210,
                        212, 214, 217, 219, 222, 224, 226, 230, 232, 236,
                        239, 242,
                        243, 246, 249, 252, 254, 260, 262, 265, 267, 269, 275,
                        280, 282, 294, 296, 299, 303, 306, 323, 327, 329, 335,
                        345,
                        347, 350, 354, 356, 360, 366, 371, 373, 375, 377, 379,
                        381, 383, 386, 389, 398, 401, 404, 407, 410, 414,
                        416, 426,
                        428, 432, 434, 439, 445, 457, 515, 546, 552, 559,
                        566, 571, 573, 646, 662, 668, 756, 867, 906, 921,
                        1027, 1046,
                        1121, 1133, 1191, 1194, 1271, 1479, 1509, 1513, 1521,
                        1536, 1574, 1579, 1585, 1587, 1626, 1639, 1643, 1652,
                        1658,
                        1671, 1786, 1805, 1884, 1991, 2019, 2094, 2119, 2213,
                        2239, 2245, 2271, 2321, 2398, 2701, 2741, 2819, 2889,
                        2907,
                        2910, 2919, 2939, 2944, 2948, 2951, 2958, 2977, 2985,
                        2988, 2991, 2993, 3002, 3008, 3014, 3027, 3029, 3036,
                        3047,
                        3049, 3053, 3058, 3064, 3069, 3087, 3093, 3098, 3105,
                        3107, 3110, 3127, 3136, 3151, 3160, 3165, 3168, 3175,
                        3178,
                        3207, 3228, 3244, 3248, 3252, 3256, 3263, 3281, 3303,
                        3309, 3312, 3322, 3339, 3375, 3378, 3411, 3438, 3440,
                        3442,
                        3444, 3446, 3448, 3450, 3452, 3454, 3458, 3467, 3476,
                        3484, 3491, 3497, 3499, 3504, 3506, 3509, 3518, 3522,
                        3527,
                        3540, 3555, 3575, 3577, 3580, 3582, 3586, 3589, 3599,
                        3653, 3658, 3661, 3943, 4032, 5130, 5368, 5371, 5379,
                        5381,
                        5383, 5397, 5399, 5401, 5403, 5405, 5455, 5480, 5483,
                        5485, 5492, 5502, 5507, 5509, 5517, 5558, 5988, 5992,
                        5994,
                        6003, 6350, 6463, 6601, 6962, 6980, 6982, 6985, 6987,
                        6989, 6991, 6993, 6995, 6997, 7267, 7269, 7424, 7426,
                        7428, 7885, 8007]
        print(">>>>>>>>>> testrunRetrieve_IASI ")
        pTrue = rmodel.project.Project()
        pBg = rmodel.project.Project()
        pTrue.setFileNameMark("TrueTestIASI")
        pBg.setFileNameMark("BgTestIASI")

        print(("pBg.profileFileName", pBg.profileFileName))
        filename = pTrue.config.ENV[
            "RTTOV_GUI_PROFILE_DIR"] + "/standard54lev_allgas.H5"
        pTrue.openProfile(filename, 1)
        pBg.openProfile(filename, 1)
        pBg.myProfile["T"] = pBg.myProfile["T"] + 4
        coefFile = pTrue.config.ENV[
            "RTTOV_GUI_COEFF_DIR"] + "/rttov7pred54L/rtcoef_metop_2_iasi.H5"
        pBg.myCoeffs.fileName["standard"] = coefFile
        pBg.loadCoefficients(channel_list)

        my_chan_list = pBg.myCoeffs.getFF_ORI_CHN()
        print("################################""getFF_ORI_CHN,", my_chan_list)

        (inst, sat, id_sat) = pBg.myCoeffs.get_instrument_and_platform_name()
        print(">>>>>>", (sat, inst, id_sat))
        self.assertEqual(sat, "metop")
        self.assertEqual(inst, "iasi")
        self.assertEqual(id_sat, 2)

        pTrue.myProfile["LATITUDE"] = -90.
        pTrue.myProfile["LONGITUDE"] = 0.
        pTrue.myProfile["ELEVATION"] = 0.0
        pTrue.myProfile["SKIN"]["SURFTYPE"] = 1
        pTrue.myProfile["SUNZENANGLE"] = 0.00
        pTrue.myProfile["ZENANGLE"] = 0.00
        pBg.myCoeffs.fileName["standard"] = coefFile
        pBg.myProfile["LATITUDE"] = -90.
        pBg.myProfile["LONGITUDE"] = 0.
        pBg.myProfile["ELEVATION"] = 0.0
        pBg.myProfile["SKIN"]["SURFTYPE"] = 1
        pBg.myProfile["SUNZENANGLE"] = 0.00
        pBg.myProfile["ZENANGLE"] = 0.00
        Bmatrix = os.environ[
            "RTTOV_GUI_PREFIX"] + "/r1Dvar/data/Sample_Bmatrices/Bmatrix_54L"
        retrieveProject = r1Dvar.r1dvar.Project1dvar(
            pTrue, pBg, matrixBfile=Bmatrix, satellite="metop-2",
            instrument="iasi", channel_list=channel_list)

        retrieveProject.setFactorB(1)
        retrieveProject.setFactorR(1)
        Xr1 = retrieveProject.retrieve1d()
        retrieveProject.plot()
        Xr2 = retrieveProject.stepRetrieve1d()
        Xr3 = retrieveProject.stepRetrieve1d()
        Xr4 = retrieveProject.stepRetrieve1d()
        retrieveProject.plot()


if __name__ == "__main__":
    unittest.main()
