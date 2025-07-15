
// Example of using the RttovScattSafe class to call RTTOV-SCATT

#include <ProfileScatt.h>
#include <RttovScattSafe.h>
#include <Atlas.h>
#include <vector>

#include "example_data_rttovscatt_cpp.h"


int main(int argc, char **argv) {

    // This example program simulates two profiles
    // The example profile data are defined in example_data_rttovscatt_cpp.h

    // ------------------------------------------------------------------------
    // Set up the profile data
    // ------------------------------------------------------------------------

    int nlevels = 61;
    int nprofiles = 2;

    // Declare the vector of ProfileScatt objects to pass to RttovScattSafe
    std::vector <rttov::ProfileScatt> profiles;

    // Declare the individual ProfileScatt objects and add them to the vector
    for (int p = 0; p < nprofiles; p++) {
        rttov::ProfileScatt myProfile(nlevels);
        profiles.push_back(myProfile);
    }

    // Populate the Profile objects in the vector profiles
    // Note that the simplecloud, icecloud and zeeman data are not mandatory and are omitted here
    int month = 8; // Used when initialising atlases later
    try {
        profiles[0].setAngles(23., 0.) ;
        profiles[1].setAngles(23., 0.);
        profiles[0].setSurfType(1);
        profiles[1].setSurfType(0);
        profiles[0].setS2m(990.138, 279.638, 8507.6938, -2.69568, 3.88115);
        profiles[1].setS2m(989.504, 278.232, 6534.2795, -6.16443, 4.91045);
        profiles[0].setSkin(280.403, 35., 0., 3.0, 5.0, 15.0, 0.1, 0.3);
        profiles[1].setSkin(280.356, 35., 0., 3.0, 5.0, 15.0, 0.1, 0.3);
        profiles[0].setSurfGeom(0., 0., 0.);
        profiles[1].setSurfGeom(10., 20., 0.);
        profiles[0].setDateTimes(2015, month, 1, 0, 0, 0);
        profiles[1].setDateTimes(2015, month, 1, 0, 0, 0);

        std::vector <double> cc_mod;
        for (int l = 0; l < nlevels; l++) {
            // Modified cloud cover to represent a rain fraction below the peak in the rain-producing cloud
            if (l <= 52) {
                cc_mod.push_back(cc_ex[l]);
            } else {
                cc_mod.push_back(cc_ex[52]);
            }
        }

        std::vector <double> ph;
        for (int i = 0; i < nlevels; i++) ph.push_back(ph_ex[i]);
        ph.push_back(0.);
        for (int i = 0; i < nprofiles; i++) {
            profiles[i].setGasUnits(rttov::ppmv_wet);
            profiles[i].setP(p_ex);
            profiles[i].setT(t_ex);
            profiles[i].setQ(q_ex);
            profiles[i].setHydroFrac(cc_mod);
            profiles[i].setClw(clw_ex);
            profiles[i].setCiw(ciw_ex);
            profiles[i].setSnow(snow_ex);
            profiles[i].setRain(rain_ex);

            std::vector <double> s2m = profiles[i].getS2m();
            ph[nlevels] = s2m[0];
            profiles[i].setPh(ph);
        }
    }
    catch (exception& e) {
        std::cerr << "Error defining the profile data " << e.what() << endl;
        exit(1);
    }

    // ------------------------------------------------------------------------
    // Set up RttovScatt instance
    // ------------------------------------------------------------------------

    bool doK = true;

    rttov::RttovScattSafe dprRttov = rttov::RttovScattSafe();

    dprRttov.setFileCoef("../rtcoef_rttov13/rttov13pred54L/rtcoef_gpm_1_dpr.dat");
    dprRttov.setFileHydrotable("../rtcoef_rttov13/hydrotable/hydrotable_gpm_dpr.dat");

    dprRttov.options.setLuserCfrac(false);
    dprRttov.options.setVerboseWrapper(true);
    dprRttov.options.setStoreRad(true);
    dprRttov.setCalcZef(true);

    // Load the instrument
    try {
        dprRttov.loadInst();
    }
    catch (exception& e) {
        std::cerr << "Error loading instrument(s) " << e.what() << std::endl;
        exit(1);
    }

    int nchan_dpr = dprRttov.getNchannels();

    // Associate the profiles: the profiles undergo some
    // checks so use a try block to catch any errors that are thrown up
    try {
      dprRttov.setTheProfiles(profiles);
    }
    catch (exception& e) {
        std::cerr << "Error setting the profiles " << e.what() << endl;
        exit(1);
    }


    // Set up the surface emissivity/reflectance arrays and associate with the RttovScatt object
    double surfemis_dpr[nprofiles][nchan_dpr];
    for (int i = 0; i < nprofiles; i++) {
        for (int c = 0; c < nchan_dpr; c++) surfemis_dpr[i][c] = -1.;
    }
    dprRttov.setSurfEmis((double *)surfemis_dpr);


    // ------------------------------------------------------------------------
    // Call RTTOV-SCATT
    // ------------------------------------------------------------------------

    // Call the RTTOV-SCATT direct model:
    // no arguments are supplied to runDirect so all loaded channels are simulated
    try {
        if ( doK ) {
            dprRttov.runK();
        } else {
            dprRttov.runDirect();
        }
    }
    catch (exception& e) {
        std::cerr << "Error running RTTOV direct model " << e.what() << std::endl;
        exit(1);
    }


    // ------------------------------------------------------------------------
    // Print out some of the output
    // ------------------------------------------------------------------------

    int prof = 0;
    int chan = 0;
    std::cout << "Zef for channel 1, profile 1" << std::endl;
    std::vector <double> zef;
    zef = dprRttov.getZef(prof, chan);
    for (int lev = 0; lev < nlevels; lev++) {
        std::cout << zef[lev] << " ";
    }
    std::cout << std::endl << std::endl;

    std::cout << "AZef for channel 1, profile 1" << std::endl;
    std::vector <double> azef;
    azef = dprRttov.getAZef(prof, chan);
    for (int lev = 0; lev < nlevels; lev++) {
        std::cout << azef[lev] << " ";
    }
    std::cout << std::endl << std::endl;


    // ------------------------------------------------------------------------
    // Deallocate memory
    // ------------------------------------------------------------------------

    // In this example nothing needs to be done: the destructors of the Rttov
    // and Profile classes will deallocate the atlases and other RTTOV data and
    // all locally defined arrays will be deallocated by the garbage collection

}
