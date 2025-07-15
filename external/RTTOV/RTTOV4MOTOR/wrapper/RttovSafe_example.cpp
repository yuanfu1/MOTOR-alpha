
// Example of using the RttovSafe class to call RTTOV for multiple instruments
// with the emissivity and BRDF atlases

// Three RttovSafe instances are created representing three instruments

#include <Profile.h>
#include <RttovSafe.h>
#include <Atlas.h>
#include <vector>

#include "example_data_cpp.h"


int main(int argc, char **argv) {

    // This example program simulates two profiles for each of three instruments
    // The example profile data are defined in example_data_cpp.h

    // ------------------------------------------------------------------------
    // Set up the profile data
    // ------------------------------------------------------------------------

    int nlevels = 101;
    int nprofiles = 2;

    // Declare the vector of Profile objects to pass to RttovSafe
    std::vector <rttov::Profile> profiles;

    // Declare the individual Profile objects and add them to the vector
    for (int p = 0; p < nprofiles; p++) {
        rttov::Profile myProfile(nlevels);
        profiles.push_back(myProfile);
    }

    std::vector <double> t[nprofiles];

    // example_data.h contains p, T, q, co2 for a single profile: copy the temperature
    // profile and modify it for the second profile
    for (int i = 0; i < nprofiles; i++) {
        for (int l = 0; l < nlevels; l++) {
            t[i].push_back(t_ex[l] + i * 2);
        }
    }

    // Populate the Profile objects in the vector profiles
    // Note that the simplecloud, clwscheme, icecloud and zeeman data are not mandatory and are omitted here
    int month = 8; // Used when initialising atlases later
    try {
        profiles[0].setAngles(0.,  0., 45., 180.) ;
        profiles[1].setAngles(60., 0., 45., 180.);
        profiles[0].setSurfType(0, 0);
        profiles[1].setSurfType(1, 0);
        profiles[0].setS2m(1013., 0.263178E+03, 0.236131E+04, 4., 2., 100000.);
        profiles[1].setS2m(1013., 0.265178E+03, 0.236131E+04, 4., 2., 100000.);

        for (int i = 0; i < nprofiles; i++) {
            profiles[i].setGasUnits(rttov::ppmv_wet);
            profiles[i].setP(p_ex);
            profiles[i].setT(t[i]);
            profiles[i].setQ(q_ex);
            profiles[i].setCO2(co2_ex);
            profiles[i].setSkin(270., 35., 0., 0., 3.0, 5.0, 15.0, 0.1, 0.3);
            profiles[i].setSurfGeom(10., 20., 0.);
            profiles[i].setDateTimes(2015, month, 1, 0, 0, 0);
        }
    }
    catch (exception& e) {
        std::cerr << "Error defining the profile data " << e.what() << endl;
        exit(1);
    }


    // ------------------------------------------------------------------------
    // Set up Rttov instances for each instrument
    // ------------------------------------------------------------------------

    // Create three RttovSafe objects for three instruments
    rttov::RttovSafe seviriRttov = rttov::RttovSafe();
    rttov::RttovSafe hirsRttov   = rttov::RttovSafe();
    rttov::RttovSafe mhsRttov    = rttov::RttovSafe();

    // For HIRS and ATMS we will read all channels, but we will read a subset for SEVIRI
    int nchan_seviri = 10;
    int nchan_hirs   = 19;
    int nchan_mhs    = 5;

    // For SEVIRI exclude ozone and hi-res vis channels (9 and 12) in this example
    std::vector <int> chan_list_seviri(nchan_seviri);
    chan_list_seviri = {1, 2, 3, 4, 5, 6, 7, 9, 10, 11};


    // Set the options for each RttovSafe instance:
    // - the path to the coefficient file must always be specified
    // - turn RTTOV interpolation on (because input pressure levels differ from coefficient file levels)
    // - set the verbose_wrapper flag to true so the wrapper provides more information
    // - enable solar simulations for SEVIRI
    // - enable CO2 simulations for HIRS (the CO2 profiles are ignored for the SEVIRI and MHS simulations)
    // - enable the store_trans wrapper option for MHS to provide access to RTTOV transmission structure

    seviriRttov.setFileCoef("../rtcoef_rttov13/rttov13pred54L/rtcoef_msg_3_seviri_o3.dat");
    seviriRttov.options.setAddInterp(true);
    seviriRttov.options.setAddSolar(true);
    seviriRttov.options.setVerboseWrapper(true);

    hirsRttov.setFileCoef("../rtcoef_rttov13/rttov13pred54L/rtcoef_metop_1_hirs_o3co2.dat");
    hirsRttov.options.setAddInterp(true);
    hirsRttov.options.setCO2Data(true);
    hirsRttov.options.setVerboseWrapper(true);

    mhsRttov.setFileCoef("../rtcoef_rttov13/rttov13pred54L/rtcoef_noaa_19_mhs.dat");
    mhsRttov.options.setAddInterp(true);
    mhsRttov.options.setStoreTrans(true);
    mhsRttov.options.setVerboseWrapper(true);


    // Load the instruments: for HIRS and MHS do not supply a channel list and so read all channels
    try {
        seviriRttov.loadInst(chan_list_seviri);
        hirsRttov.loadInst();
        mhsRttov.loadInst();
    }
    catch (exception& e) {
        std::cerr << "Error loading instrument(s) " << e.what() << std::endl;
        exit(1);
    }

    // Associate the profiles with each RttovSafe instance: the profiles undergo some
    // checks so use a try block to catch any errors that are thrown up
    try {
        seviriRttov.setTheProfiles(profiles);
        hirsRttov.setTheProfiles(profiles);
        mhsRttov.setTheProfiles(profiles);
    }
    catch (exception& e) {
        std::cerr << "Error setting the profiles " << e.what() << endl;
        exit(1);
    }


    // ------------------------------------------------------------------------
    // Load the emissivity and BRDF atlases
    // ------------------------------------------------------------------------

    // Load the emissivity and BRDF atlases:
    // - load data for the month in the profile data
    // - load the IR emissivity atlas data for multiple instruments so it can be used for SEVIRI and HIRS
    // - SEVIRI is the only VIS/NIR instrument we can use the single-instrument initialisation for the BRDF atlas

    rttov::Atlas irAtlas = rttov::Atlas();
    irAtlas.setAtlasPath("../emis_data");
    irAtlas.loadIrEmisAtlas(month, true);  // Include angular correction, but do not initialise for single-instrument

    rttov::Atlas brdfAtlas = rttov::Atlas();
    brdfAtlas.setAtlasPath("../brdf_data");
    brdfAtlas.loadBrdfAtlas(month, &seviriRttov); // Pass Rttov object to enable single-instrument initialisation
    brdfAtlas.setIncSea(false);                   // Do not use BRDF atlas for sea surface types

    rttov::Atlas mwAtlas = rttov::Atlas();
    mwAtlas.setAtlasPath("../emis_data");
    mwAtlas.loadMwEmisAtlas(month);

    // Set up the surface emissivity/reflectance arrays and associate with the Rttov objects
    double surfemisrefl_seviri[4][nprofiles][nchan_seviri];
    double surfemisrefl_hirs[4][nprofiles][nchan_hirs];
    double surfemisrefl_mhs[4][nprofiles][nchan_mhs];

    seviriRttov.setSurfEmisRefl((double *)surfemisrefl_seviri);
    hirsRttov.setSurfEmisRefl((double *)surfemisrefl_hirs);
    mhsRttov.setSurfEmisRefl((double *)surfemisrefl_mhs);


    // ------------------------------------------------------------------------
    // Call RTTOV
    // ------------------------------------------------------------------------

    // Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
    // Negative values will cause RTTOV to supply emissivity/BRDF values (i.e. equivalent to
    // calcemis/calcrefl TRUE - see RTTOV user guide)

    for (int i = 0; i < nprofiles; i++) {
        for (int j = 0; j < 4; j++) {
            for (int c = 0; c < nchan_seviri; c++) surfemisrefl_seviri[j][i][c] = -1.;
            for (int c = 0; c < nchan_hirs; c++)   surfemisrefl_hirs[j][i][c]   = -1.;
            for (int c = 0; c < nchan_mhs; c++)    surfemisrefl_mhs[j][i][c]    = -1.;
        }
    }

    // Call emissivity and BRDF atlases
    try {
        // Do not supply a channel list for SEVIRI: this returns emissivity/BRDF values for all
        // *loaded* channels which is what is required
        irAtlas.fillEmisBrdf((double *)surfemisrefl_seviri[0], &seviriRttov);
        brdfAtlas.fillEmisBrdf((double *)surfemisrefl_seviri[1], &seviriRttov);
        irAtlas.fillEmisBrdf((double *)surfemisrefl_hirs[0], &hirsRttov);
        mwAtlas.fillEmisBrdf((double *)surfemisrefl_mhs[0], &mhsRttov);
    }
    catch (exception& e) {
        // If there was an error the emissivities/BRDFs will not have been modified so it
        // is OK to continue and call RTTOV with calcemis/calcrefl set to TRUE everywhere
        std::cerr << "Error calling atlas " << e.what() << std::endl;
    }


    // Call the RTTOV direct model for each instrument:
    // no arguments are supplied to runDirect so all loaded channels are simulated
    try {
        seviriRttov.runDirect();
        hirsRttov.runDirect();
        mhsRttov.runDirect();
    }
    catch (exception& e) {
        std::cerr << "Error running RTTOV direct model " << e.what() << std::endl;
        exit(1);
    }


    // ------------------------------------------------------------------------
    // Print out some of the output
    // ------------------------------------------------------------------------
    std::cout << std::endl;
    std::cout << "SELECTED OUTPUT" << std::endl;
    std::cout << std::endl;

    std::cout << "SEVIRI visible channel reflectances, channels 1-3" << std::endl;
    std::vector <double> btrefl_seviri;
    for (int p = 0; p < nprofiles; p++) {
        btrefl_seviri = seviriRttov.getBtRefl(p);
        std::cout << "Profile " << p << " : ";
        for (int c = 0; c < 3; c++) std::cout << btrefl_seviri[c] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "HIRS radiances" << std::endl;
    std::vector <double> rads_hirs;
    for (int p = 0; p < nprofiles; p++) {
        rads_hirs = hirsRttov.getRads(p);
        std::cout << "Profile " << p << " : ";
        for (int c = 0; c < nchan_hirs; c++) std::cout << rads_hirs[c] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // We can access the RTTOV transmission structure because the store_trans option was set above for mhsRttov
    std::cout << "MHS total transmittance" << std::endl;
    std::vector <double> tautotal_mhs;
    for (int p = 0; p < nprofiles; p++) {
        tautotal_mhs = mhsRttov.getTauTotal(p);
        std::cout << "Profile " << p << " : ";
        for (int c = 0; c < nchan_mhs; c++) std::cout << tautotal_mhs[c] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;


    // ------------------------------------------------------------------------
    // Deallocate memory
    // ------------------------------------------------------------------------

    // In this example nothing needs to be done: the destructors of the Rttov
    // and Profile classes will deallocate the atlases and other RTTOV data and
    // all locally defined arrays will be deallocated by the garbage collection

}
