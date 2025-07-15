
// Example of using the Rttov class to call RTTOV for multiple instruments
// with the emissivity and BRDF atlases

// Three Rttov instances are created representing three instruments

#include <Profiles.h>
#include <Rttov.h>
#include <Atlas.h>
#include <vector>

#include "example_data.h"


int main(int argc, char **argv) {

    // This example program simulates two profiles for each of three instruments
    // The example profile data are defined in example_data.h

    // ------------------------------------------------------------------------
    // Set up the profile data
    // ------------------------------------------------------------------------

    // Declare an instance of Profiles
    int nlevels = 101;
    int nprofiles = 2;
    rttov::Profiles myProfiles = rttov::Profiles(nprofiles, nlevels);

    // Populate myProfiles with data
    double p[nprofiles][nlevels];    // Input pressure profiles
    double t[nprofiles][nlevels];    // Input temperature profiles
    double q[nprofiles][nlevels];    // Input water vapour profiles
    double co2[nprofiles][nlevels];  // Input CO2 profiles

    // example_data.h contains p, T, q, co2 for a single profile
    for (int i = 0; i < nprofiles; i++) {
        for (int l = 0; l < nlevels; l++) {
            p[i][l] = p_ex[l];
            t[i][l] = t_ex[l] + i * 2;    // Modify the second temperature profile
            q[i][l] = q_ex[l];
            co2[i][l] = co2_ex[l];
        }
    }

    // Associate the profiles and other data from example_data.h with myProfiles
    // Note that the simplecloud, clwscheme, icecloud and zeeman data are not mandatory and are omitted here
    myProfiles.setGasUnits(gas_units);
    myProfiles.setP((double *)p);
    myProfiles.setT((double *)t);
    myProfiles.setQ((double *)q);
    myProfiles.setCO2((double *)co2);
    myProfiles.setAngles((double *)angles);
    myProfiles.setS2m((double *)s2m);
    myProfiles.setSkin((double *)skin);
    myProfiles.setSurfType((int *)surftype);
    myProfiles.setSurfGeom((double *)surfgeom);
    myProfiles.setDateTimes((int *)datetimes);


    // ------------------------------------------------------------------------
    // Set up Rttov instances for each instrument
    // ------------------------------------------------------------------------

    // Create three Rttov objects for three instruments
    rttov::Rttov seviriRttov = rttov::Rttov();
    rttov::Rttov hirsRttov   = rttov::Rttov();
    rttov::Rttov mhsRttov    = rttov::Rttov();

    // For HIRS and ATMS we will read all channels, but we will read a subset for SEVIRI
    int nchan_seviri = 10;
    int nchan_hirs   = 19;
    int nchan_mhs    = 5;

    // For SEVIRI exclude ozone and hi-res vis channels (9 and 12) in this example
    std::vector <int> chan_list_seviri(nchan_seviri);
    chan_list_seviri = {1, 2, 3, 4, 5, 6, 7, 9, 10, 11};


    // Set the options for each Rttov instance:
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

    // Associate the profiles with each Rttov instance
    try {
        seviriRttov.setProfiles(&myProfiles);
        hirsRttov.setProfiles(&myProfiles);
        mhsRttov.setProfiles(&myProfiles);
    }
    catch (exception& e) {
        std::cerr << "Error setting the profiles " << e.what() << endl;
        exit(1);
    }

    // ------------------------------------------------------------------------
    // Load the emissivity and BRDF atlases
    // ------------------------------------------------------------------------

    // Load the IR emissivity and BRDF atlases:
    // - load data for the month in the profile data
    // - load the IR emissivity atlas data for multiple instruments so it can be used for SEVIRI and HIRS
    // - SEVIRI is the only VIS/NIR instrument we can use the single-instrument initialisation for the BRDF atlas

    rttov::Atlas irAtlas = rttov::Atlas();
    irAtlas.setAtlasPath("../emis_data");
    irAtlas.loadIrEmisAtlas(datetimes[0][1], true);  // Include angular correction, but do not initialise for single-instrument

    rttov::Atlas brdfAtlas = rttov::Atlas();
    brdfAtlas.setAtlasPath("../brdf_data");
    brdfAtlas.loadBrdfAtlas(datetimes[0][1], &seviriRttov); // Pass Rttov object to enable single-instrument initialisation
    brdfAtlas.setIncSea(false);                             // Do not use BRDF atlas for sea surface types

    rttov::Atlas mwAtlas = rttov::Atlas();
    mwAtlas.setAtlasPath("../emis_data");
    mwAtlas.loadMwEmisAtlas(datetimes[0][1]);

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
