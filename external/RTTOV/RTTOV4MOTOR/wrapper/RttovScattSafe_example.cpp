
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

        std::vector <double> ph;
        for (int i = 0; i < nlevels; i++) ph.push_back(ph_ex[i]);
        ph.push_back(0.);
        for (int i = 0; i < nprofiles; i++) {
            profiles[i].setGasUnits(rttov::ppmv_wet);
            profiles[i].setP(p_ex);
            profiles[i].setT(t_ex);
            profiles[i].setQ(q_ex);
            
            profiles[i].setHydroFrac(cc_ex);
            profiles[i].setClw(clw_ex);
            profiles[i].setCiw(ciw_ex);
            profiles[i].setSnow(snow_ex);
            profiles[i].setRain(rain_ex);

            // Equivalent code using the flexible hydro interface
            //profiles[i].setHydroN(rain_ex, 1);
            //profiles[i].setHydroN(snow_ex, 2);
            //profiles[i].setHydroN(clw_ex, 4);
            //profiles[i].setHydroN(ciw_ex, 5);
            //profiles[i].setHydroFracN(cc_ex, 1);

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

    rttov::RttovScattSafe amsuaRttov = rttov::RttovScattSafe();

    amsuaRttov.setFileCoef("../rtcoef_rttov13/rttov13pred54L/rtcoef_noaa_15_amsua.dat");
    amsuaRttov.setFileHydrotable("../rtcoef_rttov13/hydrotable/hydrotable_noaa_amsua.dat");

    amsuaRttov.options.setLuserCfrac(false);
    amsuaRttov.options.setVerboseWrapper(true);
    amsuaRttov.options.setStoreRad(true);

    // Load the instrument
    try {
        amsuaRttov.loadInst();
    }
    catch (exception& e) {
        std::cerr << "Error loading instrument(s) " << e.what() << std::endl;
        exit(1);
    }

    int nchan_amsua = amsuaRttov.getNchannels();

    // Associate the profiles: the profiles undergo some
    // checks so use a try block to catch any errors that are thrown up
    try {
      amsuaRttov.setTheProfiles(profiles);
    }
    catch (exception& e) {
        std::cerr << "Error setting the profiles " << e.what() << endl;
        exit(1);
    }

    // ------------------------------------------------------------------------
    // Load the emissivity atlas
    // ------------------------------------------------------------------------

    // TELSEM2 atlas does not require an RttovScatt object to initialise
    rttov::Atlas mwAtlas = rttov::Atlas();
    mwAtlas.setAtlasPath("../emis_data");
    mwAtlas.loadMwEmisAtlas(month);
    mwAtlas.setIncSea(false);

    // Set up the surface emissivity/reflectance arrays and associate with the RttovScatt object
    double surfemis_amsua[nprofiles][nchan_amsua];
    amsuaRttov.setSurfEmis((double *)surfemis_amsua);


    // ------------------------------------------------------------------------
    // Call RTTOV-SCATT
    // ------------------------------------------------------------------------

    // Surface emissivity arrays must be initialised *before every call to RTTOV*
    // Negative values will cause RTTOV to supply emissivity values (i.e. equivalent to
    // calcemis TRUE - see RTTOV user guide)

    for (int i = 0; i < nprofiles; i++) {
        for (int c = 0; c < nchan_amsua; c++) surfemis_amsua[i][c] = -1.;
    }

    // Call emissivity atlas
    try {
        mwAtlas.fillEmisBrdf((double *)surfemis_amsua, &amsuaRttov);
    }
    catch (exception& e) {
        // If there was an error the emissivities/BRDFs will not have been modified so it
        // is OK to continue and call RTTOV with calcemis/calcrefl set to TRUE everywhere
        std::cerr << "Error calling atlas " << e.what() << std::endl;
    }


    // Call the RTTOV-SCATT direct model:
    // no arguments are supplied to runDirect so all loaded channels are simulated
    try {
        //amsuaRttov.runDirect();
        amsuaRttov.runK();
    }
    catch (exception& e) {
        std::cerr << "Error running RTTOV direct model " << e.what() << std::endl;
        exit(1);
    }


    // ------------------------------------------------------------------------
    // Print out some of the output
    // ------------------------------------------------------------------------

    std::cout << "Surface emissivity used by RTTOV" << std::endl;
    for (int p = 0; p < nprofiles; p++) {
        std::cout << "Profile " << p << " : ";
        for (int c = 0; c < nchan_amsua; c++) std::cout << surfemis_amsua[p][c] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Total cloudy BT" << std::endl;
    std::vector <double> bt;
    for (int p = 0; p < nprofiles; p++) {
        bt = amsuaRttov.getBt(p);
        std::cout << "Profile " << p << " : ";
        for (int c = 0; c < nchan_amsua; c++) std::cout << bt[c] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Clear-sky BT" << std::endl;
    std::vector <double> btclear;
    for (int p = 0; p < nprofiles; p++) {
        btclear = amsuaRttov.getBtClear(p);
        std::cout << "Profile " << p << " : ";
        for (int c = 0; c < nchan_amsua; c++) std::cout << btclear[c] << " ";
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Example code for rain Jacobian using flexible hydro interface
    //   (profiles must have been specified using the same interface)
    //std::cout << "Rain Jacobian profile 1, channel 1" << std::endl;
    //std::vector <double> raink = amsuaRttov.getItemK(rttov::HYDRO1, 0, 0);
    //for (int lev = 0; lev < nlevels; lev++) std::cout << raink[lev] << " ";
    //std::cout << std::endl;


    // ------------------------------------------------------------------------
    // Deallocate memory
    // ------------------------------------------------------------------------

    // In this example nothing needs to be done: the destructors of the Rttov
    // and Profile classes will deallocate the atlases and other RTTOV data and
    // all locally defined arrays will be deallocated by the garbage collection

}
