
// Example of using the RttovScatt class to call RTTOV-SCATT

#include <ProfilesScatt.h>
#include <RttovScatt.h>
#include <Atlas.h>
#include <vector>

#include "example_data_rttovscatt.h"


int main(int argc, char **argv) {

    // This example program simulates two profiles
    // The example profile data are defined in example_data_rttovscatt.h

    // ------------------------------------------------------------------------
    // Set up the profile data
    // ------------------------------------------------------------------------

    // Declare an instance of Profiles
    int nlevels = 61;
    int nprofiles = 2;
    rttov::ProfilesScatt myProfiles = rttov::ProfilesScatt(nprofiles, nlevels);

    // Populate myProfiles with data
    double p[nprofiles][nlevels];    // Input pressure profiles
    double ph[nprofiles][nlevels+1]; // Input pressure half-level profiles
    double t[nprofiles][nlevels];    // Input temperature profiles
    double q[nprofiles][nlevels];    // Input water vapour profiles
    double cc[nprofiles][nlevels];   // Input CC profiles
    double clw[nprofiles][nlevels];  // Input CLW profiles
    double ciw[nprofiles][nlevels];  // Input CIW profiles
    double snow[nprofiles][nlevels]; // Input snow profiles
    double rain[nprofiles][nlevels]; // Input rain profiles

    // example_data.h contains p, T, q and hydrometeors for a single profile
    for (int i = 0; i < nprofiles; i++) {
        for (int l = 0; l < nlevels; l++) {
            p[i][l] = p_ex[l];
            ph[i][l] = ph_ex[l];
            t[i][l] = t_ex[l];
            q[i][l] = q_ex[l];
            cc[i][l] = cc_ex[l];
            clw[i][l] = clw_ex[l];
            ciw[i][l] = ciw_ex[l];
            snow[i][l] = snow_ex[l];
            rain[i][l] = rain_ex[l];
        }
        ph[i][nlevels] = s2m[i][0];
    }

    // Associate the profiles and other data from example_data.h with myProfiles
    // Note that the zeeman and cfrac data are not mandatory and are omitted here
    myProfiles.setGasUnits(gas_units);
    myProfiles.setP((double *)p);
    myProfiles.setPh((double *)ph);
    myProfiles.setT((double *)t);
    myProfiles.setQ((double *)q);

    myProfiles.setHydroFrac((double *)cc);
    myProfiles.setClw((double *)clw);
    myProfiles.setCiw((double *)ciw);
    myProfiles.setSnow((double *)snow);
    myProfiles.setRain((double *)rain);

    // Equivalent code using the flexible hydro interface
    //myProfiles.setGasItem((double *)rain, rttov::HYDRO1);
    //myProfiles.setGasItem((double *)snow, rttov::HYDRO2);
    //myProfiles.setGasItem((double *)clw, rttov::HYDRO4);
    //myProfiles.setGasItem((double *)ciw, rttov::HYDRO5);
    //myProfiles.setGasItem((double *)cc, rttov::HYDRO_FRAC1);

    myProfiles.setAngles((double *)angles);
    myProfiles.setS2m((double *)s2m);
    myProfiles.setSkin((double *)skin);
    myProfiles.setSurfType((int *)surftype);
    myProfiles.setSurfGeom((double *)surfgeom);
    myProfiles.setDateTimes((int *)datetimes);

    // ------------------------------------------------------------------------
    // Set up RttovScatt instance
    // ------------------------------------------------------------------------

    rttov::RttovScatt amsuaRttov = rttov::RttovScatt();

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

    // Associate the profiles

    try {
        amsuaRttov.setProfiles(&myProfiles);
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
    mwAtlas.loadMwEmisAtlas(datetimes[0][1]);
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
