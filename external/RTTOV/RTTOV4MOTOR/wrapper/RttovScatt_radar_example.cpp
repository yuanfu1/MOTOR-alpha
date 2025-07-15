
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
            clw[i][l] = clw_ex[l];
            ciw[i][l] = ciw_ex[l];
            snow[i][l] = snow_ex[l];
            rain[i][l] = rain_ex[l];

            // Modified cloud cover to represent a rain fraction below the peak in the rain-producing cloud
            cc[i][l] = cc_ex[l];
            if (l >= 53) cc[i][l] = cc_ex[52];
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
    myProfiles.setAngles((double *)angles);
    myProfiles.setS2m((double *)s2m);
    myProfiles.setSkin((double *)skin);
    myProfiles.setSurfType((int *)surftype);
    myProfiles.setSurfGeom((double *)surfgeom);
    myProfiles.setDateTimes((int *)datetimes);

    // ------------------------------------------------------------------------
    // Set up RttovScatt instance
    // ------------------------------------------------------------------------

    bool doK = true;

    rttov::RttovScatt dprRttov = rttov::RttovScatt();

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

    // Associate the profiles

    try {
        dprRttov.setProfiles(&myProfiles);
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
