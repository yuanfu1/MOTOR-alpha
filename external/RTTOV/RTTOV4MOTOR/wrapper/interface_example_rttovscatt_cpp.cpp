
// Example program calling RTTOV-SCATT from C++ directly through wrapper

// Example profile data is contained in example_data_rttovscatt.h

extern "C" {
#include <rttov_c_interface.h>
}
#include <stdio.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <cstdlib>

#include "example_data_rttovscatt.h"

using namespace std;

// Function for printing output arrays
void printout(int nprofiles, int nchannels, double (*data)) {
    for (int i = 0; i < nprofiles; i++) {
        stringstream outStream;
        outStream << "Profile " << i << " : " << std::fixed << std::setprecision(3);
        for (int c = 0; c < nchannels; c++) outStream << data[i*nchannels + c] << " ";
        std::cerr << outStream.str() << std::endl;
    }
}

int main(int argc, char **argv) {

    int err = 0;
    int calc_zef = 0;
    int multi_hydro_frac = 0;

    // =================================================================
    // Specify the profile data

    // This example demonstrates how to run a simulation for AMSUA-A for two
    // profiles.

    // example_data_rttovscatt.h contains p, T, q, and hydrometor data for a single profile
    // It also contains surface variables and other data for two profiles


    // Define dimensions: must be consistent with data in example_data_rttovscatt.h
    int nprofiles = 2;
    int nlevels = 61;
    int ngases = 6;     // Includes q and RTTOV-SCATT inputs 

    // See wrapper user guide for gas IDs
    int gas_id_q                = 1;
    int gas_id_scatt_hydro_frac = 60;
    int gas_id_scatt_clw        = 61;
    int gas_id_scatt_ciw        = 62;
    int gas_id_scatt_rain       = 63;
    int gas_id_scatt_snow       = 64;

    // The gas ID array tells the wrapper which gases, aerosol and cloud profiles are being supplied:
    // it must specify water vapour in all cases plus any other optional items;
    // the order of the elements is not important, but it must be consistent with the data in the gases array
    int gas_id[] = {gas_id_q, gas_id_scatt_hydro_frac, gas_id_scatt_clw, gas_id_scatt_rain, gas_id_scatt_ciw, gas_id_scatt_snow};

    // Define arrays for pressure, temperature and gases/clouds/aerosols;
    double p[nprofiles][nlevels];              // Input pressure profiles
    double t[nprofiles][nlevels];              // Input temperature profiles
    double gases[ngases][nprofiles][nlevels];  // Input gas profiles
    double ph[nprofiles][nlevels+1];           // Input pressure half-levels
    double cfrac[nprofiles];                   // User cloud fraction

    // Populate the pressure, temperature, q and hydrometeor arrays: these are the same for both profiles
    for (int i = 0; i < nprofiles; i++) {
        for (int l = 0; l < nlevels; l++) {
            p[i][l] = p_ex[l];
            t[i][l] = t_ex[l];
            gases[0][i][l] = q_ex[l];      // index 0 in gas_id array above is water vapour
            gases[1][i][l] = cc_ex[l];     // similarly for cloud inputs...
            gases[2][i][l] = clw_ex[l];
            gases[3][i][l] = rain_ex[l];
            gases[4][i][l] = ciw_ex[l];
            gases[5][i][l] = snow_ex[l];
            ph[i][l] = ph_ex[l];
        }
        ph[i][nlevels] = s2m[i][0];        // Bottom pressure half-level set to 2m pressure
        cfrac[i] = cfrac_ex;
    }

    // The remaining profile data is specified in example_data.h
    // =================================================================


    // =================================================================
    // Load the instrument

    // Specify RTTOV and wrapper options. In this case:
    // - provide access to the full radiance structure after calling RTTOV
    // - turn on the verbose wrapper option
    // - specify interpolation mode
    // NB the spaces in the string between option names and values are important!
    string opts;
    opts.append(" store_rad 1 ");
    opts.append(" verbose_wrapper 1 ");
    opts.append(" opts%interpolation%addinterp 1 ");
    opts.append(" nprofs_per_call 2 ");

    // Specify instrument and channel list and add coefficient files to the options string
    opts.append(" file_coef ../rtcoef_rttov13/rttov13pred54L/rtcoef_noaa_15_amsua.dat ");
    opts.append(" file_hydrotable ../rtcoef_rttov13/hydrotable/hydrotable_noaa_amsua.dat ");

    int nchannels = 15;
    int channel_list[nchannels];
    for (int i = 0; i < nchannels; i++) channel_list[i] = i + 1;

    // Build input string containing options and paths
    int len_opts_str = opts.length();

    char opts_str[len_opts_str];
    for (int i = 0; i < len_opts_str; i++) opts_str[i] = opts[i];


    // Call the wrapper subroutine to load the instrument and check we obtained a valid instrument ID;
    // note that the length of the options string is specified in the last argument
    int inst_id;
    rttov_load_inst_(&inst_id, opts_str, &nchannels, channel_list, len_opts_str);
    if (err != 0) {
        std::cerr << "Error loading instrument" << std::endl;
        exit(1);
    }
    // =================================================================


    // =================================================================
    // Initialise emissivity atlas

    char * emis_atlas_path = "../emis_data/";
    int len_atlas_path_str = strlen(emis_atlas_path);
    int month = datetimes[0][1];            // Month is taken from the profile date
    int year = datetimes[0][0];             // Only used by CNRM atlas

/*    // Call the wrapper subroutine to set up the emissivity atlas
    int version = -1;   // default atlas version
    rttov_mw_emis_atlas_setup_(&err, emis_atlas_path, &month, &version, &inst_id, &year, len_atlas_path_str);
    if (err != 0) std::cerr << "Error initialising MW emissivity atlas: atlas will not be used" << std::endl;*/
    // =================================================================


    // =================================================================
    // Declare arrays for other inputs and outputs

    // Define array for input/output surface emissivity and BRDF
    double surfemis[nprofiles][nchannels];

    // Define direct model outputs
    double bt[nprofiles][nchannels];
    // =================================================================


    // =================================================================
    // Call RTTOV

    // Initialise the surface emissivity before every call to RTTOV
    for (int i = 0; i < nprofiles; i++) {
        for (int c = 0; c < nchannels; c++) {
            surfemis[i][c] = -1.;
        }
    }

    // Call the wrapper subroutine to run RTTOV direct
    rttov_scatt_call_direct_(
        &err,
        &inst_id,
        channel_list,
        (int *)datetimes,            // profile dates/times                                     [nprofiles][6]
        (double *)angles,            // satzen, satazi angles                                   [nprofiles][2]
        (double *)surfgeom,          // lat, lon, elevation                                     [nprofiles][3]
        (int *)surftype,             // surftype                                                [nprofiles]
        (double *)skin,              // skin T, salinity, foam_frac, fastem_coefsx5             [nprofiles][8]
        (double *)s2m,               // 2m p, 2m t, 2m q, 10m wind u, v                         [nprofiles][5]
        (double *)zeeman,            // Be, cosbk                                               [nprofiles][2]
        (double *)p,                 // pressure                                                [nprofiles][nlevels]
        (double *)t,                 // temperature                                             [nprofiles][nlevels]
        &gas_units,                  // units for gas profiles
        gas_id,                      // gas ID list                                             [ngases]
        (double *)gases,             // gas profiles                                            [ngases][nprofiles][nlevels]
        (double *)ph,                // pressure half-levels                                    [nprofiles][nlevels+1]
        (double *)cfrac,             // user cloud fraction                                     [nprofiles]
        &multi_hydro_frac,           // false => single hydro_frac profile, true => one hydro_frac profile per hydrometeor
        &calc_zef,                   // enable/disable radar reflectivity calculations
        (double *)surfemis,          // input/output surface emissivities                       [nprofiles][nchannels]
        (double *)bt,                // output BTs                                              [nprofiles][nchannels]
        &nchannels, &ngases, &nlevels, &nprofiles);

    if (err != 0) {
        std::cerr << "Error running RTTOV direct" << std::endl;
        exit(1);
    }
    // =================================================================


    // =================================================================
    // Examine outputs

    // Outputs available are:
    // - surfemis array contains surface emissivities used by RTTOV
    // - rad array contains RTTOV radiance%total array
    // - bt array contains RTTOV radiance%bt arrays
    // - it is also possible to access the whole radiance structure because we set the store_rad option above

    std::cerr << "Surface emissivity used by RTTOV" << std::endl;
    printout(nprofiles, nchannels, (double *)surfemis);

    std::cerr << "Total cloudy BT" << std::endl;
    printout(nprofiles, nchannels, (double *)bt);

    // To obtain data from RTTOV output structures, declare an array and call the relevant wrapper subroutine.
    // For example for the clear-sky BTs:
    double btclear[nprofiles][nchannels];
    int nchanprof = nprofiles * nchannels;
    rttov_get_bt_clear_(&err, &inst_id, (double *)btclear, &nchanprof);
    std::cerr << "Clear-sky BT" << std::endl;
    printout(nprofiles, nchannels, (double *)btclear);
    // =================================================================


    // =================================================================
    // Deallocate memory for all instruments and atlases

    rttov_drop_all_(&err);
    if (err != 0) std::cerr << "Error deallocating wrapper" << std::endl;
    // =================================================================
}

