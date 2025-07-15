
// Example program calling RTTOV from C++ directly through wrapper

// Example profile data is contained in example_data.h

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

#include "example_data.h"

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

    // =================================================================
    // Specify the profile data

    // This example demonstrates how to run a simulation for HIRS for two
    // profiles with variable CO2 and a very simple cloud profile

    // example_data.h contains p, T, q, co2 for a single profile
    // It also contains surface variables and other data for two profiles


    // Define dimensions: must be consistent with data in example_data.h
    int nprofiles = 2;
    int nlevels = 101;
    int ngases = 5;

    // See wrapper user guide for gas IDs
    int gas_id_q = 1;
    int gas_id_o3 = 2;
    int gas_id_co2 = 3;

    int gas_id_cfrac = 20;
    int gas_id_lwc1  = 21;
    int gas_id_iwc   = 30;

    // The gas ID array tells the wrapper which gases, aerosol and cloud profiles are being supplied:
    // it must specify water vapour in all cases plus any other optional items;
    // the order of the elements is not important, but it must be consistent with the data in the gases array
    int gas_id[] = {gas_id_q, gas_id_co2, gas_id_cfrac, gas_id_lwc1, gas_id_iwc};

    // Define arrays for pressure, temperature and gases/clouds/aerosols;
    double p[nprofiles][nlevels];              // Input pressure profiles
    double t[nprofiles][nlevels];              // Input temperature profiles
    double gases[ngases][nprofiles][nlevels];  // Input gas profiles

    // Populate the pressure, temperature, q and co2 arrays: these are the same for both profiles
    for (int i = 0; i < nprofiles; i++) {
        for (int l = 0; l < nlevels; l++) {
            p[i][l] = p_ex[l];
            t[i][l] = t_ex[l];
            gases[0][i][l] = q_ex[l];       // index 0 in gas_id array above is water vapour
            gases[1][i][l] = co2_ex[l];     // index 1 in gas_id array above is co2
            for (int g = 2; g < ngases; g++) gases[g][i][l] = 0.; // Initialise cloud inputs to zero
        }
    }

    // Specify some very simple cloud inputs (in layers 50 and 60) in profile 1; profile 2 contains no cloud
    gases[2][0][49] = 0.5; // cfrac in layer 50, profile 1
    gases[4][0][49] = 1.0; // IWC in layer 50, profile 1
    gases[2][0][59] = 0.8; // cfrac in layer 60, profile 1
    gases[3][0][59] = 1.0; // LWC in layer 60, profile 1

    int mmr_cldaer = 0;
    // The remaining profile data is specified in example_data.h
    // =================================================================


    // =================================================================
    // Load the instrument

    // Specify RTTOV and wrapper options. In this case:
    // - turn interpolation on
    // - supply CO2 as a variable gas
    // - turn cloudy IR simulations on
    // - provide access to the full radiance structure after calling RTTOV
    // - turn on the verbose wrapper option
    // NB the spaces in the string between option names and values are important!
    string opts;
    opts.append(" opts%interpolation%addinterp 1 ");
    opts.append(" opts%rt_all%co2_data 1 ");
    opts.append(" opts%rt_ir%addclouds 1 ");
    opts.append(" store_rad 1 ");
    opts.append(" verbose_wrapper 1 ");

    // Specify instrument and channel list and add coefficient files to the options string
    opts.append(" file_coef ../rtcoef_rttov13/rttov13pred54L/rtcoef_metop_1_hirs_o3co2.dat ");
    opts.append(" file_sccld ../rtcoef_rttov13/cldaer_visir/sccldcoef_metop_1_hirs.dat ");

    int nchannels = 19;
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
    if (inst_id < 1) {
        std::cerr << "Error loading instrument" << std::endl;
        exit(1);
    }
    // =================================================================


    // =================================================================
    // Initialise emissivity atlas

    char * emis_atlas_path = "../emis_data/";
    int len_atlas_path_str = strlen(emis_atlas_path);
    int month = datetimes[0][1];            // Month is taken from the profile date

    // Call the wrapper subroutine to set up the IR emissivity atlas
    // NB we specify inst_id here so the atlas is initialised for this specific instrument for faster access;
    //    to initialise the atlas for use with multiple instruments pass 0 as the inst_id
    //    (see wrapper user guide for more information)
    int atlas_wrap_id;
    int atlas_id = -1;   // use default IR atlas
    int ang_corr = 0;    // no angular correction
    rttov_load_ir_emis_atlas_(&atlas_wrap_id, emis_atlas_path, &month, &atlas_id, &inst_id, &ang_corr, len_atlas_path_str);
    if (atlas_wrap_id < 1) std::cerr << "Error initialising IR emissivity atlas: atlas will not be used" << std::endl;
    // =================================================================


    // =================================================================
    // Declare arrays for other inputs and outputs

    // Define array for input/output surface emissivity and BRDF
    double surfemisrefl[4][nprofiles][nchannels];

    // Initialise to negative values: RTTOV will supply emissivity and reflectance values
    // unless positive values are provided
    for (int i = 0; i < nprofiles; i++) {
        for (int c = 0; c < nchannels; c++) {
            surfemisrefl[0][i][c] = -1.;
            surfemisrefl[1][i][c] = -1.;
            surfemisrefl[2][i][c] = -1.;
            surfemisrefl[3][i][c] = -1.;
        }
    }

    // Define direct model outputs
    double btrefl[nprofiles][nchannels];
    double rads[nprofiles][nchannels];
    // =================================================================

    // =================================================================
    // Obtain emissivities from atlas

    double latitude[nprofiles];
    double longitude[nprofiles];
    int surftypeval[nprofiles];
    int watertype[nprofiles];
    double zenangle[nprofiles];
    double azangle[nprofiles];
    double sunzenangle[nprofiles];
    double sunazangle[nprofiles];
    double snow_frac[nprofiles];
    for (int i = 0; i < nprofiles; i++) {
        latitude[i]    = surfgeom[i][0];
        longitude[i]   = surfgeom[i][1];
        surftypeval[i] = surftype[i][0];
        watertype[i]   = surftype[i][1];
        zenangle[i]    = angles[i][0];
        azangle[i]     = angles[i][1];
        sunzenangle[i] = angles[i][2];
        sunazangle[i]  = angles[i][3];
        snow_frac[i]   = skin[i][2];
    }

    if (atlas_wrap_id > 0) {
        rttov_get_emisbrdf_(
            &err,
            &atlas_wrap_id,
            latitude,         // [nprofiles]
            longitude,        // [nprofiles]
            surftypeval,      // [nprofiles]
            watertype,        // [nprofiles]
            zenangle,         // [nprofiles]
            azangle,          // [nprofiles]
            sunzenangle,      // [nprofiles]
            sunazangle,       // [nprofiles]
            snow_frac,        // [nprofiles]
            &inst_id,
            channel_list,                 // [nchannels]
            (double *)surfemisrefl[0],    // [nprofiles][nchannels]
            &nchannels, &nprofiles);

        if (err != 0) {
            std::cerr << "Error getting emissivities from atlas: atlas will not be used" << std::endl;
            for (int i = 0; i < nprofiles; i++) {
                for (int c = 0; c < nchannels; c++) {
                    surfemisrefl[0][i][c] = -1.;
                    surfemisrefl[1][i][c] = -1.;
                }
            }
        }
    }
    // =================================================================


    // =================================================================
    // Call RTTOV

    // Call the wrapper subroutine to run RTTOV direct
    rttov_call_direct_(
        &err,
        &inst_id,
        channel_list,
        (int *)datetimes,            // profile dates/times                                     [nprofiles][6]
        (double *)angles,            // satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
        (double *)surfgeom,          // lat, lon, elevation                                     [nprofiles][3]
        (int *)surftype,             // surftype, watertype                                     [nprofiles][2]
        (double *)skin,              // skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
        (double *)s2m,               // 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
        (double *)simplecloud,       // ctp, cfraction                                          [nprofiles][2]
        (int *)clwscheme,            // clw_scheme, clwde_param                                 [nprofiles][2]
        (int *)icecloud,             // ice_scheme, icede_param                                 [nprofiles][2]
        (double *)zeeman,            // Be, cosbk                                               [nprofiles][2]
        (double *)p,                 // pressure                                                [nprofiles][nlevels]
        (double *)t,                 // temperature                                             [nprofiles][nlevels]
        &gas_units,                  // units for gas profiles
        &mmr_cldaer,                 // units for cloud/aerosol profiles
        gas_id,                      // gas ID list                                             [ngases]
        (double *)gases,             // gas profiles                                            [ngases][nprofiles][nlevels]
        (double *)surfemisrefl,      // input/output surface emissivities/BRDFs                 [2][nprofiles][nchannels]
        (double *)btrefl,            // output BTs/refls (for thermal/solar chans)              [nprofiles][nchannels]
        (double *)rads,              // output radiances                                        [nprofiles][nchannels]
        &nchannels, &ngases, &nlevels, &nprofiles);

    if (err != 0) {
        std::cerr << "Error running RTTOV direct" << std::endl;
        exit(1);
    }
    // =================================================================


    // =================================================================
    // Examine outputs

    // Outputs available are:
    // - surfemisrefl array contains surface emissivities (and reflectances) used by RTTOV
    // - rad array contains RTTOV radiance%total array
    // - btrefl array contains RTTOV radiance%bt and radiance%refl arrays (depending on channel wavelength)
    // - it is also possible to access the whole radiance structure because we set the store_rad option above

    std::cerr << "Surface emissivity used by RTTOV" << std::endl;
    printout(nprofiles, nchannels, (double *)surfemisrefl);

    std::cerr << "Total cloudy BT" << std::endl; // This example has no visible/near-IR channels so this array contains BTs only
    printout(nprofiles, nchannels, (double *)btrefl);

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

