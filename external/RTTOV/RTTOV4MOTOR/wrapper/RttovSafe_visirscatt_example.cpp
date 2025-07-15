
// Example of using the RttovSafe class to do visible/IR scattering simulations
// with explicit optical parameters

#include <Profile.h>
#include <RttovSafe.h>
#include <Atlas.h>
#include <vector>

#include "example_data_cpp.h"
#include "example_data_opt_param.h"


int main(int argc, char **argv) {

    // This example program simulates two profiles
    // The example profile data are defined in example_data_cpp.h
    // The example optical properties are defined in example_data_opt_param.h

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

    // Populate the Profile objects in the vector profiles
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
            profiles[i].setT(t_ex);
            profiles[i].setQ(q_ex);
            profiles[i].setSkin(270., 35., 0., 0., 3.0, 5.0, 15.0, 0.1, 0.3);
            profiles[i].setSurfGeom(10., 20., 0.);
            profiles[i].setDateTimes(2015, 8, 1, 0, 0, 0);
        }
    }
    catch (exception& e) {
        std::cerr << "Error defining the profile data " << e.what() << endl;
        exit(1);
    }

    // ------------------------------------------------------------------------
    // Set up the Rttov instance
    // ------------------------------------------------------------------------

    rttov::RttovSafe seviriRttov = rttov::RttovSafe();

    // Set the options

    seviriRttov.setFileCoef("../rtcoef_rttov13/rttov13pred54L/rtcoef_msg_3_seviri_o3.dat");
    seviriRttov.options.setAddInterp(true);        // Use the RTTOV interpolator
    seviriRttov.options.setAddSolar(true);         // Turn on solar radiation
    seviriRttov.options.setVerboseWrapper(true);   // Turn on verbose wrapper output
    seviriRttov.options.setAddClouds(true);        // Activate cloud simulations
    seviriRttov.options.setUserCldOptParam(true);  // Use explicit cloud optical properties
    seviriRttov.options.setVisScattModel(1);       // Use DOM solver for solar radiation
    seviriRttov.options.setIrScattModel(2);        // Use Chou-scaling for thermal emitted radiation
    seviriRttov.options.setDomNstreams(16);        // Number of DOM streams to use
    seviriRttov.options.setNthreads(8);            // Take advantage of multiple threads if RTTOV was compiled with OpenMP
    seviriRttov.options.setStoreRad(true);         // Store all radiance outputs

    // Load the instrument (reads all channels)
    try {
        seviriRttov.loadInst();
    }
    catch (exception& e) {
        std::cerr << "Error loading instrument " << e.what() << std::endl;
        exit(1);
    }

    // ------------------------------------------------------------------------
    // Specify cloud optical parameters
    // ------------------------------------------------------------------------

    // The example_data_opt_param file contains data for a single cloud layer for two channels:
    //   opt_param_chan_list[nchan]  - the channels for which optical parameters are provided
    //   abscoef[nchan]              - absorption coefficient for each channel
    //   scacoef[nchan]              - scattering coefficient for each channel
    //   phangle[nphangle]           - phase function angle grid
    //   phasefn[nchan][nphangle]    - phase functions for each channel

    // In this example we define a single-layer cloud for two channels: one visible, one IR

    // For cloudy simulations the cloud fraction profile must be specified in the Profiles object
    std::vector <double> cfrac;
    for (int l = 0; l < nlevels; l++) {
      if (l == 74) {
        cfrac.push_back(1.); // cloud fraction: both profiles, layer 75
      } else {
        cfrac.push_back(0.);
      }
    }
    for (int p = 0; p < nprofiles; p++) profiles[p].setCfrac(cfrac);

    // Define optical parameters:
    // - these are defined for every *layer* for every channel being simulated for every profile
    // - absorption and scattering coefficients are always required
    // - bpr parameter is only required when using Chou-scaling, can be zero otherwise; this can
    //   be calculated from the phase function using the calcBpr method of Rttov
    // - phase functions are required only for solar-affected channels when the AddSolar option is true
    // - Legendre coefficients are only required if using the DOM solver; these can be calculated
    //   from the phase function using the calcLegcoef method of Rttov

    // In the second profile the scattering coefficient is doubled

    // Absorption coef, scattering coef, bpr parameter [3][nprofiles][nchan][nlayers]
    double asb[3][nprofiles][nchan][nlevels-1];
    for (int i = 0; i < 3; i++) {
      for (int p = 0; p < nprofiles; p++) {
        for (int c = 0; c < nchan; c++) {
          for (int l = 0; l < nlevels-1; l++) asb[i][p][c][l] = 0.;
        }
      }
    }
    for (int p = 0; p < nprofiles; p++) {
      for (int c = 0; c < nchan; c++) {
        asb[0][p][c][74] = abscoef[c];        // abs coef: both profiles, both channels, layer 75
        if (p == 0) {
          asb[1][p][c][74] = scacoef[c];      // sca coef: profile 1, both channels, layer 75
        } else {
          asb[1][p][c][74] = scacoef[c] * 2.; // sca coef: profile 2, both channels, layer 75
        }
      }
    }

    // Since Chou-scaling is used in the IR we must calculate the bpr for the IR channel.
    // Note that this is relatively slow so for performance-critical applications this
    // calculation would be done off-line for each phase function. If RTTOV was compiled
    // with OpenMP it uses the number of threads specified in the wrapper options.
    // The bpr value is left as zero for the visible channel.
    double bpr = seviriRttov.calcBpr(nphangle, phangle, phasefn[1]);
    asb[2][0][1][74] = bpr; // bpr: profile 1, channel 2, layer 75
    asb[2][1][1][74] = bpr; // bpr: profile 2, channel 2, layer 75

    // Specify phase function for the visible channel. This is left as zero for the IR channel.
    double pha[nprofiles][nchan][nlevels-1][nphangle];
    for (int p = 0; p < nprofiles; p++) {
      for (int c = 0; c < nchan; c++) {
        for (int l = 0; l < nlevels-1; l++) {
          for (int i = 0; i < nphangle; i++) {
            if (c == 0 && l == 74) {
              pha[p][c][l][i] = phasefn[c][i]; // phase fn: both profiles, channel 1, layer 75
            } else {
              pha[p][c][l][i] = 0.;
            }
          }
        }
      }
    }

    // Calculate Legendre coefficients for visible channel as DOM is being used. We only need
    // DomNstreams coefficients (note that this requires a DomNstreams+1 sized array: the zeroth
    // coefficient is always 1.) These are left as zero for the IR channel.
    int nmom = seviriRttov.options.getDomNstreams();
    double legcoef[nprofiles][nchan][nlevels-1][nmom+1];
    double lc[nmom+1];
    seviriRttov.calcLegcoef(nphangle, phangle, phasefn[0], nmom, lc, 0);
    for (int p = 0; p < nprofiles; p++) {
      for (int c = 0; c < nchan; c++) {
        for (int l = 0; l < nlevels-1; l++) {
          for (int i = 0; i < nmom+1; i++) {
            if (c == 0 && l == 74) {
              legcoef[p][c][l][i] = lc[i]; // Leg. coefs: both profiles, channel 1, layer 75
            } else {
              legcoef[p][c][l][i] = 0.;
            }
          }
        }
      }
    }

    // Now assign the optical property arrays to the Rttov object
    seviriRttov.setCldAsb((double *)asb);
    seviriRttov.setCldPha(nphangle, phangle, (double *)pha);
    seviriRttov.setCldLegcoef(nmom, (double *)legcoef);


    // Associate the profiles with the RttovSafe instance: the profiles undergo some
    // checks so use a try block to catch any errors that are thrown up
    try {
        seviriRttov.setTheProfiles(profiles);
    }
    catch (exception& e) {
        std::cerr << "Error setting the profiles " << e.what() << endl;
        exit(1);
    }

    // ------------------------------------------------------------------------
    // Call RTTOV
    // ------------------------------------------------------------------------

    // In this example we let the Rttov class automatically set calcemis/calcrefl
    // to true for simplicity. You can supply emissivity/BRDF values or use the
    // atlases in exactly the same way as for other simulation types.

    // Call the RTTOV direct model: note the channel list must be consistent with the
    // channels for which optical properties are specified.
    try {
        seviriRttov.runDirect(opt_param_chan_list);
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

    std::vector <double> bt, btclear, refl, reflclear;
    for (int p = 0; p < nprofiles; p++) {
        std::cout << "Profile " << p << " : " << std::endl;
        bt = seviriRttov.getBt(p);
        btclear = seviriRttov.getBtClear(p);
        refl = seviriRttov.getRefl(p);
        reflclear = seviriRttov.getReflClear(p);
        std::cout << "  Visible channel cloudy and clear reflectances:" << std::endl;
        std::cout << "  " << refl[0] << "  " << reflclear[0] << std::endl;
        std::cout << "  IR channel cloudy and clear BTs:" << std::endl;
        std::cout << "  " << bt[1] << "  " << btclear[1] << std::endl;
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
