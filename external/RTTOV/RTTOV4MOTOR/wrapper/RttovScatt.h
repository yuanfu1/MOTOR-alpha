/*
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2015, EUMETSAT, All Rights Reserved.
*/

/// @class RttovScatt
///  This class contains a set of functions for running RTTOV-SCATT.
///  An instance of an RttovScatt Class is defined to handle one instrument.
///  It contains an Options object which allows the Rttov options to
///  be specified and can be associated with a ProfilesScatt object.


#ifndef RTTOVSCATT_H_
#define RTTOVSCATT_H_


using namespace std;
#include <rttov_cc_interface.h>

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <Rttov_common.h>
#include <ProfilesScatt.h>
#include <ProfileScatt.h>
#include <Options.h>

namespace rttov {

class RttovScatt {

public:
    RttovScatt();
    virtual ~RttovScatt();

    rttov::Options options;

    const string& getFileCoef() const;
    const string& getFileHydrotable() const;

    void setFileCoef(const string& fileCoef);
    void setFileHydrotable(const string& fileHydrotable);

    bool isCalcZef() const;
    void setCalcZef(bool calcZef);
    bool isMultiHydroFrac() const;
    void setMultiHydroFrac(bool multiHydroFrac);

    void loadInst(); // load instrument for all channels
    int getInstId() const;

    bool isCoeffsLoaded() const;
    int getNchannels() const;
    int getCoeffsNlevels();
    double * getWaveNumbers();
    bool isProfileSet() const;
    int getNprofiles() const;

    void updateOptions();
    void printOptions();

    void setSurfEmis(double* surfemis);

    rttov::ProfilesScatt * getProfiles() const;
    void setProfiles(rttov::ProfilesScatt * profiles);
    void printGases();

    void runDirect(); // run direct for all loaded channels
    void runDirect(const vector<int>& channels);
    void runK() ;
    void runK(const vector<int>& channels);

    const double* getBt() const;
    std::vector<double> getBt(const int profile);
    const double* getSurfEmis() const;

    std::vector <double> getPK(int profile, int channel);
    std::vector <double> getPhK(int profile, int channel);
    std::vector <double> getTK(int profile, int channel);
    std::vector <double> getUserCfracK(int profile);
    std::vector <double> getSkinK(int profile, int channel);
    std::vector <double> getS2mK(int profile, int channel);
    std::vector <double> getItemK(rttov::itemIdType, int profile, int channel);
    std::vector <double> getSurfEmisK(int profile);

    std::vector <double> getBtClear(int profile);
    std::vector <int>    getRadQuality(int profile);
    bool getPlaneParallel();
    std::vector <double> getGeometricHeight(int profile, int channel);

    std::vector <double> getZef(int profile, int channel);
    std::vector <double> getAZef(int profile, int channel);

    std::vector <double> getEmisTermsCfrac(int profile);
    std::vector <double> getEmisTermsBsfc(int profile);
    std::vector <double> getEmisTermsTauCld(int profile);
    std::vector <double> getEmisTermsUpCld(int profile);
    std::vector <double> getEmisTermsDownCld(int profile);
    std::vector <double> getEmisTermsTauClr(int profile);
    std::vector <double> getEmisTermsUpClr(int profile);
    std::vector <double> getEmisTermsDownClr(int profile);

protected :

    rttov::ProfilesScatt * profiles;
    string strOptions;
    int nprofiles;
    int nlevels;
    int nchannels; // number of channels loaded from coefficients
    int nchannelsForLastRun;  // number of channels used for the last run
    int ngases;
    int gas_units;
    int inst_id;
    double * p;//[nprofiles][nlevels];                  // Input pressure profiles
    double * ph;//[nprofiles][nlevels+1];               // Input pressure half-level profiles
    double * t;//[nprofiles][nlevels];                  // Input temperature profiles
    double * usercfrac;//[nprofiles];                   // User cloud fraction
    double * gases;//[ngases][nprofiles][nlevels];      // Input gas profiles
    double * surfemis;//[nprofiles][nchannels];         // Input/output surface emissivities
    double * bt;//[nprofiles][nchannels];               // Output BTs

    double * btclear;
    int * radquality;
    int plane_parallel;
    double * geometric_height;
    double * zef;
    double * azef;

    double * emis_terms_cfrac;
    double * emis_terms_bsfc;
    double * emis_terms_tau_cld;
    double * emis_terms_up_cld;
    double * emis_terms_down_cld;
    double * emis_terms_tau_clr;
    double * emis_terms_up_clr;
    double * emis_terms_down_clr;

    int * gas_id;

    // datetimes: yy, mm, dd, hh, mm, ss
    int * datetimes; //[nprofiles][6]

    // angles: satzen, satazi
    double * angles;//[nprofiles][2]

    // surftype
    int * surftype;//[nprofiles]

    // surfgeom: lat, lon, elev
    double * surfgeom ;//[nprofiles][3]

    // s2m: 2m p, 2m t, 2m q, 10m wind u, v
    double * s2m;//[nprofiles][5]

    // skin: skin T, salinity, foam_frac, fastem_coefsx5
    double * skin; //[nprofiles][8]

    // zeeman: be, cosbk
    double * zeeman;//[nprofiles][2]

    string file_coef;
    string file_hydrotable;
    int calc_zef;
    int multi_hydro_frac;

    bool coeffsLoaded;
    bool profileSet;

    double * wavenumbers;
    bool allocatedSurfemis;

    double * skin_k; //skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][nchannels][9]
    double * s2m_k; //2m p, 2m t, 2m q, 10m wind u, v, wind-fetch    [nprofiles][nchannels][6]
    double * simplecloud_k ; // [nprofiles][nchannels][2]
    double * p_k; // [nprofiles][nchannels][nlevels]
    double * ph_k; // [nprofiles][nchannels][nlevels+1]
    double * t_k; // [nprofiles][nchannels][nlevels]
    double * usercfrac_k; // [nprofiles][nchannels]
    double * gases_k; // [ngases][nprofiles][nchannels][nlevels]
    double * surfemis_k; //[nprofiles][nchannels]
    double * bt_k; //[nprofiles][nchannels]
    double * zef_k; //[nprofiles][nchannels][nlevels]
    std::vector <double> convertPointer4D2Vector(double  ptr[],int x, int y,int z, int dim1, int dim2, int dim3, int dim4);
    std::vector <double> convertPointer3D2Vector(double  ptr[],int x, int y, int dim1, int dim2, int dim3);
    std::vector <double> convertPointer2D2Vector(double  ptr[],  int x,  int dim2, int dim3);
    std::vector <int>    convertIntPointer2D2Vector(int  ptr[],  int x,  int dim2, int dim3);
    bool debug;

    void doStoreRad();
    void doStoreZef();
    void doStoreEmisTerms();
};
} /* namespace rttov */
#endif /* RTTOVSCATT_H_ */
