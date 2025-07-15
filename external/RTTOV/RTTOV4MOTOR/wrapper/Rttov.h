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

/// @class Rttov
///  This class contains a set of functions for running Rttov.
///  An instance of an Rttov Class is defined to handle one instrument.
///  It contains an Options object which allows the Rttov options to
///  be specified and can be associated with a Profiles object.


#ifndef RTTOV_H_
#define RTTOV_H_


using namespace std;
#include <rttov_cc_interface.h>

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <exception>
#include <stdexcept>
#include <Rttov_common.h>
#include <Profiles.h>
#include <Profile.h>
#include <Options.h>

namespace rttov {

class Rttov {

public:
    Rttov();
    virtual ~Rttov();

    rttov::Options options;

    const string& getFileCoef() const;
    const string& getFileSccld() const;
    const string& getFileScaer() const;
    const string& getFileMfasisCld() const;
//     const string& getFileMfasisAer() const;

    void setFileCoef(const string& fileCoef);
    void setFileSccld(const string& fileSccld);
    void setFileScaer(const string& fileScaer);
    void setFileMfasisCld(const string& fileMfasisCld);
//     void setFileMfasisAer(const string& fileMfasisAer);

    void loadInst(); // load instrument for all channels
    void loadInst(const vector<int>& channels);
    int getInstId() const;

    bool isCoeffsLoaded() const;
    int getNchannels() const;
    double * getRefPressures();
    int getCoeffsNlevels();
    double * getWaveNumbers();
    bool isProfileSet() const;
    int getNprofiles() const;

    void updateOptions();
    void printOptions();

    void setSurfEmisRefl(double* surfemisrefl);

    void setAerAsb(double* asb);
    void setAerPha(int nphangle, double* phangle, double* pha);
    void setAerLegcoef(int nmom, double* legcoef);
    void setCldAsb(double* asb);
    void setCldPha(int nphangle, double* phangle, double* pha);
    void setCldLegcoef(int nmom, double* legcoef);

    rttov::Profiles * getProfiles() const;
    void setProfiles(rttov::Profiles * profiles);
    void printGases();

    void runDirect(); // run direct for all loaded channels
    void runDirect(const vector<int>& channels);
    void runK() ;
    void runK(const vector<int>& channels);

    const double* getBtRefl() const;
    const double* getRads() const;
    std::vector<double> getBtRefl(const int profile);
    std::vector<double> getRads(const int profile);
    const double* getSurfEmisRefl() const;

    int getAerNphangle() const;
    int getAerNmom() const;
    const double* getAerAsb() const;
    const double* getAerPhangle() const;
    const double* getAerLegcoef() const;
    const double* getAerPha() const;
    int getCldNphangle() const;
    int getCldNmom() const;
    const double* getCldAsb() const;
    const double* getCldPhangle() const;
    const double* getCldLegcoef() const;
    const double* getCldPha() const;

    double calcBpr(int nphangle, double* phangle, double* pha);
    void calcLegcoef(int nphangle, double* phangle, double* pha, int nmom, double* legcoef, int ngauss);

    std::vector <double> getPK(int profile, int channel);
    std::vector <double> getTK(int profile, int channel);
    std::vector <double> getSkinK(int profile, int channel);
    std::vector <double> getS2mK(int profile, int channel);
    std::vector <double> getSimpleCloudK(int profile, int channel);
    std::vector <double> getItemK(rttov::itemIdType, int profile, int channel);
    std::vector <double> getSurfEmisK(int profile);
    std::vector <double> getSurfReflK(int profile);
    std::vector <double> getSurfDiffuseReflK(int profile);
    std::vector <double> getSpecularityK(int profile);

    std::vector <double> getTauTotal(int profile);
    std::vector <double> getTauLevels(int profile, int channel);
    std::vector <double> getTauSunTotalPath1(int profile);
    std::vector <double> getTauSunLevelsPath1(int profile, int channel);
    std::vector <double> getTauSunTotalPath2(int profile);
    std::vector <double> getTauSunLevelsPath2(int profile, int channel);
    std::vector <double> getTauTotalCld(int profile);
    std::vector <double> getTauLevelsCld(int profile, int channel);

    std::vector <double> getRadClear(int profile);
    std::vector <double> getRadTotal(int profile);
    std::vector <double> getBtClear(int profile);
    std::vector <double> getBt(int profile);
    std::vector <double> getReflClear(int profile);
    std::vector <double> getRefl(int profile);
    std::vector <double> getRadCloudy(int profile);
    std::vector <double> getOvercast(int profile, int channel);
    std::vector <int> getRadQuality(int profile);
    bool getPlaneParallel();
    std::vector <double> getGeometricHeight(int profile, int channel);

    std::vector <double> getRad2UpClear(int profile);
    std::vector <double> getRad2DnClear(int profile);
    std::vector <double> getRad2ReflDnClear(int profile);
    std::vector <double> getRad2Up(int profile, int channel);
    std::vector <double> getRad2Down(int profile, int channel);
    std::vector <double> getRad2Surf(int profile, int channel);

protected :

    rttov::Profiles * profiles;
    string strOptions;
    int nprofiles;
    int nlevels;
    int nchannels; // number of channels loaded from coefficients
    int nchannelsForLastRun;  // number of channels used for the last run
    int ngases;
    int gas_units;
    int mmr_cldaer;
    int inst_id;
    double * p;//[nprofiles][nlevels];                  // Input pressure profiles
    double * t;//[nprofiles][nlevels];                  // Input temperature profiles
    double * gases;//[ngases][nprofiles][nlevels];      // Input gas profiles
    double * surfemisrefl;//[4][nprofiles][nchannels];  // Input/output surface emissivities and BRDFs
    double * btrefl;//[nprofiles][nchannels];           // Output BTs/refls (for thermal/solar chans)
    double * rads;//[nprofiles][nchannels];             // Output radiances

    double * tautotal;
    double * taulevels;
    double * tausuntotalpath1;
    double * tausunlevelspath1;
    double * tausuntotalpath2;
    double * tausunlevelspath2;
    double * tautotalcld;
    double * taulevelscld;

    double * radclear;
    double * radtotal;
    double * btclear;
    double * bt;
    double * reflclear;
    double * refl;
    double * radcloudy;
    double * overcast;
    int * radquality;
    int plane_parallel;
    double * geometric_height;

    double * rad2upclear;
    double * rad2dnclear;
    double * rad2refldnclear;
    double * rad2up;
    double * rad2down;
    double * rad2surf;

    int * gas_id;// {gas_id_q, gas_id_co2, gas_id_cfrac, gas_id_lwc1, gas_id_iwc};

    // datetimes: yy, mm, dd, hh, mm, ss
    int * datetimes; //[nprofiles][6]

    // angles: satzen, satazi, sunzen, sunazi
    double * angles;//[nprofiles][4]

    // surftype: surftype, watertype
    int * surftype;//[nprofiles][2]

    // surfgeom: lat, lon, elev
    double * surfgeom ;//[nprofiles][3]

    // s2m: 2m p, 2m t, 2m q, 10m wind u, v, wind fetch
    double * s2m;//[nprofiles][6]

    // skin: skin T, salinity, snow_frac, foam_frac, fastem_coefsx5
    double * skin; //[nprofiles][9]

    // simplecloud: ctp, cfraction
    double * simplecloud ;//[nprofiles][2]

    // clwscheme: clw_scheme, clwde_param
    int * clwscheme ;//[nprofiles][2]

    // icecloud: ice_scheme, icede_param
    int * icecloud ;//[nprofiles][2]

    // zeeman: be, cosbk
    double * zeeman;//[nprofiles][2]

    int aer_nphangle;
    int aer_nmom;
    double * aer_asb;
    double * aer_phangle;
    double * aer_legcoef;
    double * aer_pha;
    int cld_nphangle;
    int cld_nmom;
    double * cld_asb;
    double * cld_phangle;
    double * cld_legcoef;
    double * cld_pha;
    bool allocatedAerAsb;
    bool allocatedAerPha;
    bool allocatedAerLegcoef;
    bool allocatedCldAsb;
    bool allocatedCldPha;
    bool allocatedCldLegcoef;
    void setOptParams(int nchannels);

    string file_coef;
    string file_sccld;
    string file_scaer;
    string file_mfasis_cld;
    string file_mfasis_aer;
    string file_pccoef;
    string file_hydrotable;

    bool coeffsLoaded;
    bool profileSet;

    double * wavenumbers;
    double * refPressures;
    bool allocatedSurfemisrefl;

    double * skin_k; //skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][nchannels][9]
    double * s2m_k; //2m p, 2m t, 2m q, 10m wind u, v, wind-fetch    [nprofiles][nchannels][6]
    double * simplecloud_k ; // [nprofiles][nchannels][2]
    double * p_k; // [nprofiles][nchannels][nlevels]
    double * t_k; // [nprofiles][nchannels][nlevels]
    double * gases_k; // [ngases][nprofiles][nchannels][nlevels]
    double * surfemisrefl_k; //[4][nprofiles][nchannels]
    double * bt_k; //[nprofiles][nchannels]
    double * rad_k; //[nprofiles][nchannels]
    bool allocatedSurfemisreflk;
    bool allocatedBtk;
    bool allocatedRadk;
    std::vector <double> convertPointer4D2Vector(double  ptr[],int x, int y,int z, int dim1, int dim2, int dim3, int dim4);
    std::vector <double> convertPointer3D2Vector(double  ptr[],int x, int y, int dim1, int dim2, int dim3);
    std::vector <double> convertPointer2D2Vector(double  ptr[],  int x,  int dim2, int dim3);
    std::vector <int> convertIntPointer2D2Vector(int  ptr[],  int x,  int dim2, int dim3);
    bool debug;
    bool allocatedP;

    void doStoreTrans();
    void doStoreRad();
    void doStoreRad2();
};
} /* namespace rttov */
#endif /* RTTOV_H_ */
