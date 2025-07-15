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

/// @class RttovScattSafe
///  This class contains a set of functions for running RTTOV-SCATT.
///  An instance of a RttovScattSafe Class is defined to handle one instrument.
///  It contains an Options object which allows the Rttov options to
///  be specified and can be associated with a vector of ProfileScatt objects.
///  ProfileScatt objects are a safe way to initialise the profile values
///  therefore the RttovScattSafe presents a safe access to RTTOV-SCATT.

#ifndef RTTOVSCATTSAFE_H_
#define RTTOVSCATTSAFE_H_

#include <RttovScatt.h>
#include <ProfileScatt.h>
namespace rttov {

class RttovScattSafe: private RttovScatt {
  public:
    RttovScattSafe();
    virtual ~RttovScattSafe();

    rttov::Options &options;

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
    rttov::ProfilesScatt * getProfiles() const;

    void updateOptions();
    void printOptions();

    void setSurfEmis(double* surfemis);

    std::vector <rttov::ProfileScatt> * getTheProfiles() const;
    void setTheProfiles(std::vector <rttov::ProfileScatt>& theProfiles);
    ItemIdIndexMap gas_index;
    void printGases();

    void runDirect(); // run direct for all loaded channels
    void runDirect(const vector<int>& channels);
    void runK();
    void runK(const vector<int>& channels);

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

protected:

    std::vector<rttov::ProfileScatt> * profilesP;
};

} /* namespace rttov */

#endif /* RTTOVSCATTSAFE_H_ */
