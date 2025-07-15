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

/// @class RttovSafe
///  This class contains a set of functions for running Rttov.
///  An instance of a RttovSafe Class is defined to handle one instrument.
///  It contains an Options object which allows the Rttov options to
///  be specified and can be associated with a vector of Profile objects.
///  Profile objects are a safe way to initialise the profile values
///  therefore the RttovSafe presents a safe access to Rttov.

#ifndef RTTOVSAFE_H_
#define RTTOVSAFE_H_

#include <Rttov.h>
#include <Profile.h>
namespace rttov {

class RttovSafe: private Rttov {
  public:
    RttovSafe();
    virtual ~RttovSafe();

    rttov::Options &options;

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
    int getCoeffsNlevels();
    double * getWaveNumbers();
    bool isProfileSet() const;
    int getNprofiles() const;
    rttov::Profiles * getProfiles() const;

    void updateOptions();
    void printOptions();

    void setSurfEmisRefl(double* surfemisrefl);

    void setAerAsb(double* asb);
    void setAerPha(int nphangle, double* phangle, double* pha);
    void setAerLegcoef(int nmom, double* legcoef);
    void setCldAsb(double* asb);
    void setCldPha(int nphangle, double* phangle, double* pha);
    void setCldLegcoef(int nmom, double* legcoef);

    std::vector <rttov::Profile> * getTheProfiles() const;
    void setTheProfiles(std::vector <rttov::Profile>& theProfiles);
    ItemIdIndexMap gas_index;
    void printGases();

    void runDirect(); // run direct for all loaded channels
    void runDirect(const vector<int>& channels);
    void runK();
    void runK(const vector<int>& channels);

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

protected:

    std::vector<rttov::Profile> * profilesP;
};

} /* namespace rttov */

#endif /* RTTOVSAFE_H_ */
