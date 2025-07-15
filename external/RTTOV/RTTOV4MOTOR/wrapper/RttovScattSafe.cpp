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

///
/// @file RttovScattSafe.cpp
///
/// @brief Wrapper class for RTTOV-SCATT
///

#include <RttovScattSafe.h>

namespace rttov {

///@brief RttovScattSafe class constructor method
RttovScattSafe::RttovScattSafe() : RttovScatt::RttovScatt(), options(RttovScatt::options), profilesP(NULL) {

}

/// @brief Deallocate internal arrays and drop the coefficients
RttovScattSafe::~RttovScattSafe() {
    delete [] p;
    delete [] ph;
    delete [] t;
    delete [] usercfrac;
    delete [] datetimes;
    delete [] angles;
    delete [] surftype;
    delete [] surfgeom;
    delete [] s2m;
    delete [] skin;
    delete [] zeeman;
    delete [] gases;
    delete [] gas_id;
}

///@brief Return the coefficient filename
const string& RttovScattSafe::getFileCoef() const {
    return RttovScatt::getFileCoef();
}

///@brief Set the coefficient filename
///@param fileCoef string containing full path to "rtcoef" coefficient file
void RttovScattSafe::setFileCoef(const string& fileCoef) {
    RttovScatt::setFileCoef(fileCoef);
}

///@brief Return the hydrotable filename
const string& RttovScattSafe::getFileHydrotable() const {
    return RttovScatt::getFileHydrotable();
}

///@brief Set the hydrotable filename
///@param fileHydrotable string containing full path to hydrotable file
void RttovScattSafe::setFileHydrotable(const string& fileHydrotable) {
    RttovScatt::setFileHydrotable(fileHydrotable);
}

///@brief Return the calc_zef boolean
bool RttovScattSafe::isCalcZef() const {
    return RttovScatt::isCalcZef();
}

///@brief Set the calc_zef boolean
///@param calcZef boolean value
void RttovScattSafe::setCalcZef(bool calcZef) {
    RttovScatt::setCalcZef(calcZef);
}

///@brief Return the multi_hydro_frac boolean
bool RttovScattSafe::isMultiHydroFrac() const {
    return RttovScatt::isMultiHydroFrac();
}

///@brief Set the multi_hydro_frac boolean
///@param multiHydroFrac boolean value
void RttovScattSafe::setMultiHydroFrac(bool multiHydroFrac) {
    RttovScatt::setMultiHydroFrac(multiHydroFrac);
}

///@brief Load instrument with all channels
void RttovScattSafe::loadInst() {
    RttovScatt::loadInst();
}

///@brief Return the inst_id
int RttovScattSafe::getInstId() const {
    return RttovScatt::getInstId();
}

///@brief Return true if instrument is loaded
bool RttovScattSafe::isCoeffsLoaded() const {
    return RttovScatt::isCoeffsLoaded();
}

///@brief Return the number of loaded channels
int RttovScattSafe::getNchannels() const {
    return RttovScatt::getNchannels();
}

///@brief Return the channel central wavenumbers of the coefficient file
double* RttovScattSafe::getWaveNumbers() {
    return RttovScatt::getWaveNumbers();
}

///@brief Return the number of levels of the coefficient file
int RttovScattSafe::getCoeffsNlevels() {
    return RttovScatt::getCoeffsNlevels();
}

///@brief Return true if profiles have been associated
bool RttovScattSafe::isProfileSet() const {
    return RttovScatt::isProfileSet();
}

///@brief Return the number of associated profiles
int RttovScattSafe::getNprofiles() const {
    return RttovScatt::getNprofiles();
}

///@brief Update RTTOV options for the currently loaded instrument
void RttovScattSafe::updateOptions() {
    RttovScatt::updateOptions();
}

///@brief Print RTTOV options for the currently loaded instrument
void RttovScattSafe::printOptions() {
    RttovScatt::printOptions();
}

///@brief Set pointer to array containing input/output surface emissivity values;
///this must be previously allocated a double array of dimensions [nprofiles][nchannels]; this
///is used to pass emissivity values into RTTOV-SCATT; if this is not called the RttovScattSafe object
///will allocate an array containing the values used by RTTOV-SCATT which can be accessed by getSurfEmis
///@param surfemis pointer to a previously allocated double array [nprofiles][nchannels]
void RttovScattSafe::setSurfEmis(double* surfemis) {
    return RttovScatt::setSurfEmis(surfemis);
}

///@brief Return the associated profiles
std::vector<rttov::ProfileScatt> * RttovScattSafe::getTheProfiles() const {
    if (this->isProfileSet()) {
        return profilesP;
    } else {
        return NULL;
    }
}

/// @brief Associate a vector of ProfileScatt objects with this RttovScattSafe object;
///        carries out checks on profiles before calling RTTOV-SCATT to help prevent errors:
///        all profiles must be have the same number of levels with the same content (gases, hydrometeors)
///        and have the same gas_units
void RttovScattSafe::setTheProfiles(std::vector<rttov::ProfileScatt>& theProfiles) {
    profilesP = &theProfiles;

    std::vector<rttov::ProfileScatt>& profiles = theProfiles;
    int nbProfiles=profiles.size();
//     if (this->options.isVerboseWrapper()) std::cerr << "Number of profiles : " << nbProfiles << std::endl;
    if (nbProfiles == 0) throw logic_error("Error: no Profile instances in the vector of profiles");

    int nbLevels=profiles[0].getNlevels();
//     if (this->options.isVerboseWrapper()) std::cerr << "Number of levels : " << nbLevels << std::endl;
    std::vector < rttov::itemIdType> Profilecontents = profiles[0].getItemContents();
    this->ngases=profiles[0].getItemContents().size();
    if (this->debug) std::cerr << "ngases= " << this->ngases << std::endl;
    bool isOK;
    isOK=profiles[0].check();
    if (! isOK ) throw logic_error("Error: profile 0 is not OK");

    this->profileSet=false;

    this->gas_units=profiles[0].getGasUnits();
    for (int i=1; i< nbProfiles;i++) {
        if (profiles[i].getNlevels() != nbLevels)
            throw logic_error("Error: Profiles must all have the same number of levels");
        if (profiles[i].getItemContents() != Profilecontents )
            throw logic_error("Error: Profiles must all have the same gas/cloud/aerosol data defined");
        if ( ! profiles[i].check() )
            throw logic_error("Error: Profile check failed");
        if (  profiles[i].getGasUnits() != gas_units )
            throw logic_error("Error: Profiles must all have the same gas_units");
    }
    if ( this->gas_units == unknown ) {
        std::cerr << "Profile warning: gas_units not specified, assuming ppmv wet" << std::endl;
        this->gas_units = ppmv_wet;
    }

    this->nprofiles = nbProfiles;
    this->nlevels = nbLevels;
    // convert to pointers to doubles for the wrapper
    std::vector <double> vOfd;
    std::vector <int> vOfi;
    delete [] p;
    p=new double[nbProfiles*nbLevels]();
    delete [] ph;
    ph=new double[nbProfiles*(nbLevels+1)]();
    delete [] t;
    t=new double[nbProfiles*nbLevels]();
    delete [] usercfrac;
    usercfrac=new double[nbProfiles]();
    delete [] datetimes;
    datetimes=new int[nbProfiles*6]();
    delete [] angles;
    angles=new double[nbProfiles*2]();
    delete [] surftype;
    surftype=new int[nbProfiles]();
    delete [] surfgeom;
    surfgeom=new double[nbProfiles*3]();
    delete [] s2m;
    s2m=new double[nbProfiles*5]();
    delete [] skin;
    skin=new double[nbProfiles*8]();
    delete [] zeeman;
    zeeman=new double[nbProfiles*2]();
    delete [] gases;
    gases=new double[this->ngases*nbProfiles*nbLevels]();
    delete [] gas_id;
    gas_id=new int[this->ngases]();

    for (int i=0; i< nbProfiles; i++) {

        vOfd=profiles[i].getP();
        for (int j=0; j <vOfd.size() ; j++) {
            p[nbLevels*i+j]=vOfd[j];
        }

        vOfd=profiles[i].getPh();
        for (int j=0; j <vOfd.size() ; j++) {
            ph[(nbLevels+1)*i+j]=vOfd[j];
        }

        vOfd=profiles[i].getT();
        for (int j=0; j <vOfd.size() ; j++) {
            t[nbLevels*i+j]=vOfd[j];
        }

        usercfrac[i]=profiles[i].getUserCfrac();

        vOfi=(profiles[i].getDateTimes());
        for (int j=0; j <vOfi.size() ; j++) {
            datetimes[6*i+j]=vOfi[j];
        }

        vOfd=profiles[i].getAngles();
        for (int j=0; j <vOfd.size() ; j++) {
            angles[2*i+j]=vOfd[j];
        }

        surftype[i]=profiles[i].getSurfType();

        vOfd=profiles[i].getSurfGeom();
        for (int j=0; j <vOfd.size() ; j++) {
            surfgeom[3*i+j]=vOfd[j];
        }

        vOfd=profiles[i].getS2m();
        for (int j=0; j <vOfd.size() ; j++) {
            s2m[5*i+j]=vOfd[j];
        }

        vOfd=profiles[i].getSkin();
        for (int j=0; j <vOfd.size() ; j++) {
            skin[8*i+j]=vOfd[j];
        }

        vOfd=profiles[i].getZeeman();
        for (int j=0; j <vOfd.size() ; j++) {
            zeeman[2*i+j]=vOfd[j];
        }
    }
    rttov::ItemIdPointerMap myItems;
    rttov::itemIdType key;
    int i=0;
    for (int prof=0; prof< nbProfiles; prof++) {
        i=0;
        if (this->debug) std::cerr << "Profile " << prof << std::endl;
        myItems=profiles[prof].getItems();
        for (rttov::ItemIdPointerMap::iterator iter = myItems.begin(); iter != myItems.end(); ++iter) {
            key = iter->first;
            if ( myItems.count(key) > 0 ) {
                if (prof==0) this->gas_id[i]=key;
                for (int lev=0;lev<nbLevels;lev++) {
                    gases[i*nbProfiles*nbLevels+nbLevels*prof+lev]=myItems[key][lev];
                    gas_index[key] = i;
                }
            }
            i++;
        }
    }

    this->profileSet=true;
    //this->profiles=profiles;
    if (this->debug) {
        for (int i=0;i< this->ngases;i++) std::cerr << "gas_id : " << this->gas_id[i] << std::endl;
    }
}

///@brief Print gases array contents on standard output
void rttov::RttovScattSafe::printGases() {
    RttovScatt::printGases();
}

///@brief Run the RTTOV-SCATT direct model for a list of channels
///@param channels vector containing the list of channels to simulate
void RttovScattSafe::runDirect(const vector<int>& channels) {
    RttovScatt::runDirect(channels);
}

///@brief Run the RTTOV-SCATT direct model for all channels
void RttovScattSafe::runDirect() {
    RttovScatt::runDirect();
}

///@brief Run the RTTOV-SCATT K model for a list of channels
///@param channels vector containing the list of channels to simulate
void RttovScattSafe::runK(const vector<int>& channels) {
    RttovScatt::runK(channels);
}

///@brief Run the RTTOV-SCATT K model for all channels
void RttovScattSafe::runK() {
    RttovScatt::runK();
}

///@brief Return a pointer to an array of dimensions [nprofiles][nchannels]
/// containing output values of surface emissivity;
/// this array can be initialised by the user and set by calling the setSurfEmis method;
/// alternatively if the emissivity array is allocated by the RttovScattSafe object
/// it is deleted at the next run or when the RttovScattSafe instance is destroyed.
const double* RttovScattSafe::getSurfEmis() const {
    return RttovScatt::getSurfEmis();
}

///@brief Return vector of brightness temperatures
///       computed by the previous run for the given profile number
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovScattSafe::getBt(const int profile) {
    return RttovScatt::getBt(profile);
}

///@brief Return the computed pressure Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScattSafe::getPK(int profile, int channel) {
    return RttovScatt::getPK(profile,channel);
}

///@brief Return the computed pressure half-level Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScattSafe::getPhK(int profile, int channel) {
    return RttovScatt::getPhK(profile,channel);
}

///@brief Return computed temperature Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1
std::vector<double> RttovScattSafe::getTK(int profile, int channel) {
    return RttovScatt::getTK(profile,channel);
}

///@brief Return vector of user cfrac Jacobians
///       computed by the previous run for the given profile number
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovScattSafe::getUserCfracK(const int profile) {
    return RttovScatt::getUserCfracK(profile);
}

///@brief Return computed skin variable Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScattSafe::getSkinK(int profile, int channel) {
    return RttovScatt::getSkinK(profile,channel);
}

///@brief Return computed 2m variable Jacobian for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScattSafe::getS2mK(int profile, int channel) {
    return RttovScatt::getS2mK(profile,channel);
}

///@brief Return computed gas, cloud and aerosol Jacobian values for a given profile and channel
///@param itemIdType item type
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScattSafe::getItemK(rttov::itemIdType itemIdType,
        int profile, int channel) {
    return this->convertPointer4D2Vector(this->gases_k,
            gas_index[itemIdType],
            profile,
            channel,
            this->ngases,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

///@brief Return computed surface emissivity Jacobians for a given profile
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovScattSafe::getSurfEmisK(int profile) {
    return RttovScatt::getSurfEmisK(profile);
}

///@brief Return RTTOV radiance bt_clear output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> RttovScattSafe::getBtClear(int profile) {
    return RttovScatt::getBtClear(profile);}

///@brief Return RTTOV radiance radiance quality output array of size [nchannels] for given profile, requires store_rad true
std::vector <int> RttovScattSafe::getRadQuality(int profile) {
    return RttovScatt::getRadQuality(profile);}

///@brief Return RTTOV radiance plane_parallel flag, requires store_rad true
bool RttovScattSafe::getPlaneParallel() {
    return RttovScatt::getPlaneParallel();}

///@brief Return RTTOV radiance geometric_height output array of size [nlevels] for given profile and channel, requires store_rad true
std::vector <double> RttovScattSafe::getGeometricHeight(int profile, int channel) {
    return RttovScatt::getGeometricHeight(profile, channel);}

///@brief Return RTTOV-SCATT radar reflectivity zef output array of size [nlevels] for given profile and channel
std::vector <double> RttovScattSafe::getZef(int profile, int channel) {
    return RttovScatt::getZef(profile, channel);
}

///@brief Return RTTOV-SCATT radar reflectivity azef output array of size [nlevels] for given profile and channel
std::vector <double> RttovScattSafe::getAZef(int profile, int channel) {
    return RttovScatt::getAZef(profile, channel);
}

///@brief Return RTTOV-SCATT emis retrieval cfrac output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScattSafe::getEmisTermsCfrac(int profile) {
    return RttovScatt::getEmisTermsCfrac(profile);}

///@brief Return RTTOV-SCATT emis retrieval bsfc output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScattSafe::getEmisTermsBsfc(int profile) {
    return RttovScatt::getEmisTermsBsfc(profile);}

///@brief Return RTTOV-SCATT emis retrieval tau_cld output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScattSafe::getEmisTermsTauCld(int profile) {
    return RttovScatt::getEmisTermsTauCld(profile);}

///@brief Return RTTOV-SCATT emis retrieval up_cld output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScattSafe::getEmisTermsUpCld(int profile) {
    return RttovScatt::getEmisTermsUpCld(profile);}

///@brief Return RTTOV-SCATT emis retrieval down_cld output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScattSafe::getEmisTermsDownCld(int profile) {
    return RttovScatt::getEmisTermsDownCld(profile);}

///@brief Return RTTOV-SCATT emis retrieval tau_clr output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScattSafe::getEmisTermsTauClr(int profile) {
    return RttovScatt::getEmisTermsTauClr(profile);}

///@brief Return RTTOV-SCATT emis retrieval up_clr output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScattSafe::getEmisTermsUpClr(int profile) {
    return RttovScatt::getEmisTermsUpClr(profile);}

///@brief Return RTTOV-SCATT emis retrieval down_clr output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScattSafe::getEmisTermsDownClr(int profile) {
    return RttovScatt::getEmisTermsDownClr(profile);}

} /* namespace rttov */
