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
/// @file RttovSafe.cpp
///
/// @brief Wrapper class for rttov
///

#include <RttovSafe.h>

namespace rttov {

///@brief RttovSafe class constructor method
RttovSafe::RttovSafe() : Rttov::Rttov(), options(Rttov::options),profilesP(NULL) {

}

/// @brief Deallocate internal arrays and drop the coefficients
RttovSafe::~RttovSafe() {
    delete [] p;
    this->allocatedP = false;
    delete [] t;
    delete [] datetimes;
    delete [] angles;
    delete [] surftype;
    delete [] surfgeom;
    delete [] s2m;
    delete [] skin;
    delete [] simplecloud;
    delete [] clwscheme;
    delete [] icecloud;
    delete [] zeeman;
    delete [] gases;
    delete [] gas_id;
}

///@brief Return the coefficient filename
const string& RttovSafe::getFileCoef() const {
    return Rttov::getFileCoef();
}

///@brief Set the coefficient filename
///@param fileCoef string containing full path to "rtcoef" coefficient file
void RttovSafe::setFileCoef(const string& fileCoef) {
    Rttov::setFileCoef(fileCoef);
}

///@brief Return the cloud coefficient filename
const string& RttovSafe::getFileSccld() const {
    return Rttov::getFileSccld();
}

///@brief Set the cloud coefficient filename
///@param fileSccld string containing full path to cloud coefficient file
void RttovSafe::setFileSccld(const string& fileSccld) {
    Rttov::setFileSccld(fileSccld);
}

///@brief Return the aerosol coefficient filename
const string& RttovSafe::getFileScaer() const {
    return Rttov::getFileScaer();
}

///@brief Set the aerosol coefficient filename
///@param fileScaer string containing full path to aerosol coefficient file
void RttovSafe::setFileScaer(const string& fileScaer) {
    Rttov::setFileScaer(fileScaer);
}

///@brief Return the MFASIS cloud LUT filename
const string& RttovSafe::getFileMfasisCld() const {
    return Rttov::getFileMfasisCld();
}

///@brief Set the MFASIS cloud LUT filename
///@param fileMfasisCld string containing full path to cloud coefficient file
void RttovSafe::setFileMfasisCld(const string& fileMfasisCld) {
    Rttov::setFileMfasisCld(fileMfasisCld);
}

/////@brief Return the MFASIS aerosol LUT filename
// const string& RttovSafe::getFileMfasisAer() const {
//     return Rttov::getFileMfasisAer();
// }

/////@brief Set the MFASIS aerosol LUT filename
/////@param fileMfasisAer string containing full path to aerosol coefficient file
// void RttovSafe::setFileMfasisAer(const string& fileMfasisAer) {
//     Rttov::setFileMfasisAer(fileMfasisAer);
// }



///@brief Load instrument with all channels
void RttovSafe::loadInst() {
    Rttov::loadInst();
}

///@brief Load instrument for a list of channels;
///       the method setFileCoef() must have been called previously
///@param channels vector containing the list of channels to be loaded
void RttovSafe::loadInst(const vector<int>& channels) {
    Rttov::loadInst(channels);
}

///@brief Return the inst_id
int RttovSafe::getInstId() const {
    return Rttov::getInstId();
}

///@brief Return true if instrument is loaded
bool RttovSafe::isCoeffsLoaded() const {
    return Rttov::isCoeffsLoaded();
}

///@brief Return the number of loaded channels
int RttovSafe::getNchannels() const {
    return Rttov::getNchannels();
}

///@brief Return the channel central wavenumbers of the coefficient file
double* RttovSafe::getWaveNumbers() {
    return Rttov::getWaveNumbers();
}

///@brief Return the number of levels of the coefficient file
int RttovSafe::getCoeffsNlevels() {
    return Rttov::getCoeffsNlevels();
}

///@brief Return true if profiles have been associated
bool RttovSafe::isProfileSet() const {
    return Rttov::isProfileSet();
}

///@brief Return the number of associated profiles
int RttovSafe::getNprofiles() const {
    return Rttov::getNprofiles();
}


///@brief Update RTTOV options for the currently loaded instrument
void RttovSafe::updateOptions() {
    Rttov::updateOptions();
}

///@brief Print RTTOV options for the currently loaded instrument
void RttovSafe::printOptions() {
    Rttov::printOptions();
}



///@brief Set pointer to array containing input/output surface emissivity and reflectance values;
///this must be previously allocated a double array of dimensions [4][nprofiles][nchannels]; this
///is used to pass emissivity/reflectance values into RTTOV; if this is not called the RttovSafe object
///will allocate an array containing the values used by RTTOV which can be accessed by getSurfEmisRefl
///@param surfemisrefl pointer to a previously allocated double array [4][nprofiles][nchannels]
void RttovSafe::setSurfEmisRefl(double* surfemisrefl) {
    return Rttov::setSurfEmisRefl(surfemisrefl);
}


///@brief Return the associated profiles
std::vector<rttov::Profile> * RttovSafe::getTheProfiles() const {
    if (this->isProfileSet()) {
        return profilesP;
    } else {
        return NULL;
    }
}

/// @brief Associate a vector of Profile objects with this RttovSafe object; 
///        carries out checks on profiles before calling RTTOV to help prevent errors:
///        all profiles must be have the same number of levels with the same content (gases, clouds, aerosols)
///        and have the same gas_units
void RttovSafe::setTheProfiles(std::vector<rttov::Profile>& theProfiles) {
    profilesP = &theProfiles;

    std::vector<rttov::Profile>& profiles = theProfiles;
    int nbProfiles=profiles.size();
//     if (this->options.isVerboseWrapper()) std::cerr << "Number of profiles : " << nbProfiles << std::endl;
    if (nbProfiles == 0) throw logic_error("Error: no Profile instances in the vector of profiles");

    if (allocatedP) delete[] p;
    // check if pressure is set for the first profile
    if (profiles[0].isDefaultPressureLevels() )
        if (this->coeffsLoaded ) {
            int nlevcoeff=this->getCoeffsNlevels();
            if (nlevcoeff == profiles[0].getNlevels()) {
                if (this->options.isVerboseWrapper()) std::cerr << "Using pressure levels from coefficient file" << std::endl;
                double * pressure=this->getRefPressures();
                delete [] p;
                p= new double[nlevcoeff*nbProfiles];
                for (int prof=0;prof<nbProfiles;prof++)
                    for (int l=0;l<nlevcoeff;l++) p[prof*nlevcoeff+l]=pressure[l];
                allocatedP=true;
            }
            else throw logic_error("Error: number of levels differs to coefficient file, cannot use coefficient file pressure levels");
        }
        else throw logic_error("Error: instrument not loaded, cannot use coefficient file pressure levels");

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
    this->mmr_cldaer=profiles[0].isMmrCldAer();
    for (int i=1; i< nbProfiles;i++) {
        if (profiles[i].getNlevels() != nbLevels)
            throw logic_error("Error: Profiles must all have the same number of levels");
        if (profiles[i].getItemContents() != Profilecontents )
            throw logic_error("Error: Profiles must all have the same gas/cloud/aerosol data defined");
        if ( ! profiles[i].check() )
            throw logic_error("Error: Profile check failed");
        if (  profiles[i].getGasUnits() != gas_units )
            throw logic_error("Error: Profiles must all have the same gas_units");
        if (  profiles[i].isMmrCldAer() != mmr_cldaer )
            throw logic_error("Error: Profiles must all have the same mmr_cldaer");
    }
    if ( this->gas_units == unknown ) {
        std::cerr << "Profile warning: gas_units not specified, assuming ppmv wet" << std::endl;
        this->gas_units=ppmv_wet;
    }

    this->nprofiles= nbProfiles;
    this->nlevels=nbLevels;
    // convert to pointers to doubles for the wrapper
    std::vector <double> vOfd;
    std::vector <int>  vOfi;
    if ( ! allocatedP) { // default pressure level is used and previously set
       delete [] p;
       p=new double[nbProfiles*nbLevels]();
       for (int p=0; p< nbProfiles; p++) {

               vOfd=profiles[p].getP();
               for (int j=0; j <vOfd.size() ; j++) {
                   this->p[nbLevels*p+j]=vOfd[j];
               }
       }
    }
    delete [] t;
    t=new double[nbProfiles*nbLevels]();
    delete [] datetimes;
    datetimes=new int[nbProfiles*6]();
    delete [] angles;
    angles=new double[nbProfiles*4]();
    delete [] surftype;
    surftype=new int[nbProfiles*2]();
    delete [] surfgeom;
    surfgeom=new double[nbProfiles*3]();
    delete [] s2m;
    s2m=new double[nbProfiles*6]();
    delete [] skin;
    skin=new double[nbProfiles*9]();
    delete [] simplecloud;
    simplecloud=new double[nbProfiles*2]();
    delete [] clwscheme;
    clwscheme=new int[nbProfiles*2]();
    delete [] icecloud;
    icecloud=new int[nbProfiles*2]();
    delete [] zeeman;
    zeeman=new double[nbProfiles*2]();
    delete [] gases;
    gases=new double[this->ngases*nbProfiles*nbLevels]();
    delete [] gas_id;
    gas_id=new int[this->ngases]();

    for (int p=0; p< nbProfiles; p++) {

        vOfd=profiles[p].getT();
        for (int j=0; j <vOfd.size() ; j++) {
            t[nbLevels*p+j]=vOfd[j];
        }

        vOfi=(profiles[p].getDateTimes());
        for (int j=0; j <vOfi.size() ; j++) {
            datetimes[6*p+j]=vOfi[j];
        }

        vOfd=profiles[p].getAngles();
        for (int j=0; j <vOfd.size() ; j++) {
            angles[4*p+j]=vOfd[j];
        }

        vOfi=profiles[p].getSurfType();
        for (int j=0; j <vOfi.size() ; j++) {
            surftype[2*p+j]=vOfi[j];
        }

        vOfd=profiles[p].getSurfGeom();
        for (int j=0; j <vOfd.size() ; j++) {
            surfgeom[3*p+j]=vOfd[j];
        }

        vOfd=profiles[p].getS2m();
        for (int j=0; j <vOfd.size() ; j++) {
            s2m[6*p+j]=vOfd[j];
        }

        vOfd=profiles[p].getSkin();
        for (int j=0; j <vOfd.size() ; j++) {
            skin[9*p+j]=vOfd[j];
        }

        vOfd=profiles[p].getSimpleCloud();
        for (int j=0; j <vOfd.size() ; j++) {
            simplecloud[2*p+j]=vOfd[j];
        }

        vOfi=profiles[p].getClwScheme();
        for (int j=0; j <vOfi.size() ; j++) {
            clwscheme[2*p+j]=vOfi[j];
        }

        vOfi=profiles[p].getIceCloud();
        for (int j=0; j <vOfi.size() ; j++) {
            icecloud[2*p+j]=vOfi[j];
        }

        vOfd=profiles[p].getZeeman();
        for (int j=0; j <vOfd.size() ; j++) {
            zeeman[2*p+j]=vOfd[j];
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
void rttov::RttovSafe::printGases() {
    Rttov::printGases();
}


///@brief Run the RTTOV direct model for a list of channels
///@param channels vector containing the list of channels to simulate
void RttovSafe::runDirect(const vector<int>& channels) {
    Rttov::runDirect(channels);
}

///@brief Run the RTTOV direct model for all channels
void RttovSafe::runDirect() {
    Rttov::runDirect();
}


///@brief Run the RTTOV K model for a list of channels
///@param channels vector containing the list of channels to simulate
void RttovSafe::runK(const vector<int>& channels) {
    Rttov::runK(channels);
}

///@brief Run the RTTOV K model for all channels
void RttovSafe::runK() {
    Rttov::runK();
}

///@brief Set the aerosol absorption coefs, scattering coefs and bpr parameters
///@param asb Abs coefs, sca coefs, bpr parameters [3][nprofiles][nchannels][nlayers]
void RttovSafe::setAerAsb(double* asb) {
    Rttov::setAerAsb(asb);
}

///@brief Set the aerosol phase functions
///@param nphangle Number of angles on which phase functions are specified
///@param phangle  Angle grid spanning 0-180 degrees on which phase functions are specified [nphangle]
///@param pha      Phase functions [nprofiles][nchannels][nlayers][nphangle]
void RttovSafe::setAerPha(int nphangle, double* phangle, double* pha) {
    Rttov::setAerPha(nphangle, phangle, pha);
}

///@brief Set the aerosol phase function Legendre coefficients
///@param nmom    Number of Legendre coefficients per phase function
///@param legcoef Legendre coefficients [nprofiles][nchannels][nlayers][nmom+1]
void RttovSafe::setAerLegcoef(int nmom, double* legcoef) {
    Rttov::setAerLegcoef(nmom, legcoef);
}

///@brief Set the cloud absorption coefs, scattering coefs and bpr parameters
///@param asb Abs coefs, sca coefs, bpr parameters [3][nprofiles][nchannels][nlayers]
void RttovSafe::setCldAsb(double* asb) {
    Rttov::setCldAsb(asb);
}

///@brief Set the cloud phase functions
///@param nphangle Number of angles on which phase functions are specified
///@param phangle  Angle grid spanning 0-180 degrees on which phase functions are specified [nphangle]
///@param pha      Phase functions [nprofiles][nchannels][nlayers][nphangle]
void RttovSafe::setCldPha(int nphangle, double* phangle, double* pha) {
    Rttov::setCldPha(nphangle, phangle, pha);
}

///@brief Set the cloud phase function Legendre coefficients
///@param nmom    Number of Legendre coefficients per phase function
///@param legcoef Legendre coefficients [nprofiles][nchannels][nlayers][nmom+1]
void RttovSafe::setCldLegcoef(int nmom, double* legcoef) {
    Rttov::setCldLegcoef(nmom, legcoef);
}

///@brief Return the number of aerosol phase function angles
int RttovSafe::getAerNphangle() const {
    return Rttov::getAerNphangle();
}

///@brief Return the number of aerosol phase function Legendre coefficients
int RttovSafe::getAerNmom() const {
    return Rttov::getAerNmom();
}

///@brief Return the aerosol absorption coefs, scattering coefs and bpr parameters
const double* RttovSafe::getAerAsb() const {
    return Rttov::getAerAsb();
}

///@brief Return the aerosol phase function angles
const double* RttovSafe::getAerPhangle() const {
    return Rttov::getAerPhangle();
}

///@brief Return the aerosol phase functions
const double* RttovSafe::getAerPha() const {
    return Rttov::getAerPha();
}

///@brief Return the aerosol phase function Legendre coefficients
const double* RttovSafe::getAerLegcoef() const {
    return Rttov::getAerLegcoef();
}

///@brief Return the number of cloud phase function angles
int RttovSafe::getCldNphangle() const {
    return Rttov::getCldNphangle();
}

///@brief Return the number of cloud phase function Legendre coefficients
int RttovSafe::getCldNmom() const {
    return Rttov::getCldNmom();
}

///@brief Return the cloud absorption coefs, scattering coefs and bpr parameters
const double* RttovSafe::getCldAsb() const {
    return Rttov::getCldAsb();
}

///@brief Return the cloud phase function angles
const double* RttovSafe::getCldPhangle() const {
    return Rttov::getCldPhangle();
}

///@brief Return the cloud phase functions
const double* RttovSafe::getCldPha() const {
    return Rttov::getCldPha();
}

///@brief Return the cloudphase function Legendre coefficients
const double* RttovSafe::getCldLegcoef() const {
    return Rttov::getCldLegcoef();
}

///@brief Calculate bpr parameter for given phase function
///@param nphangle Number of angles on which phase function is specified
///@param phangle  Angle grid spanning 0-180 degrees on which phase function is specified [nphangle]
///@param pha      Phase function [nphangle]
double RttovSafe::calcBpr(int nphangle, double* phangle, double* pha) {
    return Rttov::calcBpr(nphangle, phangle, pha);
}

///@brief Calculate Legendre coefficients for given phase function
///@param nphangle Number of angles on which phase function is specified
///@param phangle  Angle grid spanning 0-180 degrees on which phase function is specified [nphangle]
///@param pha      Phase function [nphangle]
///@param nmom     Number of Legendre coefficients to calculate
///@param legcoef  Pointer to pre-allocated array of size [nmom+1] where Legendre coefficients are to be written
///@param ngauss   Size of Gaussian quadrature to use, ignored if less than nmom, default value is 1000
void RttovSafe::calcLegcoef(int nphangle, double* phangle, double* pha, int nmom, double* legcoef, int ngauss) {
    Rttov::calcLegcoef(nphangle, phangle, pha, nmom, legcoef, ngauss);
}


///@brief Return a pointer to an array of dimensions [4][nprofiles][nchannels]
/// containing output values of surface emissivity and reflectance;
/// this array can be initialised by the user and set by calling the setSurfEmisRefl method;
/// alternatively if the emissivity/reflectance array is allocated by the RttovSafe object
/// it is deleted at the next run or when the RttovSafe instance is destroyed.
const double* RttovSafe::getSurfEmisRefl() const {
    return Rttov::getSurfEmisRefl();
}

///@brief Return vector of brightness temperatures/reflectances
///       computed by the previous run for the given profile number
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovSafe::getBtRefl(const int profile) {
    return Rttov::getBtRefl(profile);
}

///@brief Return a vector of radiances computed by the previous run
///       for the given profile number
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovSafe::getRads(const int profile) {
    return Rttov::getRads(profile);
}


///@brief Return the computed pressure Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovSafe::getPK(int profile, int channel) {
    return Rttov::getPK(profile,channel);
}

///@brief Return computed temperature Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1
std::vector<double> RttovSafe::getTK(int profile, int channel) {
    return Rttov::getTK(profile,channel);
}

///@brief Return computed skin variable Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovSafe::getSkinK(int profile, int channel) {
    return Rttov::getSkinK(profile,channel);
}

///@brief Return computed 2m variable Jacobian for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovSafe::getS2mK(int profile, int channel) {
    return Rttov::getS2mK(profile,channel);
}

///@brief Return computed simple cloud variable Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovSafe::getSimpleCloudK(int profile, int channel) {
    return Rttov::getSimpleCloudK(profile,channel);
}

///@brief Return computed gas, cloud and aerosol Jacobian values for a given profile and channel
///@param itemIdType item type
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovSafe::getItemK(rttov::itemIdType itemIdType,
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
std::vector<double> RttovSafe::getSurfEmisK(int profile) {
    return Rttov::getSurfEmisK(profile);
}

///@brief Return computed surface BRDF Jacobians for a given profile
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovSafe::getSurfReflK(int profile) {
    return Rttov::getSurfReflK(profile);
}

///@brief Return computed surface diffuse reflectance Jacobians for a given profile
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovSafe::getSurfDiffuseReflK(int profile) {
    return Rttov::getSurfDiffuseReflK(profile);
}

///@brief Return computed surface specularity Jacobians for a given profile
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovSafe::getSpecularityK(int profile) {
    return Rttov::getSpecularityK(profile);
}


///@brief Return RTTOV transmission tau_total output array of size [nchannels] for given profile, requires store_trans true
std::vector <double> RttovSafe::getTauTotal(int profile) {
    return Rttov::getTauTotal(profile);
}

///@brief Return RTTOV transmission tau_levels output array of size [nlevels] for given profile and channel, requires store_trans true
std::vector <double> RttovSafe::getTauLevels(int profile, int channel) {
    return Rttov::getTauLevels(profile, channel);
}

///@brief Return RTTOV transmission tausun_total_path1 output array of size [nchannels] for given profile, requires store_trans true
std::vector <double> RttovSafe::getTauSunTotalPath1(int profile) {
    return Rttov::getTauSunTotalPath1(profile);
}

///@brief Return RTTOV transmission tausun_levels_path1 output array of size [nlevels] for given profile and channel, requires store_trans true
std::vector <double> RttovSafe::getTauSunLevelsPath1(int profile, int channel) {
    return Rttov::getTauSunLevelsPath1(profile, channel);
}

///@brief Return RTTOV transmission tausun_total_path2 output array of size [nchannels] for given profile, requires store_trans true
std::vector <double> RttovSafe::getTauSunTotalPath2(int profile) {
    return Rttov::getTauSunTotalPath2(profile);
}

///@brief Return RTTOV transmission tausun_levels_path2 output array of size [nlevels] for given profile and channel, requires store_trans true
std::vector <double> RttovSafe::getTauSunLevelsPath2(int profile, int channel) {
    return Rttov::getTauSunLevelsPath2(profile, channel);
}

///@brief Return RTTOV transmission tau_total_cld output array of size [nchannels] for given profile, requires store_trans true
std::vector <double> RttovSafe::getTauTotalCld(int profile) {
    return Rttov::getTauTotalCld(profile);
}

///@brief Return RTTOV transmission tau_levels_cld output array of size [nlevels] for given profile and channel, requires store_trans true
std::vector <double> RttovSafe::getTauLevelsCld(int profile, int channel) {
    return Rttov::getTauLevelsCld(profile, channel);
}


///@brief Return RTTOV radiance clear output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> RttovSafe::getRadClear(int profile) {
    return Rttov::getRadClear(profile);
}

///@brief Return RTTOV radiance total output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> RttovSafe::getRadTotal(int profile) {
    return Rttov::getRadTotal(profile);}

///@brief Return RTTOV radiance bt_clear output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> RttovSafe::getBtClear(int profile) {
    return Rttov::getBtClear(profile);}

///@brief Return RTTOV radiance bt output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> RttovSafe::getBt(int profile) {
    return Rttov::getBt(profile);}

///@brief Return RTTOV radiance refl_clear output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> RttovSafe::getReflClear(int profile) {
    return Rttov::getReflClear(profile);}

///@brief Return RTTOV radiance refl output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> RttovSafe::getRefl(int profile) {
    return Rttov::getRefl(profile);}

///@brief Return RTTOV radiance cloudy output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> RttovSafe::getRadCloudy(int profile) {
    return Rttov::getRadCloudy(profile);}

///@brief Return RTTOV radiance overcast output array of size [nlayers] for given profile and channel, requires store_rad true
std::vector <double> RttovSafe::getOvercast(int profile, int channel) {
    return Rttov::getOvercast(profile, channel);}

///@brief Return RTTOV radiance quality flag array of size [nchannels] for given profile, requires store_rad true
std::vector <int> RttovSafe::getRadQuality(int profile) {
    return Rttov::getRadQuality(profile);}

///@brief Return RTTOV radiance plane_parallel flag, requires store_rad true
bool RttovSafe::getPlaneParallel() {
    return Rttov::getPlaneParallel();}

///@brief Return RTTOV radiance geometric_height output array of size [nlevels] for given profile and channel, requires store_rad true
std::vector <double> RttovSafe::getGeometricHeight(int profile, int channel) {
    return Rttov::getGeometricHeight(profile, channel);}


///@brief Return RTTOV radiance2 upclear output array of size [nchannels] for given profile, requires store_rad2 true
std::vector <double> RttovSafe::getRad2UpClear(int profile) {
    return Rttov::getRad2UpClear(profile);
}

///@brief Return RTTOV radiance2 dnclear output array of size [nchannels] for given profile, requires store_rad2 true
std::vector <double> RttovSafe::getRad2DnClear(int profile) {
    return Rttov::getRad2DnClear(profile);
}

///@brief Return RTTOV radiance2 refldnclear output array of size [nchannels] for given profile, requires store_rad2 true
std::vector <double> RttovSafe::getRad2ReflDnClear(int profile) {
    return Rttov::getRad2ReflDnClear(profile);
}

///@brief Return RTTOV radiance2 up output array of size [nlayers] for given profile and channel, requires store_rad2 true
std::vector <double> RttovSafe::getRad2Up(int profile, int channel) {
    return Rttov::getRad2Up(profile, channel);
}

///@brief Return RTTOV radiance2 down output array of size [nlayers] for given profile and channel, requires store_rad2 true
std::vector <double> RttovSafe::getRad2Down(int profile, int channel) {
    return Rttov::getRad2Down(profile, channel);
}

///@brief Return RTTOV radiance2 surf output array of size [nlayers] for given profile and channel, requires store_rad2 true
std::vector <double> RttovSafe::getRad2Surf(int profile, int channel) {
    return Rttov::getRad2Surf(profile, channel);
}

} /* namespace rttov */
