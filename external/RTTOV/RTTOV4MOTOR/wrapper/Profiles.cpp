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

///@file Profiles.cpp
///
///@brief Method definition for Profiles class

#include <Profiles.h>

namespace rttov {

///@brief Constructor method for manual gas array creation
///@param nbprofiles number of profiles
///@param nblevels number of vertical levels
///@param nbgases number of gases
Profiles::Profiles(int nbprofiles, const int nblevels, int nbgases)
: nprofiles(nbprofiles), nlevels(nblevels), ngases(nbgases),
    p(NULL),t(NULL),gases(NULL),angles(NULL),
    datetimes(NULL),gas_id(NULL),gas_units(unknown),mmr_cldaer(true),clwscheme(NULL),icecloud(NULL),
    s2m(NULL),simplecloud(NULL),skin(NULL),surfgeom(NULL),surftype(NULL),
    zeeman(NULL),allocatedGases(false) {
    this->initialisePointers();
}

///@brief Constructor method for individual gas specification
///@param nbprofiles number of profiles
///@param nblevels number of vertical levels
Profiles::Profiles(int nbprofiles, const int nblevels)
: nprofiles(nbprofiles), nlevels(nblevels), ngases(0),
    p(NULL),t(NULL),gases(NULL),angles(NULL),
    datetimes(NULL),gas_id(NULL),gas_units(unknown),mmr_cldaer(true),clwscheme(NULL),icecloud(NULL),
    s2m(NULL),simplecloud(NULL),skin(NULL),surfgeom(NULL),surftype(NULL),
    zeeman(NULL),allocatedGases(false) {
    this->initialisePointers();
}

///@brief Constructor method for one profiles and 54 levels
Profiles::Profiles(): nprofiles(1), nlevels(54), ngases(0),
        p(NULL),t(NULL),gases(NULL),angles(NULL),
        datetimes(NULL),gas_id(NULL),gas_units(unknown),mmr_cldaer(true),clwscheme(NULL),icecloud(NULL),
        s2m(NULL),simplecloud(NULL),skin(NULL),surfgeom(NULL),surftype(NULL),
        zeeman(NULL), allocatedGases(false) {
    this->initialisePointers();
}

///@brief Destructor method
Profiles::~Profiles() {
    // deallocate allocated arrays if necessary
    if (this->allocatedPointers["zeeman"]) delete [] this->zeeman;
    if (this->allocatedPointers["clwscheme"])delete [] this->clwscheme;
    if (this->allocatedPointers["icecloud"])delete [] this->icecloud;
    if (this->allocatedPointers["simplecloud"]) delete [] this->simplecloud;
    if (this->allocatedGases) {
        delete[] this-> gases;
        delete[] this-> gas_id;
    }
}

///@brief Set the gas_units
///@param gasUnits 0=ppmv_dry,1=kg_per_kg,2=ppmv_wet
void Profiles::setGasUnits(int gasUnits) {
    gas_units = gasUnits;
}

///@brief Set the mmr_cldaer flag
///@param mmrCldAer false => cm-3 (aer), g/m3 (cld); true => kg/kg (cld, aer)
void Profiles::setMmrCldAer(bool mmrCldAer) {
    mmr_cldaer = mmrCldAer;
}

///@brief Set the pointer to the p array of size [nprofiles][nlevels]
///@param p pointer to an array of doubles [nprofiles][nlevels]
void Profiles::setP(double* p) {
    this->p = p;
}

///@brief Set the pointer to the t array of size [nprofiles][nlevels]
///@param t pointer to an array of doubles [nprofiles][nlevels]
void Profiles::setT(double* t) {
    this->t = t;
}

///@brief Set the pointer to the q array of size [nprofiles][nlevels]
///@param q pointer to an array of double [nprofiles][nlevels]
void Profiles::setQ(double* q) {
    myGases[Q] = q;
}

///@brief Set the pointer to the o3 array of size [nprofiles][nlevels]
///@param o3 pointer to an array of double [nprofiles][nlevels]
void Profiles::setO3(double* o3) {
    myGases[O3] = o3;
}

///@brief Set the pointer to the co2 array of size [nprofiles][nlevels]
///@param co2 pointer to an array of double [nprofiles][nlevels]
void Profiles::setCO2(double* co2) {
    this->myGases[CO2] = co2;
}

///@brief Set the pointer to the n2o array of size [nprofiles][nlevels]
///@param n2o pointer to an array of double [nprofiles][nlevels]
void Profiles::setN2O(double* n2o) {
    this->myGases[N2O] = n2o;
}

///@brief Set the pointer to the co array of size [nprofiles][nlevels]
///@param co pointer to an array of double [nprofiles][nlevels]
void Profiles::setCO(double* co) {
    this->myGases[CO] = co;
}

///@brief Set the pointer to the ch4 array of size [nprofiles][nlevels]
///@param ch4 pointer to an array of double [nprofiles][nlevels]
void Profiles::setCH4(double* ch4) {
    myGases[CH4] = ch4;
}

///@brief Set the pointer to the so2 array of size [nprofiles][nlevels]
///@param so2 pointer to an array of double [nprofiles][nlevels]
void Profiles::setSO2(double* so2) {
    myGases[SO2] = so2;
}

///@brief Set the pointer to the clw array of size [nprofiles][nlevels]
///@param clw pointer to an array of double [nprofiles][nlevels]
void Profiles::setCLW(double* clw) {
    myGases[CLW] = clw;
}

///@brief Set the pointer to the angles array of size [nprofiles][4]
///       containing satzen, satazi, sunzen, sunazi for each profile
///@param angles: pointer to an array of double [nprofiles][4]
///containing angles: satzen, satazi, sunzen, sunazi
///for all the profiles
void Profiles::setAngles(double* angles) {
    this->angles = angles;
}

///@brief Set the pointer to the s2m array of size [nprofiles][6]
///       containing 2m p, 2m t, 2m q, 10m wind u, v, wind fetch for each profile
///@param s2m pointer to an array of double [nprofiles][6]
/// containing the s2m parameters :
/// 2m p, 2m t, 2m q, 10m wind u, v, wind fetch
void Profiles::setS2m(double* s2m) {
    this->s2m = s2m;
}

///@brief Set the pointer to the skin array of size [nprofiles][9]
///       containing skin T, salinity, snow_fraction, foam_fraction, fastem_coefs(1:5) for each profile
///@param skin pointer to an array of double [nprofiles][9]
/// skin T, salinity, snow_fraction, foam_fraction, fastem_coefs(1:5)
void Profiles::setSkin(double* skin) {
    this->skin = skin;
}

///@brief Set the pointer to the surftype array of size [nprofiles][2]
///       containing surftype, watertype for each profile
///@param surftype pointer to an array of int [nprofiles][2]
/// surftype: surftype, watertype
void Profiles::setSurfType(int* surftype) {
    this->surftype = surftype;
}

///@brief Set the pointer to the surfgeom array of size [nprofiles][3]
///       containing latitude, longitude, elevation for each profile
///@param surfgeom pointer to an array of double [nprofiles][3]
/// surfgeom: latitude, longitude, elevation
void Profiles::setSurfGeom(double* surfgeom) {
    this->surfgeom = surfgeom;
}

///@brief Set the pointer to the datetimes array of size [nprofiles][6]
///       containing yy, mm, dd, hh, mm, ss for each profile
///@param datetimes  pointer to an array of int [nprofiles][6]
///containing  yy, mm, dd, hh, mm, ss
///for all the profiles
void Profiles::setDateTimes(int* datetimes) {
    this->datetimes = datetimes;
}

///@brief Set the pointer to the simplecloud array of size [nprofiles][2]
///       containing ctp, cfraction for each profile
///@param simplecloud pointer to an array of double [nprofiles][2]
/// containing ctp, cfraction
void Profiles::setSimpleCloud(double* simplecloud) {
    this->simplecloud = simplecloud;
}

///@brief Set the pointer to the clwscheme array of size [nprofiles][2]
///       containing clw_scheme, clwde_param for each profile
///@param clwscheme pointer to an array of int [nprofiles][2]
/// containing clw_scheme
void Profiles::setClwScheme(int* clwscheme) {
    this->clwscheme = clwscheme;
}

///@brief Set the pointer to the icecloud array of size [nprofiles][2]
///       containing ice_scheme, icede_param for each profile
///@param icecloud pointer to an array of int [nprofiles][2]
/// containing ice_scheme, icede_param
void Profiles::setIceCloud(int* icecloud) {
    this->icecloud = icecloud;
}

///@brief Set the pointer to the zeeman array of size [nprofiles][2]
///       containing be, cosbk for each profile
///@param zeeman pointer to an array of double [nprofiles][2]
///zeeman: be, cosbk
void Profiles::setZeeman(double* zeeman) {
    this->zeeman = zeeman;
}

///@brief Set a gas, cloud or aerosol profile variable; item likes clouds, cfrac or aerosols
/// must have the same dimensions as temperature or water vapour [nprofiles][nlevels]
void Profiles::setGasItem(double* gasItem, rttov::itemIdType item_id) {
    myGases[item_id] = gasItem;
}

///@brief Set the pointer to the gas_id array
///@param gasId pointer to an array of int [nbgases]
/// this method must be called if all the gases, clouds and aerosols are
/// set by the setGases method.
void Profiles::setGasId(int* gasId) {
    gas_id = gasId;
}

///@brief Set the pointer to the gases array
///@param gases pointer to an array of double
///dimensions [ngases][nprofiles][nlevels]
/// the gases can be set 2 different methods :
/// the setGases method for setting all the gases contents
/// (gases, clouds and aerosols)
/// if this method is used, the user must explicitly tell
/// which component are set by calling the setGasId methods
/// the values for the GasId are :
/// Q=1, O3, CO2, N2O, CO, CH4, SO2, CLW=15, CFRAC=20,STCO,STMA,CUCC,CUCP,CUMA,CIRR=30,ICEDE,CLWDE,
/// INSO=41,WASO,SOOT,SSAM,SSCM,MINM,MIAM,MICM,MITR,SUSO,VOLA,VAPO,ASDU
/// BCAR=81,DUS1,DUS2,DUS3,SULP,SSA1,SSA2,SSA3,OMAT
/// AER1=101, AER2=102, ...
void Profiles::setGases(double* gases) {
    if (this->allocatedGases ) {
        delete[] this->gases;
        delete[] this->gas_id;
    }
    this->gases = gases;
}


///@brief Return the number of levels
int Profiles::getNlevels() const {
    return nlevels;
}

///@brief Return the number of profiles
int Profiles::getNprofiles() const {
    return nprofiles;
}

///@brief Return the gas_units
int Profiles::getGasUnits() const {
    return gas_units;
}

///@brief Return the mmr_cldaer flag
bool Profiles::isMmrCldAer() const {
    return mmr_cldaer;
}

///@brief Return the pointer to the P array
double* Profiles::getP() const {
    return p;
}

///@brief Return the pointer to the T array
double* Profiles::getT() const {
    return t;
}

///@brief Return the q pointer
double* Profiles::getQ()  {
    return myGases[Q];
}

///@brief Return the o3 pointer
double* Profiles::getO3()  {
    return myGases[O3];
}

///@brief Return the co2 pointer
double* Profiles::getCO2()  {
    return this->myGases[CO2];
}

///@brief Return the n2o pointer
double* Profiles::getN2O()  {
    return this->myGases[N2O];
}

///@brief Return the co pointer
double* Profiles::getCO()  {
    return this->myGases[CO];
}

///@brief Return the ch4 pointer
double* Profiles::getCH4()  {
    return myGases[CH4];
}

///@brief Return the so2 pointer
double* Profiles::getSO2()  {
    return myGases[SO2];
}

///@brief Return the clw pointer
double* Profiles::getCLW()  {
    return myGases[CLW];
}

///@brief Return a pointer to the angles array
double* Profiles::getAngles() const {
    return angles;
}

///@brief Return the pointer to the s2m array
double* Profiles::getS2m() const {
    return s2m;
}

///@brief Return the pointer to the skin array :
double* Profiles::getSkin() const {
    return skin;
}

///@brief Return the pointer to the surftype array
int* Profiles::getSurfType() const {
    return surftype;
}

///@brief Return the pointer to the surfgeom array
double* Profiles::getSurfGeom() const {
    return surfgeom;
}

///@brief Return a pointer to the datetimes array
int* Profiles::getDateTimes() const {
    return datetimes;
}

///@brief Return the pointer to the simplecloud array
double* Profiles::getSimpleCloud() const {
    return simplecloud;
}

///@brief Return the pointer to the clwscheme array
int* Profiles::getClwScheme() const {
    return clwscheme;
}

///@brief Return the pointer to the icecloud array
int* Profiles::getIceCloud() const {
    return icecloud;
}

///@brief Return the pointer to the zeeman array
double* Profiles::getZeeman() const {
    return zeeman;
}

///@brief Return the number of gas, cloud and aerosol profiles defined per profile
int Profiles::getNgases() const {
    return ngases;
}

///@brief Return the pointer to the gas_id
int* Profiles::getGasId() const {
    return gas_id;
}

///@brief Return a pointer to the gases array
/// dimension [ngases][nprofiles][nlevels]
double* Profiles::getGases() const {
    return gases;
}



///@brief Return true if coefficient file pressure levels will be used
///       (i.e. RTTOV interpolator is not being used)
bool Profiles::isDefaultPressureLevels() const {
    return ( this->p == NULL );
}


///@brief Check all mandatory profile data has been defined and provide
///       zero-filled arrays for optional fields.
bool Profiles::check() {

    if ( this->t == NULL ) {
        std::cerr << "Profiles error: missing T" << std::endl;
        return false;
    }
    if ( this->gases == NULL ) {
        bool ok=this->buildGasesArray();
        if ( ! ok ) {
            std::cerr << "Profiles error: missing gases" << std::endl;
            return false;
        }
    }
    if ( this->gas_id == NULL ) {
        std::cerr << "Profiles error: missing gas_id" << std::endl;
        return false;
    }
    if ( this->gas_units == unknown ) {
        std::cerr << "Profiles warning: gas_units not specified, assuming kg/kg wet" << std::endl;
        this->gas_units=kg_per_kg;
    }
    if ( this->datetimes == NULL ) {
        std::cerr << "Profiles error: missing datetimes" << std::endl;
        return false;
    }
    if ( this->angles == NULL ) {
        std::cerr << "Profiles error: missing angles" << std::endl;
        return false;
    }
    if ( this->surftype == NULL ) {
        std::cerr << "Profiles error: missing surftype" << std::endl;
        return false;
    }
    if ( this->surfgeom == NULL ) {
        std::cerr << "Profiles error: missing surfgeom" << std::endl;
        return false;
    }
    if ( this->zeeman == NULL ) {
//                 std::cerr << "warning profile missing zeeman, assigning zero values" << std::endl;
                this->zeeman = new double[nprofiles*2];
                for (int i=0; i<nprofiles*2; i++){this->zeeman[i]=0.;}
                this->allocatedPointers["zeeman"]=true;
    }
    if ( this->clwscheme == NULL ) {
//                 std::cerr << "warning profile missing clwscheme, assigning default values" << std::endl;
                this->clwscheme = new int[nprofiles*2];
                for (int i=0; i<nprofiles*2; i++){this->clwscheme[i]=1;}
                this->allocatedPointers["clwscheme"]=true;
    }
    if ( this->icecloud == NULL ) {
//                 std::cerr << "warning profile missing icecloud, assigning default values" << std::endl;
                this->icecloud = new int[nprofiles*2];
                for (int i=0; i<nprofiles*2; i++){this->icecloud[i]=1;}
                this->allocatedPointers["icecloud"]=true;
    }
    if ( this->simplecloud == NULL ) {
//                 std::cerr << "warning profile missing simplecloud, assigning zero values" << std::endl;
                this->simplecloud= new double[nprofiles*2];
                for (int i=0; i<nprofiles*2; i++){this->simplecloud[i]=0.;}
                this->allocatedPointers["simplecloud"]=true;
    }
    return (true);
}

/// @brief Build the gases array from individually allocated gases
bool Profiles::buildGasesArray() {
    //
    if (myGases[Q]==NULL) {
        return false;
    }
    // determine dimension
    ngases=0;
    rttov:itemIdType key;
    for(GasIdPointerMap::iterator iter = myGases.begin(); iter != myGases.end(); ++iter) {
        key = iter->first;
        if ( myGases[key] != NULL ) {
            ngases++;
//             std::cerr << "found allocated gas_id: " << key << std::endl;
        }
    }
    // allocate gases and fill it
//     std::cerr << "found " << ngases << " gases" << std::endl;
    gases=new double[ngases*nprofiles*nlevels];
    gas_id=new int[ngases];
    int i=0;
    for(GasIdPointerMap::iterator iter = myGases.begin(); iter != myGases.end(); ++iter) {
         key= iter->first;
         if ( myGases[key] != NULL ) {
            gas_id[i]=key;
            gas_index[key]=i;
            for (int prof=0;prof< nprofiles;prof++) {
               for (int lev=0;lev<nlevels;lev++) {
                 gases[i*nprofiles*nlevels+nlevels*prof+lev]=myGases[key][nlevels*prof+lev];
               }
            }
         i++;
         }
    }
    this->allocatedGases=true;
    return true;
}

///@brief Initialise internal pointers to NULL
void Profiles::initialisePointers() {
    this->allocatedPointers["zeeman"]=false;
    this->allocatedPointers["clwscheme"]=false;
    this->allocatedPointers["icecloud"]=false;
    this->allocatedPointers["simplecloud"]=false;
    for(unsigned i = 0; i < itemIds.size(); ++i) {
        myGases[itemIds[i]]=NULL;
    }
}
} /* namespace rttov */
