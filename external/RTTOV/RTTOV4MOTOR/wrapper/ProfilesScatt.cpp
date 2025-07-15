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

///@file ProfilesScatt.cpp
///
///@brief Method definition for ProfilesScatt class

#include <ProfilesScatt.h>

namespace rttov {

///@brief Constructor method for manual gas array creation
///@param nbprofiles number of profiles
///@param nblevels number of vertical levels
///@param nbgases number of gases
ProfilesScatt::ProfilesScatt(int nbprofiles, const int nblevels, int nbgases)
: nprofiles(nbprofiles), nlevels(nblevels), ngases(nbgases),
    p(NULL),ph(NULL),t(NULL),usercfrac(NULL),gases(NULL),angles(NULL),
    datetimes(NULL),gas_id(NULL),gas_units(unknown),
    s2m(NULL),skin(NULL),surfgeom(NULL),surftype(NULL),
    zeeman(NULL),allocatedGases(false) {
    this->initialisePointers();
}

///@brief Constructor method for individual gas specification
///@param nbprofiles number of profiles
///@param nblevels number of vertical levels
ProfilesScatt::ProfilesScatt(int nbprofiles, const int nblevels)
: nprofiles(nbprofiles), nlevels(nblevels), ngases(0),
    p(NULL),ph(NULL),t(NULL),usercfrac(NULL),gases(NULL),angles(NULL),
    datetimes(NULL),gas_id(NULL),gas_units(unknown),
    s2m(NULL),skin(NULL),surfgeom(NULL),surftype(NULL),
    zeeman(NULL),allocatedGases(false) {
    this->initialisePointers();
}

///@brief Constructor method for one profiles and 54 levels
ProfilesScatt::ProfilesScatt(): nprofiles(1), nlevels(54), ngases(0),
        p(NULL),ph(NULL),t(NULL),usercfrac(NULL),gases(NULL),angles(NULL),
        datetimes(NULL),gas_id(NULL),gas_units(unknown),
        s2m(NULL),skin(NULL),surfgeom(NULL),surftype(NULL),
        zeeman(NULL), allocatedGases(false) {
    this->initialisePointers();
}

///@brief Destructor method
ProfilesScatt::~ProfilesScatt() {
    // deallocate allocated arrays if necessary
    if (this->allocatedPointers["zeeman"]) delete [] this->zeeman;
    if (this->allocatedPointers["usercfrac"]) delete [] this->usercfrac;
    if (this->allocatedGases) {
        delete[] this-> gases;
        delete[] this-> gas_id;
    }
}

///@brief Set the gas_units
///@param gasUnits 0=ppmv_dry,1=kg_per_kg,2=ppmv_wet
void ProfilesScatt::setGasUnits(int gasUnits) {
    gas_units = gasUnits;
}

///@brief Set the pointer to the p array of size [nprofiles][nlevels]
///@param p pointer to an array of doubles [nprofiles][nlevels]
void ProfilesScatt::setP(double* p) {
    this->p = p;
}

///@brief Set the pointer to the ph array of size [nprofiles][nlevels+1]
///@param ph pointer to an array of doubles [nprofiles][nlevels+1]
void ProfilesScatt::setPh(double* ph) {
    this->ph = ph;
}

///@brief Set the pointer to the t array of size [nprofiles][nlevels]
///@param t pointer to an array of doubles [nprofiles][nlevels]
void ProfilesScatt::setT(double* t) {
    this->t = t;
}

///@brief Set the pointer to the q array of size [nprofiles][nlevels]
///@param q pointer to an array of double [nprofiles][nlevels]
void ProfilesScatt::setQ(double* q) {
    this->myGases[Q] = q;
}

///@brief Set the pointer to the o3 array of size [nprofiles][nlevels]
///@param o3 pointer to an array of double [nprofiles][nlevels]
void ProfilesScatt::setO3(double* o3) {
    this->myGases[O3] = o3;
}

///@brief Set the pointer to the hydro_frac array of size [nprofiles][nlevels]
///@param hydro_frac pointer to an array of double [nprofiles][nlevels]
void ProfilesScatt::setHydroFrac(double* hydro_frac) {
    this->myGases[SCATT_HYDRO_FRAC] = hydro_frac;
}

///@brief Set the pointer to the clw array of size [nprofiles][nlevels]
///@param clw pointer to an array of double [nprofiles][nlevels]
void ProfilesScatt::setClw(double* clw) {
    this->myGases[SCATT_CLW] = clw;
}

///@brief Set the pointer to the ciw array of size [nprofiles][nlevels]
///@param ciw pointer to an array of double [nprofiles][nlevels]
void ProfilesScatt::setCiw(double* ciw) {
    this->myGases[SCATT_CIW] = ciw;
}

///@brief Set the pointer to the snow array of size [nprofiles][nlevels]
///@param snow pointer to an array of double [nprofiles][nlevels]
void ProfilesScatt::setSnow(double* snow) {
    this->myGases[SCATT_SNOW] = snow;
}

///@brief Set the pointer to the rain array of size [nprofiles][nlevels]
///@param rain pointer to an array of double [nprofiles][nlevels]
void ProfilesScatt::setRain(double* rain) {
    this->myGases[SCATT_RAIN] = rain;
}

///@brief Set the pointer to the graupel array of size [nprofiles][nlevels]
///@param graupel pointer to an array of double [nprofiles][nlevels]
void ProfilesScatt::setGraupel(double* graupel) {
    this->myGases[SCATT_GRAUPEL] = graupel;
}

///@brief Set the pointer to the user cfrac array of size [nprofiles]
///@param usercfrac pointer to an array of double [nprofiles]
void ProfilesScatt::setUserCfrac(double* usercfrac) {
    this->usercfrac = usercfrac;
}

///@brief Set the pointer to the angles array of size [nprofiles][2]
///       containing satzen, satazi for each profile
///@param angles: pointer to an array [nprofiles][2]
///containing angles: satzen, satazi
///for all the profiles
void ProfilesScatt::setAngles(double* angles) {
    this->angles = angles;
}

///@brief Set the pointer to the s2m array of size [nprofiles][5]
///       containing 2m p, 2m t, 2m q, 10m wind u, v for each profile
///@param s2m pointer to an array of double [nprofiles][5]
/// containing the s2m parameters :
/// 2m p, 2m t, 2m q, 10m wind u, v
void ProfilesScatt::setS2m(double* s2m) {
    this->s2m = s2m;
}

///@brief Set the pointer to the skin array of size [nprofiles][8]
///       containing skin T, salinity, foam_fraction, fastem_coefs(1:5) for each profile
///@param skin pointer to an array of doubles [nprofiles][8]
/// skin T, salinity, foam_fraction, fastem_coefs(1:5)
void ProfilesScatt::setSkin(double* skin) {
    this->skin = skin;
}

///@brief Set the pointer to the surftype array of size [nprofiles]
///       containing surftype for each profile
///@param surftype pointer to an array of int [nprofiles]
void ProfilesScatt::setSurfType(int* surftype) {
    this->surftype = surftype;
}

///@brief Set the pointer to the surfgeom array of size [nprofiles][3]
///       containing latitude, longitude, elevation for each profile
///@param surfgeom pointer to an array of doubles [nprofiles][3]
/// surfgeom: latitude, longitude, elevation
void ProfilesScatt::setSurfGeom(double* surfgeom) {
    this->surfgeom = surfgeom;
}

///@brief Set the pointer to the datetimes array of size [nprofiles][6]
///       containing yy, mm, dd, hh, mm, ss for each profile
///@param datetimes  pointer to an array [nprofiles][6]
///containing  yy, mm, dd, hh, mm, ss
///for all the profiles
void ProfilesScatt::setDateTimes(int* datetimes) {
    this->datetimes = datetimes;
}

///@brief Set the pointer to the zeeman array of size [nprofiles][2]
///       containing be, cosbk for each profile
///@param zeeman pointer to an array of int [nprofiles][2]
///zeeman: be, cosbk
void ProfilesScatt::setZeeman(double* zeeman) {
    this->zeeman = zeeman;
}

///@brief Set a gas or hydrometeor profile variable
/// must have the same dimensions as temperature or water vapour [nprofiles][nlevels]
void ProfilesScatt::setGasItem(double* gasItem, rttov::itemIdType item_id) {
    myGases[item_id] = gasItem;
}

///@brief Set the pointer to the gas_id array
///@param gasId pointer to an array of int [nbgases]
/// this method must be called if all the gases, clouds and aerosols are
/// set by the setGases method.
void ProfilesScatt::setGasId(int* gasId) {
    gas_id = gasId;
}

///@brief Set the pointer to the gases array
///@param gases pointer to an array of double
///dimensions [ngases][nprofiles][nlevels]
/// the gases can be set 2 different methods :
/// the setGases method for setting all the gases contents
/// (gases and hydrometeors)
/// if this method is used, the user must explicitly tell
/// which component are set by calling the setGasId methods
/// the values for the GasId are :
/// Q=1,O3=2,SCATT_HYDRO_FRAC=60,SCATT_CLW,SCATT_CIW,SCATT_RAIN,SCATT_SNOW,SCATT_GRAUPEL
/// HYDRO1=201,..., HYDRO30=230,HYDRO_FRAC1=301,...HYDRO_FRAC30=330
///
void ProfilesScatt::setGases(double* gases) {
    if (this->allocatedGases ) {
        delete[] this->gases;
        delete[] this->gas_id;
    }
    this->gases = gases;
}


///@brief Return the number of levels
int ProfilesScatt::getNlevels() const {
    return nlevels;
}

///@brief Return the number of profiles
int ProfilesScatt::getNprofiles() const {
    return nprofiles;
}

///@brief Return the gas_units
int ProfilesScatt::getGasUnits() const {
    return gas_units;
}

///@brief Return the pointer to the P array
double* ProfilesScatt::getP() const {
    return p;
}

///@brief Return the pointer to the Ph array
double* ProfilesScatt::getPh() const {
    return ph;
}

///@brief Return the pointer to the T array
double* ProfilesScatt::getT() const {
    return t;
}

///@brief Return the pointer to the user cfrac array
double* ProfilesScatt::getUserCfrac() const {
    return usercfrac;
}

///@brief Return the q pointer
double* ProfilesScatt::getQ()  {
    return myGases[Q];
}

///@brief Return the o3 pointer
double* ProfilesScatt::getO3()  {
    return myGases[O3];
}

///@brief Return the hydro_frac pointer
double* ProfilesScatt::getHydroFrac()  {
    return myGases[SCATT_HYDRO_FRAC];
}

///@brief Return the clw pointer
double* ProfilesScatt::getClw()  {
    return this->myGases[SCATT_CLW];
}

///@brief Return the ciw pointer
double* ProfilesScatt::getCiw()  {
    return this->myGases[SCATT_CIW];
}

///@brief Return the snow pointer
double* ProfilesScatt::getSnow()  {
    return this->myGases[SCATT_SNOW];
}

///@brief Return the rain pointer
double* ProfilesScatt::getRain()  {
    return myGases[SCATT_RAIN];
}

///@brief Return the graupel pointer
double* ProfilesScatt::getGraupel()  {
    return myGases[SCATT_GRAUPEL];
}

///@brief Return a pointer to the angles array
double* ProfilesScatt::getAngles() const {
    return angles;
}

///@brief Return the pointer to the s2m array
double* ProfilesScatt::getS2m() const {
    return s2m;
}

///@brief Return the pointer to the skin array :
double* ProfilesScatt::getSkin() const {
    return skin;
}

///@brief Return the pointer to the surftype array
int* ProfilesScatt::getSurfType() const {
    return surftype;
}

///@brief Return the pointer to the surfgeom array
double* ProfilesScatt::getSurfGeom() const {
    return surfgeom;
}

///@brief Return a pointer to the datetimes array
int* ProfilesScatt::getDateTimes() const {
    return datetimes;
}

///@brief Return the pointer to the zeeman array
double* ProfilesScatt::getZeeman() const {
    return zeeman;
}

///@brief Return the number of gas and hydrometeor profiles defined per profile
int ProfilesScatt::getNgases() const {
    return ngases;
}

///@brief Return the pointer to the gas_id
int* ProfilesScatt::getGasId() const {
    return gas_id;
}

///@brief Return a pointer to the gases array
/// dimension [ngases][nprofiles][nlevels]
double* ProfilesScatt::getGases() const {
    return gases;
}


///@brief Check all mandatory profile data has been defined and provide
///       zero-filled arrays for optional fields.
bool ProfilesScatt::check() {

    if ( this->p == NULL ) {
        std::cerr << "Profiles error: missing P" << std::endl;
        return false;
    }
    if ( this->ph == NULL ) {
        std::cerr << "Profiles error: missing Ph" << std::endl;
        return false;
    }
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
    if ( this->usercfrac == NULL ) {
//                 std::cerr << "warning profile missing usercfrac, assigning default values" << std::endl;
                this->usercfrac = new double[nprofiles];
                for (int i=0; i<nprofiles; i++){this->usercfrac[i]=0.;}
                this->allocatedPointers["usercfrac"]=true;
    }
    return (true);
}

/// @brief Build the gases array from individually allocated gases
bool ProfilesScatt::buildGasesArray() {
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
void ProfilesScatt::initialisePointers() {
    this->allocatedPointers["zeeman"]=false;
    this->allocatedPointers["usercfrac"]=false;
    for(unsigned i = 0; i < itemIdsScatt.size(); ++i) {
        myGases[itemIdsScatt[i]]=NULL;
    }
}
} /* namespace rttov */
