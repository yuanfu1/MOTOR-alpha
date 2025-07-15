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
/// @file RttovScatt.cpp
///
/// @brief Wrapper class for RTTOV-SCATT
///

#include <RttovScatt.h>

namespace rttov {

/// @brief RttovScatt class constructor method
RttovScatt::RttovScatt() : coeffsLoaded(false),profileSet(false),calc_zef(0),multi_hydro_frac(0),
        p(NULL),ph(NULL),t(NULL),usercfrac(NULL),gases(NULL),angles(NULL),
        datetimes(NULL),gas_id(NULL),gas_units(NULL),inst_id(-1),
        nchannels(0),s2m(NULL),skin(NULL),surfgeom(NULL),surftype(NULL),
        zeeman(NULL),ngases(0),nprofiles(0),nlevels(0),
        bt(NULL),surfemis(NULL),wavenumbers(NULL),allocatedSurfemis(false),
        p_k(NULL),ph_k(NULL),t_k(NULL),usercfrac_k(NULL),gases_k(NULL),bt_k(NULL),zef_k(NULL),
        skin_k(NULL),s2m_k(NULL),surfemis_k(NULL),profiles(NULL),
        nchannelsForLastRun(0),debug(false),
        btclear(NULL),radquality(NULL),plane_parallel(false),geometric_height(NULL),zef(NULL),azef(NULL),
        emis_terms_cfrac(NULL),emis_terms_bsfc(NULL),emis_terms_tau_cld(NULL),emis_terms_up_cld(NULL),
        emis_terms_down_cld(NULL),emis_terms_tau_clr(NULL),emis_terms_up_clr(NULL),emis_terms_down_clr(NULL)
{

}

/// @brief Deallocate RTTOV structures and drop the coefficients
RttovScatt::~RttovScatt() {
    // Deallocate RTTOV structures
    delete [] bt;
    delete [] bt_k;
    delete [] zef_k;
    delete [] surfemis_k;
    if ( wavenumbers != NULL ) delete [] wavenumbers;
    if (this->allocatedSurfemis) delete [] this->surfemis;
    this->allocatedSurfemis = false;
    delete [] skin_k;
    delete [] s2m_k;
    delete [] p_k;
    delete [] ph_k;
    delete [] t_k;
    delete [] usercfrac_k;
    delete [] gases_k;

    delete [] btclear;
    delete [] radquality;
    delete [] geometric_height;
    delete [] zef;
    delete [] azef;

    delete [] emis_terms_cfrac;
    delete [] emis_terms_bsfc;
    delete [] emis_terms_tau_cld;
    delete [] emis_terms_up_cld;
    delete [] emis_terms_down_cld;
    delete [] emis_terms_tau_clr;
    delete [] emis_terms_up_clr;
    delete [] emis_terms_down_clr;

    int err;
    if ( inst_id > 0 ) {
        rttov_drop_inst_(&err, &inst_id);
        inst_id = -1;
        if ( err!=0) {
            throw runtime_error("Error in rttov_drop_inst");
        }
    }
}


///@brief Return the coefficient filename
const string& RttovScatt::getFileCoef() const {
    return file_coef;
}

///@brief Set the coefficient filename
///@param fileCoef string containing full path to "rtcoef" coefficient file
void RttovScatt::setFileCoef(const string& fileCoef) {
    file_coef = fileCoef;
}

///@brief Return the hydrotable filename
const string& RttovScatt::getFileHydrotable() const {
    return file_hydrotable;
}

///@brief Set the hydrotable filename
///@param fileHydrotable string containing full path to hydrotable file
void RttovScatt::setFileHydrotable(const string& fileHydrotable) {
    file_hydrotable = fileHydrotable;
}

///@brief Return the calc_zef boolean
bool RttovScatt::isCalcZef() const {
    return (calc_zef != 0);
}

///@brief Set the calc_zef boolean
///@param calcZef boolean value
void RttovScatt::setCalcZef(bool calcZef) {
    calc_zef = 0;
    if (calcZef) calc_zef = 1;
}

///@brief Return the multi_hydro_frac boolean
bool RttovScatt::isMultiHydroFrac() const {
    return (multi_hydro_frac != 0);
}

///@brief Set the multi_hydro_frac boolean
///@param multiHydroFrac boolean value
void RttovScatt::setMultiHydroFrac(bool multiHydroFrac) {
    multi_hydro_frac = 0;
    if (multiHydroFrac) multi_hydro_frac = 1;
}

///@brief Load instrument with all channels
///       the methods setFileCoef() and setFileHydrotable() must have been called previously
void RttovScatt::loadInst() {

    string filesPaths;
    if (file_coef.empty()) {
        throw logic_error("Error: call setFileCoef() method first");
    }
    if (is_readable(file_coef)) {
        filesPaths.append(" file_coef ");
        filesPaths.append(file_coef + " ");
    } else {
        throw logic_error("Error: " + file_coef + " does not exist or is not readable");
    }

    if (file_hydrotable.empty()) {
        throw logic_error("Error: call setFileHydrotable() method first");
    }
    if (is_readable(file_hydrotable))  {
        filesPaths.append(" file_hydrotable ");
        filesPaths.append(file_hydrotable + " ");
    } else {
        throw logic_error("Error: " + file_hydrotable + " does not exist or is not readable");
    }

    this->coeffsLoaded = false;

    string optionsAndPaths;
    strOptions=this->options.defineStrOptions();
    optionsAndPaths.append(strOptions);
    optionsAndPaths.append(" ");
    optionsAndPaths.append(filesPaths);

    int len_opts_str=optionsAndPaths.length();
    char opts_str[len_opts_str];
    for (int c = 0; c < len_opts_str; c++ ) {
        opts_str[c]=optionsAndPaths[c];
    }
    if (this->debug) {
        std::cerr << optionsAndPaths << std::endl;
    }

    int nchan = 1;
    int c_channel_list[1];
    c_channel_list[0] = 0;

    rttov_load_inst_(&inst_id, opts_str, &nchan, c_channel_list, len_opts_str);
    if (inst_id > 0) {
        int err=0;
        rttov_get_coef_val_i0_(&err,&inst_id,"NCHANNELS",&nchannels,9);
        if ( err != 0 ) throw runtime_error("Error getting number of channels");
        this->coeffsLoaded=true;
    }
    else {
        throw runtime_error("Error in rttov_load_inst");
    }
    if (this->options.isVerboseWrapper()) {
        std::cerr << "Load successful >>>>> inst_id : " << inst_id << ", nchannels : " << this->nchannels << std::endl;
    }
}

///@brief Return the inst_id
int RttovScatt::getInstId() const {
    return inst_id;
}

///@brief Return true if instrument is loaded
bool RttovScatt::isCoeffsLoaded() const {
    return coeffsLoaded;
}

///@brief Return the number of loaded channels
int RttovScatt::getNchannels() const {
    return nchannels;
}

///@brief Return true if profiles have been associated
bool RttovScatt::isProfileSet() const {
    return profileSet;
}

///@brief Return the number of associated profiles
int RttovScatt::getNprofiles() const {
    return nprofiles;
}

///@brief Return the number of levels of the coefficient file
int RttovScatt::getCoeffsNlevels() {
    if ( this->isCoeffsLoaded() ){
        int coeffNlevels;
        int err=0;
        rttov_get_coef_val_i0_(&err,&inst_id,"SIZE_REF_PRFL_P",&coeffNlevels,15);
        if ( err != 0 ) throw runtime_error("Error getting number of levels");
        return coeffNlevels;
    }
    else {
        return 0;
    }
}

///@brief Return the channel central wavenumbers of the coefficient file
double* RttovScatt::getWaveNumbers() {
    if ( this->isCoeffsLoaded() ){
        int err=0;
        if ( wavenumbers != NULL ) delete [] wavenumbers;
        wavenumbers=new double[nchannels]();
        if (this->debug) std::cerr << "Get wave numbers" << nchannels << std::endl;
        rttov_get_coef_val_r1_(&err,&inst_id,"WAVENUMBERS",&nchannels,wavenumbers,11);
        if ( err != 0 ) throw runtime_error("Error getting channel central wave numbers");
        return wavenumbers;
    }
    else {
        return NULL;
    }
}

///@brief Update RTTOV options for the currently loaded instrument
void RttovScatt::updateOptions() {
    if (! this-> coeffsLoaded  ) {
        throw logic_error( "Error: coefficients not loaded");
    }

    int err=0;
    strOptions=this->options.defineStrOptions();
    int len_opts_str = strOptions.length();
    char opts_str[len_opts_str];
    for (int c = 0; c < len_opts_str; c++ ) {
        opts_str[c]=strOptions[c];
    }
    if (this->debug) std::cerr << "Set RTTOV options " << strOptions << std::endl;
    rttov_set_options_(&err, &inst_id, opts_str, len_opts_str);
    if ( err != 0 ) {
        throw runtime_error("Error in rttov_set_options");
    }
}

///@brief Print RTTOV options for the currently loaded instrument
void RttovScatt::printOptions() {
    if (! this-> coeffsLoaded  ) {
        throw logic_error( "Error: coefficients not loaded");
    }
    int err=0;
    rttov_print_options_(&err, &inst_id);
    if ( err != 0 ) {
        throw runtime_error("Error in rttov_print_options");
    }
}


///@brief Set pointer to array containing input/output surface emissivity values;
///this must be previously allocated a double array of dimensions [nprofiles][nchannels]; this
///is used to pass emissivity values into RTTOV-SCATT; if this is not called the RttovScatt object
///will allocate an array containing the values used by RTTOV-SCATT which can be accessed by getSurfEmis
///@param surfemis pointer to a previously allocated double array [nprofiles][nchannels]
void RttovScatt::setSurfEmis(double* surfemis) {
    if (this->allocatedSurfemis) {
        delete [] this->surfemis;
        this->allocatedSurfemis = false;
    }
    this->surfemis = surfemis;
}

///@brief Return the associated profiles
rttov::ProfilesScatt * RttovScatt::getProfiles() const {
    if (this->isProfileSet()) {
        return profiles;
    } else {
        return NULL;
    }
}

/// @brief Associate a ProfilesScatt object with this RttovScatt object; 
///        this is fast, but does not carry out any checks on profiles before calling RTTOV-SCATT
void  RttovScatt::setProfiles(rttov::ProfilesScatt * profiles) {
    this->profiles = profiles;
    bool profileOK = this->profiles->check();
    if (! profileOK )  {
        this->profileSet = false;
        throw logic_error( "Error: some mandatory profile fields are missing");
    }

    // Make some pointers to simplify code
    this->nlevels=this->profiles->getNlevels();
    this->ngases=this->profiles->getNgases();
    this->nprofiles=this->profiles->getNprofiles();

    this->zeeman=this->profiles->getZeeman();
    this->skin=this->profiles->getSkin();
    this->angles=this->profiles->getAngles();
    this->gas_id=this->profiles->getGasId();
    if (this->debug) for (int i=0; i < this->ngases; i++) std::cerr << "gasid" << " " << i << " : " << this->gas_id[i] << std::endl;
    this->datetimes=this->profiles->getDateTimes();
    this->surfgeom=this->profiles->getSurfGeom();
    this->s2m=this->profiles->getS2m();
    this->surftype=this->profiles->getSurfType();
    this->surfgeom=this->profiles->getSurfGeom();
    this->p=this->profiles->getP();
    this->ph=this->profiles->getPh();
    this->gases=this->profiles->getGases();
    this->t=this->profiles->getT();
    this->usercfrac=this->profiles->getUserCfrac();
    this->gas_units=this->profiles->getGasUnits();

    this->profileSet=true;
}

///@brief Print gases array contents on standard output
void RttovScatt::printGases() {
    if (this->gas_id != NULL) {
        std::cout << "gas_id : " << std::endl;
        for (int g=0; g < ngases; g++) std::cout << "gases " << g << " gas_id " << gas_id[g] << std::endl;
    }
    if (this->p != NULL && this->t != NULL) {
        std::cout << "profile pressure t" << std::endl;
        for  (int prof=0; prof < nprofiles; prof++)
            for (int lev=0; lev < nlevels; lev++) {
                std::cout << prof << " " << p[prof*nlevels+lev] << " " << t[prof*nlevels+lev] << std::endl;
            }
    }
    if (this->gases != NULL) {
        std::cout << "profile gas level value" << std::endl;
        for  (int prof=0; prof < nprofiles; prof++)
            for (int g=0; g < ngases; g++)
                for (int lev=0; lev < nlevels; lev++) {
                    std::cout << prof << " " << g << " " << lev << " " << gases[g*nprofiles*nlevels+nlevels*prof+lev] << std::endl;
                }
    }
    if (this->surfemis != NULL) {
        for  (int prof=0; prof < nprofiles; prof++){
            std::cout << "surfemis prof " << prof << " " << this->surfemis << " " << std::endl;
        }
    }
}


///@brief Run the RTTOV-SCATT direct model for all channels
void RttovScatt::runDirect() {
    vector<int> all_channels;
    for (int c = 0; c < this->getNchannels(); c++) all_channels.push_back(c+1);
    this->runDirect(all_channels);
}

///@brief Run the RTTOV-SCATT direct model for a list of channels
///@param channels vector containing the list of channels to simulate
void RttovScatt::runDirect(const vector<int>& channels) {
    if (! this-> coeffsLoaded  ) {
        throw logic_error( "Error: coefficients not loaded");
    }
    if (! this->profileSet) {
        throw logic_error( "Error: profiles not set");
    }

    int err = 0;

    // update the options which are locally stored in the class instance for RTTOV
    this->updateOptions();

    // when doing the run "chanprof" channel list must be now 1...nchan2
    int nchannels = channels.size();
    int channel_list[nchannels];
    for (int i = 0; i < nchannels; i++) {
        channel_list[i] = channels[i];
    }

    if (this->surfemis == NULL) {
        if (this->options.isVerboseWrapper())
            std::cerr << "No surface emissivity supplied: setting calcemis to true" << std::endl;
        this->allocatedSurfemis = true;
        this->surfemis = new double [nprofiles*nchannels]();
        for (int i = 0; i < nprofiles * nchannels; i++) this->surfemis[i] = -1.;
    }
    else if (allocatedSurfemis) {
        // Reallocate in case nprofiles or nchannels changed since the previous call
        delete [] this->surfemis;
        surfemis = new double [nprofiles*nchannels]();
        for (int i = 0; i < nprofiles * nchannels; i++) surfemis[i] = -1;
    }

    // alloc memory for bt
    delete [] bt;
    if (this->debug)
        std::cerr << "Allocating bt for " << nprofiles << " profiles and " << nchannels << " channels" << std::endl;
    bt = new double[nprofiles*nchannels]();
    if (this->debug) {
        std::cerr << "gas_units= " << this->gas_units << std::endl;
        for (int p=0; p< nprofiles;p++) {
            std::cerr << "Profile " << p << std::endl;
            std::cerr << "datetimes ";
            for (int i=0;i<6;i++) std::cerr << datetimes[6*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "angles ";
            for (int i=0;i<2;i++) std::cerr << angles[2*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "surfgeom ";
            for (int i=0;i<3;i++) std::cerr << surfgeom[3*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "surftype ";
            std::cerr << surftype[p]<<" ";
            std::cerr <<std::endl;
            std::cerr << "skin ";
            for (int i=0;i<8;i++) std::cerr << skin[8*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "s2m ";
            for (int i=0;i<5;i++) std::cerr << s2m[5*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr << "zeeman ";
            for (int i=0;i<2;i++) std::cerr << zeeman[2*p+i]<<" ";
            std::cerr <<std::endl;
            std::cerr <<std::endl;
        }
    }

    rttov_scatt_call_direct_(
            &err,
            &inst_id,
            channel_list,
            datetimes,         // profile dates/times                                     [nprofiles][6]
            angles,            // satzen, satazi angles                                   [nprofiles][2]
            surfgeom,          // lat, lon, elevation                                     [nprofiles][3]
            surftype,          // surftype                                                [nprofiles]
            skin,              // skin T, salinity, foam_frac, fastem_coefsx5             [nprofiles][8]
            s2m,               // 2m p, 2m t, 2m q, 10m wind u, v                         [nprofiles][5]
            zeeman,            // Be, cosbk                                               [nprofiles][2]
            p,                 // pressure                                                [nprofiles][nlevels]
            t,                 // temperature                                             [nprofiles][nlevels]
            &gas_units,        // units for gas profiles
            gas_id,            // gas ID list                                             [ngases]
            gases,             // gas profiles                                            [ngases][nprofiles][nlevels]
            ph,                // pressure half-levels                                    [nprofiles][nlevels+1]
            usercfrac,         // user cloud fraction                                     [nprofiles]
            &multi_hydro_frac, // false => single hydro_frac profile, true => one hydro_frac profile per hydrometeor
            &calc_zef,         // enable/disable radar reflectivity calculations
            surfemis,          // input/output surface emissivities                       [nprofiles][nchannels]
            bt,                // output BTs                                              [nprofiles][nchannels]
            &nchannels, &ngases, &nlevels, &nprofiles);
    if ( err != 0 ) throw runtime_error ("Error in rttov_scatt_call_direct");
    this->nchannelsForLastRun=nchannels;

    if (options.isStoreRad())       doStoreRad();
    if (options.isStoreEmisTerms()) doStoreEmisTerms();
    if (calc_zef != 0)              doStoreZef();
}


///@brief Run the RTTOV-SCATT K model for all channels
void RttovScatt::runK() {
    vector<int> all_channels;
    for (int c=0; c < this->getNchannels(); c++) all_channels.push_back(c+1);
    this->runK(all_channels);
}

///@brief Run the RTTOV-SCATT K model for a list of channels
///@param channels vector containing the list of channels to simulate
void RttovScatt::runK(const vector<int>& channels) {
    int err=0;

    if (! this-> coeffsLoaded  ) {
        throw logic_error( "Error: coefficients not loaded");
    }
    if (! this->profileSet) {
        throw logic_error( "Error: profiles not set");
    }

    // update the options which are locally stored in the class instance for RTTOV
    this->updateOptions();

    // when doing the run "chanprof" channel list must be now 1...nchan2
    int nchannels = channels.size();
    int channel_list[nchannels];
    for (int i = 0; i < nchannels; i++) {
        channel_list[i] = channels[i];
    }

    // allocate arrays for rttov_scatt_call_k output
    delete [] skin_k;
    skin_k= new double [nprofiles*nchannels*8]();
    delete [] s2m_k;
    s2m_k= new double [nprofiles*nchannels*5]();
    delete [] p_k;
    p_k= new double [nprofiles*nchannels*nlevels]();
    delete [] ph_k;
    ph_k= new double [nprofiles*nchannels*(nlevels+1)]();
    delete [] t_k;
    t_k = new double [nprofiles*nchannels*nlevels]();
    delete [] usercfrac_k;
    usercfrac_k = new double [nprofiles*nchannels]();
    delete [] gases_k;
    gases_k= new double [ngases*nprofiles*nchannels*nlevels]();

    if (surfemis == NULL) {
        if (this->options.isVerboseWrapper())
            std::cerr << "No surface emissivity supplied: setting calcemis to true" << std::endl;
        surfemis = new double [nprofiles*nchannels]();
        for (int i = 0; i < nprofiles * nchannels; i++) surfemis[i]=-1;
        allocatedSurfemis=true;
    }
    else if (allocatedSurfemis) {
        // Reallocate in case nprofiles or nchannels changed since the previous call
        delete [] this->surfemis;
        surfemis = new double [nprofiles*nchannels]();
        for (int i = 0; i < nprofiles * nchannels; i++) surfemis[i]=-1;
    }

    delete [] this->surfemis_k;
    surfemis_k = new double [nprofiles*nchannels]();
    for (int i = 0; i < nprofiles * nchannels; i++) surfemis_k[i]=0;

    delete [] this->bt_k;
    bt_k = new double [nprofiles*nchannels]();
    delete [] this->zef_k;
    zef_k = new double [nprofiles*nchannels*nlevels]();
    if (calc_zef != 0) {
        for (int i = 0; i < nprofiles * nchannels; i++) bt_k[i]=0;
        for (int i = 0; i < nprofiles * nchannels * nlevels; i++) zef_k[i]=1;
    }
    else {
        for (int i = 0; i < nprofiles * nchannels; i++) bt_k[i]=1;
        for (int i = 0; i < nprofiles * nchannels * nlevels; i++) zef_k[i]=0;
    }

    if (this->debug) {
        std::cerr << "inst_id " << inst_id << std::endl;
        std::cerr << "channels " << nchannels << std::endl;
        std::cerr << "ngases " << ngases << std::endl;
        std::cerr << "nlevels " << nlevels << std::endl;
        std::cerr << "nprofiles " << nprofiles << std::endl;
    }

    delete [] bt;
    if (this->debug)
        std::cerr << "Allocating bt for " << nprofiles << " profiles and " << nchannels << " channels" << std::endl;
    bt = new double[nprofiles*nchannels]();

    rttov_scatt_call_k_(
            &err,
            &inst_id,
            channel_list,
            datetimes,         // profile dates/times                                     [nprofiles][6]
            angles,            // satzen, satazi, sunzen, sunazi angles                   [nprofiles][4]
            surfgeom,          // lat, lon, elevation                                     [nprofiles][3]
            surftype,          // surftype, watertype                                     [nprofiles][2]
            skin,              // skin T, salinity, snow_frac, foam_frac, fastem_coefsx5  [nprofiles][9]
            skin_k,            // output skin K                                           [nprofiles][nchannels][9]
            s2m,               // 2m p, 2m t, 2m q, 10m wind u, v, wind-fetch             [nprofiles][6]
            s2m_k,             // output 2m K                                             [nprofiles][nchannels][6]
            zeeman,            // Be, cosbk                                               [nprofiles][2]
            p,                 // pressure                                                [nprofiles][nlevels]
            p_k,               // output pressure K                                       [nprofiles][nchannels][nlevels]
            t,                 // temperature                                             [nprofiles][nlevels]
            t_k,               // output temperature K                                    [nprofiles][nchannels][nlevels]
            &gas_units,        // units for gas profiles
            gas_id,            // gas ID list                                             [ngases]
            gases,             // gas profiles                                            [ngases][nprofiles][nlevels]
            gases_k,           // output gas profiles K                                   [ngases][nprofiles][nchannels][nlevels]
            ph,                // pressure half-levels                                    [nprofiles][nlevels+1]
            ph_k,              // output pressure half-levels K                           [nprofiles][nchannels][nlevels+1]
            usercfrac,         // user cloud fraction                                     [nprofiles]
            usercfrac_k,       // user cloud fraction K                                   [nprofiles][nchannels]
            &multi_hydro_frac, // false => single hydro_frac profile, true => one hydro_frac profile per hydrometeor
            &calc_zef,         // enable/disable radar reflectivity calculations
            surfemis,          // input/output surface emissivities                       [nprofiles][nchannels]
            surfemis_k,        // input/output surface emissivities K                     [nprofiles][nchannels]
            bt,                // output BTs                                              [nprofiles][nchannels]
            bt_k,              // input BT perturbations                                  [nprofiles][nchannels]
            zef_k,             // input radar reflectivity perturbations                  [nprofiles][nchannels][nlevels]
            &nchannels, &ngases, &nlevels, &nprofiles);

    if ( err != 0 ) throw runtime_error ("Error in rttov_scatt_call_k");
    this->nchannelsForLastRun=nchannels;

    if (options.isStoreRad())   doStoreRad();
    if (calc_zef != 0)          doStoreZef();
}

///@brief Return a pointer to an array of dimensions [nprofiles][nchannels]
/// containing output values of surface emissivity;
/// this array can be initialised by the user and set by calling the setSurfEmis method;
/// alternatively if the emissivity array is allocated by the RttovScatt object
/// it is deleted at the next run or when the RttovScatt instance is destroyed.
const double* RttovScatt::getSurfEmis() const {
    return surfemis;
}

///@brief Return a pointer to an array of dimensions [nprofiles][nchannels]
/// filled with computed brightness temperatures by the previous run;
/// this array is allocated by the RttovScatt object and is destroyed when
/// a new run is performed or if the instance is destroyed
const double* RttovScatt::getBt() const {
    return bt;
}

///@brief Return a vector for x y z from a pointer to [dim1][dim2][dim3][dim4]
std::vector<double> RttovScatt::convertPointer4D2Vector(double ptr[], int x, int y,
        int z, int dim1, int dim2, int dim3, int dim4) {
    std::vector <double>myvector{};
    if (ptr == NULL) throw logic_error("Error: data not present");
    if (x < 0 || x >= dim1) throw logic_error("Error: wrong call to convertPointer4D2Vector");
    if (y < 0 || y >= dim2) throw logic_error("Error: wrong call to convertPointer4D2Vector");
    if (z < 0 || z >= dim3) throw logic_error("Error: wrong call to convertPointer4D2Vector");
    if (ptr == NULL) throw logic_error("Error: wrong call to convertPointer4D2Vector");
    for (int w=0; w < dim4; w++) myvector.push_back(ptr[x*dim4*dim3*dim2 + dim4*dim3*y + dim4*z+w]);
    return myvector;
}

///@brief Return a vector for x and y from a pointer to [dim1][dim2][dim3]
std::vector<double> RttovScatt::convertPointer3D2Vector(double ptr[], const int x,
        const int y, const int dim1, const int dim2, const int dim3) {
    std::vector <double>myvector{};
    if (ptr == NULL) throw logic_error("Error: data not present");
    if (x < 0 || x >= dim1) throw logic_error("Error: wrong call to convertPointer3D2Vector");
    if (y < 0 || y >= dim2) throw logic_error("Error: wrong call to convertPointer3D2Vector");
    if (ptr == NULL) throw logic_error("Error: wrong call to convertPointer3D2Vector");
    for (int z=0; z < dim3; z++) myvector.push_back(ptr[x*dim3*dim2 + dim3*y + z]);
    return myvector;
}

///@brief Return a vector for x from a pointer to [dim1][dim2]
std::vector<double> RttovScatt::convertPointer2D2Vector(double ptr[], int x, int dim1, int dim2) {
    if (ptr == NULL) throw logic_error("Error: data not present");
    if (x < 0 || x >= dim1) throw logic_error("Error: wrong call to convertPointer2D2Vector");
    if (ptr == NULL) throw logic_error("Error: wrong call to convertPointer2D2Vector");
    std::vector <double>myvector{};
    for (int y=0; y < dim2; y++) {
        myvector.push_back(ptr[x*dim2 + y]);
    }
    return myvector;
}

///@brief Return an integer vector for x from a pointer to integer [dim1][dim2]
std::vector<int> RttovScatt::convertIntPointer2D2Vector(int ptr[], int x, int dim1, int dim2) {
    if (ptr == NULL) throw logic_error("Error: data not present");
    if (x < 0 || x >= dim1) throw logic_error("Error: wrong call to convertIntPointer2D2Vector");
    if (ptr == NULL) throw logic_error("Error: wrong call to convertIntPointer2D2Vector");
    std::vector <int>myvector{};
    for (int y=0; y < dim2; y++) {
        myvector.push_back(ptr[x*dim2 + y]);
    }
    return myvector;
}

///@brief Return vector of brightness temperatures
///       computed by the previous run for the given profile number
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovScatt::getBt(const int profile)  {
    return this->convertPointer2D2Vector(this->bt,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return the computed pressure Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScatt::getPK(int profile, int channel) {
    return this->convertPointer3D2Vector(this->p_k,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

///@brief Return the computed pressure half-level Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScatt::getPhK(int profile, int channel) {
    return this->convertPointer3D2Vector(this->ph_k,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels+1);
}

///@brief Return vector of user cloud fraction Jacobians for a given profile
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovScatt::getUserCfracK(const int profile)  {
    return this->convertPointer2D2Vector(this->usercfrac_k,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}


///@brief Return computed temperature Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScatt::getTK(int profile, int channel) {
    return this->convertPointer3D2Vector(this->t_k,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

///@brief Return computed skin variable Jacobians for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScatt::getSkinK(int profile, int channel) {
    return this->convertPointer3D2Vector(this->skin_k,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            8);
}

///@brief Return computed 2m variable Jacobian for a given profile and channel
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScatt::getS2mK(int profile, int channel) {
    return this->convertPointer3D2Vector(this->s2m_k,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            5);
}

///@brief Return computed gas and hydrometeor Jacobian values for a given profile and channel
///@param itemIdType item type
///@param profile index of the profile (0..nprofiles-1)
///@param channel index of the channel (0..nchannels-1)
std::vector<double> RttovScatt::getItemK(rttov::itemIdType itemIdType, int profile,
        int channel) {
    return this->convertPointer4D2Vector(this->gases_k,
            this->profiles->gas_index[itemIdType],
            profile,
            channel,
            this->ngases,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

///@brief Return computed surface emissivity Jacobians for a given profile
///@param profile index of the profile (0..nprofiles-1)
std::vector<double> RttovScatt::getSurfEmisK(int profile) {
    return this->convertPointer2D2Vector(this->surfemis_k,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

void RttovScatt::doStoreRad() {
    int err = 0;
    int nchanprof = nprofiles * nchannelsForLastRun;
    int nlayers = nlevels-1;

    delete [] btclear;
    btclear = new double [nchanprof]();
    rttov_get_bt_clear_(&err, &inst_id, btclear, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting btclear");

    delete [] radquality;
    radquality = new int [nchanprof]();
    rttov_get_rad_quality_(&err, &inst_id, radquality, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting radquality");

    rttov_get_plane_parallel_(&err, &inst_id, &plane_parallel);
    if ( err != 0 ) throw runtime_error("Error getting plane_parallel");

    delete [] geometric_height;
    geometric_height = new double [nchanprof*nlevels]();
    rttov_get_geometric_height_(&err, &inst_id, geometric_height, &nchanprof, &nlevels);
    if ( err != 0 ) throw runtime_error("Error getting geometric_height");
}

///@brief Return RTTOV radiance bt_clear output array of size [nchannels] for given profile, requires store_rad true
std::vector <double> RttovScatt::getBtClear(int profile) {
    return this->convertPointer2D2Vector(this->btclear,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance quality flag array of size [nchannels] for given profile, requires store_rad true
std::vector <int> RttovScatt::getRadQuality(int profile) {
    return this->convertIntPointer2D2Vector(this->radquality,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV radiance plane_parallel flag, requires store_rad true
bool RttovScatt::getPlaneParallel() {
    return (this->plane_parallel != 0);
}

///@brief Return RTTOV radiance geometric_height output array of size [nlevels] for given profile and channel, requires store_rad true
std::vector <double> RttovScatt::getGeometricHeight(int profile, int channel) {
    return this->convertPointer3D2Vector(this->geometric_height,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

void RttovScatt::doStoreZef() {
    int err = 0;
    int nchanprof = nprofiles * nchannelsForLastRun;

    delete [] zef;
    zef = new double [nchanprof*nlevels]();
    rttov_get_zef_(&err, &inst_id, zef, &nchanprof, &nlevels);
    if ( err != 0 ) throw runtime_error("Error getting zef");

    delete [] azef;
    azef = new double [nchanprof*nlevels]();
    rttov_get_azef_(&err, &inst_id, azef, &nchanprof, &nlevels);
    if ( err != 0 ) throw runtime_error("Error getting azef");
}

///@brief Return RTTOV-SCATT radar reflectivity zef output array of size [nlevels] for given profile and channel
std::vector <double> RttovScatt::getZef(int profile, int channel) {
    return this->convertPointer3D2Vector(this->zef,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

///@brief Return RTTOV-SCATT radar reflectivity azef output array of size [nlevels] for given profile and channel
std::vector <double> RttovScatt::getAZef(int profile, int channel) {
    return this->convertPointer3D2Vector(this->azef,
            profile,
            channel,
            this->nprofiles,
            this->nchannelsForLastRun,
            this->nlevels);
}

void RttovScatt::doStoreEmisTerms() {
    int err = 0;
    int nchanprof = nprofiles * nchannelsForLastRun;

    delete [] emis_terms_cfrac;
    emis_terms_cfrac = new double [nchanprof]();
    rttov_get_emis_terms_cfrac_(&err, &inst_id, emis_terms_cfrac, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting emis_terms_cfrac");

    delete [] emis_terms_bsfc;
    emis_terms_bsfc = new double [nchanprof]();
    rttov_get_emis_terms_bsfc_(&err, &inst_id, emis_terms_bsfc, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting emis_terms_bsfc");

    delete [] emis_terms_tau_cld;
    emis_terms_tau_cld = new double [nchanprof]();
    rttov_get_emis_terms_tau_cld_(&err, &inst_id, emis_terms_tau_cld, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting emis_terms_tau_cld");

    delete [] emis_terms_up_cld;
    emis_terms_up_cld = new double [nchanprof]();
    rttov_get_emis_terms_up_cld_(&err, &inst_id, emis_terms_up_cld, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting emis_terms_up_cld");

    delete [] emis_terms_down_cld;
    emis_terms_down_cld = new double [nchanprof]();
    rttov_get_emis_terms_down_cld_(&err, &inst_id, emis_terms_down_cld, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting emis_terms_down_cld");

    delete [] emis_terms_tau_clr;
    emis_terms_tau_clr = new double [nchanprof]();
    rttov_get_emis_terms_tau_clr_(&err, &inst_id, emis_terms_tau_clr, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting emis_terms_tau_clr");

    delete [] emis_terms_up_clr;
    emis_terms_up_clr = new double [nchanprof]();
    rttov_get_emis_terms_up_clr_(&err, &inst_id, emis_terms_up_clr, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting emis_terms_up_clr");

    delete [] emis_terms_down_clr;
    emis_terms_down_clr = new double [nchanprof]();
    rttov_get_emis_terms_down_clr_(&err, &inst_id, emis_terms_down_clr, &nchanprof);
    if ( err != 0 ) throw runtime_error("Error getting emis_terms_down_clr");
}

///@brief Return RTTOV-SCATT emis retrieval cfrac output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScatt::getEmisTermsCfrac(int profile) {
    return this->convertPointer2D2Vector(this->emis_terms_cfrac,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV-SCATT emis retrieval bsfc output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScatt::getEmisTermsBsfc(int profile) {
    return this->convertPointer2D2Vector(this->emis_terms_bsfc,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV-SCATT emis retrieval tau_cld output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScatt::getEmisTermsTauCld(int profile) {
    return this->convertPointer2D2Vector(this->emis_terms_tau_cld,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV-SCATT emis retrieval up_cld output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScatt::getEmisTermsUpCld(int profile) {
    return this->convertPointer2D2Vector(this->emis_terms_up_cld,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV-SCATT emis retrieval down_cld output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScatt::getEmisTermsDownCld(int profile) {
    return this->convertPointer2D2Vector(this->emis_terms_down_cld,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV-SCATT emis retrieval tau_clr output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScatt::getEmisTermsTauClr(int profile) {
    return this->convertPointer2D2Vector(this->emis_terms_tau_clr,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV-SCATT emis retrieval up_clr output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScatt::getEmisTermsUpClr(int profile) {
    return this->convertPointer2D2Vector(this->emis_terms_up_clr,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

///@brief Return RTTOV-SCATT emis retrieval down_clr output array of size [nchannels] for given profile, requires store_emis_terms true
std::vector <double> RttovScatt::getEmisTermsDownClr(int profile) {
    return this->convertPointer2D2Vector(this->emis_terms_down_clr,
            profile,
            this->nprofiles,
            this->nchannelsForLastRun);
}

} /* namespace rttov */
