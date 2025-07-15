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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
*/

///
/// @file Atlas.cpp
///
/// @brief Wrapper class for RTTOV emissivity and BRDF atlases
///

#include <Atlas.h>

namespace rttov {

///@brief Atlas class constructor method
Atlas::Atlas() : verbose(true),
        atlas_wrap_id(-1),inc_land(true),inc_seaice(true),inc_sea(true)
{

}

///@brief Atlas class constructor method
///@param[in] is_verbose   verbosity flag
Atlas::Atlas(bool is_verbose) :
        atlas_wrap_id(-1),inc_land(true),inc_seaice(true),inc_sea(true)
{
    this->verbose = is_verbose;
}

/// @brief Deallocate Atlas structures
Atlas::~Atlas() {
    dropAtlas();
}

///@brief Return true if atlas has been loaded
bool Atlas::isAtlasLoaded() const {
    return (this->atlas_wrap_id > 0);
}

///@brief Return the path for the atlas files
const string& Atlas::getAtlasPath() const {
    return atlas_path;
}

///@brief Set the path for the atlas files
///@param[in] atlasPath   string containing path to directory containing atlas files
void Atlas::setAtlasPath(const string& atlasPath) {
    atlas_path = atlasPath;
}

///@brief Deallocate memory for the atlas
void Atlas::dropAtlas() {
    int err = 0;
    if (this->isAtlasLoaded()) rttov_drop_atlas_(&err, &atlas_wrap_id);
    this->atlas_wrap_id = -1;
    if (err != 0) throw runtime_error ("Error deallocating atlas");
    if (this->verbose) std::cerr << "Atlas deallocated" << std::endl;
}


bool Atlas::_loadBrdfAtlas(int month, int inst_id, int atlas_id) {
    if (this->isAtlasLoaded()) {
        if (this->verbose) std::cerr << "Warning: this atlas object already initialised" << std::endl;
        return false;
    }

    int len_path = this->atlas_path.length();
    if (len_path == 0) {
        if (this->verbose) std::cerr << "Warning: atlas_path not set, BRDF atlas not initialised" << std::endl;
        return false;
    }
    char path_str[len_path];
    for (int c = 0; c < len_path; c++) {
        path_str[c]=this->atlas_path[c];
    }

    rttov_load_brdf_atlas_ (
            &atlas_wrap_id,
            path_str,
            &month,
            &atlas_id,
            &inst_id,
            len_path);
    if (this->isAtlasLoaded()) {
        if (this->verbose) std::cerr << "BRDF atlas loaded successfully" << std::endl;
        return true;
    }
    else {
        if (this->verbose) std::cerr << "Warning: BRDF atlas not initialised" << std::endl;
        return false;
    }
}

///@brief Initialise the BRDF atlas for use with any instrument
///@param[in] month month (1-12 => Jan-Dec) for which atlas should be initialised
///@param[in] atlas_id ID for BRDF atlas, optional (default value : 1)
bool Atlas::loadBrdfAtlas(int month, int atlas_id) {
    return _loadBrdfAtlas(month, -1, atlas_id);
}

///@brief Initialise the BRDF atlas for a specific instrument
///@param[in] month month (1-12 => Jan-Dec) for which atlas should be initialised
///@param[in] rttov Rttov instance for which to initialise atlas, must be loaded
///@param[in] atlas_id ID for BRDF atlas, optional (default value : 1)
bool Atlas::loadBrdfAtlas(int month, rttov::Rttov * rttov, int atlas_id) {
    if (! rttov->isCoeffsLoaded()) {
        if (this->verbose) std::cerr << "Warning: instrument must be loaded, BRDF atlas not initialised" << std::endl;
        return false;
    }
    return _loadBrdfAtlas(month, rttov->getInstId(), atlas_id);
}

///@brief Initialise the BRDF atlas for a specific instrument
///@param[in] month month (1-12 => Jan-Dec) for which atlas should be initialised
///@param[in] rttov RttovSafe instance for which to initialise atlas, must be loaded
///@param[in] atlas_id ID for BRDF atlas, optional (default value : 1)
bool Atlas::loadBrdfAtlas(int month, rttov::RttovSafe * rttov, int atlas_id) {
    if (! rttov->isCoeffsLoaded()) {
        if (this->verbose) std::cerr << "Warning: instrument must be loaded, BRDF atlas not initialised" << std::endl;
        return false;
    }
    return _loadBrdfAtlas(month, rttov->getInstId(), atlas_id);
}

bool Atlas::_loadIrEmisAtlas(int month, int inst_id, bool ang_corr, int atlas_id) {
    if (this->isAtlasLoaded()) {
        if (this->verbose) std::cerr << "Warning: this atlas object already initialised" << std::endl;
        return false;
    }

    int len_path = this->atlas_path.length();
    if (len_path == 0) {
        if (this->verbose) std::cerr << "Warning: atlas_path not set, IR emissivity atlas not initialised" << std::endl;
        return false;
    }
    char path_str[len_path];
    for (int c = 0; c < len_path; c++) {
        path_str[c]=this->atlas_path[c];
    }

    int int_ang_corr = 0;
    if (ang_corr) int_ang_corr = 1;
    rttov_load_ir_emis_atlas_ (
            &atlas_wrap_id,
            path_str,
            &month,
            &atlas_id,
            &inst_id,
            &int_ang_corr,
            len_path);
    if (this->isAtlasLoaded()) {
        if (this->verbose) std::cerr << "IR emissivity atlas loaded successfully" << std::endl;
        return true;
    }
    else {
        if (this->verbose) std::cerr<< "Warning: IR emissivity atlas not initialised" << std::endl;
        return false;
    }
}

///@brief Initialise the IR emissivity atlas for use with any instrument
///@param[in] month month (1-12 => Jan-Dec) for which atlas should be initialised
///@param[in] ang_corr bool apply zenith angle correction to IR emissivities, optional
///@param[in] atlas_id ID for IR atlas, optional, (default value : 1)
bool Atlas::loadIrEmisAtlas(int month, bool ang_corr, int atlas_id) {
    return _loadIrEmisAtlas(month, -1, ang_corr, atlas_id);
}

///@brief Initialise the IR emissivity atlas for a specific instrument
///@param[in] month month (1-12 => Jan-Dec) for which atlas should be initialised
///@param[in] rttov Rttov instance for which to initialise atlas, must be loaded
///@param[in] ang_corr bool apply zenith angle correction to IR emissivities, optional
///@param[in] atlas_id ID for IR atlas, optional, (default value : 1)
bool Atlas::loadIrEmisAtlas(int month, rttov::Rttov * rttov, bool ang_corr, int atlas_id) {
    if (! rttov->isCoeffsLoaded()) {
        if (this->verbose) std::cerr << "Warning: instrument must be loaded, IR emissivity atlas not initialised" << std::endl;
        return false;
    }
    return _loadIrEmisAtlas(month, rttov->getInstId(), ang_corr, atlas_id);
}

///@brief Initialise the IR emissivity atlas for a specific instrument
///@param[in] month month (1-12 => Jan-Dec) for which atlas should be initialised
///@param[in] rttov RttovSafe instance for which to initialise atlas, must be loaded
///@param[in] ang_corr bool apply zenith angle correction to IR emissivities, optional
///@param[in] atlas_id ID for IR atlas, optional, (default value : 1)
bool Atlas::loadIrEmisAtlas(int month, rttov::RttovSafe * rttov, bool ang_corr, int atlas_id) {
    if (! rttov->isCoeffsLoaded()) {
        if (this->verbose) std::cerr << "Warning: instrument must be loaded, IR emissivity atlas not initialised" << std::endl;
        return false;
    }
    return _loadIrEmisAtlas(month, rttov->getInstId(), ang_corr, atlas_id);
}


bool Atlas::_loadMwEmisAtlas(int month, int inst_id, int year, int atlas_id) {
    if (this->isAtlasLoaded()) {
        if (this->verbose) std::cerr << "Warning: this atlas object already initialised" << std::endl;
        return false;
    }

    int len_path = this->atlas_path.length();
    if (len_path == 0) {
        if (this->verbose) std::cerr << "Warning: atlas_path not set, MW emissivity atlas not initialised" << std::endl;
        return false;
    }
    char path_str[len_path];
    for (int c = 0; c < len_path; c++) {
        path_str[c]=this->atlas_path[c];
    }

    rttov_load_mw_emis_atlas_ (
            &atlas_wrap_id,
            path_str,
            &month,
            &atlas_id,
            &inst_id,
            &year,
            len_path);
    if (this->isAtlasLoaded()) {
        if (this->verbose) std::cerr << "MW emissivity atlas loaded successfully" << std::endl;
        return true;
    }
    else {
        if (this->verbose) std::cerr << "Warning: MW emissivity atlas not initialised" << std::endl;
        return false;
    }
}

///@brief Initialise the MW emissivity atlas for use with any instrument (TELSEM2)
///@param[in] month month (1-12) of the atlas
///@param[in] atlas_id ID for MW atlas, optional (default value : 1)
bool Atlas::loadMwEmisAtlas(int month, int atlas_id) {
    return _loadMwEmisAtlas(month, -1, 0, atlas_id);
}

///@brief Initialise the MW emissivity atlas for a specific instrument (CNRM MW atlas)
///@param[in] month month (1-12) of the atlas
///@param[in] rttov Rttov instance for which to initialise atlas, must be loaded, required for CNRM atlas, ignored for TELSEM2
///@param[in] year year of the atlas (CNRM atlas only, uses default if value less than 1970)
///@param[in] atlas_id ID for MW atlas, optional (default value : 1)
bool Atlas::loadMwEmisAtlas(int month, rttov::Rttov * rttov, int year, int atlas_id) {
    if (atlas_id == 2 && ! rttov->isCoeffsLoaded()) {
        if (this->verbose) std::cerr << "Warning: instrument must be loaded, CNRM MW emissivity atlas not initialised" << std::endl;
        return false;
    }
    return _loadMwEmisAtlas(month, rttov->getInstId(), year, atlas_id);
}

///@brief Initialise the MW emissivity atlas for a specific instrument (CNRM MW atlas)
///@param[in] month month (1-12) of the atlas
///@param[in] rttov RttovSafe instance for which to initialise atlas, must be loaded, required for CNRM atlas, ignored for TELSEM2
///@param[in] year year of the atlas (CNRM atlas only, uses default if value less than 1970)
///@param[in] atlas_id ID for MW atlas, optional (default value : 1)
bool Atlas::loadMwEmisAtlas(int month, rttov::RttovSafe * rttov, int year, int atlas_id) {
    if (atlas_id == 2 && ! rttov->isCoeffsLoaded()) {
        if (this->verbose) std::cerr << "Warning: instrument must be loaded, CNRM MW emissivity atlas not initialised" << std::endl;
        return false;
    }
    return _loadMwEmisAtlas(month, rttov->getInstId(), year, atlas_id);
}

///@brief Initialise the MW emissivity atlas for a specific instrument (CNRM MW atlas)
///@param[in] month month (1-12) of the atlas
///@param[in] rttov RttovScatt instance for which to initialise atlas, must be loaded, required for CNRM atlas, ignored for TELSEM2
///@param[in] year year of the atlas (CNRM atlas only, uses default if value less than 1970)
///@param[in] atlas_id ID for MW atlas, optional (default value : 1)
bool Atlas::loadMwEmisAtlas(int month, rttov::RttovScatt * rttov, int year, int atlas_id) {
    if (atlas_id == 2 && ! rttov->isCoeffsLoaded()) {
        if (this->verbose) std::cerr << "Warning: instrument must be loaded, CNRM MW emissivity atlas not initialised" << std::endl;
        return false;
    }
    return _loadMwEmisAtlas(month, rttov->getInstId(), year, atlas_id);
}

///@brief Initialise the MW emissivity atlas for a specific instrument (CNRM MW atlas)
///@param[in] month month (1-12) of the atlas
///@param[in] rttov RttovScattSafe instance for which to initialise atlas, must be loaded, required for CNRM atlas, ignored for TELSEM2
///@param[in] year year of the atlas (CNRM atlas only, uses default if value less than 1970)
///@param[in] atlas_id ID for MW atlas, optional (default value : 1)
bool Atlas::loadMwEmisAtlas(int month, rttov::RttovScattSafe * rttov, int year, int atlas_id) {
    if (atlas_id == 2 && ! rttov->isCoeffsLoaded()) {
        if (this->verbose) std::cerr << "Warning: instrument must be loaded, CNRM MW emissivity atlas not initialised" << std::endl;
        return false;
    }
    return _loadMwEmisAtlas(month, rttov->getInstId(), year, atlas_id);
}


///@brief Return emissivities/BRDFs
///@param[in,out] emisBrdf   pointer to returned emissivity/BRDF array
///@param[in]     rttov      Rttov instance for which to return emissivities
///@param[in]     channels   vector containing the list of channels to simulate (indices with respect to channel list loaded from coef file), optional
void Atlas::fillEmisBrdf(double* emisBrdf, rttov::Rttov * rttov, const vector<int>& channels) {

    if (! this->isAtlasLoaded()) throw logic_error( "Error: atlas data not loaded");
    if (! rttov->isCoeffsLoaded()) throw logic_error( "Error: coefficients not loaded in instrument");
    if (! rttov->isProfileSet()) throw logic_error( "Error: profiles not set for instrument");

    int inst_id = rttov->getInstId();
    int nprofiles = rttov->getNprofiles();

    int nchannels = channels.size();
    if (nchannels == 0) nchannels = rttov->getNchannels();
    int channel_list[nchannels];
    if (channels.size() == 0) {
        for (int i = 0; i < nchannels; i++) channel_list[i] = i+1;
    } else {
        for (int i = 0; i < nchannels; i++) {
            if (channels[i] < 1 || channels[i] > nchannels) throw runtime_error( "Error: channel index out of range" );
            channel_list[i] = channels[i];
        }
    }

    double latitude[nprofiles];
    double longitude[nprofiles];
    int surftypeval[nprofiles];
    int watertype[nprofiles];
    double zenangle[nprofiles];
    double azangle[nprofiles];
    double sunzenangle[nprofiles];
    double sunazangle[nprofiles];
    double snow_fraction[nprofiles];

    rttov::Profiles * profiles = rttov->getProfiles();
    int * surftype    = profiles->getSurfType();
    double * surfgeom = profiles->getSurfGeom();
    double * angles   = profiles->getAngles();
    double * skin     = profiles->getSkin();
    for (int p = 0; p < nprofiles; p++) {
        latitude[p]      = surfgeom[p*3];
        longitude[p]     = surfgeom[p*3+1];
        surftypeval[p]   = surftype[p*2];
        watertype[p]     = surftype[p*2+1];
        zenangle[p]      = angles[p*4];
        azangle[p]       = angles[p*4+1];
        sunzenangle[p]   = angles[p*4+2];
        sunazangle[p]    = angles[p*4+3];
        snow_fraction[p] = skin[p*9+2];
    }
    double thisEmisBrdf[nprofiles][nchannels];

    int err;
    rttov_get_emisbrdf_(
        &err,
        &atlas_wrap_id,
        latitude,
        longitude,
        surftypeval,
        watertype,
        zenangle,
        azangle,
        sunzenangle,
        sunazangle,
        snow_fraction,
        &inst_id,
        channel_list,
        (double *)thisEmisBrdf,
        &nchannels,
        &nprofiles);

    // For the requested surface types only copy data into output array if err==0, leave data untouched on error
    if (err == 0) {
        for (int p = 0; p < nprofiles; p++) {
            if ((surftypeval[p] == 0 && this->inc_land) ||
                (surftypeval[p] == 1 && this->inc_sea)  ||
                (surftypeval[p] == 2 && this->inc_seaice)) {
                for (int c = 0; c < nchannels; c++) emisBrdf[p*nchannels+c] = thisEmisBrdf[p][c];
            }
        }
    }
    if (err != 0) throw runtime_error ("Error in rttov_get_emisbrdf");
}

///@brief Return emissivities/BRDFs
///@param[in,out] emisBrdf   pointer to returned emissivity/BRDF array
///@param[in]     rttov      RttovSafe instance for which to return emissivities
///@param[in]     channels   vector containing the list of channels to simulate (indices with respect to channel list loaded from coef file), optional
void Atlas::fillEmisBrdf(double* emisBrdf, rttov::RttovSafe * rttov, const vector<int>& channels) {

    if (! this->isAtlasLoaded()) throw logic_error( "Error: atlas data not loaded");
    if (! rttov->isCoeffsLoaded()) throw logic_error( "Error: coefficients not loaded in instrument");
    if (! rttov->isProfileSet()) throw logic_error( "Error: profiles not set for instrument");

    int inst_id = rttov->getInstId();
    int nprofiles = rttov->getNprofiles();

    int nchannels = channels.size();
    if (nchannels == 0) nchannels = rttov->getNchannels();
    int channel_list[nchannels];
    if (channels.size() == 0) {
        for (int i = 0; i < nchannels; i++) channel_list[i] = i+1;
    } else {
        for (int i = 0; i < nchannels; i++) {
            if (channels[i] < 1 || channels[i] > nchannels) throw runtime_error( "Error: channel index out of range" );
            channel_list[i] = channels[i];
        }
    }

    double latitude[nprofiles];
    double longitude[nprofiles];
    int surftypeval[nprofiles];
    int watertype[nprofiles];
    double zenangle[nprofiles];
    double azangle[nprofiles];
    double sunzenangle[nprofiles];
    double sunazangle[nprofiles];
    double snow_fraction[nprofiles];

    std::vector <rttov::Profile> * profiles = rttov->getTheProfiles();
    for (int p = 0; p < nprofiles; p++) {
        const std::vector <int>&    surftype = (*profiles)[p].getSurfType();
        const std::vector <double>& surfgeom = (*profiles)[p].getSurfGeom();
        const std::vector <double>& angles   = (*profiles)[p].getAngles();
        const std::vector <double>& skin     = (*profiles)[p].getSkin();
        latitude[p]      = surfgeom[0];
        longitude[p]     = surfgeom[1];
        surftypeval[p]   = surftype[0];
        watertype[p]     = surftype[1];
        zenangle[p]      = angles[0];
        azangle[p]       = angles[1];
        sunzenangle[p]   = angles[2];
        sunazangle[p]    = angles[3];
        snow_fraction[p] = skin[2];
    }
    double thisEmisBrdf[nprofiles][nchannels];

    int err;
    rttov_get_emisbrdf_(
        &err,
        &atlas_wrap_id,
        latitude,
        longitude,
        surftypeval,
        watertype,
        zenangle,
        azangle,
        sunzenangle,
        sunazangle,
        snow_fraction,
        &inst_id,
        channel_list,
        (double *)thisEmisBrdf,
        &nchannels,
        &nprofiles);

    // For the requested surface types only copy data into output array if err==0, leave data untouched on error
    if (err == 0) {
        for (int p = 0; p < nprofiles; p++) {
            if ((surftypeval[p] == 0 && this->inc_land) ||
                (surftypeval[p] == 1 && this->inc_sea)  ||
                (surftypeval[p] == 2 && this->inc_seaice)) {
                for (int c = 0; c < nchannels; c++) emisBrdf[p*nchannels+c] = thisEmisBrdf[p][c];
            }
        }
    }
    if (err != 0) throw runtime_error ("Error in rttov_get_emisbrdf");
}

///@brief Return emissivities
///@param[in,out] emisBrdf   pointer to returned emissivity/BRDF array
///@param[in]     rttov      RttovScatt instance for which to return emissivities
///@param[in]     channels   vector containing the list of channels to simulate (indices with respect to channel list loaded from coef file), optional
void Atlas::fillEmisBrdf(double* emisBrdf, rttov::RttovScatt * rttov, const vector<int>& channels) {

    if (! this->isAtlasLoaded()) throw logic_error( "Error: atlas data not loaded");
    if (! rttov->isCoeffsLoaded()) throw logic_error( "Error: coefficients not loaded in instrument");
    if (! rttov->isProfileSet()) throw logic_error( "Error: profiles not set for instrument");

    int inst_id = rttov->getInstId();
    int nprofiles = rttov->getNprofiles();

    int nchannels = channels.size();
    if (nchannels == 0) nchannels = rttov->getNchannels();
    int channel_list[nchannels];
    if (channels.size() == 0) {
        for (int i = 0; i < nchannels; i++) channel_list[i] = i+1;
    } else {
        for (int i = 0; i < nchannels; i++) {
            if (channels[i] < 1 || channels[i] > nchannels) throw runtime_error( "Error: channel index out of range" );
            channel_list[i] = channels[i];
        }
    }

    double latitude[nprofiles];
    double longitude[nprofiles];
    int surftypeval[nprofiles];
    int watertype[nprofiles];
    double zenangle[nprofiles];
    double azangle[nprofiles];
    double sunzenangle[nprofiles];
    double sunazangle[nprofiles];
    double snow_fraction[nprofiles];

    rttov::ProfilesScatt * profiles = rttov->getProfiles();
    int * surftype    = profiles->getSurfType();
    double * surfgeom = profiles->getSurfGeom();
    double * angles   = profiles->getAngles();
    for (int p = 0; p < nprofiles; p++) {
        latitude[p]      = surfgeom[p*3];
        longitude[p]     = surfgeom[p*3+1];
        surftypeval[p]   = surftype[p];
        watertype[p]     = 1;
        zenangle[p]      = angles[p*2];
        azangle[p]       = angles[p*2+1];
        sunzenangle[p]   = 0.;
        sunazangle[p]    = 0.;
        snow_fraction[p] = 0.;
    }
    double thisEmisBrdf[nprofiles][nchannels];

    int err;
    rttov_get_emisbrdf_(
        &err,
        &atlas_wrap_id,
        latitude,
        longitude,
        surftypeval,
        watertype,
        zenangle,
        azangle,
        sunzenangle,
        sunazangle,
        snow_fraction,
        &inst_id,
        channel_list,
        (double *)thisEmisBrdf,
        &nchannels,
        &nprofiles);

    // For the requested surface types only copy data into output array if err==0, leave data untouched on error
    if (err == 0) {
        for (int p = 0; p < nprofiles; p++) {
            if ((surftypeval[p] == 0 && this->inc_land) ||
                (surftypeval[p] == 1 && this->inc_sea)  ||
                (surftypeval[p] == 2 && this->inc_seaice)) {
                for (int c = 0; c < nchannels; c++) emisBrdf[p*nchannels+c] = thisEmisBrdf[p][c];
            }
        }
    }
    if (err != 0) throw runtime_error ("Error in rttov_get_emisbrdf");
}

///@brief Return emissivities
///@param[in,out] emisBrdf   pointer to returned emissivity/BRDF array
///@param[in]     rttov      RttovScattSafe instance for which to return emissivities
///@param[in]     channels   vector containing the list of channels to simulate (indices with respect to channel list loaded from coef file), optional
void Atlas::fillEmisBrdf(double* emisBrdf, rttov::RttovScattSafe * rttov, const vector<int>& channels) {

    if (! this->isAtlasLoaded()) throw logic_error( "Error: atlas data not loaded");
    if (! rttov->isCoeffsLoaded()) throw logic_error( "Error: coefficients not loaded in instrument");
    if (! rttov->isProfileSet()) throw logic_error( "Error: profiles not set for instrument");

    int inst_id = rttov->getInstId();
    int nprofiles = rttov->getNprofiles();

    int nchannels = channels.size();
    if (nchannels == 0) nchannels = rttov->getNchannels();
    int channel_list[nchannels];
    if (channels.size() == 0) {
        for (int i = 0; i < nchannels; i++) channel_list[i] = i+1;
    } else {
        for (int i = 0; i < nchannels; i++) {
            if (channels[i] < 1 || channels[i] > nchannels) throw runtime_error( "Error: channel index out of range" );
            channel_list[i] = channels[i];
        }
    }

    double latitude[nprofiles];
    double longitude[nprofiles];
    int surftypeval[nprofiles];
    int watertype[nprofiles];
    double zenangle[nprofiles];
    double azangle[nprofiles];
    double sunzenangle[nprofiles];
    double sunazangle[nprofiles];
    double snow_fraction[nprofiles];

    std::vector <rttov::ProfileScatt> * profiles = rttov->getTheProfiles();
    for (int p = 0; p < nprofiles; p++) {
        int surftype = (*profiles)[p].getSurfType();
        const std::vector <double>& surfgeom = (*profiles)[p].getSurfGeom();
        const std::vector <double>& angles   = (*profiles)[p].getAngles();
        latitude[p]      = surfgeom[0];
        longitude[p]     = surfgeom[1];
        surftypeval[p]   = surftype;
        watertype[p]     = 1;
        zenangle[p]      = angles[0];
        azangle[p]       = angles[1];
        sunzenangle[p]   = 0.;
        sunazangle[p]    = 0.;
        snow_fraction[p] = 0.;
    }
    double thisEmisBrdf[nprofiles][nchannels];

    int err;
    rttov_get_emisbrdf_(
        &err,
        &atlas_wrap_id,
        latitude,
        longitude,
        surftypeval,
        watertype,
        zenangle,
        azangle,
        sunzenangle,
        sunazangle,
        snow_fraction,
        &inst_id,
        channel_list,
        (double *)thisEmisBrdf,
        &nchannels,
        &nprofiles);

    // For the requested surface types only copy data into output array if err==0, leave data untouched on error
    if (err == 0) {
        for (int p = 0; p < nprofiles; p++) {
            if ((surftypeval[p] == 0 && this->inc_land) ||
                (surftypeval[p] == 1 && this->inc_sea)  ||
                (surftypeval[p] == 2 && this->inc_seaice)) {
                for (int c = 0; c < nchannels; c++) emisBrdf[p*nchannels+c] = thisEmisBrdf[p][c];
            }
        }
    }
    if (err != 0) throw runtime_error ("Error in rttov_get_emisbrdf");
}


///@brief Return the inc_land boolean
bool Atlas::getIncLand() const {
    return inc_land;
}

///@brief Return the inc_seaice boolean
bool Atlas::getIncSeaIce() const {
    return inc_seaice;
}

///@brief Return the inc_sea boolean
bool Atlas::getIncSea() const {
    return inc_sea;
}

///@brief Set the inc_land boolean
///@param[in] incLand   value to assign to inc_land
void Atlas::setIncLand(bool incLand) {
    inc_land = incLand;
}

///@brief Set the inc_seaice boolean
///@param[in] incSeaice   value to assign to inc_seaice
void Atlas::setIncSeaIce(bool incSeaice) {
    inc_seaice = incSeaice;
}

///@brief Set the inc_sea boolean
///@param[in] incSea   value to assign to inc_sea
void Atlas::setIncSea(bool incSea) {
    inc_sea = incSea;
}

///@brief Set the verbose boolean
///@param[in] is_verbose   value to assign to verbose
void Atlas::setVerbose(bool is_verbose) {
    verbose = is_verbose;
}


} /* namespace rttov */
