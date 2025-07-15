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

///@class Profile
///  This class represents emissivity or BRDF atlas data for one month and, where relevant, one instrument

#ifndef ATLAS_H_
#define ATLAS_H_

using namespace std;
#include <rttov_cc_interface.h>

#include <string>
#include <vector>
#include <iostream>
#include <string>
#include <exception>
#include <stdexcept>

#include <Rttov.h>
#include <RttovSafe.h>
#include <RttovScatt.h>
#include <RttovScattSafe.h>
namespace rttov{

class Atlas {
public:
    Atlas();
    Atlas(bool verbose);
    virtual ~Atlas();

    const string& getAtlasPath() const;
    void setAtlasPath(const string& atlasPath);

    bool isAtlasLoaded() const;

    void setVerbose(bool verbose);

    void setIncLand(bool incLand);
    void setIncSeaIce(bool incSeaIce);
    void setIncSea(bool incSea);

    bool getIncLand() const;
    bool getIncSeaIce() const;
    bool getIncSea() const;

    bool loadBrdfAtlas(int month, int atlas_id=-1);
    bool loadBrdfAtlas(int month, rttov::Rttov * rttov, int atlas_id=-1);
    bool loadBrdfAtlas(int month, rttov::RttovSafe * rttov, int atlas_id=-1);

    bool loadIrEmisAtlas(int month, bool ang_corr=false, int atlas_id=-1);
    bool loadIrEmisAtlas(int month, rttov::Rttov * rttov, bool ang_corr=false, int atlas_id=-1);
    bool loadIrEmisAtlas(int month, rttov::RttovSafe * rttov, bool ang_corr=false, int atlas_id=-1);

    bool loadMwEmisAtlas(int month, int atlas_id=-1);
    bool loadMwEmisAtlas(int month, rttov::Rttov * rttov, int year=0, int atlas_id=-1);
    bool loadMwEmisAtlas(int month, rttov::RttovSafe * rttov, int year=0, int atlas_id=-1);
    bool loadMwEmisAtlas(int month, rttov::RttovScatt * rttov, int year=0, int atlas_id=-1);
    bool loadMwEmisAtlas(int month, rttov::RttovScattSafe * rttov, int year=0, int atlas_id=-1);

    void fillEmisBrdf(double * emisBrdf, rttov::Rttov * rttov, const vector<int>& channels=vector<int>{});
    void fillEmisBrdf(double * emisBrdf, rttov::RttovSafe * rttov, const vector<int>& channels=vector<int>{});
    void fillEmisBrdf(double * emisBrdf, rttov::RttovScatt * rttov, const vector<int>& channels=vector<int>{});
    void fillEmisBrdf(double * emisBrdf, rttov::RttovScattSafe * rttov, const vector<int>& channels=vector<int>{});

    void dropAtlas();

protected:

    bool _loadBrdfAtlas(int month, int inst_id, int atlas_id);
    bool _loadIrEmisAtlas(int month, int inst_id, bool ang_corr, int atlas_id);
    bool _loadMwEmisAtlas(int month, int inst_id, int year, int atlas_id);

    string atlas_path;
    int atlas_wrap_id;
    bool verbose;

    bool inc_land;
    bool inc_seaice;
    bool inc_sea;
    };

} /* namespace rttov */

#endif /* ATLAS_H_ */
