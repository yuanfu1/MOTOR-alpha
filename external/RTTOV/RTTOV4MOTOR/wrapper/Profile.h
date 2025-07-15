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

///@class Profile
///  This class represents one atmospheric profile

#ifndef PROFILE_H_
#define PROFILE_H_
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <string>
#include <exception>
#include <stdexcept>

#include <Rttov_common.h>
namespace rttov{

class Profile {
public:
    Profile(int nlevels);
    virtual ~Profile();

    void setGasUnits(rttov::gasUnitType gasUnits);
    void setMmrCldAer(const bool mmrCldAer);

    void setP(const std::vector<double>& p);
    void setT(const std::vector<double>& t);
    void setQ (const std::vector <double>& q);
    void setO3 (const std::vector <double>& o3);
    void setCO2 (const std::vector <double>& co2);
    void setN2O (const std::vector <double>& n2o);
    void setCO (const std::vector <double>& co);
    void setCH4 (const std::vector <double>& ch4);
    void setSO2 (const std::vector <double>& so2);
    void setCLW (const std::vector <double>& clw);
    void setCfrac (const std::vector <double>& cfrac);
    void setStco (const std::vector <double>& stco);
    void setStma (const std::vector <double>& stma);
    void setCucc (const std::vector <double>& cucc);
    void setCucp (const std::vector <double>& cucp);
    void setCuma (const std::vector <double>& cuma);
    void setCirr (const std::vector <double>& cirr);
    void setClwde (const std::vector <double>& clwde);
    void setIcede (const std::vector <double>& icede);
    void setInso (const std::vector <double>& inso);
    void setWaso (const std::vector <double>& waso);
    void setSoot (const std::vector <double>& soot);
    void setSsam (const std::vector <double>& ssam);
    void setSscm (const std::vector <double>& sscm);
    void setMinm (const std::vector <double>& minm);
    void setMiam (const std::vector <double>& miam);
    void setMicm (const std::vector <double>& micm);
    void setMitr (const std::vector <double>& mitr);
    void setSuso (const std::vector <double>& suso);
    void setVola (const std::vector <double>& vola);
    void setVapo (const std::vector <double>& vapo);
    void setAsdu (const std::vector <double>& asdu);
    void setBcar (const std::vector <double>& bcar);
    void setDus1 (const std::vector <double>& dus1);
    void setDus2 (const std::vector <double>& dus2);
    void setDus3 (const std::vector <double>& dus3);
    void setSulp (const std::vector <double>& sulp);
    void setSsa1 (const std::vector <double>& ssa1);
    void setSsa2 (const std::vector <double>& ssa2);
    void setSsa3 (const std::vector <double>& ssa3);
    void setOmat (const std::vector <double>& omat);
    void setUserAerN (const std::vector <double>& aer, const int n);

    void setAngles(const double satzen, const double satazi, const double sunzen,
        const double sunazi);
    void setS2m(const double p_2m, const double t_2m, const double q_2m,
        const double u_10m, const double v_10m, const double wind_fetch );
    void setSkin(const double t, const double salinity, const double snow_fraction,
        const double foam_fraction, const double fastem_coef_1,
        const double fastem_coef_2, const double fastem_coef_3,
        const double fastem_coef_4, const double fastem_coef_5);
    void setSurfType(const int surftype, const int watertype);
    void setSurfGeom(const double lat, const double lon, const double elevation);
    void setDateTimes(const int yy, const int mm, const int dd, const int hh,
        const int mn, const int ss);
    void setSimpleCloud(const double ctp, const double cfraction);
    void setClwScheme(const int clw_scheme, const int clwde_param);
    void setIceCloud(const int ice_scheme, const int icede_param);
    void setZeeman(const double Be, const double cosbk);


    int getNlevels() const;
    int getGasUnits() const;
    bool isMmrCldAer() const;
    const std::vector<double>& getP() const;
    const std::vector<double>& getT() const;
    const std::vector<double>& getAngles() const;
    const std::vector<double>& getSurfGeom() const;
    const std::vector<int>& getSurfType() const;
    const std::vector<double>& getS2m() const;
    const std::vector<double>& getSkin() const;
    const std::vector<int>& getDateTimes() const;
    const std::vector<double>& getSimpleCloud() const;
    const std::vector<int>& getClwScheme() const;
    const std::vector<int>& getIceCloud() const;
    const std::vector<double>& getZeeman() const;
    const std::vector<itemIdType> getItemContents() const;
    const ItemIdPointerMap& getItems() const;

    bool isDefaultPressureLevels() const;
    bool check();

private :
    bool verbose;
    int nlevels;
    int gas_units;
    int mmr_cldaer;
    std::vector <double> P;
    std::vector <double> T;
    ItemIdPointerMap items;

    // datetimes: yy, mm, dd, hh, mm, ss
    std::vector <int> datetimes;

    // surfgeom: lat, lon, elev
    std::vector <double> surfgeom;

    // angles: satzen, satazi, sunzen, sunazi
    std::vector <double> angles;

    // surftype: surftype, watertype
    std::vector <int> surftype;

    // skin: skin T, salinity, snow_frac, foam_frac, fastem_coefsx5
    std::vector <double> skin;

    // s2m: 2m p, 2m t, 2m q, 10m wind u, v, wind fetch
    std::vector <double> s2m;

    // simplecloud: ctp, cfraction
    std::vector <double> simplecloud;

    // clwscheme: clw_scheme, clwde_param
    std::vector <int> clwscheme;

    // icecloud: ice_scheme, icede_param
    std::vector <int> icecloud;

    // zeeman: be, cosbk
    std::vector <double> zeeman;
};

} /* namespace rttov */

#endif /* PROFILE_H_ */
