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

///@file Profile.cpp
///
///@brief Method definition for Profile class

#include <Profile.h>

namespace rttov {

///@brief Constructor method
Profile::Profile(int nlevels):gas_units(unknown),verbose(true),angles{},datetimes{},s2m{},
        skin{},surfgeom{},surftype{},simplecloud{},zeeman{},clwscheme{},icecloud{},mmr_cldaer(1),
        P{},T{},items{} {
            // constructor
            this->nlevels=nlevels;
        }

        Profile::~Profile() {

        }

        ///@brief Set the gas_units
        ///@param gasUnits unknown=-1,ppmv_dry=0,kg_per_kg=1,ppmv_wet=2
        void Profile::setGasUnits(rttov::gasUnitType gasUnits) {
            gas_units = gasUnits;
        }

        ///@brief Set the mmr_cldaer flag
        ///@param mmrCldAer false => cm-3 (aer), g/m3 (cld); true => kg/kg (cld, aer)
        void Profile::setMmrCldAer(bool mmrCldAer) {
            mmr_cldaer = mmrCldAer;
        }

        ///@brief Set the p (pressure) vector
        ///@param p reference to a vector containing the pressures (size nlevels)
        void Profile::setP(const std::vector<double>& p) {
            if ( p.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for p");
            P = p;
        }

        ///@brief Set the temperatures vector
        ///@param t reference to a vector containing the temperatures (size nlevels)
        void Profile::setT(const std::vector<double>& t) {
            if ( t.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for t");
            T = t;
        }

        ///@brief Set satellite an solar angles
        ///@param satzen, satazi, sunzen, sunazi
        void Profile::setAngles(const double satzen, const double satazi,
                const double sunzen, const double sunazi) {

            if ( ! angles.empty()) angles.erase(angles.begin(),angles.end());
            angles.push_back(satzen);
            angles.push_back(satazi);
            angles.push_back(sunzen);
            angles.push_back(sunazi);
        }

        ///@brief Set date and time
        ///@param yy, mm, dd, hh, mn, ss
        void Profile::setDateTimes(const int yy, const int mm, const int dd,
                const int hh, const int mn, const int ss) {

            if ( ! datetimes.empty()) datetimes.erase(datetimes.begin(),datetimes.end());
            datetimes.push_back(yy);
            datetimes.push_back(mm);
            datetimes.push_back(dd);
            datetimes.push_back(hh);
            datetimes.push_back(mn);
            datetimes.push_back(ss);
        }

        ///@brief Set clw scheme parameters
        ///@param clw_scheme, clwde_param
        void Profile::setClwScheme(const int clw_scheme, const int clwde_param) {

            if ( !clwscheme.empty() ) clwscheme.erase(clwscheme.begin(),clwscheme.end());
            clwscheme.push_back(clw_scheme);
            clwscheme.push_back(clwde_param);
        }

        ///@brief Set ice cloud parameters
        ///@param ice_scheme, icede_param
        void Profile::setIceCloud(const int ice_scheme, const int icede_param) {

            if ( !icecloud.empty() ) icecloud.erase(icecloud.begin(),icecloud.end());
            icecloud.push_back(ice_scheme);
            icecloud.push_back(icede_param);
        }

        ///@brief Set simple cloud parameters
        ///@param ctp, cfraction
        void Profile::setSimpleCloud(const double ctp, const double cfraction) {

            if ( !simplecloud.empty()) simplecloud.erase(simplecloud.begin(),simplecloud.end());
            simplecloud.push_back(ctp);
            simplecloud.push_back(cfraction);
        }

        ///@brief Set surface geometry parameters
        ///@param lat, lon, elevation
        void Profile::setSurfGeom(const double lat, const double lon,
                const double elevation) {

            if ( ! surfgeom.empty()) surfgeom.erase(surfgeom.begin(),surfgeom.end());
            surfgeom.push_back(lat);
            surfgeom.push_back(lon);
            surfgeom.push_back(elevation);
        }

        ///@brief Set surface type parameters
        ///@param surftype, watertype
        void Profile::setSurfType(const int surftype, const int watertype) {

            if ( ! this->surftype.empty())  this->surftype.erase( this->surftype.begin(), this->surftype.end());
            this->surftype.push_back(surftype);
            this->surftype.push_back(watertype);
        }

        ///@brief Set zeeman parameters
        ///@param Be, cosbk
        void Profile::setZeeman(const double Be, const double cosbk) {

            if ( ! zeeman.empty()) zeeman.erase(zeeman.begin(),zeeman.end());
            zeeman.push_back(Be);
            zeeman.push_back(cosbk);
        }

        ///@brief Set surface 2m and 10m parameters
        ///@param p_2m, t_2m, q_2m, u_10m, v_10m, wind_fetch
        void Profile::setS2m(const double p_2m, const double t_2m, const double q_2m,
                const double u_10m, const double v_10m, const double wind_fetch) {

            if (! s2m.empty()) s2m.erase(s2m.begin(),s2m.end());
            s2m.push_back(p_2m);
            s2m.push_back(t_2m);
            s2m.push_back(q_2m);
            s2m.push_back(u_10m);
            s2m.push_back(v_10m);
            s2m.push_back(wind_fetch);
        }

        ///@brief Set skin parameters
        ///@param t, salinity, snow_fraction, foam_fraction, fastem_coef_1
        ///@param fastem_coef_2, fastem_coef_3, fastem_coef_4, fastem_coef_5
        void Profile::setSkin(const double t, const double salinity, const double snow_fraction,
                const double foam_fraction, const double fastem_coef_1,
                const double fastem_coef_2, const double fastem_coef_3,
                const double fastem_coef_4, const double fastem_coef_5){
            if (! skin.empty()) skin.erase(skin.begin(),skin.end());
            skin.push_back(t);
            skin.push_back(salinity);
            skin.push_back(snow_fraction);
            skin.push_back(foam_fraction);
            skin.push_back(fastem_coef_1);
            skin.push_back(fastem_coef_2);
            skin.push_back(fastem_coef_3);
            skin.push_back(fastem_coef_4);
            skin.push_back(fastem_coef_5);
        }

        ///@brief Set item q for the profile (vector size must equal nlevels)
        void Profile::setQ (const std::vector <double>& q) {
            if ( q.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for q");
            this->items[Q]=q;
        }

        ///@brief Set item o3 for the profile (vector size must equal nlevels)
        void Profile::setO3 (const std::vector <double>& o3) {
            if ( o3.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for o3");
            this->items[O3]=o3;
        }

        ///@brief Set item co2 for the profile (vector size must equal nlevels)
        void Profile::setCO2 (const std::vector <double>& co2) {
            if ( co2.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for co2");
            this->items[CO2]=co2;
        }

        ///@brief Set item n2o for the profile (vector size must equal nlevels)
        void Profile::setN2O (const std::vector <double>& n2o) {
            if ( n2o.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for n2o");
            this->items[N2O]=n2o;
        }

        ///@brief Set item co for the profile (vector size must equal nlevels)
        void Profile::setCO (const std::vector <double>& co) {
            if ( co.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for co");
            this->items[CO]=co;
        }

        ///@brief Set item ch4 for the profile (vector size must equal nlevels)
        void Profile::setCH4 (const std::vector <double>& ch4) {
            if ( ch4.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for ch4");
            this->items[CH4]=ch4;
        }

        ///@brief Set item so2 for the profile (vector size must equal nlevels)
        void Profile::setSO2 (const std::vector <double>& so2) {
            if ( so2.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for so2");
            this->items[SO2]=so2;
        }

        ///@brief Set item clw for the profile (vector size must equal nlevels)
        void Profile::setCLW (const std::vector <double>& clw) {
            if ( clw.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for clw");
            this->items[CLW]=clw;
        }

        ///@brief Set item cfrac for the profile (vector size must equal nlevels)
        void Profile::setCfrac (const std::vector <double>& cfrac) {
            if ( cfrac.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for cfrac");
            this->items[CFRAC]=cfrac;
        }

        ///@brief Set item stco for the profile (vector size must equal nlevels)
        void Profile::setStco (const std::vector <double>& stco) {
            if ( stco.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for stco");
            this->items[STCO]=stco;
        }

        ///@brief Set item stma for the profile (vector size must equal nlevels)
        void Profile::setStma (const std::vector <double>& stma) {
            if ( stma.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for stma");
            this->items[STMA]=stma;
        }

        ///@brief Set item cucc for the profile (vector size must equal nlevels)
        void Profile::setCucc (const std::vector <double>& cucc) {
            if ( cucc.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for cucc");
            this->items[CUCC]=cucc;
        }

        ///@brief Set item cucp for the profile (vector size must equal nlevels)
        void Profile::setCucp (const std::vector <double>& cucp) {
            if ( cucp.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for cucp");
            this->items[CUCP]=cucp;
        }

        ///@brief Set item cuma for the profile (vector size must equal nlevels)
        void Profile::setCuma (const std::vector <double>& cuma) {
            if ( cuma.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for cuma");
            this->items[CUMA]=cuma;
        }

        ///@brief Set item cirr for the profile (vector size must equal nlevels)
        void Profile::setCirr (const std::vector <double>& cirr) {
            if ( cirr.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for cirr");
            this->items[CIRR]=cirr;
        }

        ///@brief Set item clwde for the profile (vector size must equal nlevels)
        void Profile::setClwde (const std::vector <double>& clwde) {
            if ( clwde.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for clwde");
            this->items[CLWDE]=clwde;
        }

        ///@brief Set item icede for the profile (vector size must equal nlevels)
        void Profile::setIcede (const std::vector <double>& icede) {
            if ( icede.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for icede");
            this->items[ICEDE]=icede;
        }

        ///@brief Set item inso for the profile (vector size must equal nlevels)
        void Profile::setInso (const std::vector <double>& inso) {
            if ( inso.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for inso");
            this->items[INSO]=inso;
        }

        ///@brief Set item waso for the profile (vector size must equal nlevels)
        void Profile::setWaso (const std::vector <double>& waso) {
            if ( waso.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for waso");
            this->items[WASO]=waso;
        }

        ///@brief Set item soot for the profile (vector size must equal nlevels)
        void Profile::setSoot (const std::vector <double>& soot) {
            if ( soot.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for soot");
            this->items[SOOT]=soot;
        }

        ///@brief Set item ssam for the profile (vector size must equal nlevels)
        void Profile::setSsam (const std::vector <double>& ssam) {
            if ( ssam.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for ssam");
            this->items[SSAM]=ssam;
        }

        ///@brief Set item sscm for the profile (vector size must equal nlevels)
        void Profile::setSscm (const std::vector <double>& sscm) {
            if ( sscm.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for sscm");
            this->items[SSCM]=sscm;
        }

        ///@brief Set item minm for the profile (vector size must equal nlevels)
        void Profile::setMinm (const std::vector <double>& minm) {
            if ( minm.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for minm");
            this->items[MINM]=minm;
        }

        ///@brief Set item miam for the profile (vector size must equal nlevels)
        void Profile::setMiam (const std::vector <double>& miam) {
            if ( miam.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for miam");
            this->items[MIAM]=miam;
        }

        ///@brief Set item micm for the profile (vector size must equal nlevels)
        void Profile::setMicm (const std::vector <double>& micm) {
            if ( micm.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for micm");
            this->items[MICM]=micm;
        }

        ///@brief Set item mitr for the profile (vector size must equal nlevels)
        void Profile::setMitr (const std::vector <double>& mitr) {
            if ( mitr.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for mitr");
            this->items[MITR]=mitr;
        }

        ///@brief Set item suso for the profile (vector size must equal nlevels)
        void Profile::setSuso (const std::vector <double>& suso) {
            if ( suso.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for suso");
            this->items[SUSO]=suso;
        }

        ///@brief Set item vola for the profile (vector size must equal nlevels)
        void Profile::setVola (const std::vector <double>& vola) {
            if ( vola.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for vola");
            this->items[VOLA]=vola;
        }

        ///@brief Set item vapo for the profile (vector size must equal nlevels)
        void Profile::setVapo (const std::vector <double>& vapo) {
            if ( vapo.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for vapo");
            this->items[VAPO]=vapo;
        }

        ///@brief Set item asdu for the profile (vector size must equal nlevels)
        void Profile::setAsdu (const std::vector <double>& asdu) {
            if ( asdu.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for asdu");
            this->items[ASDU]=asdu;
        }

        ///@brief Set item bcar for the profile (vector size must equal nlevels)
        void Profile::setBcar (const std::vector <double>& bcar) {
            if ( bcar.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for bcar");
            this->items[BCAR]=bcar;
        }

        ///@brief Set item dus1 for the profile (vector size must equal nlevels)
        void Profile::setDus1 (const std::vector <double>& dus1) {
            if ( dus1.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for dus1");
            this->items[DUS1]=dus1;
        }

        ///@brief Set item dus2 for the profile (vector size must equal nlevels)
        void Profile::setDus2 (const std::vector <double>& dus2) {
            if ( dus2.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for dus2");
            this->items[DUS2]=dus2;
        }

        ///@brief Set item dus3 for the profile (vector size must equal nlevels)
        void Profile::setDus3 (const std::vector <double>& dus3) {
            if ( dus3.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for dus3");
            this->items[DUS3]=dus3;
        }

        ///@brief Set item sulp for the profile (vector size must equal nlevels)
        void Profile::setSulp (const std::vector <double>& sulp) {
            if ( sulp.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for sulp");
            this->items[SULP]=sulp;
        }

        ///@brief Set item ssa1 for the profile (vector size must equal nlevels)
        void Profile::setSsa1 (const std::vector <double>& ssa1) {
            if ( ssa1.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for ssa1");
            this->items[SSA1]=ssa1;
        }

        ///@brief Set item ssa2 for the profile (vector size must equal nlevels)
        void Profile::setSsa2 (const std::vector <double>& ssa2) {
            if ( ssa2.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for ssa2");
            this->items[SSA2]=ssa2;
        }

        ///@brief Set item ssa3 for the profile (vector size must equal nlevels)
        void Profile::setSsa3 (const std::vector <double>& ssa3) {
            if ( ssa3.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for ssa3");
            this->items[SSA3]=ssa3;
        }

        ///@brief Set item omat for the profile (vector size must equal nlevels)
        void Profile::setOmat (const std::vector <double>& omat) {
            if ( omat.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for omat");
            this->items[OMAT]=omat;
        }

        ///@brief Set item user aerosol n for the profile (vector size must equal nlevels)
        void Profile::setUserAerN (const std::vector <double>& aer, const int n) {
            if ( aer.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for aer");
            if ( n < 0 || n > MAXUSERAER ) throw std::logic_error("Profile error: user aerosol index outside bounds");
            this->items[(itemIdType) (USERAER1+n-1)]=aer;
        }


        ///@brief Return the number of levels
        int Profile::getNlevels() const {
            return nlevels;
        }

        ///@brief Return the gas_units
        int Profile::getGasUnits() const {
            return gas_units;
        }

        ///@brief Return the mmr_cldaer flag
        bool Profile::isMmrCldAer() const {
            return mmr_cldaer;
        }

        ///@brief Return true if profile will use the coefficient pressure levels
        ///       (i.e. the RTTOV interpolator is not being used)
        bool Profile::isDefaultPressureLevels() const {
            return (this->P.size() == 0);
        }

        ///@brief Return a reference to a vector containing the pressures
        const std::vector<double>& Profile::getP() const {
            return P;
        }

        ///@brief Return a reference to a vector containing the temperatures
        const std::vector<double>& Profile::getT() const {
            return T;
        }

        ///@brief Return the angles vector :
        /// satzen, satazi, sunzen, sunazi
        const std::vector<double>& Profile::getAngles() const {
            return angles;
        }

        ///@brief Return the datetimes vector :
        ///  yy, mm, dd, hh, mm, ss
        const std::vector<int>& Profile::getDateTimes() const {
            return datetimes;
        }

        ///@brief Return the clwscheme vector :
        ///  clw_scheme, clwde_param
        const std::vector<int>& Profile::getClwScheme() const {
            return clwscheme;
        }

        ///@brief Return the icecloud vector :
        ///  ice_scheme, icede_param
        const std::vector<int>& Profile::getIceCloud() const {
            return icecloud;
        }

        ///@brief Return the simplecloud vector :
        /// ctp, cfraction
        const std::vector<double>& Profile::getSimpleCloud() const {
            return simplecloud;
        }

        ///@brief Return the surfgeom vector :
        /// lat, lon, elev
        const std::vector<double>& Profile::getSurfGeom() const {
            return surfgeom;
        }

        ///@brief Return the surftype vector :
        /// surftype, watertype
        const std::vector<int>& Profile::getSurfType() const {
            return surftype;
        }

        ///@brief Return the zeeman vector :
        ///  be, cosbk
        const std::vector<double>& Profile::getZeeman() const {
            return zeeman;
        }

        ///@brief Return the s2m vector :
        ///  2m p, 2m t, 2m q, 10m wind u, v, wind fetch
        const std::vector<double>& Profile::getS2m() const {
            return s2m;
        }

        ///@brief Return the skin vector :
        ///  Tskin, salinity, snow_fraction, foam_fraction, fastem_coefs(1:5)
        const std::vector<double>& Profile::getSkin() const {
            return skin;
        }

        ///@brief Return a vector containing item contents
        const std::vector<itemIdType> Profile::getItemContents() const {
            std::vector <itemIdType> contents;
            for (int i=0; i< itemIds.size(); i++ ) {
                if ( this->items.count(itemIds[i]) == 1  ) contents.push_back(itemIds[i]);
            }
            return contents;
        }

        ///@brief Return items : the internal ItemIdPointerMap where gas, cloud and aerosol
        /// vertical profiles are stored (not for general use)
        const ItemIdPointerMap& Profile::getItems() const {
            return items;
        }

        ///@brief Check the profile
        /// @details Checks if all mandatory fields have been provided,
        /// but does not perform a check upon the values (this is done within RTTOV itself);
        /// if simplecoud, clwscheme, icecloud or zeeman have not been set initialise them with default values
        bool Profile::check() {
            bool ok=true;
            if ( this->P.size() > 0 ) {
                if (this->P.size() != this->nlevels ) {
                    ok=false;
                    if (this->verbose ) std::cerr << "Profile error: wrong size for P"<< std::endl;
                }
            }
            if  (this->T.size() != this->nlevels)  {
                ok=false;
                if (this->verbose )std::cerr << "Profile error: wrong size for T" << std::endl;
            }
            if  (this->items[Q].size() != this->nlevels){
                ok=false;
                if (this->verbose )std::cerr << "Profile error: wrong size for gas/cloud/aerosol items" << std::endl;
            }
            if   (this->datetimes.size() != 6 ) {
                ok=false;
                if (this->verbose )std::cerr << "Profile error: wrong size for datetimes" << std::endl;
            }
            if  (this->angles.size() != 4 ) {
                ok=false;
                if (this->verbose )std::cerr <<"Profile error: wrong size for angles" << std::endl;
            }
            if     (this->surftype.size() != 2 ){
                ok=false;
                if (this->verbose )std::cerr << "Profile error: wrong size for surftype" << std::endl;
            }
            if   (this->surfgeom.size() != 3 ){
                ok=false;
                if (this->verbose )std::cerr << "Profile error: wrong size for surfgeom" << std::endl;
            }
            if   (this->s2m.size() != 6 ){
                ok=false;
                if (this->verbose )std::cerr <<"Profile error: wrong size for s2m" << std::endl;
            }
            if   (this->skin.size() != 9 ) {
                ok=false;
                if (this->verbose )std::cerr <<"Profile error: wrong size for skin" << std::endl;
            }
            // assume some fields maybe missing .. initialize them to zeros in this case
            if   (this->simplecloud.size() == 0 ) {
//                 std::cerr << "allocate missing fields for simplecloud " << std::endl;
                this->simplecloud.push_back(0.);
                this->simplecloud.push_back(0.);
            }
            else if   (this->simplecloud.size() != 2 ){
                ok=false;
                if (this->verbose )std::cerr << "Profile error: wrong size for simplecloud" << std::endl;
            }
            if (this->clwscheme.size() == 0 ){
//                 std::cerr << "allocate missing fields for clwscheme " << std::endl;
                this->clwscheme.push_back(1);
                this->clwscheme.push_back(1);
            }
            else if  (this->clwscheme.size()!= 2 ) {
                ok=false;
                if (this->verbose )std::cerr << "Profile error: wrong size for clwscheme" << std::endl;
            }
            if (this->icecloud.size() == 0 ){
//                 std::cerr << "allocate missing fields for icecloud " << std::endl;
                this->icecloud.push_back(1);
                this->icecloud.push_back(1);
            }
            else if  (this->icecloud.size()!= 2 ) {
                ok=false;
                if (this->verbose )std::cerr << "Profile error: wrong size for icecloud" << std::endl;
            }
            if (this->zeeman.size() == 0. ) {
//                 std::cerr << "allocate missing fields zeeman " << std::endl;
                this->zeeman.push_back(0.);
                this->zeeman.push_back(0.);
            }
            else if (this->zeeman.size() != 2 ) {
                ok=false;
                if (this->verbose )std::cerr << "Profile error: wrong size for zeeman" << std::endl;
            }

            return ok;
        }
} /* namespace rttov */
