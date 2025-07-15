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

///@file ProfileScatt.cpp
///
///@brief Method definition for ProfileScatt class

#include <ProfileScatt.h>

namespace rttov {

///@brief Constructor method
ProfileScatt::ProfileScatt(int nlevels):gas_units(unknown),verbose(true),angles{},datetimes{},s2m{},
        skin{},surfgeom{},surftype(0),zeeman{},usercfrac(0.),
        P{},Ph{},T{},items{} {
            // constructor
            this->nlevels=nlevels;
        }

        ProfileScatt::~ProfileScatt() {

        }

        ///@brief Set the gas_units
        ///@param gasUnits unknown=-1,ppmv_dry=0,kg_per_kg=1,ppmv_wet=2
        void ProfileScatt::setGasUnits(rttov::gasUnitType gasUnits) {
            gas_units = gasUnits;
        }

        ///@brief Set the p (pressure) vector
        ///@param p reference to a vector containing the pressures (size nlevels)
        void ProfileScatt::setP(const std::vector<double>& p) {
            if ( p.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for p");
            P = p;
        }

        ///@brief Set the ph (pressure half-levels) vector
        ///@param ph reference to a vector containing the pressure half-levels (size nlevels+1)
        void ProfileScatt::setPh(const std::vector<double>& ph) {
            if ( ph.size() != this->nlevels+1 ) throw std::logic_error("Profile error: wrong size for ph");
            Ph = ph;
        }

        ///@brief Set the temperatures vector
        ///@param t reference to a vector containing the temperatures (size nlevels)
        void ProfileScatt::setT(const std::vector<double>& t) {
            if ( t.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for t");
            T = t;
        }

        ///@brief Set satellite angles
        ///@param satzen, satazi
        void ProfileScatt::setAngles(const double satzen, const double satazi) {
            if ( ! angles.empty()) angles.erase(angles.begin(),angles.end());
            angles.push_back(satzen);
            angles.push_back(satazi);
        }

        ///@brief Set date and time
        ///@param yy, mm, dd, hh, mn, ss
        void ProfileScatt::setDateTimes(const int yy, const int mm, const int dd,
                const int hh, const int mn, const int ss) {
            if ( ! datetimes.empty()) datetimes.erase(datetimes.begin(),datetimes.end());
            datetimes.push_back(yy);
            datetimes.push_back(mm);
            datetimes.push_back(dd);
            datetimes.push_back(hh);
            datetimes.push_back(mn);
            datetimes.push_back(ss);
        }

        ///@brief Set surface geometry parameters
        ///@param lat, lon, elevation
        void ProfileScatt::setSurfGeom(const double lat, const double lon,
                const double elevation) {
            if ( ! surfgeom.empty()) surfgeom.erase(surfgeom.begin(),surfgeom.end());
            surfgeom.push_back(lat);
            surfgeom.push_back(lon);
            surfgeom.push_back(elevation);
        }

        ///@brief Set surface type
        ///@param surftype_in
        void ProfileScatt::setSurfType(const int surftype_in) {
            surftype = surftype_in;
        }

        ///@brief Set zeeman parameters
        ///@param Be, cosbk
        void ProfileScatt::setZeeman(const double Be, const double cosbk) {
            if ( ! zeeman.empty()) zeeman.erase(zeeman.begin(),zeeman.end());
            zeeman.push_back(Be);
            zeeman.push_back(cosbk);
        }

        ///@brief Set surface 2m and 10m parameters
        ///@param p_2m, t_2m, q_2m, u_10m, v_10m
        void ProfileScatt::setS2m(const double p_2m, const double t_2m, const double q_2m,
                const double u_10m, const double v_10m) {
            if (! s2m.empty()) s2m.erase(s2m.begin(),s2m.end());
            s2m.push_back(p_2m);
            s2m.push_back(t_2m);
            s2m.push_back(q_2m);
            s2m.push_back(u_10m);
            s2m.push_back(v_10m);
        }

        ///@brief Set skin parameters
        ///@param t, salinity, foam_fraction, fastem_coef_1
        ///@param fastem_coef_2, fastem_coef_3, fastem_coef_4, fastem_coef_5
        void ProfileScatt::setSkin(const double t, const double salinity,
                const double foam_fraction, const double fastem_coef_1,
                const double fastem_coef_2, const double fastem_coef_3,
                const double fastem_coef_4, const double fastem_coef_5){
            if (! skin.empty()) skin.erase(skin.begin(),skin.end());
            skin.push_back(t);
            skin.push_back(salinity);
            skin.push_back(foam_fraction);
            skin.push_back(fastem_coef_1);
            skin.push_back(fastem_coef_2);
            skin.push_back(fastem_coef_3);
            skin.push_back(fastem_coef_4);
            skin.push_back(fastem_coef_5);
        }

        ///@brief Set item q for the profile (vector size must equal nlevels)
        void ProfileScatt::setQ (const std::vector <double>& q) {
            if ( q.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for q");
            this->items[Q]=q;
        }

        ///@brief Set item o3 for the profile (vector size must equal nlevels)
        void ProfileScatt::setO3 (const std::vector <double>& o3) {
            if ( o3.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for o3");
            this->items[O3]=o3;
        }

        ///@brief Set item hydro_frac for the profile (vector size must equal nlevels)
        void ProfileScatt::setHydroFrac (const std::vector <double>& hydro_frac) {
            if ( hydro_frac.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for hydro_frac");
            this->items[SCATT_HYDRO_FRAC]=hydro_frac;
        }

        ///@brief Set item clw for the profile (vector size must equal nlevels)
        void ProfileScatt::setClw (const std::vector <double>& clw) {
            if ( clw.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for clw");
            this->items[SCATT_CLW]=clw;
        }

        ///@brief Set item ciw for the profile (vector size must equal nlevels)
        void ProfileScatt::setCiw (const std::vector <double>& ciw) {
            if ( ciw.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for ciw");
            this->items[SCATT_CIW]=ciw;
        }

        ///@brief Set item snow for the profile (vector size must equal nlevels)
        void ProfileScatt::setSnow (const std::vector <double>& snow) {
            if ( snow.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for snow");
            this->items[SCATT_SNOW]=snow;
        }

        ///@brief Set item rain for the profile (vector size must equal nlevels)
        void ProfileScatt::setRain (const std::vector <double>& rain) {
            if ( rain.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for rain");
            this->items[SCATT_RAIN]=rain;
        }

        ///@brief Set item graupel for the profile (vector size must equal nlevels)
        void ProfileScatt::setGraupel (const std::vector <double>& graupel) {
            if ( graupel.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for graupel");
            this->items[SCATT_GRAUPEL]=graupel;
        }

        ///@brief Set user cfrac for the profile
        void ProfileScatt::setUserCfrac (const double usercfrac_in) {
            usercfrac = usercfrac_in;
        }

        ///@brief Set item hydrometeor n for the profile (vector size must equal nlevels)
        void ProfileScatt::setHydroN (const std::vector <double>& hydro, const int n) {
            if ( hydro.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for hydro");
            if ( n < 0 || n > MAXSCATTHYDRO ) throw std::logic_error("Profile error: hydro index outside bounds");
            this->items[(itemIdType) (SCATTHYDRO1+n-1)]=hydro;
        }

        ///@brief Set item hydrometeor fraction n for the profile (vector size must equal nlevels)
        void ProfileScatt::setHydroFracN (const std::vector <double>& hydro_frac, const int n) {
            if ( hydro_frac.size() != this->nlevels ) throw std::logic_error("Profile error: wrong size for hydro_frac");
            if ( n < 0 || n > MAXSCATTHYDRO ) throw std::logic_error("Profile error: hydro_frac index outside bounds");
            this->items[(itemIdType) (SCATTHYDROFRAC1+n-1)]=hydro_frac;
        }

        ///@brief Return the number of levels
        int ProfileScatt::getNlevels() const {
            return nlevels;
        }

        ///@brief Return the gas_units
        int ProfileScatt::getGasUnits() const {
            return gas_units;
        }

        ///@brief Return a reference to a vector containing the pressures
        const std::vector<double>& ProfileScatt::getP() const {
            return P;
        }

        ///@brief Return a reference to a vector containing the pressures
        const std::vector<double>& ProfileScatt::getPh() const {
            return Ph;
        }

        ///@brief Return a reference to a vector containing the temperatures
        const std::vector<double>& ProfileScatt::getT() const {
            return T;
        }

        ///@brief Return the user cfrac
        double ProfileScatt::getUserCfrac() const {
            return usercfrac;
        }

        ///@brief Return the angles vector :
        /// satzen, satazi
        const std::vector<double>& ProfileScatt::getAngles() const {
            return angles;
        }

        ///@brief Return the datetimes vector :
        ///  yy, mm, dd, hh, mm, ss
        const std::vector<int>& ProfileScatt::getDateTimes() const {
            return datetimes;
        }

        ///@brief Return the surfgeom vector :
        /// lat, lon, elev
        const std::vector<double>& ProfileScatt::getSurfGeom() const {
            return surfgeom;
        }

        ///@brief Return the surftype
        int ProfileScatt::getSurfType() const {
            return surftype;
        }

        ///@brief Return the zeeman vector :
        ///  be, cosbk
        const std::vector<double>& ProfileScatt::getZeeman() const {
            return zeeman;
        }

        ///@brief Return the s2m vector :
        ///  2m p, 2m t, 2m q, 10m wind u, v
        const std::vector<double>& ProfileScatt::getS2m() const {
            return s2m;
        }

        ///@brief Return the skin vector :
        ///  Tskin, salinity, foam_fraction, fastem_coefs(1:5)
        const std::vector<double>& ProfileScatt::getSkin() const {
            return skin;
        }

        ///@brief Return a vector containing item contents
        const std::vector<itemIdType> ProfileScatt::getItemContents() const {
            std::vector <itemIdType> contents;
            for (int i = 0; i < itemIdsScatt.size(); i++) {
                if (this->items.count(itemIdsScatt[i]) == 1) contents.push_back(itemIdsScatt[i]);
            }
            return contents;
        }

        ///@brief Return items : the internal ItemIdPointerMap where gas and hydrometeor
        /// vertical profiles are stored (not for general use)
        const ItemIdPointerMap& ProfileScatt::getItems() const {
            return items;
        }

        ///@brief Check the profile
        /// @details Checks if all mandatory fields have been provided,
        /// but does not perform a check upon the values (this is done within RTTOV itself);
        /// if zeeman or userCfrac have not been set initialise them with default values
        bool ProfileScatt::check() {
            bool ok=true;
            if  (this->P.size() != this->nlevels) {
                ok=false;
                if (this->verbose) std::cerr << "Profile error: wrong size for P" << std::endl;
            }
            if  (this->Ph.size() != this->nlevels+1) {
                ok=false;
                if (this->verbose) std::cerr << "Profile error: wrong size for Ph" << std::endl;
            }
            if  (this->T.size() != this->nlevels) {
                ok=false;
                if (this->verbose) std::cerr << "Profile error: wrong size for T" << std::endl;
            }
            if  (this->items[Q].size() != this->nlevels) {
                ok=false;
                if (this->verbose) std::cerr << "Profile error: wrong size for gas/hydrometeor items (water vapour is mandatory)" << std::endl;
            }
            if   (this->datetimes.size() != 6) {
                ok=false;
                if (this->verbose) std::cerr << "Profile error: wrong size for datetimes" << std::endl;
            }
            if  (this->angles.size() != 2) {
                ok=false;
                if (this->verbose) std::cerr << "Profile error: wrong size for angles" << std::endl;
            }
            if   (this->surfgeom.size() != 3) {
                ok=false;
                if (this->verbose) std::cerr << "Profile error: wrong size for surfgeom" << std::endl;
            }
            if   (this->s2m.size() != 5) {
                ok=false;
                if (this->verbose) std::cerr << "Profile error: wrong size for s2m" << std::endl;
            }
            if   (this->skin.size() != 8) {
                ok=false;
                if (this->verbose) std::cerr << "Profile error: wrong size for skin" << std::endl;
            }
            // assume some fields maybe missing .. initialize them to zeros in this case
            if (this->zeeman.size() == 0.) {
//                 std::cerr << "allocate missing fields zeeman " << std::endl;
                this->zeeman.push_back(0.);
                this->zeeman.push_back(0.);
            }
            else if (this->zeeman.size() != 2) {
                ok=false;
                if (this->verbose) std::cerr << "Profile error: wrong size for zeeman" << std::endl;
            }
            return ok;
        }
} /* namespace rttov */
