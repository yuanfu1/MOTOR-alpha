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

///@class Options

#include <Options.h>
#include <sstream>
/**
 * @file Options.cpp
 * Options class definitions
 *
 * @brief Options class definitions for managing rttov options.
 *
 * @author Pascale

 * @version
 */

namespace rttov {

///@brief Constructor method
Options::Options(): nprofs_per_call(1),
                    nthreads(1),
                    cloud_overlap(1),
                    cc_low_cloud_top(750.),
                    cldcol_threshold(-1.),
                    ir_scatt_model(chou),
                    vis_scatt_model(dom_vis),
                    rayleigh_max_wavelength(2.),
                    rayleigh_min_pressure(0.),
                    dom_nstreams(8),
                    dom_accuracy(0.),
                    dom_opdep_threshold(0.),
                    solar_sea_brdf_model(elfouhaily),
                    ir_sea_emis_model(iremis),
                    clw_scheme(rosenkranz),
                    clw_cloud_top(322.),
                    fastem_version(fastem6),
                    interp_mode(1),
                    ipcbnd(band1),
                    ipcreg(index1),
                    cc_threshold(0.001),
                    ice_polarisation(1.40){
    // General configuration options
    this->config["apply_reg_limits"]=false;
    this->config["verbose"]=true;
    this->config["do_checkinput"]=true;
    this->config["fix_hgpl"]=true;
    // Interpolation options
    this->interp["addinterp"]=true ; // to prevent crashing if user forget to set it if profile number of level different from coeffs
    this->interp["spacetop"]=true;
    this->interp["lgradp"]=false;
    this->interp["reg_limit_extrap"]=true;
    // General RT options
    this->rt_all["ozone_data"]=false;
    this->rt_all["co2_data"]=false;
    this->rt_all["n2o_data"]=false;
    this->rt_all["co_data"]=false;
    this->rt_all["ch4_data"]=false;
    this->rt_all["so2_data"]=false;
    this->rt_all["switchrad"]=false;
    this->rt_all["use_q2m"]=true;
    this->rt_all["use_t2m_opdep"]=true;
    this->rt_all["addrefrac"]=true;
    this->rt_all["plane_parallel"]=false;
    this->rt_all["do_lambertian"]=false;
    this->rt_all["lambertian_fixed_angle"]=true;
    this->rt_all["rad_down_lin_tau"]=true;
    this->rt_all["dtau_test"]=false;
    // MW-only RT options
    this->mw["clw_data"]=false;
    this->mw["supply_foam_fraction"]=false;
    // VIS/IR-only RT options
    this->ir["addsolar"]=false;
    this->ir["rayleigh_single_scatt"]=true;
    this->ir["dom_rayleigh"]=false;
    this->ir["do_nlte_correction"]=false;
    this->ir["addaerosl"]=false;
    this->ir["user_aer_opt_param"]=false;
    this->ir["addclouds"]=false;
    this->ir["user_cld_opt_param"]=false;
    this->ir["grid_box_avg_cloud"]=true;
    // PC-RTTOV options
    this->pc["addpc"]=false;
//     this->pc["addradrec"]=false;
    // RTTOV-SCATT options
    this->scatt["lusercfrac"]=false;
    this->scatt["hydro_cfrac_tlad"]=true;
    this->scatt["zero_hydro_tlad"]=false;
    // Developer options
    this->dev["do_opdep_calc"]=true;
    // Wrapper options
    this->verbose_wrapper=false;
    this->check_opts=true;
    this->store_rad=false;
    this->store_rad2=false;
    this->store_trans=false;
    this->store_emis_terms=false;
}

///@brief Destructor method
Options::~Options() {

}

// RTTOV opts%config%

///@brief Return the opts\%config\%apply_reg_limits and opts_scatt\%config\%apply_reg_limits options
bool Options::isApplyRegLimits() {
    return this->config["apply_reg_limits"] ;
}

///@brief Set the opts\%config\%apply_reg_limits and opts_scatt\%config\%apply_reg_limits options
void Options::setApplyRegLimits(bool applyRegLimits) {
    this->config["apply_reg_limits"] = applyRegLimits;
}

///@brief Return the opts\%config\%do_checkinput and opts_scatt\%config\%do_checkinput options
bool Options::isDoCheckinput() {
    return this->config["do_checkinput"];
}

///@brief Set the opts\%config\%do_checkinput and opts_scatt\%config\%do_checkinput options
void Options::setDoCheckinput(bool doCheckinput) {
    this->config["do_checkinput"] = doCheckinput;
}

///@brief Return the opts\%config\%verbose and opts_scatt\%config\%verbose options
bool Options::isVerbose() {
    return this->config["verbose"];
}

///@brief Set the opts\%config\%verbose and opts_scatt\%config\%verbose options
void Options::setVerbose(bool verbose) {
    this->config["verbose"] = verbose;
}

///@brief Return the opts\%config\%fix_hgpl and opts_scatt\%config\%fix_hgpl options
bool Options::isFixHgpl() {
    return this->config["fix_hgpl"];
}

///@brief Set the opts\%config\%verbose and opts_scatt\%config\%verbose options
void Options::setFixHgpl(bool fixHgpl) {
    this->config["fix_hgpl"] = fixHgpl;
}


// RTTOV opts%interpolation%

///@brief Return the opts\%interpolation\%addinterp option
bool Options::isAddInterp() {
    return this->interp["addinterp"];
}

///@brief Set the opts\%interpolation\%addinterp option
void Options::setAddInterp(bool addinterp) {
    this->interp["addinterp"]= addinterp;
}

///@brief Return the opts\%interpolation\%interp_mode and opts_scatt\%interp_mode options
int Options::getInterpMode() const {
    return interp_mode;
}

///@brief Set the opts\%interpolation\%interp_mode and opts_scatt\%interp_mode options
void Options::setInterpMode(int interpMode) {
    interp_mode = interpMode;
}

///@brief Return the opts\%interpolation\%reg_limit_extrap and opts_scatt\%reg_limit_extrap options
bool Options::isRegLimitExtrap() {
    return this->interp["reg_limit_extrap"];
}

///@brief Set the opts\%interpolation\%reg_limit_extrap and opts_scatt\%reg_limit_extrap options
void Options::setRegLimitExtrap(bool regLimitExtrap) {
    this->interp["reg_limit_extrap"] = regLimitExtrap;
}

///@brief Return the opts\%interpolation\%spacetop option
bool Options::isSpacetop() {
    return this->interp["spacetop"];
}

///@brief Set the opts\%interpolation\%spacetop option
void Options::setSpacetop(bool spacetop) {
    this->interp["spacetop"] = spacetop;
}

///@brief Return the opts\%interpolation\%lgradp and opts_scatt\%lgradp options
bool Options::isLgradp() {
    return this->interp["lgradp"];
}

///@brief Set the opts\%interpolation\%lgradp and opts_scatt\%lgradp options
void Options::setLgradp(bool lgradp) {
    this->interp["lgradp"]= lgradp;
}


// RTTOV opts%rt_all%

///@brief Return the opts\%rt_all\%ozone_data option
bool Options::isOzoneData() {
    return this->rt_all["ozone_data"];
}

///@brief Set  the opts\%rt_all\%ozone_data option
void Options::setOzoneData(bool ozoneData) {
    this->rt_all["ozone_data"] = ozoneData;
}

///@brief Return the opts\%rt_all\%co2_data option
bool Options::isCO2Data() {
    return this->rt_all["co2_data"];
}

///@brief Set the opts\%rt_all\%co2_data option
void Options::setCO2Data(bool co2Data) {
    this->rt_all["co2_data"] = co2Data;
}

///@brief Return the opts\%rt_all\%ch4_data option
bool Options::isCH4Data() {
    return this->rt_all["ch4_data"];
}

///@brief Set the opts\%rt_all\%ch4_data option
void Options::setCH4Data(bool ch4Data) {
    this->rt_all["ch4_data"] = ch4Data;
}

///@brief Return the opts\%rt_all\%co_data option
bool Options::isCOData() {
    return this->rt_all["co_data"];
}

///@brief Set the opts\%rt_all\%co_data option
void Options::setCOData(bool coData) {
    this->rt_all["co_data"] = coData;
}

///@brief Return the opts\%rt_all\%n2o_data option
bool Options::isN2OData() {
    return this->rt_all["n2o_data"];
}

///@brief Set the opts\%rt_all\%n2o_data option
void Options::setN2OData(bool n2oData) {
    this->rt_all["n2o_data"] = n2oData;
}

///@brief Return the opts\%rt_all\%so2_data option
bool Options::isSO2Data() {
    return this->rt_all["so2_data"];
}

///@brief Set the opts\%rt_all\%so2_data option
void Options::setSO2Data(bool so2Data) {
    this->rt_all["so2_data"] = so2Data;
}

///@brief Return the opts\%rt_all\%switchrad option
bool Options::isSwitchrad() {
    return this->rt_all["switchrad"];
}

///@brief Set the opts\%rt_all\%switchrad option
void Options::setSwitchrad(bool switchrad) {
    this->rt_all["switchrad"] = switchrad;
}

///@brief Return the opts\%rt_all\%use_t2m_opdep and opts_scatt\%use_t2m_opdep options
bool Options::isUseT2mOpdep() {
    return this->rt_all["use_t2m_opdep"];
}

///@brief Set the opts\%rt_all\%use_t2m_opdep and opts_scatt\%use_t2m_opdep options
void Options::setUseT2mOpdep(bool useT2mOpdep) {
    this->rt_all["use_t2m_opdep"] = useT2mOpdep;
}

///@brief Return the opts\%rt_all\%use_q2m and opts_scatt\%use_q2m options
bool Options::isUseQ2m() {
    return this->rt_all["use_q2m"];
}

///@brief Set the opts\%rt_all\%use_q2m and opts_scatt\%use_q2m options
void Options::setUseQ2m(bool useQ2m) {
    this->rt_all["use_q2m"] = useQ2m;
}

///@brief Return the opts\%rt_all\%addrefrac option
bool Options::isAddRefrac() {
    return this->rt_all["addrefrac"];
}

///@brief Set the opts\%rt_all\%addrefrac option
void Options::setAddRefrac(bool addRefrac) {
    this->rt_all["addrefrac"] = addRefrac;
}

///@brief Return the opts\%rt_all\%plane_parallel option
bool Options::isPlaneParallel() {
    return this->rt_all["plane_parallel"];
}

///@brief Set the opts\%rt_all\%plane_parallel option
void Options::setPlaneParallel(bool planeParallel) {
    this->rt_all["plane_parallel"] = planeParallel;
}

///@brief Return the opts\%rt_all\%do_lambertian option
bool Options::isDoLambertian() {
    return this->rt_all["do_lambertian"];
}

///@brief Set the opts\%rt_all\%do_lambertian option
void Options::setDoLambertian(bool doLambertian) {
    this->rt_all["do_lambertian"] = doLambertian;
}

///@brief Return the opts\%rt_all\%lambertian_fixed_angle option
bool Options::isLambertianFixedAngle() {
    return this->rt_all["lambertian_fixed_angle"];
}

///@brief Set the opts\%rt_all\%lambertian_fixed_angle option
void Options::setLambertianFixedAngle(bool lambertianFixedAngle) {
    this->rt_all["lambertian_fixed_angle"] = lambertianFixedAngle;
}

///@brief Return the opts\%rt_all\%rad_down_lin_tau and opts_scatt\%rad_down_lin_tau options
bool Options::isRadDownLinTau() {
    return this->rt_all["rad_down_lin_tau"];
}

///@brief Set the opts\%rt_all\%rad_down_lin_tau and opts_scatt\%rad_down_lin_tau options
void Options::setRadDownLinTau(bool radDownLinTau) {
    this->rt_all["rad_down_lin_tau"] = radDownLinTau;
}

///@brief Return the opts\%rt_all\%dtau_test and opts_scatt\%dtau_test options
bool Options::isDtauTest() {
    return this->rt_all["dtau_test"];
}

///@brief Set the opts\%rt_all\%dtau_test and opts_scatt\%dtau_test options
void Options::setDtauTest(bool dtauTest) {
    this->rt_all["dtau_test"] = dtauTest;
}


// RTTOV opts%rt_mw%

///@brief Return the opts\%rt_mw\%supply_foam_fraction and opts_scatt\%supply_foam_fraction options
bool Options::isSupplyFoamFraction() {
    return this->mw["supply_foam_fraction"];
}

///@brief Set the opts\%rt_mw\%supply_foam_fraction and opts_scatt\%supply_foam_fraction options
void Options::setSupplyFoamFraction(bool supplyFoamFraction) {
    this->mw["supply_foam_fraction"] = supplyFoamFraction;
}

///@brief Return the opts\%rt_mw\%clw_data option
bool Options::isCLWData() {
    return this->mw["clw_data"];
}

///@brief Set the opts\%rt_mw\%clw_data option
void Options::setCLWData(bool clwData) {
    this->mw["clw_data"] = clwData;
}

///@brief Return the opts\%rt_mw\%clw_scheme option
int Options::getCLWScheme() const {
    return clw_scheme;
}

///@brief Set the opts\%rt_mw\%clw_scheme option
void Options::setCLWScheme(int clwScheme) {
    clw_scheme = clwScheme;
}

///@brief Return the opts\%rt_mw\%clw_cloud_top option
double Options::getCLWCloudTop() const {
    return clw_cloud_top;
}

///@brief Set the opts\%rt_mw\%clw_cloud_top option
void Options::setCLWCloudTop(double clwCloudTop) {
    clw_cloud_top = clwCloudTop;
}

///@brief Return the opts\%rt_mw\%fastem_version and opts_scatt\%fastem_version options
int Options::getFastemVersion() const {
    return fastem_version;
}

///@brief Set the opts\%rt_mw\%fastem_version and opts_scatt\%fastem_version options
void Options::setFastemVersion(int fastemVersion) {
    fastem_version = fastemVersion;
}


// RTTOV opts%rt_ir%

///@brief Return the opts\%rt_ir\%solar_sea_brdf_model option
int Options::getSolarSeaBrdfModel() const {
    return solar_sea_brdf_model;
}

///@brief Set the opts\%rt_ir\%solar_sea_brdf_model option
void Options::setSolarSeaBrdfModel(int solarSeaBrdfModel) {
    solar_sea_brdf_model = solarSeaBrdfModel;
}

///@brief Return the opts\%rt_ir\%ir_sea_emis_model option
int Options::getIrSeaEmisModel() const {
    return ir_sea_emis_model;
}

///@brief Set the opts\%rt_ir\%ir_sea_emis_model option
void Options::setIrSeaEmisModel(int irSeaEmisModel) {
    ir_sea_emis_model = irSeaEmisModel;
}

///@brief Return the opts\%rt_ir\%addsolar option
bool Options::isAddSolar() {
    return this->ir["addsolar"];
}

///@brief Set the opts\%rt_ir\%addsolar option
void Options::setAddSolar(bool addsolar) {
    this->ir["addsolar"] = addsolar;
}

///@brief Return the opts\%rt_ir\%rayleigh_max_wavelength option
double Options::getRayleighMaxWavelength() const {
    return rayleigh_max_wavelength;
}

///@brief Set the opts\%rt_ir\%rayleigh_max_wavelength option
void Options::setRayleighMaxWavelength(double rayleighMaxWavelength) {
    rayleigh_max_wavelength = rayleighMaxWavelength;
}

///@brief Return the opts\%rt_ir\%rayleigh_min_pressure option
double Options::getRayleighMinPressure() const {
    return rayleigh_min_pressure;
}

///@brief Set the opts\%rt_ir\%rayleigh_min_pressure option
void Options::setRayleighMinPressure(double rayleighMinPressure) {
    rayleigh_min_pressure = rayleighMinPressure;
}

///@brief Return the opts\%rt_ir\%rayleigh_single_scatt option
bool Options::isRayleighSingleScatt() {
    return this->ir["rayleigh_single_scatt"];
}

///@brief Set the opts\%rt_ir\%rayleigh_single_scatt option
void Options::setRayleighSingleScatt(bool rayleighSingleScatt) {
    this->ir["rayleigh_single_scatt"] = rayleighSingleScatt;
}

///@brief Return the opts\%rt_ir\%do_nlte_correction option
bool Options::isDoNlteCorrection() {
    return this->ir["do_nlte_correction"];
}

///@brief Set the opts\%rt_ir\%do_nlte_correction option
void Options::setDoNlteCorrection(bool doNlteCorrection) {
    this->ir["do_nlte_correction"] = doNlteCorrection;
}

///@brief Return the opts\%rt_ir\%addaerosl option
bool Options::isAddAerosl() {
    return this->ir["addaerosl"];
}

///@brief Set the opts\%rt_ir\%addaerosl option
void Options::setAddAerosl(bool addaerosl) {
    this->ir["addaerosl"] = addaerosl;
}

///@brief Return the opts\%rt_ir\%user_aer_opt_param option
bool Options::isUserAerOptParam() {
    return this->ir["user_aer_opt_param"];
}

///@brief Set the opts\%rt_ir\%user_aer_opt_param option
void Options::setUserAerOptParam(bool userAerOptParam) {
    this->ir["user_aer_opt_param"] = userAerOptParam;
}

///@brief Return the opts\%rt_ir\%addclouds option
bool Options::isAddClouds() {
    return this->ir["addclouds"];
}

///@brief Set the opts\%rt_ir\%addclouds option
void Options::setAddClouds(bool addclouds) {
    this->ir["addclouds"]= addclouds;
}

///@brief Return the opts\%rt_ir\%user_cld_opt_param option
bool Options::isUserCldOptParam() {
    return this->ir["user_cld_opt_param"];
}

///@brief Set the opts\%rt_ir\%user_cld_opt_param option
void Options::setUserCldOptParam(bool userCldOptParam) {
    this->ir["user_cld_opt_param"] = userCldOptParam;
}

///@brief Return the opts\%rt_ir\%grid_box_avg_cloud option
bool Options::isGridBoxAvgCloud() {
    return this->ir["grid_box_avg_cloud"];
}

///@brief Set the opts\%rt_ir\%grid_box_avg_cloud option
void Options::setGridBoxAvgCloud(bool gridBoxAvgCloud) {
    this->ir["grid_box_avg_cloud"] = gridBoxAvgCloud;
}

///@brief Return the opts\%rt_ir\%cldcol_threshold option
double Options::getCldcolThreshold() const {
    return cldcol_threshold;
}

///@brief Set the opts\%rt_ir\%cldcol_threshold option
void Options::setCldcolThreshold(double cldcolThreshold) {
    cldcol_threshold = cldcolThreshold;
}

///@brief Return the opts\%rt_ir\%cloud_overlap option
int Options::getCloudOverlap() const {
    return cloud_overlap;
}

///@brief Set the opts\%rt_ir\%cloud_overlap option
void Options::setCloudOverlap(int cloudOverlap) {
    cloud_overlap = cloudOverlap;
}

///@brief Return the opts\%rt_ir\%cc_low_cloud_top option
double Options::getCCLowCloudTop() const {
    return cc_low_cloud_top;
}

///@brief Set the opts\%rt_ir\%cc_low_cloud_top option
void Options::setCCLowCloudTop(double ccLowCloudTop) {
    cc_low_cloud_top = ccLowCloudTop;
}

///@brief Return the opts\%rt_ir\%ir_scatt_model option
int Options::getIrScattModel() const {
    return ir_scatt_model;
}

///@brief Set the opts\%rt_ir\%ir_scatt_model option
void Options::setIrScattModel(int irScattModel) {
    ir_scatt_model = irScattModel;
}

///@brief Return the opts\%rt_ir\%vis_scatt_model option
int Options::getVisScattModel() const {
    return vis_scatt_model;
}

///@brief Set the opts\%rt_ir\%vis_scatt_model option
void Options::setVisScattModel(int visScattModel) {
    vis_scatt_model = visScattModel;
}

///@brief Return the opts\%rt_ir\%dom_nstreams option
int Options::getDomNstreams() const {
    return dom_nstreams;
}

///@brief Set the opts\%rt_ir\%dom_nstreams option
void Options::setDomNstreams(int domNstreams) {
    dom_nstreams = domNstreams;
}

///@brief Return the opts\%rt_ir\%dom_accuracy option
double Options::getDomAccuracy() const {
    return dom_accuracy;
}

///@brief Set the opts\%rt_ir\%dom_accuracy option
void Options::setDomAccuracy(double domAccuracy) {
    dom_accuracy = domAccuracy;
}

///@brief Return the opts\%rt_ir\%dom_opdep_threshold option
double Options::getDomOpdepThreshold() const {
    return dom_opdep_threshold;
}

///@brief Set the opts\%rt_ir\%dom_opdep_threshold option
void Options::setDomOpdepThreshold(double domOpdepThreshold) {
    dom_opdep_threshold = domOpdepThreshold;
}

///@brief Return the opts\%rt_ir\%dom_rayleigh option
bool Options::isDomRayleigh() {
    return this->ir["dom_rayleigh"];
}

///@brief Set the opts\%rt_ir\%dom_rayleigh option
void Options::setDomRayleigh(bool domRayleigh) {
    this->ir["dom_rayleigh"] = domRayleigh;
}


// RTTOV opts%rt_ir%pc%

///@brief Return the opts\%rt_ir\%pc\%addpc option
bool Options::isAddPC() {
    return this->pc["addpc"];
}

/////@brief Set the opts\%rt_ir\%pc\%addpc option
// void Options::setAddPC(bool addpc) {
//     this->pc["addpc"]= addpc;
// }

///@brief Return the opts\%rt_ir\%pc\%addradrec option
bool Options::isAddRadrec() {
    return this->pc["addradrec"];
}

/////@brief Set the opts\%rt_ir\%pc\%addradrec option
// void Options::setAddRadrec(bool addradrec) {
//     this->pc["addradrec"]= addradrec;
// }

///@brief Return the opts\%rt_ir\%pc\%ipcreg option
int Options::getIpcreg() const {
    return ipcreg;
}

/////@brief Set the opts\%rt_ir\%pc\%ipcreg option
// void Options::setIpcreg(int ipcregx) {
//     ipcreg = ipcregx;
// }

///@brief Return the opts\%rt_ir\%pc\%ipcbnd option
int Options::getIpcbnd() const {
    return ipcbnd;
}

/////@brief Set the opts\%rt_ir\%pc\%ipcbnd option
// void Options::setIpcbnd(int ipcbndx) {
//     ipcbnd = ipcbndx;
// }


// RTTOV opts_scatt%

///@brief Return the opts_scatt\%lusercfrac option
bool Options::isLuserCfrac() {
    return this->scatt["lusercfrac"];
}

///@brief Set the opts_scatt\%lusercfrac option
void Options::setLuserCfrac(bool lusercfrac) {
    this->scatt["lusercfrac"]= lusercfrac;
}

///@brief Return the opts_scatt\%cc_threshold option
double Options::getCCThreshold() const {
    return cc_threshold;
}

///@brief Set the opts_scatt\%cc_threshold option
void Options::setCCThreshold(double ccThreshold) {
    cc_threshold = ccThreshold;
}

///@brief Return the opts_scatt\%ice_polarisation option
double Options::getIcePolarisation() const {
    return ice_polarisation;
}

///@brief Set the opts_scatt\%ice_polarisation option
void Options::setIcePolarisation(double icePolarisation) {
    ice_polarisation = icePolarisation;
}

///@brief Return the opts_scatt\%hydro_cfrac_tlad option
bool Options::isHydroCfracTLAD() {
    return this->scatt["hydro_cfrac_tlad"];
}

///@brief Set the opts_scatt\%hydro_cfrac_tlad option
void Options::setHydroCfracTLAD(bool hydroCfracTLAD) {
    this->scatt["hydro_cfrac_tlad"]= hydroCfracTLAD;
}

///@brief Return the opts_scatt\%zero_hydro_tlad option
bool Options::isZeroHydroTLAD() {
    return this->scatt["zero_hydro_tlad"];
}

///@brief Set the opts_scatt\%zero_hydro_tlad option
void Options::setZeroHydroTLAD(bool zeroHydroTLAD) {
    this->scatt["zero_hydro_tlad"]= zeroHydroTLAD;
}


// RTTOV opts%dev%

///@brief Return the opts\%dev\%do_opdep_calc option
bool Options::isDoOpdepCalc() {
    return this->dev["do_opdep_calc"];
}

///@brief Set the opts\%dev\%do_opdep_calc option
void Options::setDoOpdepCalc(bool doOpdepCalc) {
    this->dev["do_opdep_calc"] = doOpdepCalc;
}


// RTTOV wrapper options

///@brief Return the number of profiles passed into rttov_direct or rttov_k per call
int Options::getNprofsPerCall() const {
    return nprofs_per_call;
}

///@brief Set the number of profiles passed into rttov_direct or rttov_k per call
void Options::setNprofsPerCall(int nprofsPerCall) {
    nprofs_per_call = nprofsPerCall;
}

///@brief Return the number of threads RTTOV will use (compile RTTOV with OpenMP to make use of this)
int Options::getNthreads() const {
    return nthreads;
}

///@brief Set the number of threads RTTOV will use (compile RTTOV with OpenMP to make use of this)
void Options::setNthreads(int nthreads) {
    this->nthreads = nthreads;
}

///@brief Return set the verbose_wrapper option
bool Options::isVerboseWrapper() const {
    return verbose_wrapper;
}

///@brief Set the verbose_wrapper option
void Options::setVerboseWrapper(bool verboseWrapper) {
    verbose_wrapper = verboseWrapper;
}

///@brief Return set the check_opts option
bool Options::isCheckOpts() const {
    return check_opts;
}

///@brief Set the check_opts option
void Options::setCheckOpts(bool checkOpts) {
    check_opts = checkOpts;
}

///@brief Return the store_rad wrapper option
bool Options::isStoreRad() const {
    return store_rad;
}

///@brief Set the store_rad wrapper option
void Options::setStoreRad(bool storeRad) {
    store_rad = storeRad;
}

///@brief Return the store_rad2 wrapper option
bool Options::isStoreRad2() const {
    return store_rad2;
}

///@brief Set the store_rad2 wrapper option
void Options::setStoreRad2(bool storeRad2) {
    store_rad2 = storeRad2;
}

///@brief Return the store_trans wrapper option
bool Options::isStoreTrans() const {
    return store_trans;
}

///@brief Set the store_trans wrapper option
void Options::setStoreTrans(bool storeTrans) {
    store_trans = storeTrans;
}

///@brief Return the store_emis_terms wrapper option
bool Options::isStoreEmisTerms() const {
    return store_emis_terms;
}

///@brief Set the store_emis_terms wrapper option
void Options::setStoreEmisTerms(bool storeEmisTerms) {
    store_emis_terms = storeEmisTerms;
}


///@internal generate the options string needed by the wrapper rttov_set_options subroutine
///does not perform any control upon the options
///return a string option example :
///... opts\%interpolation\%addinterp 1 opts\%rt_ir\%co2_data 1 opts\%rt_ir\%addclouds 1 nprofs_per_call 2 nthreads 1 verbose_wrapper 0 ..."
std::string Options::defineStrOptions() {

    std::string strOptions="";
    std::string key;
    bool val;

    for(StrBoolMap::iterator iter = this->config.begin(); iter != this->config.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->config[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%config%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->interp.begin(); iter != this->interp.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->interp[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%interpolation%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->rt_all.begin(); iter != this->rt_all.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->rt_all[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%rt_all%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->mw.begin(); iter != this->mw.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->mw[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%rt_mw%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->ir.begin(); iter != this->ir.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->ir[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%rt_ir%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->pc.begin(); iter != this->pc.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->pc[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts%rt_ir%pc%");
        strOptions.append(key);
        strOptions.append(strval);
    }
    for(StrBoolMap::iterator iter = this->scatt.begin(); iter != this->scatt.end(); ++iter)
    {
        std::string strval;
        key =  iter->first;
        val=this->scatt[key];
        if (val) strval.append(" 1"); else strval.append(" 0");
        strOptions.append(" opts_scatt%");
        strOptions.append(key);
        strOptions.append(strval);
    }

    std::ostringstream oss;
    if (this->isAddClouds()) {
        oss << " opts%rt_ir%cc_low_cloud_top "<< this->cc_low_cloud_top;
        oss << " opts%rt_ir%cldcol_threshold "<< this->cldcol_threshold;
    }
    if (this->isAddClouds() || this->isAddAerosl() ) {
      oss << " opts%rt_ir%ir_scatt_model "<< this->ir_scatt_model;
      oss << " opts%rt_ir%vis_scatt_model "<< this->vis_scatt_model;
      oss << " opts%rt_ir%dom_nstreams "<< this->dom_nstreams;
      oss << " opts%rt_ir%dom_accuracy "<< this->dom_accuracy;
      oss << " opts%rt_ir%dom_opdep_threshold "<< this->dom_opdep_threshold;
    }
    if (this->isAddPC()) {
        oss << " opts%rt_ir%pc%ipcbnd "<< this->ipcbnd;
        oss << " opts%rt_ir%pc%ipcreg "<< this->ipcreg;
    }
    oss << " opts%interpolation%interp_mode "<< this->interp_mode;
    oss << " opts%rt_mw%fastem_version "<< this->fastem_version;
    oss << " opts%rt_mw%clw_scheme "<< this->clw_scheme;
    oss << " opts%rt_mw%clw_cloud_top "<< this->clw_cloud_top;
    oss << " opts%rt_ir%ir_sea_emis_model "<< this->ir_sea_emis_model;
    oss << " opts%rt_ir%solar_sea_brdf_model "<< this->solar_sea_brdf_model;
    oss << " opts_scatt%cc_threshold "<< this->cc_threshold;
    oss << " opts_scatt%ice_polarisation "<< this->ice_polarisation;
    oss << " nthreads "<< this->nthreads;
    oss << " nprofs_per_call "<< this->nprofs_per_call;
    strOptions.append(oss.str());

    strOptions.append(" verbose_wrapper");
    val=this->verbose_wrapper;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    strOptions.append(" check_opts");
    val=this->check_opts;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    strOptions.append(" store_trans");
    val=this->store_trans;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    strOptions.append(" store_rad");
    val=this->store_rad;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    strOptions.append(" store_rad2");
    val=this->store_rad2;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    strOptions.append(" store_emis_terms");
    val=this->store_emis_terms;
    if (val) strOptions.append(" 1"); else strOptions.append(" 0");

    return strOptions;
}
} /* namespace rttov */
