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
///@brief This class allows the RTTOV and wrapper options to be specified

#ifndef OPTIONS_H_
#define OPTIONS_H_
#include <map>
#include <string>

namespace rttov {

typedef std::map<std::string,bool>StrBoolMap;
enum ipcbndType {band1=1,band2,band3};
enum ipcregType {index1=1,index2,index3,index4};
enum fastemType {tessem2=0,fastem1,fastem2,fastem3,fastem4,fastem5,fastem6};
enum solarseabrdfType{jonswap=1,elfouhaily};
enum irseaemisType{isem=1,iremis};
enum cloudoverlapType{maxrandom_overlap=1,simple_overlap};
enum irscattmodelType{dom_ir=1,chou};
enum visscattmodelType{dom_vis=1,singlescatt,mfasis};
enum mwclwschemeType{liebe=1,rosenkranz,tkc};

class Options {
public:
    Options();
    virtual ~Options();

    // RTTOV opts%config%
    void setApplyRegLimits(bool applyRegLimts);
    void setDoCheckinput(bool doCheckinput);
    void setVerbose(bool verbose);
    void setFixHgpl(bool fixHgpl);

    //RTTOV opts%interpolation%
    void setAddInterp(bool addinterp);
    void setInterpMode(int interpMode);
    void setRegLimitExtrap(bool regLimitExtrap);
    void setSpacetop(bool spacetop);
    void setLgradp(bool lgradp);

    //RTTOV opts%rt_all%
    void setOzoneData(bool ozoneData);
    void setCO2Data(bool co2Data);
    void setCH4Data(bool ch4Data);
    void setCOData(bool coData);
    void setN2OData(bool n2oData);
    void setSO2Data(bool so2Data);

    void setDoLambertian(bool doLambertian);
    void setLambertianFixedAngle(bool lambertianFixedAngle);
    void setUseT2mOpdep(bool useT2mOpdep);
    void setUseQ2m(bool useQ2m);
    void setSwitchrad(bool switchrad);
    void setAddRefrac(bool addRefrac);
    void setPlaneParallel(bool planeParallel);
    void setRadDownLinTau(bool radDownLinTau);
    void setDtauTest(bool dtauTest);

    //RTTOV opts%rt_mw%
    void setCLWData(bool clwData);
    void setCLWScheme(int clwScheme);
    void setCLWCloudTop(double clwCloudTop);
    void setFastemVersion(int fastemVersion);
    void setSupplyFoamFraction(bool supplyFoamFraction);

    //RTTOV opts%rt_ir%
    void setSolarSeaBrdfModel(int solarSeaBrdfModel);
    void setIrSeaEmisModel(int irSeaEmisModel);
    void setAddSolar(bool addsolar);
    void setRayleighMaxWavelength(double rayleighMaxWavelength);
    void setRayleighMinPressure(double rayleighMinPressure);
    void setRayleighSingleScatt(bool rayleighSingleScatt);
    void setDoNlteCorrection(bool doNlteCorrection);

    void setAddAerosl(bool addaerosl);
    void setUserAerOptParam(bool userAerOptParam);
    void setAddClouds(bool addclouds);
    void setUserCldOptParam(bool userCldOptParam);
    void setGridBoxAvgCloud(bool gridBoxAvgCloud);
    void setCloudOverlap(int cloudOverlap);
    void setCCLowCloudTop(double ccLowCloudTop);
    void setCldcolThreshold(double cldcolThreshold);
    void setIrScattModel(int irScattModel);
    void setVisScattModel(int visScattModel);
    void setDomNstreams(int domNstreams);
    void setDomAccuracy(double domAccuracy);
    void setDomOpdepThreshold(double domOpdepThreshold);
    void setDomRayleigh(bool domRayleigh);

//     void setAddPC(bool addpc);
//     void setAddRadrec(bool addradrec);
//     void setIpcreg(int ipcreg);
//     void setIpcbnd(int ipcbnd);

    //RTTOV-SCATT
    void setLuserCfrac(bool lusercfrac);
    void setCCThreshold(double ccThreshold);
    void setIcePolarisation(double icePolarisation);
    void setHydroCfracTLAD(bool hydroCfracTLAD);
    void setZeroHydroTLAD(bool zeroHydroTLAD);

    //RTTOV opts%dev%
    void setDoOpdepCalc(bool doOpdepCalc);

    // Wrapper options
    void setNthreads(int nthreads);
    void setNprofsPerCall(int nprofsPerCall);
    void setVerboseWrapper(bool verboseWrapper);
    void setCheckOpts(bool checkOpts);
    void setStoreRad(bool storeRad);
    void setStoreRad2(bool storeRad2);
    void setStoreTrans(bool storeTrans);
    void setStoreEmisTerms(bool storeEmisTerms);


    // RTTOV opts%config%
    bool isApplyRegLimits();
    bool isDoCheckinput();
    bool isVerbose();
    bool isFixHgpl();

    //RTTOV opts%interpolation%
    bool isAddInterp();
    int getInterpMode() const;
    bool isRegLimitExtrap();
    bool isSpacetop();
    bool isLgradp();

    //RTTOV opts%rt_all%
    bool isOzoneData();
    bool isCO2Data();
    bool isCH4Data();
    bool isCOData();
    bool isN2OData();
    bool isSO2Data();

    bool isDoLambertian();
    bool isLambertianFixedAngle();
    bool isUseT2mOpdep();
    bool isUseQ2m();
    bool isSwitchrad();
    bool isAddRefrac();
    bool isPlaneParallel();
    bool isRadDownLinTau();
    bool isDtauTest();

    //RTTOV opts%rt_mw%
    bool isCLWData();
    int getCLWScheme() const;
    double getCLWCloudTop() const;
    int getFastemVersion() const;
    bool isSupplyFoamFraction();

    //RTTOV opts%rt_ir%
    int getSolarSeaBrdfModel() const;
    int getIrSeaEmisModel() const;
    bool isAddSolar();
    double getRayleighMaxWavelength() const;
    double getRayleighMinPressure() const;
    bool isRayleighSingleScatt();
    bool isDoNlteCorrection();

    bool isAddAerosl();
    bool isUserAerOptParam();
    bool isAddClouds();
    bool isUserCldOptParam();
    bool isGridBoxAvgCloud();
    int getCloudOverlap() const;
    double getCCLowCloudTop() const;
    double getCldcolThreshold() const;
    int getIrScattModel() const;
    int getVisScattModel() const;
    int getDomNstreams() const;
    double getDomAccuracy() const;
    double getDomOpdepThreshold() const;
    bool isDomRayleigh();

    bool isAddPC();
    bool isAddRadrec();
    int getIpcreg() const;
    int getIpcbnd() const;

    //RTTOV-SCATT
    bool isLuserCfrac();
    double getCCThreshold() const;
    double getIcePolarisation() const;
    bool isHydroCfracTLAD();
    bool isZeroHydroTLAD();

    //RTTOV opts%dev%
    bool isDoOpdepCalc();

    // Wrapper options
    int getNthreads() const;
    int getNprofsPerCall() const;
    bool isVerboseWrapper() const;
    bool isCheckOpts() const;
    bool isStoreRad() const;
    bool isStoreRad2() const;
    bool isStoreTrans() const;
    bool isStoreEmisTerms() const;

    std::string defineStrOptions();

private :

    StrBoolMap config;
    StrBoolMap ir;
    StrBoolMap mw;
    StrBoolMap pc;
    StrBoolMap interp;
    StrBoolMap rt_all;
    StrBoolMap scatt;
    StrBoolMap dev;
    bool verbose_wrapper;
    bool check_opts;
    bool store_trans;
    bool store_rad;
    bool store_rad2;
    bool store_emis_terms;

    int nthreads;
    int nprofs_per_call;
    int cloud_overlap;
    double cc_low_cloud_top;
    double cldcol_threshold;
    int ir_scatt_model;
    int vis_scatt_model;
    int dom_nstreams;
    double rayleigh_max_wavelength;
    double rayleigh_min_pressure;
    double dom_accuracy;
    double dom_opdep_threshold;
    int clw_scheme;
    double clw_cloud_top;
    int fastem_version;
    int ir_sea_emis_model;
    int solar_sea_brdf_model;
    int interp_mode;
    int ipcbnd;
    int ipcreg;
    double cc_threshold;
    double ice_polarisation;
};
} /* namespace rttov */
#endif /* OPTIONS_H_ */
