  subroutine rttov_setpressure (p_sfc, p, ph)

! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------

  Use parkind1, Only : jpim     ,jprb
!INTF_OFF
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
!INTF_ON
  implicit none
  
  Integer (Kind=jpim), parameter :: nlev = 60
  
  Real (Kind=jprb) :: p_sfc
  Real (Kind=jprb) :: p   (nlev)  , ph  (nlev+1)


!INTF_END
  Real (Kind=jprb) :: vah (nlev+1), vbh (nlev+1)
  Integer (Kind=jpim)            :: ilev

REAL(KIND=JPRB) :: ZHOOK_HANDLE
        
      data vah / &
          & 0.000000_JPRB,    20.000000_JPRB,    38.425343_JPRB, &
         & 63.647804_JPRB,    95.636963_JPRB,   134.483307_JPRB, &
        & 180.584351_JPRB,   234.779053_JPRB,   298.495789_JPRB, &
        & 373.971924_JPRB,   464.618134_JPRB,   575.651001_JPRB, &
        & 713.218079_JPRB,   883.660522_JPRB,  1094.834717_JPRB, &
       & 1356.474609_JPRB,  1680.640259_JPRB,  2082.273926_JPRB, &
       & 2579.888672_JPRB,  3196.421631_JPRB,  3960.291504_JPRB, &
       & 4906.708496_JPRB,  6018.019531_JPRB,  7306.631348_JPRB, &
       & 8765.053711_JPRB, 10376.126953_JPRB, 12077.446289_JPRB, &
      & 13775.325195_JPRB, 15379.805664_JPRB, 16819.474609_JPRB, &
      & 18045.183594_JPRB, 19027.695313_JPRB, 19755.109375_JPRB, &
      & 20222.205078_JPRB, 20429.863281_JPRB, 20384.480469_JPRB, &
      & 20097.402344_JPRB, 19584.330078_JPRB, 18864.750000_JPRB, &
      & 17961.357422_JPRB, 16899.468750_JPRB, 15706.447266_JPRB, &
      & 14411.124023_JPRB, 13043.218750_JPRB, 11632.758789_JPRB, &
      & 10209.500977_JPRB,  8802.356445_JPRB,  7438.803223_JPRB, &
       & 6144.314941_JPRB,  4941.778320_JPRB,  3850.913330_JPRB, &
       & 2887.696533_JPRB,  2063.779785_JPRB,  1385.912598_JPRB, &
        & 855.361755_JPRB,   467.333588_JPRB,   210.393890_JPRB, &
         & 65.889244_JPRB,     7.367743_JPRB,     0.000000_JPRB, &
          & 0.000000_JPRB &
              & / 
      data vbh / &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000000000_JPRB, 0.0000000000_JPRB, 0.0000000000_JPRB, &
      & 0.0000758235_JPRB, 0.0004613950_JPRB, 0.0018151561_JPRB, &
      & 0.0050811190_JPRB, 0.0111429105_JPRB, 0.0206778757_JPRB, &
      & 0.0341211632_JPRB, 0.0516904071_JPRB, 0.0735338330_JPRB, &
      & 0.0996746942_JPRB, 0.1300225109_JPRB, 0.1643843204_JPRB, &
      & 0.2024759352_JPRB, 0.2439331412_JPRB, 0.2883229554_JPRB, &
      & 0.3351548910_JPRB, 0.3838921487_JPRB, 0.4339629412_JPRB, &
      & 0.4847715795_JPRB, 0.5357099175_JPRB, 0.5861684084_JPRB, &
      & 0.6355474591_JPRB, 0.6832686067_JPRB, 0.7287858129_JPRB, &
      & 0.7715966105_JPRB, 0.8112534285_JPRB, 0.8473749161_JPRB, &
      & 0.8796569109_JPRB, 0.9078838825_JPRB, 0.9319403172_JPRB, &
      & 0.9518215060_JPRB, 0.9676452279_JPRB, 0.9796627164_JPRB, &
      & 0.9882701039_JPRB, 0.9940194488_JPRB, 0.9976301193_JPRB, &
      & 1.0000000000_JPRB &
              & / 

  !- End of header ------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RTTOV_SETPRESSURE',0_jpim,ZHOOK_HANDLE)
      ph (:) = vah (:) + vbh (:) * p_sfc
      
      do ilev = 1, nlev
         p (ilev) = 0.5_JPRB * (ph (ilev) + ph (ilev+1))
      end do
        
IF (LHOOK) CALL DR_HOOK('RTTOV_SETPRESSURE',1_jpim,ZHOOK_HANDLE)
      end subroutine rttov_setpressure
