! Description:
!> @file
!!   Data and subroutines for visible/near-IR water reflectance
!
!> @brief
!!   Data and subroutines for visible/near-IR water reflectance
!!
!! @details
!!   This data is from the USGS Digital Spectral Library 06:
!!   http://speclab.cr.usgs.gov/spectral.lib06/
!!
!!   This is used as the reflectance for water surfaces when computing
!!   the contribution from downward-scattered radiation in visible/near-IR
!!   channels. For the direct solar beam the sunglint BRDF model is used.
!!
!!   The wavenumber spectrum must match that in the BRDF atlas.
!
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
!
MODULE rttov_solar_refl_mod

  USE parkind1, ONLY : jprb, jpim

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: coastal_waters_ref, ocean_waters_ref, rttov_refl_water_interp

  INTEGER(jpim), PARAMETER :: numwave = 2101  !< Number of wavenumbers
  REAL(jprb),    PROTECTED :: coastal_waters_ref(numwave)   !< Coastal water reflectances
  REAL(jprb),    PROTECTED :: ocean_waters_ref(numwave)     !< Ocean water reflectances

  DATA coastal_waters_ref / &
    0.01199_jprb, 0.01232_jprb, 0.01261_jprb, 0.01290_jprb, 0.01319_jprb, 0.01348_jprb, 0.01371_jprb, 0.01388_jprb, 0.01405_jprb, &
    0.01421_jprb, 0.01436_jprb, 0.01443_jprb, 0.01449_jprb, 0.01455_jprb, 0.01460_jprb, 0.01464_jprb, 0.01468_jprb, 0.01472_jprb, &
    0.01478_jprb, 0.01484_jprb, 0.01490_jprb, 0.01497_jprb, 0.01504_jprb, 0.01511_jprb, 0.01518_jprb, 0.01524_jprb, 0.01530_jprb, &
    0.01535_jprb, 0.01539_jprb, 0.01543_jprb, 0.01547_jprb, 0.01551_jprb, 0.01555_jprb, 0.01558_jprb, 0.01562_jprb, 0.01566_jprb, &
    0.01570_jprb, 0.01575_jprb, 0.01579_jprb, 0.01584_jprb, 0.01588_jprb, 0.01592_jprb, 0.01596_jprb, 0.01599_jprb, 0.01602_jprb, &
    0.01605_jprb, 0.01607_jprb, 0.01609_jprb, 0.01611_jprb, 0.01613_jprb, 0.01615_jprb, 0.01618_jprb, 0.01621_jprb, 0.01626_jprb, &
    0.01630_jprb, 0.01635_jprb, 0.01638_jprb, 0.01642_jprb, 0.01645_jprb, 0.01647_jprb, 0.01650_jprb, 0.01652_jprb, 0.01654_jprb, &
    0.01657_jprb, 0.01659_jprb, 0.01662_jprb, 0.01665_jprb, 0.01668_jprb, 0.01672_jprb, 0.01675_jprb, 0.01679_jprb, 0.01683_jprb, &
    0.01687_jprb, 0.01691_jprb, 0.01694_jprb, 0.01698_jprb, 0.01701_jprb, 0.01703_jprb, 0.01706_jprb, 0.01708_jprb, 0.01710_jprb, &
    0.01712_jprb, 0.01714_jprb, 0.01715_jprb, 0.01717_jprb, 0.01718_jprb, 0.01719_jprb, 0.01720_jprb, 0.01721_jprb, 0.01722_jprb, &
    0.01724_jprb, 0.01725_jprb, 0.01726_jprb, 0.01727_jprb, 0.01729_jprb, 0.01730_jprb, 0.01732_jprb, 0.01734_jprb, 0.01736_jprb, &
    0.01738_jprb, 0.01741_jprb, 0.01744_jprb, 0.01747_jprb, 0.01750_jprb, 0.01753_jprb, 0.01756_jprb, 0.01759_jprb, 0.01761_jprb, &
    0.01764_jprb, 0.01766_jprb, 0.01768_jprb, 0.01769_jprb, 0.01770_jprb, 0.01771_jprb, 0.01771_jprb, 0.01771_jprb, 0.01771_jprb, &
    0.01772_jprb, 0.01773_jprb, 0.01775_jprb, 0.01776_jprb, 0.01779_jprb, 0.01781_jprb, 0.01783_jprb, 0.01785_jprb, 0.01788_jprb, &
    0.01790_jprb, 0.01792_jprb, 0.01794_jprb, 0.01796_jprb, 0.01798_jprb, 0.01800_jprb, 0.01801_jprb, 0.01803_jprb, 0.01804_jprb, &
    0.01805_jprb, 0.01807_jprb, 0.01808_jprb, 0.01809_jprb, 0.01810_jprb, 0.01811_jprb, 0.01811_jprb, 0.01812_jprb, 0.01813_jprb, &
    0.01813_jprb, 0.01814_jprb, 0.01815_jprb, 0.01815_jprb, 0.01816_jprb, 0.01816_jprb, 0.01817_jprb, 0.01818_jprb, 0.01818_jprb, &
    0.01819_jprb, 0.01820_jprb, 0.01821_jprb, 0.01822_jprb, 0.01823_jprb, 0.01824_jprb, 0.01825_jprb, 0.01826_jprb, 0.01827_jprb, &
    0.01828_jprb, 0.01830_jprb, 0.01831_jprb, 0.01833_jprb, 0.01834_jprb, 0.01835_jprb, 0.01837_jprb, 0.01838_jprb, 0.01839_jprb, &
    0.01841_jprb, 0.01842_jprb, 0.01843_jprb, 0.01845_jprb, 0.01846_jprb, 0.01847_jprb, 0.01848_jprb, 0.01849_jprb, 0.01850_jprb, &
    0.01851_jprb, 0.01851_jprb, 0.01852_jprb, 0.01852_jprb, 0.01852_jprb, 0.01852_jprb, 0.01852_jprb, 0.01852_jprb, 0.01852_jprb, &
    0.01851_jprb, 0.01850_jprb, 0.01849_jprb, 0.01848_jprb, 0.01848_jprb, 0.01847_jprb, 0.01846_jprb, 0.01846_jprb, 0.01846_jprb, &
    0.01846_jprb, 0.01846_jprb, 0.01848_jprb, 0.01849_jprb, 0.01851_jprb, 0.01853_jprb, 0.01856_jprb, 0.01859_jprb, 0.01861_jprb, &
    0.01864_jprb, 0.01865_jprb, 0.01867_jprb, 0.01868_jprb, 0.01868_jprb, 0.01869_jprb, 0.01869_jprb, 0.01869_jprb, 0.01868_jprb, &
    0.01868_jprb, 0.01867_jprb, 0.01867_jprb, 0.01866_jprb, 0.01865_jprb, 0.01864_jprb, 0.01863_jprb, 0.01863_jprb, 0.01862_jprb, &
    0.01862_jprb, 0.01861_jprb, 0.01861_jprb, 0.01861_jprb, 0.01861_jprb, 0.01862_jprb, 0.01862_jprb, 0.01862_jprb, 0.01863_jprb, &
    0.01863_jprb, 0.01864_jprb, 0.01865_jprb, 0.01865_jprb, 0.01866_jprb, 0.01867_jprb, 0.01868_jprb, 0.01869_jprb, 0.01869_jprb, &
    0.01870_jprb, 0.01871_jprb, 0.01872_jprb, 0.01872_jprb, 0.01873_jprb, 0.01874_jprb, 0.01874_jprb, 0.01875_jprb, 0.01875_jprb, &
    0.01876_jprb, 0.01876_jprb, 0.01877_jprb, 0.01877_jprb, 0.01878_jprb, 0.01878_jprb, 0.01879_jprb, 0.01879_jprb, 0.01880_jprb, &
    0.01880_jprb, 0.01880_jprb, 0.01881_jprb, 0.01881_jprb, 0.01881_jprb, 0.01882_jprb, 0.01882_jprb, 0.01883_jprb, 0.01883_jprb, &
    0.01883_jprb, 0.01884_jprb, 0.01884_jprb, 0.01885_jprb, 0.01885_jprb, 0.01885_jprb, 0.01886_jprb, 0.01886_jprb, 0.01887_jprb, &
    0.01887_jprb, 0.01888_jprb, 0.01888_jprb, 0.01888_jprb, 0.01889_jprb, 0.01889_jprb, 0.01890_jprb, 0.01890_jprb, 0.01891_jprb, &
    0.01891_jprb, 0.01892_jprb, 0.01892_jprb, 0.01893_jprb, 0.01893_jprb, 0.01894_jprb, 0.01895_jprb, 0.01895_jprb, 0.01896_jprb, &
    0.01896_jprb, 0.01897_jprb, 0.01897_jprb, 0.01898_jprb, 0.01898_jprb, 0.01898_jprb, 0.01899_jprb, 0.01899_jprb, 0.01900_jprb, &
    0.01900_jprb, 0.01901_jprb, 0.01901_jprb, 0.01901_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, &
    0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, &
    0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, &
    0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01903_jprb, &
    0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01904_jprb, 0.01904_jprb, 0.01904_jprb, 0.01904_jprb, &
    0.01905_jprb, 0.01905_jprb, 0.01905_jprb, 0.01906_jprb, 0.01906_jprb, 0.01906_jprb, 0.01907_jprb, 0.01907_jprb, 0.01908_jprb, &
    0.01908_jprb, 0.01909_jprb, 0.01909_jprb, 0.01910_jprb, 0.01910_jprb, 0.01910_jprb, 0.01911_jprb, 0.01912_jprb, 0.01912_jprb, &
    0.01913_jprb, 0.01913_jprb, 0.01914_jprb, 0.01914_jprb, 0.01915_jprb, 0.01915_jprb, 0.01916_jprb, 0.01917_jprb, 0.01917_jprb, &
    0.01918_jprb, 0.01918_jprb, 0.01919_jprb, 0.01919_jprb, 0.01920_jprb, 0.01920_jprb, 0.01920_jprb, 0.01921_jprb, 0.01921_jprb, &
    0.01921_jprb, 0.01922_jprb, 0.01922_jprb, 0.01922_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, &
    0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01924_jprb, &
    0.01924_jprb, 0.01924_jprb, 0.01924_jprb, 0.01925_jprb, 0.01925_jprb, 0.01925_jprb, 0.01926_jprb, 0.01926_jprb, 0.01927_jprb, &
    0.01927_jprb, 0.01928_jprb, 0.01928_jprb, 0.01928_jprb, 0.01929_jprb, 0.01929_jprb, 0.01930_jprb, 0.01930_jprb, 0.01930_jprb, &
    0.01931_jprb, 0.01931_jprb, 0.01932_jprb, 0.01932_jprb, 0.01932_jprb, 0.01932_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, &
    0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, &
    0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01934_jprb, 0.01934_jprb, 0.01934_jprb, 0.01934_jprb, 0.01934_jprb, 0.01934_jprb, &
    0.01934_jprb, 0.01935_jprb, 0.01935_jprb, 0.01935_jprb, 0.01935_jprb, 0.01935_jprb, 0.01936_jprb, 0.01936_jprb, 0.01936_jprb, &
    0.01936_jprb, 0.01937_jprb, 0.01937_jprb, 0.01937_jprb, 0.01937_jprb, 0.01937_jprb, 0.01938_jprb, 0.01938_jprb, 0.01938_jprb, &
    0.01938_jprb, 0.01938_jprb, 0.01939_jprb, 0.01939_jprb, 0.01939_jprb, 0.01939_jprb, 0.01939_jprb, 0.01939_jprb, 0.01939_jprb, &
    0.01939_jprb, 0.01939_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, &
    0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01941_jprb, 0.01941_jprb, 0.01941_jprb, &
    0.01941_jprb, 0.01941_jprb, 0.01941_jprb, 0.01941_jprb, 0.01941_jprb, 0.01942_jprb, 0.01942_jprb, 0.01942_jprb, 0.01942_jprb, &
    0.01942_jprb, 0.01943_jprb, 0.01943_jprb, 0.01943_jprb, 0.01943_jprb, 0.01944_jprb, 0.01944_jprb, 0.01944_jprb, 0.01945_jprb, &
    0.01945_jprb, 0.01945_jprb, 0.01946_jprb, 0.01946_jprb, 0.01946_jprb, 0.01947_jprb, 0.01947_jprb, 0.01947_jprb, 0.01948_jprb, &
    0.01948_jprb, 0.01949_jprb, 0.01949_jprb, 0.01949_jprb, 0.01950_jprb, 0.01950_jprb, 0.01950_jprb, 0.01951_jprb, 0.01951_jprb, &
    0.01951_jprb, 0.01952_jprb, 0.01952_jprb, 0.01952_jprb, 0.01953_jprb, 0.01953_jprb, 0.01953_jprb, 0.01953_jprb, 0.01953_jprb, &
    0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, &
    0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, &
    0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01955_jprb, 0.01955_jprb, 0.01955_jprb, 0.01955_jprb, &
    0.01955_jprb, 0.01955_jprb, 0.01955_jprb, 0.01956_jprb, 0.01956_jprb, 0.01956_jprb, 0.01956_jprb, 0.01956_jprb, 0.01956_jprb, &
    0.01957_jprb, 0.01957_jprb, 0.01957_jprb, 0.01957_jprb, 0.01958_jprb, 0.01958_jprb, 0.01958_jprb, 0.01958_jprb, 0.01959_jprb, &
    0.01959_jprb, 0.01959_jprb, 0.01960_jprb, 0.01960_jprb, 0.01960_jprb, 0.01961_jprb, 0.01961_jprb, 0.01961_jprb, 0.01962_jprb, &
    0.01962_jprb, 0.01962_jprb, 0.01963_jprb, 0.01963_jprb, 0.01964_jprb, 0.01964_jprb, 0.01964_jprb, 0.01965_jprb, 0.01965_jprb, &
    0.01966_jprb, 0.01966_jprb, 0.01966_jprb, 0.01967_jprb, 0.01967_jprb, 0.01968_jprb, 0.01968_jprb, 0.01969_jprb, 0.01969_jprb, &
    0.01969_jprb, 0.01970_jprb, 0.01970_jprb, 0.01971_jprb, 0.01971_jprb, 0.01971_jprb, 0.01972_jprb, 0.01972_jprb, 0.01972_jprb, &
    0.01973_jprb, 0.01973_jprb, 0.01973_jprb, 0.01973_jprb, 0.01974_jprb, 0.01974_jprb, 0.01974_jprb, 0.01974_jprb, 0.01974_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, &
    0.01977_jprb, 0.01977_jprb, 0.01977_jprb, 0.01977_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01979_jprb, &
    0.01979_jprb, 0.01979_jprb, 0.01979_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01981_jprb, 0.01981_jprb, &
    0.01981_jprb, 0.01981_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01983_jprb, 0.01983_jprb, &
    0.01983_jprb, 0.01983_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, &
    0.01984_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, &
    0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01984_jprb, &
    0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01983_jprb, 0.01983_jprb, 0.01983_jprb, 0.01983_jprb, &
    0.01983_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01981_jprb, 0.01981_jprb, 0.01981_jprb, &
    0.01981_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01979_jprb, 0.01979_jprb, 0.01979_jprb, &
    0.01979_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01977_jprb, 0.01977_jprb, 0.01977_jprb, &
    0.01977_jprb, 0.01977_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, &
    0.01977_jprb, 0.01977_jprb, 0.01977_jprb, 0.01977_jprb, 0.01977_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, &
    0.01978_jprb, 0.01979_jprb, 0.01979_jprb, 0.01979_jprb, 0.01979_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, &
    0.01981_jprb, 0.01981_jprb, 0.01981_jprb, 0.01981_jprb, 0.01981_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, &
    0.01982_jprb, 0.01983_jprb, 0.01983_jprb, 0.01983_jprb, 0.01983_jprb, 0.01983_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, &
    0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, &
    0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, &
    0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, &
    0.01985_jprb, 0.01985_jprb, 0.01986_jprb, 0.01986_jprb, 0.01986_jprb, 0.01986_jprb, 0.01986_jprb, 0.01986_jprb, 0.01986_jprb, &
    0.01987_jprb, 0.01987_jprb, 0.01987_jprb, 0.01987_jprb, 0.01987_jprb, 0.01987_jprb, 0.01988_jprb, 0.01988_jprb, 0.01988_jprb, &
    0.01988_jprb, 0.01988_jprb, 0.01989_jprb, 0.01989_jprb, 0.01989_jprb, 0.01989_jprb, 0.01989_jprb, 0.01990_jprb, 0.01990_jprb, &
    0.01990_jprb, 0.01990_jprb, 0.01991_jprb, 0.01991_jprb, 0.01991_jprb, 0.01991_jprb, 0.01991_jprb, 0.01991_jprb, 0.01991_jprb, &
    0.01992_jprb, 0.01993_jprb, 0.01995_jprb, 0.01997_jprb, 0.01998_jprb, 0.02000_jprb, 0.02002_jprb, 0.02004_jprb, 0.02006_jprb, &
    0.02008_jprb, 0.02009_jprb, 0.02011_jprb, 0.02013_jprb, 0.02015_jprb, 0.02017_jprb, 0.02019_jprb, 0.02021_jprb, 0.02022_jprb, &
    0.02024_jprb, 0.02026_jprb, 0.02028_jprb, 0.02030_jprb, 0.02031_jprb, 0.02033_jprb, 0.02035_jprb, 0.02037_jprb, 0.02038_jprb, &
    0.02040_jprb, 0.02042_jprb, 0.02043_jprb, 0.02045_jprb, 0.02047_jprb, 0.02048_jprb, 0.02050_jprb, 0.02052_jprb, 0.02053_jprb, &
    0.02055_jprb, 0.02057_jprb, 0.02058_jprb, 0.02060_jprb, 0.02062_jprb, 0.02063_jprb, 0.02064_jprb, 0.02065_jprb, 0.02065_jprb, &
    0.02065_jprb, 0.02066_jprb, 0.02065_jprb, 0.02065_jprb, 0.02064_jprb, 0.02064_jprb, 0.02065_jprb, 0.02066_jprb, 0.02066_jprb, &
    0.02068_jprb, 0.02069_jprb, 0.02070_jprb, 0.02072_jprb, 0.02075_jprb, 0.02077_jprb, 0.02080_jprb, 0.02082_jprb, 0.02084_jprb, &
    0.02085_jprb, 0.02087_jprb, 0.02088_jprb, 0.02089_jprb, 0.02090_jprb, 0.02092_jprb, 0.02095_jprb, 0.02098_jprb, 0.02101_jprb, &
    0.02104_jprb, 0.02107_jprb, 0.02110_jprb, 0.02113_jprb, 0.02116_jprb, 0.02119_jprb, 0.02121_jprb, 0.02124_jprb, 0.02128_jprb, &
    0.02133_jprb, 0.02137_jprb, 0.02141_jprb, 0.02145_jprb, 0.02149_jprb, 0.02154_jprb, 0.02157_jprb, 0.02160_jprb, 0.02164_jprb, &
    0.02167_jprb, 0.02170_jprb, 0.02174_jprb, 0.02177_jprb, 0.02181_jprb, 0.02184_jprb, 0.02187_jprb, 0.02191_jprb, 0.02194_jprb, &
    0.02198_jprb, 0.02202_jprb, 0.02206_jprb, 0.02210_jprb, 0.02215_jprb, 0.02219_jprb, 0.02223_jprb, 0.02228_jprb, 0.02233_jprb, &
    0.02238_jprb, 0.02243_jprb, 0.02248_jprb, 0.02253_jprb, 0.02258_jprb, 0.02263_jprb, 0.02268_jprb, 0.02273_jprb, 0.02278_jprb, &
    0.02283_jprb, 0.02288_jprb, 0.02295_jprb, 0.02302_jprb, 0.02308_jprb, 0.02315_jprb, 0.02321_jprb, 0.02328_jprb, 0.02335_jprb, &
    0.02341_jprb, 0.02348_jprb, 0.02355_jprb, 0.02361_jprb, 0.02369_jprb, 0.02379_jprb, 0.02388_jprb, 0.02397_jprb, 0.02406_jprb, &
    0.02416_jprb, 0.02425_jprb, 0.02434_jprb, 0.02443_jprb, 0.02452_jprb, 0.02460_jprb, 0.02469_jprb, 0.02477_jprb, 0.02486_jprb, &
    0.02494_jprb, 0.02502_jprb, 0.02510_jprb, 0.02518_jprb, 0.02526_jprb, 0.02534_jprb, 0.02542_jprb, 0.02548_jprb, 0.02554_jprb, &
    0.02560_jprb, 0.02566_jprb, 0.02572_jprb, 0.02578_jprb, 0.02583_jprb, 0.02588_jprb, 0.02593_jprb, 0.02598_jprb, 0.02602_jprb, &
    0.02607_jprb, 0.02607_jprb, 0.02606_jprb, 0.02606_jprb, 0.02605_jprb, 0.02604_jprb, 0.02604_jprb, 0.02602_jprb, 0.02600_jprb, &
    0.02598_jprb, 0.02597_jprb, 0.02595_jprb, 0.02593_jprb, 0.02591_jprb, 0.02586_jprb, 0.02582_jprb, 0.02578_jprb, 0.02574_jprb, &
    0.02570_jprb, 0.02566_jprb, 0.02561_jprb, 0.02557_jprb, 0.02552_jprb, 0.02548_jprb, 0.02543_jprb, 0.02539_jprb, 0.02534_jprb, &
    0.02528_jprb, 0.02522_jprb, 0.02516_jprb, 0.02510_jprb, 0.02503_jprb, 0.02497_jprb, 0.02491_jprb, 0.02486_jprb, 0.02480_jprb, &
    0.02474_jprb, 0.02468_jprb, 0.02462_jprb, 0.02456_jprb, 0.02452_jprb, 0.02448_jprb, 0.02444_jprb, 0.02439_jprb, 0.02435_jprb, &
    0.02431_jprb, 0.02429_jprb, 0.02426_jprb, 0.02423_jprb, 0.02420_jprb, 0.02419_jprb, 0.02417_jprb, 0.02416_jprb, 0.02415_jprb, &
    0.02414_jprb, 0.02414_jprb, 0.02413_jprb, 0.02412_jprb, 0.02412_jprb, 0.02412_jprb, 0.02412_jprb, 0.02412_jprb, 0.02412_jprb, &
    0.02413_jprb, 0.02413_jprb, 0.02414_jprb, 0.02415_jprb, 0.02416_jprb, 0.02417_jprb, 0.02418_jprb, 0.02419_jprb, 0.02421_jprb, &
    0.02422_jprb, 0.02423_jprb, 0.02424_jprb, 0.02425_jprb, 0.02426_jprb, 0.02428_jprb, 0.02430_jprb, 0.02432_jprb, 0.02434_jprb, &
    0.02437_jprb, 0.02439_jprb, 0.02442_jprb, 0.02444_jprb, 0.02447_jprb, 0.02450_jprb, 0.02451_jprb, 0.02452_jprb, 0.02453_jprb, &
    0.02454_jprb, 0.02455_jprb, 0.02457_jprb, 0.02459_jprb, 0.02460_jprb, 0.02462_jprb, 0.02464_jprb, 0.02465_jprb, 0.02467_jprb, &
    0.02469_jprb, 0.02471_jprb, 0.02473_jprb, 0.02474_jprb, 0.02476_jprb, 0.02477_jprb, 0.02479_jprb, 0.02480_jprb, 0.02482_jprb, &
    0.02483_jprb, 0.02485_jprb, 0.02486_jprb, 0.02488_jprb, 0.02490_jprb, 0.02491_jprb, 0.02493_jprb, 0.02495_jprb, 0.02496_jprb, &
    0.02497_jprb, 0.02498_jprb, 0.02500_jprb, 0.02501_jprb, 0.02502_jprb, 0.02503_jprb, 0.02504_jprb, 0.02505_jprb, 0.02506_jprb, &
    0.02508_jprb, 0.02509_jprb, 0.02510_jprb, 0.02512_jprb, 0.02513_jprb, 0.02514_jprb, 0.02515_jprb, 0.02517_jprb, 0.02518_jprb, &
    0.02519_jprb, 0.02520_jprb, 0.02521_jprb, 0.02523_jprb, 0.02524_jprb, 0.02526_jprb, 0.02528_jprb, 0.02530_jprb, 0.02532_jprb, &
    0.02535_jprb, 0.02537_jprb, 0.02539_jprb, 0.02542_jprb, 0.02544_jprb, 0.02546_jprb, 0.02548_jprb, 0.02550_jprb, 0.02553_jprb, &
    0.02555_jprb, 0.02557_jprb, 0.02559_jprb, 0.02562_jprb, 0.02564_jprb, 0.02566_jprb, 0.02569_jprb, 0.02571_jprb, 0.02573_jprb, &
    0.02576_jprb, 0.02578_jprb, 0.02580_jprb, 0.02582_jprb, 0.02584_jprb, 0.02585_jprb, 0.02587_jprb, 0.02590_jprb, 0.02592_jprb, &
    0.02594_jprb, 0.02597_jprb, 0.02600_jprb, 0.02604_jprb, 0.02607_jprb, 0.02610_jprb, 0.02613_jprb, 0.02616_jprb, 0.02619_jprb, &
    0.02622_jprb, 0.02624_jprb, 0.02626_jprb, 0.02629_jprb, 0.02631_jprb, 0.02634_jprb, 0.02637_jprb, 0.02641_jprb, 0.02645_jprb, &
    0.02649_jprb, 0.02654_jprb, 0.02659_jprb, 0.02665_jprb, 0.02670_jprb, 0.02675_jprb, 0.02680_jprb, 0.02686_jprb, 0.02691_jprb, &
    0.02696_jprb, 0.02702_jprb, 0.02708_jprb, 0.02713_jprb, 0.02719_jprb, 0.02727_jprb, 0.02735_jprb, 0.02742_jprb, 0.02750_jprb, &
    0.02758_jprb, 0.02767_jprb, 0.02775_jprb, 0.02784_jprb, 0.02792_jprb, 0.02800_jprb, 0.02808_jprb, 0.02816_jprb, 0.02824_jprb, &
    0.02833_jprb, 0.02842_jprb, 0.02851_jprb, 0.02860_jprb, 0.02868_jprb, 0.02877_jprb, 0.02888_jprb, 0.02899_jprb, 0.02910_jprb, &
    0.02921_jprb, 0.02932_jprb, 0.02943_jprb, 0.02955_jprb, 0.02967_jprb, 0.02979_jprb, 0.02991_jprb, 0.03003_jprb, 0.03015_jprb, &
    0.03029_jprb, 0.03042_jprb, 0.03056_jprb, 0.03070_jprb, 0.03083_jprb, 0.03097_jprb, 0.03111_jprb, 0.03126_jprb, 0.03140_jprb, &
    0.03154_jprb, 0.03169_jprb, 0.03182_jprb, 0.03195_jprb, 0.03208_jprb, 0.03222_jprb, 0.03235_jprb, 0.03248_jprb, 0.03260_jprb, &
    0.03273_jprb, 0.03286_jprb, 0.03298_jprb, 0.03311_jprb, 0.03323_jprb, 0.03336_jprb, 0.03348_jprb, 0.03361_jprb, 0.03374_jprb, &
    0.03386_jprb, 0.03399_jprb, 0.03412_jprb, 0.03425_jprb, 0.03437_jprb, 0.03450_jprb, 0.03463_jprb, 0.03476_jprb, 0.03489_jprb, &
    0.03501_jprb, 0.03514_jprb, 0.03527_jprb, 0.03540_jprb, 0.03553_jprb, 0.03566_jprb, 0.03578_jprb, 0.03591_jprb, 0.03604_jprb, &
    0.03616_jprb, 0.03628_jprb, 0.03638_jprb, 0.03649_jprb, 0.03660_jprb, 0.03671_jprb, 0.03682_jprb, 0.03691_jprb, 0.03700_jprb, &
    0.03709_jprb, 0.03718_jprb, 0.03727_jprb, 0.03736_jprb, 0.03744_jprb, 0.03752_jprb, 0.03760_jprb, 0.03768_jprb, 0.03776_jprb, &
    0.03784_jprb, 0.03791_jprb, 0.03797_jprb, 0.03803_jprb, 0.03809_jprb, 0.03815_jprb, 0.03821_jprb, 0.03827_jprb, 0.03831_jprb, &
    0.03836_jprb, 0.03841_jprb, 0.03845_jprb, 0.03850_jprb, 0.03855_jprb, 0.03859_jprb, 0.03863_jprb, 0.03866_jprb, 0.03870_jprb, &
    0.03874_jprb, 0.03878_jprb, 0.03881_jprb, 0.03883_jprb, 0.03885_jprb, 0.03887_jprb, 0.03889_jprb, 0.03892_jprb, 0.03892_jprb, &
    0.03892_jprb, 0.03892_jprb, 0.03892_jprb, 0.03891_jprb, 0.03891_jprb, 0.03890_jprb, 0.03888_jprb, 0.03886_jprb, 0.03884_jprb, &
    0.03882_jprb, 0.03880_jprb, 0.03878_jprb, 0.03876_jprb, 0.03874_jprb, 0.03872_jprb, 0.03869_jprb, 0.03867_jprb, 0.03865_jprb, &
    0.03862_jprb, 0.03860_jprb, 0.03858_jprb, 0.03856_jprb, 0.03854_jprb, 0.03851_jprb, 0.03849_jprb, 0.03847_jprb, 0.03845_jprb, &
    0.03843_jprb, 0.03841_jprb, 0.03838_jprb, 0.03836_jprb, 0.03833_jprb, 0.03829_jprb, 0.03826_jprb, 0.03822_jprb, 0.03819_jprb, &
    0.03815_jprb, 0.03811_jprb, 0.03807_jprb, 0.03803_jprb, 0.03798_jprb, 0.03794_jprb, 0.03790_jprb, 0.03785_jprb, 0.03781_jprb, &
    0.03777_jprb, 0.03772_jprb, 0.03768_jprb, 0.03763_jprb, 0.03759_jprb, 0.03754_jprb, 0.03750_jprb, 0.03745_jprb, 0.03741_jprb, &
    0.03736_jprb, 0.03732_jprb, 0.03727_jprb, 0.03723_jprb, 0.03719_jprb, 0.03715_jprb, 0.03711_jprb, 0.03707_jprb, 0.03703_jprb, &
    0.03699_jprb, 0.03695_jprb, 0.03692_jprb, 0.03688_jprb, 0.03685_jprb, 0.03681_jprb, 0.03677_jprb, 0.03674_jprb, 0.03669_jprb, &
    0.03665_jprb, 0.03660_jprb, 0.03655_jprb, 0.03650_jprb, 0.03645_jprb, 0.03641_jprb, 0.03637_jprb, 0.03633_jprb, 0.03629_jprb, &
    0.03625_jprb, 0.03621_jprb, 0.03617_jprb, 0.03613_jprb, 0.03609_jprb, 0.03606_jprb, 0.03603_jprb, 0.03600_jprb, 0.03597_jprb, &
    0.03594_jprb, 0.03591_jprb, 0.03588_jprb, 0.03584_jprb, 0.03581_jprb, 0.03577_jprb, 0.03574_jprb, 0.03570_jprb, 0.03567_jprb, &
    0.03564_jprb, 0.03561_jprb, 0.03558_jprb, 0.03555_jprb, 0.03552_jprb, 0.03549_jprb, 0.03546_jprb, 0.03543_jprb, 0.03541_jprb, &
    0.03538_jprb, 0.03536_jprb, 0.03533_jprb, 0.03530_jprb, 0.03528_jprb, 0.03525_jprb, 0.03522_jprb, 0.03519_jprb, 0.03516_jprb, &
    0.03513_jprb, 0.03510_jprb, 0.03507_jprb, 0.03503_jprb, 0.03500_jprb, 0.03497_jprb, 0.03494_jprb, 0.03491_jprb, 0.03488_jprb, &
    0.03485_jprb, 0.03482_jprb, 0.03479_jprb, 0.03475_jprb, 0.03472_jprb, 0.03469_jprb, 0.03466_jprb, 0.03463_jprb, 0.03460_jprb, &
    0.03457_jprb, 0.03454_jprb, 0.03451_jprb, 0.03448_jprb, 0.03445_jprb, 0.03443_jprb, 0.03440_jprb, 0.03437_jprb, 0.03435_jprb, &
    0.03432_jprb, 0.03430_jprb, 0.03428_jprb, 0.03425_jprb, 0.03423_jprb, 0.03420_jprb, 0.03418_jprb, 0.03415_jprb, 0.03412_jprb, &
    0.03410_jprb, 0.03407_jprb, 0.03405_jprb, 0.03402_jprb, 0.03399_jprb, 0.03396_jprb, 0.03393_jprb, 0.03390_jprb, 0.03387_jprb, &
    0.03384_jprb, 0.03381_jprb, 0.03378_jprb, 0.03375_jprb, 0.03373_jprb, 0.03371_jprb, 0.03368_jprb, 0.03366_jprb, 0.03363_jprb, &
    0.03361_jprb, 0.03358_jprb, 0.03356_jprb, 0.03354_jprb, 0.03352_jprb, 0.03350_jprb, 0.03348_jprb, 0.03346_jprb, 0.03344_jprb, &
    0.03343_jprb, 0.03341_jprb, 0.03339_jprb, 0.03338_jprb, 0.03336_jprb, 0.03335_jprb, 0.03333_jprb, 0.03332_jprb, 0.03331_jprb, &
    0.03329_jprb, 0.03328_jprb, 0.03327_jprb, 0.03326_jprb, 0.03324_jprb, 0.03323_jprb, 0.03322_jprb, 0.03321_jprb, 0.03321_jprb, &
    0.03320_jprb, 0.03319_jprb, 0.03318_jprb, 0.03317_jprb, 0.03317_jprb, 0.03316_jprb, 0.03315_jprb, 0.03314_jprb, 0.03313_jprb, &
    0.03312_jprb, 0.03312_jprb, 0.03311_jprb, 0.03310_jprb, 0.03309_jprb, 0.03307_jprb, 0.03306_jprb, 0.03304_jprb, 0.03303_jprb, &
    0.03302_jprb, 0.03300_jprb, 0.03299_jprb, 0.03298_jprb, 0.03296_jprb, 0.03295_jprb, 0.03293_jprb, 0.03292_jprb, 0.03290_jprb, &
    0.03289_jprb, 0.03287_jprb, 0.03286_jprb, 0.03284_jprb, 0.03283_jprb, 0.03282_jprb, 0.03280_jprb, 0.03279_jprb, 0.03278_jprb, &
    0.03276_jprb, 0.03275_jprb, 0.03274_jprb, 0.03272_jprb, 0.03271_jprb, 0.03270_jprb, 0.03269_jprb, 0.03267_jprb, 0.03266_jprb, &
    0.03264_jprb, 0.03262_jprb, 0.03261_jprb, 0.03259_jprb, 0.03257_jprb, 0.03255_jprb, 0.03253_jprb, 0.03251_jprb, 0.03250_jprb, &
    0.03248_jprb, 0.03246_jprb, 0.03244_jprb, 0.03243_jprb, 0.03241_jprb, 0.03239_jprb, 0.03238_jprb, 0.03236_jprb, 0.03235_jprb, &
    0.03234_jprb, 0.03233_jprb, 0.03232_jprb, 0.03231_jprb, 0.03230_jprb, 0.03229_jprb, 0.03228_jprb, 0.03228_jprb, 0.03227_jprb, &
    0.03227_jprb, 0.03227_jprb, 0.03227_jprb, 0.03227_jprb, 0.03226_jprb, 0.03226_jprb, 0.03226_jprb, 0.03226_jprb, 0.03225_jprb, &
    0.03225_jprb, 0.03225_jprb, 0.03225_jprb, 0.03224_jprb, 0.03224_jprb, 0.03224_jprb, 0.03223_jprb, 0.03222_jprb, 0.03221_jprb, &
    0.03220_jprb, 0.03220_jprb, 0.03219_jprb, 0.03218_jprb, 0.03217_jprb, 0.03216_jprb, 0.03214_jprb, 0.03212_jprb, 0.03211_jprb, &
    0.03209_jprb, 0.03208_jprb, 0.03206_jprb, 0.03204_jprb, 0.03203_jprb, 0.03201_jprb, 0.03199_jprb, 0.03198_jprb, 0.03196_jprb, &
    0.03195_jprb, 0.03193_jprb, 0.03192_jprb, 0.03190_jprb, 0.03189_jprb, 0.03187_jprb, 0.03186_jprb, 0.03185_jprb, 0.03183_jprb, &
    0.03182_jprb, 0.03181_jprb, 0.03179_jprb, 0.03178_jprb, 0.03177_jprb, 0.03175_jprb, 0.03174_jprb, 0.03173_jprb, 0.03172_jprb, &
    0.03171_jprb, 0.03170_jprb, 0.03168_jprb, 0.03167_jprb, 0.03165_jprb, 0.03164_jprb, 0.03162_jprb, 0.03160_jprb, 0.03159_jprb, &
    0.03157_jprb, 0.03155_jprb, 0.03153_jprb, 0.03152_jprb, 0.03150_jprb, 0.03148_jprb, 0.03146_jprb, 0.03144_jprb, 0.03142_jprb, &
    0.03140_jprb, 0.03139_jprb, 0.03137_jprb, 0.03135_jprb, 0.03133_jprb, 0.03130_jprb, 0.03128_jprb, 0.03126_jprb, 0.03124_jprb, &
    0.03121_jprb, 0.03119_jprb, 0.03117_jprb, 0.03115_jprb, 0.03112_jprb, 0.03110_jprb, 0.03108_jprb, 0.03105_jprb, 0.03103_jprb, &
    0.03100_jprb, 0.03098_jprb, 0.03096_jprb, 0.03093_jprb, 0.03091_jprb, 0.03089_jprb, 0.03087_jprb, 0.03085_jprb, 0.03083_jprb, &
    0.03080_jprb, 0.03078_jprb, 0.03076_jprb, 0.03074_jprb, 0.03072_jprb, 0.03070_jprb, 0.03069_jprb, 0.03067_jprb, 0.03065_jprb, &
    0.03063_jprb, 0.03061_jprb, 0.03059_jprb, 0.03058_jprb, 0.03056_jprb, 0.03054_jprb, 0.03053_jprb, 0.03051_jprb, 0.03049_jprb, &
    0.03048_jprb, 0.03046_jprb, 0.03044_jprb, 0.03043_jprb, 0.03041_jprb, 0.03039_jprb, 0.03038_jprb, 0.03036_jprb, 0.03034_jprb, &
    0.03033_jprb, 0.03031_jprb, 0.03029_jprb, 0.03028_jprb, 0.03026_jprb, 0.03025_jprb, 0.03024_jprb, 0.03024_jprb, 0.03023_jprb, &
    0.03022_jprb, 0.03022_jprb, 0.03021_jprb, 0.03020_jprb, 0.03019_jprb, 0.03019_jprb, 0.03018_jprb, 0.03018_jprb, 0.03017_jprb, &
    0.03017_jprb, 0.03016_jprb, 0.03016_jprb, 0.03015_jprb, 0.03015_jprb, 0.03015_jprb, 0.03014_jprb, 0.03013_jprb, 0.03012_jprb, &
    0.03012_jprb, 0.03011_jprb, 0.03010_jprb, 0.03010_jprb, 0.03009_jprb, 0.03008_jprb, 0.03007_jprb, 0.03007_jprb, 0.03006_jprb, &
    0.03005_jprb, 0.03005_jprb, 0.03004_jprb, 0.03003_jprb, 0.03002_jprb, 0.03002_jprb, 0.03001_jprb, 0.03000_jprb, 0.02999_jprb, &
    0.02998_jprb, 0.02997_jprb, 0.02996_jprb, 0.02995_jprb, 0.02994_jprb, 0.02993_jprb, 0.02992_jprb, 0.02990_jprb, 0.02989_jprb, &
    0.02988_jprb, 0.02987_jprb, 0.02986_jprb, 0.02986_jprb, 0.02985_jprb, 0.02985_jprb, 0.02984_jprb, 0.02983_jprb, 0.02983_jprb, &
    0.02982_jprb, 0.02982_jprb, 0.02981_jprb, 0.02980_jprb, 0.02980_jprb, 0.02979_jprb, 0.02979_jprb, 0.02978_jprb, 0.02977_jprb, &
    0.02975_jprb, 0.02974_jprb, 0.02973_jprb, 0.02971_jprb, 0.02970_jprb, 0.02968_jprb, 0.02967_jprb, 0.02966_jprb, 0.02964_jprb, &
    0.02963_jprb, 0.02961_jprb, 0.02960_jprb, 0.02959_jprb, 0.02957_jprb, 0.02957_jprb, 0.02956_jprb, 0.02956_jprb, 0.02956_jprb, &
    0.02955_jprb, 0.02955_jprb, 0.02955_jprb, 0.02955_jprb, 0.02954_jprb, 0.02954_jprb, 0.02954_jprb, 0.02954_jprb, 0.02953_jprb, &
    0.02953_jprb, 0.02953_jprb, 0.02952_jprb, 0.02951_jprb, 0.02949_jprb, 0.02948_jprb, 0.02947_jprb, 0.02945_jprb, 0.02944_jprb, &
    0.02942_jprb, 0.02941_jprb, 0.02940_jprb, 0.02938_jprb, 0.02937_jprb, 0.02936_jprb, 0.02934_jprb, 0.02933_jprb, 0.02931_jprb, &
    0.02930_jprb, 0.02929_jprb, 0.02929_jprb, 0.02928_jprb, 0.02927_jprb, 0.02926_jprb, 0.02925_jprb, 0.02924_jprb, 0.02923_jprb, &
    0.02922_jprb, 0.02922_jprb, 0.02921_jprb, 0.02920_jprb, 0.02919_jprb, 0.02918_jprb, 0.02917_jprb, 0.02917_jprb, 0.02917_jprb, &
    0.02917_jprb, 0.02918_jprb, 0.02918_jprb, 0.02918_jprb, 0.02918_jprb, 0.02918_jprb, 0.02918_jprb, 0.02918_jprb, 0.02919_jprb, &
    0.02919_jprb, 0.02919_jprb, 0.02919_jprb, 0.02919_jprb, 0.02919_jprb, 0.02919_jprb, 0.02920_jprb, 0.02920_jprb, 0.02920_jprb, &
    0.02920_jprb, 0.02920_jprb, 0.02921_jprb, 0.02921_jprb, 0.02921_jprb, 0.02921_jprb, 0.02922_jprb, 0.02922_jprb, 0.02922_jprb, &
    0.02922_jprb, 0.02922_jprb, 0.02923_jprb, 0.02923_jprb, 0.02924_jprb, 0.02925_jprb, 0.02925_jprb, 0.02926_jprb, 0.02927_jprb, &
    0.02928_jprb, 0.02929_jprb, 0.02930_jprb, 0.02930_jprb, 0.02931_jprb, 0.02932_jprb, 0.02933_jprb, 0.02934_jprb, 0.02934_jprb, &
    0.02935_jprb, 0.02936_jprb, 0.02938_jprb, 0.02939_jprb, 0.02941_jprb, 0.02943_jprb, 0.02945_jprb, 0.02946_jprb, 0.02948_jprb, &
    0.02950_jprb, 0.02952_jprb, 0.02953_jprb, 0.02955_jprb, 0.02957_jprb, 0.02959_jprb, 0.02960_jprb, 0.02962_jprb, 0.02964_jprb, &
    0.02966_jprb, 0.02967_jprb, 0.02969_jprb, 0.02970_jprb, 0.02972_jprb, 0.02973_jprb, 0.02975_jprb, 0.02976_jprb, 0.02978_jprb, &
    0.02979_jprb, 0.02981_jprb, 0.02982_jprb, 0.02984_jprb, 0.02985_jprb, 0.02987_jprb, 0.02988_jprb, 0.02990_jprb, 0.02991_jprb, &
    0.02992_jprb, 0.02993_jprb, 0.02993_jprb, 0.02994_jprb, 0.02995_jprb, 0.02996_jprb, 0.02996_jprb, 0.02997_jprb, 0.02998_jprb, &
    0.02999_jprb, 0.03000_jprb, 0.03000_jprb, 0.03001_jprb, 0.03002_jprb, 0.03003_jprb, 0.03003_jprb, 0.03004_jprb, 0.03005_jprb, &
    0.03006_jprb, 0.03007_jprb, 0.03008_jprb, 0.03008_jprb, 0.03009_jprb, 0.03010_jprb, 0.03011_jprb, 0.03012_jprb, 0.03013_jprb, &
    0.03013_jprb, 0.03014_jprb, 0.03015_jprb, 0.03016_jprb, 0.03017_jprb, 0.03018_jprb, 0.03019_jprb, 0.03020_jprb, 0.03021_jprb, &
    0.03023_jprb, 0.03024_jprb, 0.03025_jprb, 0.03027_jprb, 0.03028_jprb, 0.03030_jprb, 0.03031_jprb, 0.03032_jprb, 0.03034_jprb, &
    0.03036_jprb, 0.03038_jprb, 0.03040_jprb, 0.03042_jprb, 0.03044_jprb, 0.03046_jprb, 0.03048_jprb, 0.03050_jprb, 0.03051_jprb, &
    0.03053_jprb, 0.03055_jprb, 0.03057_jprb, 0.03059_jprb, 0.03061_jprb, 0.03063_jprb, 0.03065_jprb, 0.03067_jprb, 0.03069_jprb, &
    0.03071_jprb, 0.03073_jprb, 0.03074_jprb, 0.03076_jprb, 0.03078_jprb, 0.03080_jprb, 0.03081_jprb, 0.03083_jprb, 0.03084_jprb, &
    0.03085_jprb, 0.03087_jprb, 0.03088_jprb, 0.03090_jprb, 0.03091_jprb, 0.03092_jprb, 0.03094_jprb, 0.03095_jprb, 0.03097_jprb, &
    0.03098_jprb, 0.03099_jprb, 0.03100_jprb, 0.03101_jprb, 0.03102_jprb, 0.03103_jprb, 0.03104_jprb, 0.03106_jprb, 0.03107_jprb, &
    0.03108_jprb, 0.03109_jprb, 0.03110_jprb, 0.03112_jprb, 0.03113_jprb, 0.03114_jprb, 0.03116_jprb, 0.03117_jprb, 0.03118_jprb, &
    0.03120_jprb, 0.03121_jprb, 0.03122_jprb, 0.03124_jprb, 0.03125_jprb, 0.03126_jprb, 0.03127_jprb, 0.03128_jprb, 0.03129_jprb, &
    0.03129_jprb, 0.03130_jprb, 0.03131_jprb, 0.03132_jprb /

  DATA ocean_waters_ref / &
    0.01199_jprb, 0.01232_jprb, 0.01261_jprb, 0.01290_jprb, 0.01319_jprb, 0.01348_jprb, 0.01371_jprb, 0.01388_jprb, 0.01405_jprb, &
    0.01421_jprb, 0.01436_jprb, 0.01443_jprb, 0.01449_jprb, 0.01455_jprb, 0.01460_jprb, 0.01464_jprb, 0.01468_jprb, 0.01472_jprb, &
    0.01478_jprb, 0.01484_jprb, 0.01490_jprb, 0.01497_jprb, 0.01504_jprb, 0.01511_jprb, 0.01518_jprb, 0.01524_jprb, 0.01530_jprb, &
    0.01535_jprb, 0.01539_jprb, 0.01543_jprb, 0.01547_jprb, 0.01551_jprb, 0.01555_jprb, 0.01558_jprb, 0.01562_jprb, 0.01566_jprb, &
    0.01570_jprb, 0.01575_jprb, 0.01579_jprb, 0.01584_jprb, 0.01588_jprb, 0.01592_jprb, 0.01596_jprb, 0.01599_jprb, 0.01602_jprb, &
    0.01605_jprb, 0.01607_jprb, 0.01609_jprb, 0.01611_jprb, 0.01613_jprb, 0.01615_jprb, 0.01618_jprb, 0.01621_jprb, 0.01626_jprb, &
    0.01630_jprb, 0.01635_jprb, 0.01638_jprb, 0.01642_jprb, 0.01645_jprb, 0.01647_jprb, 0.01650_jprb, 0.01652_jprb, 0.01654_jprb, &
    0.01657_jprb, 0.01659_jprb, 0.01662_jprb, 0.01665_jprb, 0.01668_jprb, 0.01672_jprb, 0.01675_jprb, 0.01679_jprb, 0.01683_jprb, &
    0.01687_jprb, 0.01691_jprb, 0.01694_jprb, 0.01698_jprb, 0.01701_jprb, 0.01703_jprb, 0.01706_jprb, 0.01708_jprb, 0.01710_jprb, &
    0.01712_jprb, 0.01714_jprb, 0.01715_jprb, 0.01717_jprb, 0.01718_jprb, 0.01719_jprb, 0.01720_jprb, 0.01721_jprb, 0.01722_jprb, &
    0.01724_jprb, 0.01725_jprb, 0.01726_jprb, 0.01727_jprb, 0.01729_jprb, 0.01730_jprb, 0.01732_jprb, 0.01734_jprb, 0.01736_jprb, &
    0.01738_jprb, 0.01741_jprb, 0.01744_jprb, 0.01747_jprb, 0.01750_jprb, 0.01753_jprb, 0.01756_jprb, 0.01759_jprb, 0.01761_jprb, &
    0.01764_jprb, 0.01766_jprb, 0.01768_jprb, 0.01769_jprb, 0.01770_jprb, 0.01771_jprb, 0.01771_jprb, 0.01771_jprb, 0.01771_jprb, &
    0.01772_jprb, 0.01773_jprb, 0.01775_jprb, 0.01776_jprb, 0.01779_jprb, 0.01781_jprb, 0.01783_jprb, 0.01785_jprb, 0.01788_jprb, &
    0.01790_jprb, 0.01792_jprb, 0.01794_jprb, 0.01796_jprb, 0.01798_jprb, 0.01800_jprb, 0.01801_jprb, 0.01803_jprb, 0.01804_jprb, &
    0.01805_jprb, 0.01807_jprb, 0.01808_jprb, 0.01809_jprb, 0.01810_jprb, 0.01811_jprb, 0.01811_jprb, 0.01812_jprb, 0.01813_jprb, &
    0.01813_jprb, 0.01814_jprb, 0.01815_jprb, 0.01815_jprb, 0.01816_jprb, 0.01816_jprb, 0.01817_jprb, 0.01818_jprb, 0.01818_jprb, &
    0.01819_jprb, 0.01820_jprb, 0.01821_jprb, 0.01822_jprb, 0.01823_jprb, 0.01824_jprb, 0.01825_jprb, 0.01826_jprb, 0.01827_jprb, &
    0.01828_jprb, 0.01830_jprb, 0.01831_jprb, 0.01833_jprb, 0.01834_jprb, 0.01835_jprb, 0.01837_jprb, 0.01838_jprb, 0.01839_jprb, &
    0.01841_jprb, 0.01842_jprb, 0.01843_jprb, 0.01845_jprb, 0.01846_jprb, 0.01847_jprb, 0.01848_jprb, 0.01849_jprb, 0.01850_jprb, &
    0.01851_jprb, 0.01851_jprb, 0.01852_jprb, 0.01852_jprb, 0.01852_jprb, 0.01852_jprb, 0.01852_jprb, 0.01852_jprb, 0.01852_jprb, &
    0.01851_jprb, 0.01850_jprb, 0.01849_jprb, 0.01848_jprb, 0.01848_jprb, 0.01847_jprb, 0.01846_jprb, 0.01846_jprb, 0.01846_jprb, &
    0.01846_jprb, 0.01846_jprb, 0.01848_jprb, 0.01849_jprb, 0.01851_jprb, 0.01853_jprb, 0.01856_jprb, 0.01859_jprb, 0.01861_jprb, &
    0.01864_jprb, 0.01865_jprb, 0.01867_jprb, 0.01868_jprb, 0.01868_jprb, 0.01869_jprb, 0.01869_jprb, 0.01869_jprb, 0.01868_jprb, &
    0.01868_jprb, 0.01867_jprb, 0.01867_jprb, 0.01866_jprb, 0.01865_jprb, 0.01864_jprb, 0.01863_jprb, 0.01863_jprb, 0.01862_jprb, &
    0.01862_jprb, 0.01861_jprb, 0.01861_jprb, 0.01861_jprb, 0.01861_jprb, 0.01862_jprb, 0.01862_jprb, 0.01862_jprb, 0.01863_jprb, &
    0.01863_jprb, 0.01864_jprb, 0.01865_jprb, 0.01865_jprb, 0.01866_jprb, 0.01867_jprb, 0.01868_jprb, 0.01869_jprb, 0.01869_jprb, &
    0.01870_jprb, 0.01871_jprb, 0.01872_jprb, 0.01872_jprb, 0.01873_jprb, 0.01874_jprb, 0.01874_jprb, 0.01875_jprb, 0.01875_jprb, &
    0.01876_jprb, 0.01876_jprb, 0.01877_jprb, 0.01877_jprb, 0.01878_jprb, 0.01878_jprb, 0.01879_jprb, 0.01879_jprb, 0.01880_jprb, &
    0.01880_jprb, 0.01880_jprb, 0.01881_jprb, 0.01881_jprb, 0.01881_jprb, 0.01882_jprb, 0.01882_jprb, 0.01883_jprb, 0.01883_jprb, &
    0.01883_jprb, 0.01884_jprb, 0.01884_jprb, 0.01885_jprb, 0.01885_jprb, 0.01885_jprb, 0.01886_jprb, 0.01886_jprb, 0.01887_jprb, &
    0.01887_jprb, 0.01888_jprb, 0.01888_jprb, 0.01888_jprb, 0.01889_jprb, 0.01889_jprb, 0.01890_jprb, 0.01890_jprb, 0.01891_jprb, &
    0.01891_jprb, 0.01892_jprb, 0.01892_jprb, 0.01893_jprb, 0.01893_jprb, 0.01894_jprb, 0.01895_jprb, 0.01895_jprb, 0.01896_jprb, &
    0.01896_jprb, 0.01897_jprb, 0.01897_jprb, 0.01898_jprb, 0.01898_jprb, 0.01898_jprb, 0.01899_jprb, 0.01899_jprb, 0.01900_jprb, &
    0.01900_jprb, 0.01901_jprb, 0.01901_jprb, 0.01901_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, &
    0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, &
    0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, &
    0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01902_jprb, 0.01903_jprb, &
    0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01903_jprb, 0.01904_jprb, 0.01904_jprb, 0.01904_jprb, 0.01904_jprb, &
    0.01905_jprb, 0.01905_jprb, 0.01905_jprb, 0.01906_jprb, 0.01906_jprb, 0.01906_jprb, 0.01907_jprb, 0.01907_jprb, 0.01908_jprb, &
    0.01908_jprb, 0.01909_jprb, 0.01909_jprb, 0.01910_jprb, 0.01910_jprb, 0.01910_jprb, 0.01911_jprb, 0.01912_jprb, 0.01912_jprb, &
    0.01913_jprb, 0.01913_jprb, 0.01914_jprb, 0.01914_jprb, 0.01915_jprb, 0.01915_jprb, 0.01916_jprb, 0.01917_jprb, 0.01917_jprb, &
    0.01918_jprb, 0.01918_jprb, 0.01919_jprb, 0.01919_jprb, 0.01920_jprb, 0.01920_jprb, 0.01920_jprb, 0.01921_jprb, 0.01921_jprb, &
    0.01921_jprb, 0.01922_jprb, 0.01922_jprb, 0.01922_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, &
    0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01923_jprb, 0.01924_jprb, &
    0.01924_jprb, 0.01924_jprb, 0.01924_jprb, 0.01925_jprb, 0.01925_jprb, 0.01925_jprb, 0.01926_jprb, 0.01926_jprb, 0.01927_jprb, &
    0.01927_jprb, 0.01928_jprb, 0.01928_jprb, 0.01928_jprb, 0.01929_jprb, 0.01929_jprb, 0.01930_jprb, 0.01930_jprb, 0.01930_jprb, &
    0.01931_jprb, 0.01931_jprb, 0.01932_jprb, 0.01932_jprb, 0.01932_jprb, 0.01932_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, &
    0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01933_jprb, &
    0.01933_jprb, 0.01933_jprb, 0.01933_jprb, 0.01934_jprb, 0.01934_jprb, 0.01934_jprb, 0.01934_jprb, 0.01934_jprb, 0.01934_jprb, &
    0.01934_jprb, 0.01935_jprb, 0.01935_jprb, 0.01935_jprb, 0.01935_jprb, 0.01935_jprb, 0.01936_jprb, 0.01936_jprb, 0.01936_jprb, &
    0.01936_jprb, 0.01937_jprb, 0.01937_jprb, 0.01937_jprb, 0.01937_jprb, 0.01937_jprb, 0.01938_jprb, 0.01938_jprb, 0.01938_jprb, &
    0.01938_jprb, 0.01938_jprb, 0.01939_jprb, 0.01939_jprb, 0.01939_jprb, 0.01939_jprb, 0.01939_jprb, 0.01939_jprb, 0.01939_jprb, &
    0.01939_jprb, 0.01939_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, &
    0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01940_jprb, 0.01941_jprb, 0.01941_jprb, 0.01941_jprb, &
    0.01941_jprb, 0.01941_jprb, 0.01941_jprb, 0.01941_jprb, 0.01941_jprb, 0.01942_jprb, 0.01942_jprb, 0.01942_jprb, 0.01942_jprb, &
    0.01942_jprb, 0.01943_jprb, 0.01943_jprb, 0.01943_jprb, 0.01943_jprb, 0.01944_jprb, 0.01944_jprb, 0.01944_jprb, 0.01945_jprb, &
    0.01945_jprb, 0.01945_jprb, 0.01946_jprb, 0.01946_jprb, 0.01946_jprb, 0.01947_jprb, 0.01947_jprb, 0.01947_jprb, 0.01948_jprb, &
    0.01948_jprb, 0.01949_jprb, 0.01949_jprb, 0.01949_jprb, 0.01950_jprb, 0.01950_jprb, 0.01950_jprb, 0.01951_jprb, 0.01951_jprb, &
    0.01951_jprb, 0.01952_jprb, 0.01952_jprb, 0.01952_jprb, 0.01953_jprb, 0.01953_jprb, 0.01953_jprb, 0.01953_jprb, 0.01953_jprb, &
    0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, &
    0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, &
    0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01954_jprb, 0.01955_jprb, 0.01955_jprb, 0.01955_jprb, 0.01955_jprb, &
    0.01955_jprb, 0.01955_jprb, 0.01955_jprb, 0.01956_jprb, 0.01956_jprb, 0.01956_jprb, 0.01956_jprb, 0.01956_jprb, 0.01956_jprb, &
    0.01957_jprb, 0.01957_jprb, 0.01957_jprb, 0.01957_jprb, 0.01958_jprb, 0.01958_jprb, 0.01958_jprb, 0.01958_jprb, 0.01959_jprb, &
    0.01959_jprb, 0.01959_jprb, 0.01960_jprb, 0.01960_jprb, 0.01960_jprb, 0.01961_jprb, 0.01961_jprb, 0.01961_jprb, 0.01962_jprb, &
    0.01962_jprb, 0.01962_jprb, 0.01963_jprb, 0.01963_jprb, 0.01964_jprb, 0.01964_jprb, 0.01964_jprb, 0.01965_jprb, 0.01965_jprb, &
    0.01966_jprb, 0.01966_jprb, 0.01966_jprb, 0.01967_jprb, 0.01967_jprb, 0.01968_jprb, 0.01968_jprb, 0.01969_jprb, 0.01969_jprb, &
    0.01969_jprb, 0.01970_jprb, 0.01970_jprb, 0.01971_jprb, 0.01971_jprb, 0.01971_jprb, 0.01972_jprb, 0.01972_jprb, 0.01972_jprb, &
    0.01973_jprb, 0.01973_jprb, 0.01973_jprb, 0.01973_jprb, 0.01974_jprb, 0.01974_jprb, 0.01974_jprb, 0.01974_jprb, 0.01974_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, &
    0.01977_jprb, 0.01977_jprb, 0.01977_jprb, 0.01977_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01979_jprb, &
    0.01979_jprb, 0.01979_jprb, 0.01979_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01981_jprb, 0.01981_jprb, &
    0.01981_jprb, 0.01981_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01983_jprb, 0.01983_jprb, &
    0.01983_jprb, 0.01983_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, &
    0.01984_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, &
    0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01984_jprb, &
    0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01983_jprb, 0.01983_jprb, 0.01983_jprb, 0.01983_jprb, &
    0.01983_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01981_jprb, 0.01981_jprb, 0.01981_jprb, &
    0.01981_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01979_jprb, 0.01979_jprb, 0.01979_jprb, &
    0.01979_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01977_jprb, 0.01977_jprb, 0.01977_jprb, &
    0.01977_jprb, 0.01977_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01975_jprb, &
    0.01975_jprb, 0.01975_jprb, 0.01975_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, 0.01976_jprb, &
    0.01977_jprb, 0.01977_jprb, 0.01977_jprb, 0.01977_jprb, 0.01977_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, 0.01978_jprb, &
    0.01978_jprb, 0.01979_jprb, 0.01979_jprb, 0.01979_jprb, 0.01979_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, 0.01980_jprb, &
    0.01981_jprb, 0.01981_jprb, 0.01981_jprb, 0.01981_jprb, 0.01981_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, 0.01982_jprb, &
    0.01982_jprb, 0.01983_jprb, 0.01983_jprb, 0.01983_jprb, 0.01983_jprb, 0.01983_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, &
    0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01984_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, &
    0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, &
    0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, 0.01985_jprb, &
    0.01985_jprb, 0.01985_jprb, 0.01986_jprb, 0.01986_jprb, 0.01986_jprb, 0.01986_jprb, 0.01986_jprb, 0.01986_jprb, 0.01986_jprb, &
    0.01987_jprb, 0.01987_jprb, 0.01987_jprb, 0.01987_jprb, 0.01987_jprb, 0.01987_jprb, 0.01988_jprb, 0.01988_jprb, 0.01988_jprb, &
    0.01988_jprb, 0.01988_jprb, 0.01989_jprb, 0.01989_jprb, 0.01989_jprb, 0.01989_jprb, 0.01989_jprb, 0.01990_jprb, 0.01990_jprb, &
    0.01990_jprb, 0.01990_jprb, 0.01991_jprb, 0.01991_jprb, 0.01991_jprb, 0.01991_jprb, 0.01991_jprb, 0.01992_jprb, 0.01992_jprb, &
    0.01992_jprb, 0.01992_jprb, 0.01992_jprb, 0.01993_jprb, 0.01993_jprb, 0.01993_jprb, 0.01993_jprb, 0.01993_jprb, 0.01994_jprb, &
    0.01994_jprb, 0.01994_jprb, 0.01994_jprb, 0.01994_jprb, 0.01994_jprb, 0.01994_jprb, 0.01995_jprb, 0.01995_jprb, 0.01995_jprb, &
    0.01995_jprb, 0.01995_jprb, 0.01995_jprb, 0.01995_jprb, 0.01995_jprb, 0.01995_jprb, 0.01995_jprb, 0.01995_jprb, 0.01996_jprb, &
    0.01996_jprb, 0.01996_jprb, 0.01996_jprb, 0.01996_jprb, 0.01996_jprb, 0.01996_jprb, 0.01996_jprb, 0.01996_jprb, 0.01996_jprb, &
    0.01996_jprb, 0.01996_jprb, 0.01995_jprb, 0.01995_jprb, 0.01995_jprb, 0.01995_jprb, 0.01996_jprb, 0.01997_jprb, 0.01998_jprb, &
    0.01999_jprb, 0.02000_jprb, 0.02000_jprb, 0.02000_jprb, 0.02000_jprb, 0.02000_jprb, 0.02000_jprb, 0.02001_jprb, 0.02001_jprb, &
    0.02001_jprb, 0.02001_jprb, 0.02001_jprb, 0.02001_jprb, 0.02001_jprb, 0.02001_jprb, 0.02001_jprb, 0.02001_jprb, 0.02001_jprb, &
    0.02002_jprb, 0.02002_jprb, 0.02002_jprb, 0.02002_jprb, 0.02002_jprb, 0.02002_jprb, 0.02002_jprb, 0.02003_jprb, 0.02003_jprb, &
    0.02003_jprb, 0.02003_jprb, 0.02003_jprb, 0.02003_jprb, 0.02004_jprb, 0.02004_jprb, 0.02005_jprb, 0.02006_jprb, 0.02006_jprb, &
    0.02007_jprb, 0.02007_jprb, 0.02008_jprb, 0.02009_jprb, 0.02009_jprb, 0.02010_jprb, 0.02011_jprb, 0.02011_jprb, 0.02012_jprb, &
    0.02013_jprb, 0.02013_jprb, 0.02014_jprb, 0.02014_jprb, 0.02015_jprb, 0.02016_jprb, 0.02016_jprb, 0.02016_jprb, 0.02017_jprb, &
    0.02018_jprb, 0.02018_jprb, 0.02019_jprb, 0.02019_jprb, 0.02020_jprb, 0.02021_jprb, 0.02021_jprb, 0.02022_jprb, 0.02022_jprb, &
    0.02022_jprb, 0.02023_jprb, 0.02023_jprb, 0.02023_jprb, 0.02024_jprb, 0.02024_jprb, 0.02025_jprb, 0.02025_jprb, 0.02025_jprb, &
    0.02026_jprb, 0.02026_jprb, 0.02027_jprb, 0.02027_jprb, 0.02027_jprb, 0.02028_jprb, 0.02028_jprb, 0.02029_jprb, 0.02029_jprb, &
    0.02030_jprb, 0.02030_jprb, 0.02031_jprb, 0.02031_jprb, 0.02032_jprb, 0.02032_jprb, 0.02032_jprb, 0.02033_jprb, 0.02033_jprb, &
    0.02033_jprb, 0.02034_jprb, 0.02034_jprb, 0.02035_jprb, 0.02035_jprb, 0.02036_jprb, 0.02037_jprb, 0.02037_jprb, 0.02038_jprb, &
    0.02038_jprb, 0.02039_jprb, 0.02040_jprb, 0.02040_jprb, 0.02041_jprb, 0.02041_jprb, 0.02042_jprb, 0.02042_jprb, 0.02042_jprb, &
    0.02043_jprb, 0.02043_jprb, 0.02044_jprb, 0.02044_jprb, 0.02044_jprb, 0.02045_jprb, 0.02045_jprb, 0.02046_jprb, 0.02046_jprb, &
    0.02047_jprb, 0.02047_jprb, 0.02047_jprb, 0.02047_jprb, 0.02047_jprb, 0.02048_jprb, 0.02048_jprb, 0.02048_jprb, 0.02049_jprb, &
    0.02049_jprb, 0.02049_jprb, 0.02050_jprb, 0.02050_jprb, 0.02051_jprb, 0.02051_jprb, 0.02052_jprb, 0.02052_jprb, 0.02052_jprb, &
    0.02053_jprb, 0.02053_jprb, 0.02054_jprb, 0.02054_jprb, 0.02054_jprb, 0.02055_jprb, 0.02055_jprb, 0.02055_jprb, 0.02056_jprb, &
    0.02056_jprb, 0.02056_jprb, 0.02057_jprb, 0.02057_jprb, 0.02057_jprb, 0.02058_jprb, 0.02058_jprb, 0.02058_jprb, 0.02059_jprb, &
    0.02059_jprb, 0.02059_jprb, 0.02060_jprb, 0.02060_jprb, 0.02060_jprb, 0.02060_jprb, 0.02060_jprb, 0.02060_jprb, 0.02060_jprb, &
    0.02061_jprb, 0.02061_jprb, 0.02061_jprb, 0.02061_jprb, 0.02061_jprb, 0.02061_jprb, 0.02061_jprb, 0.02061_jprb, 0.02062_jprb, &
    0.02062_jprb, 0.02062_jprb, 0.02063_jprb, 0.02063_jprb, 0.02063_jprb, 0.02063_jprb, 0.02064_jprb, 0.02064_jprb, 0.02065_jprb, &
    0.02065_jprb, 0.02065_jprb, 0.02065_jprb, 0.02065_jprb, 0.02065_jprb, 0.02065_jprb, 0.02066_jprb, 0.02066_jprb, 0.02066_jprb, &
    0.02066_jprb, 0.02067_jprb, 0.02068_jprb, 0.02068_jprb, 0.02069_jprb, 0.02070_jprb, 0.02070_jprb, 0.02071_jprb, 0.02071_jprb, &
    0.02072_jprb, 0.02072_jprb, 0.02073_jprb, 0.02073_jprb, 0.02073_jprb, 0.02074_jprb, 0.02074_jprb, 0.02075_jprb, 0.02075_jprb, &
    0.02076_jprb, 0.02077_jprb, 0.02077_jprb, 0.02077_jprb, 0.02078_jprb, 0.02078_jprb, 0.02078_jprb, 0.02078_jprb, 0.02078_jprb, &
    0.02079_jprb, 0.02079_jprb, 0.02079_jprb, 0.02079_jprb, 0.02080_jprb, 0.02080_jprb, 0.02080_jprb, 0.02081_jprb, 0.02081_jprb, &
    0.02081_jprb, 0.02082_jprb, 0.02082_jprb, 0.02082_jprb, 0.02082_jprb, 0.02082_jprb, 0.02082_jprb, 0.02082_jprb, 0.02083_jprb, &
    0.02083_jprb, 0.02083_jprb, 0.02083_jprb, 0.02084_jprb, 0.02084_jprb, 0.02084_jprb, 0.02085_jprb, 0.02085_jprb, 0.02085_jprb, &
    0.02086_jprb, 0.02087_jprb, 0.02087_jprb, 0.02088_jprb, 0.02088_jprb, 0.02088_jprb, 0.02089_jprb, 0.02089_jprb, 0.02089_jprb, &
    0.02090_jprb, 0.02090_jprb, 0.02090_jprb, 0.02090_jprb, 0.02091_jprb, 0.02091_jprb, 0.02091_jprb, 0.02092_jprb, 0.02092_jprb, &
    0.02093_jprb, 0.02093_jprb, 0.02093_jprb, 0.02094_jprb, 0.02094_jprb, 0.02094_jprb, 0.02095_jprb, 0.02095_jprb, 0.02096_jprb, &
    0.02096_jprb, 0.02096_jprb, 0.02097_jprb, 0.02097_jprb, 0.02098_jprb, 0.02098_jprb, 0.02098_jprb, 0.02099_jprb, 0.02099_jprb, &
    0.02099_jprb, 0.02100_jprb, 0.02100_jprb, 0.02101_jprb, 0.02101_jprb, 0.02102_jprb, 0.02102_jprb, 0.02102_jprb, 0.02103_jprb, &
    0.02103_jprb, 0.02104_jprb, 0.02105_jprb, 0.02105_jprb, 0.02106_jprb, 0.02107_jprb, 0.02107_jprb, 0.02108_jprb, 0.02108_jprb, &
    0.02109_jprb, 0.02109_jprb, 0.02109_jprb, 0.02110_jprb, 0.02110_jprb, 0.02110_jprb, 0.02111_jprb, 0.02112_jprb, 0.02112_jprb, &
    0.02113_jprb, 0.02114_jprb, 0.02115_jprb, 0.02117_jprb, 0.02118_jprb, 0.02119_jprb, 0.02120_jprb, 0.02121_jprb, 0.02122_jprb, &
    0.02123_jprb, 0.02124_jprb, 0.02125_jprb, 0.02126_jprb, 0.02127_jprb, 0.02128_jprb, 0.02129_jprb, 0.02131_jprb, 0.02132_jprb, &
    0.02134_jprb, 0.02135_jprb, 0.02137_jprb, 0.02139_jprb, 0.02140_jprb, 0.02142_jprb, 0.02144_jprb, 0.02146_jprb, 0.02148_jprb, &
    0.02150_jprb, 0.02152_jprb, 0.02154_jprb, 0.02156_jprb, 0.02158_jprb, 0.02160_jprb, 0.02163_jprb, 0.02165_jprb, 0.02168_jprb, &
    0.02170_jprb, 0.02173_jprb, 0.02175_jprb, 0.02178_jprb, 0.02180_jprb, 0.02183_jprb, 0.02186_jprb, 0.02188_jprb, 0.02191_jprb, &
    0.02194_jprb, 0.02197_jprb, 0.02200_jprb, 0.02202_jprb, 0.02205_jprb, 0.02208_jprb, 0.02212_jprb, 0.02215_jprb, 0.02218_jprb, &
    0.02221_jprb, 0.02224_jprb, 0.02228_jprb, 0.02231_jprb, 0.02234_jprb, 0.02238_jprb, 0.02241_jprb, 0.02244_jprb, 0.02247_jprb, &
    0.02251_jprb, 0.02254_jprb, 0.02257_jprb, 0.02260_jprb, 0.02264_jprb, 0.02267_jprb, 0.02270_jprb, 0.02273_jprb, 0.02276_jprb, &
    0.02279_jprb, 0.02283_jprb, 0.02286_jprb, 0.02289_jprb, 0.02293_jprb, 0.02296_jprb, 0.02300_jprb, 0.02303_jprb, 0.02307_jprb, &
    0.02310_jprb, 0.02314_jprb, 0.02318_jprb, 0.02321_jprb, 0.02325_jprb, 0.02328_jprb, 0.02332_jprb, 0.02336_jprb, 0.02339_jprb, &
    0.02343_jprb, 0.02346_jprb, 0.02349_jprb, 0.02353_jprb, 0.02356_jprb, 0.02359_jprb, 0.02363_jprb, 0.02366_jprb, 0.02369_jprb, &
    0.02372_jprb, 0.02376_jprb, 0.02379_jprb, 0.02382_jprb, 0.02385_jprb, 0.02389_jprb, 0.02392_jprb, 0.02395_jprb, 0.02399_jprb, &
    0.02402_jprb, 0.02405_jprb, 0.02408_jprb, 0.02411_jprb, 0.02414_jprb, 0.02417_jprb, 0.02420_jprb, 0.02423_jprb, 0.02426_jprb, &
    0.02429_jprb, 0.02431_jprb, 0.02434_jprb, 0.02437_jprb, 0.02440_jprb, 0.02443_jprb, 0.02446_jprb, 0.02449_jprb, 0.02452_jprb, &
    0.02455_jprb, 0.02457_jprb, 0.02460_jprb, 0.02463_jprb, 0.02465_jprb, 0.02468_jprb, 0.02470_jprb, 0.02473_jprb, 0.02475_jprb, &
    0.02477_jprb, 0.02480_jprb, 0.02482_jprb, 0.02484_jprb, 0.02486_jprb, 0.02489_jprb, 0.02491_jprb, 0.02494_jprb, 0.02496_jprb, &
    0.02499_jprb, 0.02501_jprb, 0.02504_jprb, 0.02506_jprb, 0.02509_jprb, 0.02511_jprb, 0.02514_jprb, 0.02516_jprb, 0.02519_jprb, &
    0.02521_jprb, 0.02524_jprb, 0.02527_jprb, 0.02529_jprb, 0.02532_jprb, 0.02535_jprb, 0.02537_jprb, 0.02540_jprb, 0.02543_jprb, &
    0.02546_jprb, 0.02549_jprb, 0.02552_jprb, 0.02555_jprb, 0.02558_jprb, 0.02561_jprb, 0.02564_jprb, 0.02567_jprb, 0.02570_jprb, &
    0.02573_jprb, 0.02576_jprb, 0.02579_jprb, 0.02582_jprb, 0.02585_jprb, 0.02588_jprb, 0.02591_jprb, 0.02594_jprb, 0.02596_jprb, &
    0.02599_jprb, 0.02602_jprb, 0.02605_jprb, 0.02608_jprb, 0.02610_jprb, 0.02613_jprb, 0.02616_jprb, 0.02619_jprb, 0.02623_jprb, &
    0.02626_jprb, 0.02629_jprb, 0.02632_jprb, 0.02635_jprb, 0.02639_jprb, 0.02642_jprb, 0.02645_jprb, 0.02649_jprb, 0.02652_jprb, &
    0.02656_jprb, 0.02659_jprb, 0.02663_jprb, 0.02666_jprb, 0.02669_jprb, 0.02673_jprb, 0.02676_jprb, 0.02680_jprb, 0.02683_jprb, &
    0.02686_jprb, 0.02689_jprb, 0.02691_jprb, 0.02694_jprb, 0.02697_jprb, 0.02700_jprb, 0.02703_jprb, 0.02706_jprb, 0.02709_jprb, &
    0.02711_jprb, 0.02714_jprb, 0.02717_jprb, 0.02720_jprb, 0.02723_jprb, 0.02726_jprb, 0.02729_jprb, 0.02732_jprb, 0.02735_jprb, &
    0.02738_jprb, 0.02741_jprb, 0.02745_jprb, 0.02748_jprb, 0.02751_jprb, 0.02755_jprb, 0.02758_jprb, 0.02761_jprb, 0.02764_jprb, &
    0.02767_jprb, 0.02771_jprb, 0.02774_jprb, 0.02777_jprb, 0.02780_jprb, 0.02783_jprb, 0.02786_jprb, 0.02789_jprb, 0.02791_jprb, &
    0.02794_jprb, 0.02797_jprb, 0.02800_jprb, 0.02803_jprb, 0.02806_jprb, 0.02808_jprb, 0.02811_jprb, 0.02814_jprb, 0.02817_jprb, &
    0.02820_jprb, 0.02823_jprb, 0.02826_jprb, 0.02829_jprb, 0.02831_jprb, 0.02834_jprb, 0.02837_jprb, 0.02840_jprb, 0.02843_jprb, &
    0.02846_jprb, 0.02849_jprb, 0.02851_jprb, 0.02854_jprb, 0.02857_jprb, 0.02860_jprb, 0.02863_jprb, 0.02866_jprb, 0.02869_jprb, &
    0.02872_jprb, 0.02875_jprb, 0.02878_jprb, 0.02881_jprb, 0.02884_jprb, 0.02887_jprb, 0.02890_jprb, 0.02893_jprb, 0.02896_jprb, &
    0.02900_jprb, 0.02903_jprb, 0.02906_jprb, 0.02909_jprb, 0.02912_jprb, 0.02915_jprb, 0.02919_jprb, 0.02923_jprb, 0.02927_jprb, &
    0.02930_jprb, 0.02934_jprb, 0.02938_jprb, 0.02942_jprb, 0.02946_jprb, 0.02951_jprb, 0.02956_jprb, 0.02961_jprb, 0.02966_jprb, &
    0.02971_jprb, 0.02976_jprb, 0.02981_jprb, 0.02986_jprb, 0.02992_jprb, 0.02998_jprb, 0.03003_jprb, 0.03009_jprb, 0.03014_jprb, &
    0.03020_jprb, 0.03026_jprb, 0.03033_jprb, 0.03040_jprb, 0.03046_jprb, 0.03053_jprb, 0.03059_jprb, 0.03066_jprb, 0.03072_jprb, &
    0.03080_jprb, 0.03088_jprb, 0.03095_jprb, 0.03103_jprb, 0.03111_jprb, 0.03119_jprb, 0.03126_jprb, 0.03134_jprb, 0.03143_jprb, &
    0.03152_jprb, 0.03160_jprb, 0.03169_jprb, 0.03177_jprb, 0.03186_jprb, 0.03194_jprb, 0.03203_jprb, 0.03212_jprb, 0.03221_jprb, &
    0.03230_jprb, 0.03239_jprb, 0.03248_jprb, 0.03257_jprb, 0.03266_jprb, 0.03276_jprb, 0.03285_jprb, 0.03295_jprb, 0.03305_jprb, &
    0.03314_jprb, 0.03324_jprb, 0.03334_jprb, 0.03344_jprb, 0.03354_jprb, 0.03365_jprb, 0.03375_jprb, 0.03386_jprb, 0.03396_jprb, &
    0.03407_jprb, 0.03418_jprb, 0.03428_jprb, 0.03438_jprb, 0.03448_jprb, 0.03458_jprb, 0.03468_jprb, 0.03479_jprb, 0.03489_jprb, &
    0.03499_jprb, 0.03509_jprb, 0.03518_jprb, 0.03528_jprb, 0.03538_jprb, 0.03547_jprb, 0.03557_jprb, 0.03566_jprb, 0.03576_jprb, &
    0.03586_jprb, 0.03595_jprb, 0.03603_jprb, 0.03612_jprb, 0.03621_jprb, 0.03630_jprb, 0.03639_jprb, 0.03648_jprb, 0.03657_jprb, &
    0.03666_jprb, 0.03675_jprb, 0.03683_jprb, 0.03692_jprb, 0.03701_jprb, 0.03709_jprb, 0.03718_jprb, 0.03727_jprb, 0.03735_jprb, &
    0.03744_jprb, 0.03752_jprb, 0.03760_jprb, 0.03768_jprb, 0.03777_jprb, 0.03785_jprb, 0.03793_jprb, 0.03801_jprb, 0.03809_jprb, &
    0.03817_jprb, 0.03824_jprb, 0.03832_jprb, 0.03840_jprb, 0.03848_jprb, 0.03855_jprb, 0.03863_jprb, 0.03870_jprb, 0.03878_jprb, &
    0.03885_jprb, 0.03892_jprb, 0.03899_jprb, 0.03907_jprb, 0.03914_jprb, 0.03921_jprb, 0.03928_jprb, 0.03935_jprb, 0.03942_jprb, &
    0.03948_jprb, 0.03955_jprb, 0.03962_jprb, 0.03968_jprb, 0.03975_jprb, 0.03982_jprb, 0.03987_jprb, 0.03993_jprb, 0.03999_jprb, &
    0.04005_jprb, 0.04011_jprb, 0.04016_jprb, 0.04022_jprb, 0.04028_jprb, 0.04033_jprb, 0.04037_jprb, 0.04042_jprb, 0.04047_jprb, &
    0.04051_jprb, 0.04056_jprb, 0.04060_jprb, 0.04065_jprb, 0.04069_jprb, 0.04073_jprb, 0.04077_jprb, 0.04081_jprb, 0.04085_jprb, &
    0.04088_jprb, 0.04092_jprb, 0.04096_jprb, 0.04100_jprb, 0.04104_jprb, 0.04107_jprb, 0.04110_jprb, 0.04114_jprb, 0.04117_jprb, &
    0.04121_jprb, 0.04124_jprb, 0.04128_jprb, 0.04131_jprb, 0.04135_jprb, 0.04138_jprb, 0.04142_jprb, 0.04145_jprb, 0.04149_jprb, &
    0.04152_jprb, 0.04156_jprb, 0.04159_jprb, 0.04163_jprb, 0.04167_jprb, 0.04170_jprb, 0.04174_jprb, 0.04178_jprb, 0.04182_jprb, &
    0.04185_jprb, 0.04189_jprb, 0.04193_jprb, 0.04197_jprb, 0.04200_jprb, 0.04203_jprb, 0.04207_jprb, 0.04210_jprb, 0.04214_jprb, &
    0.04217_jprb, 0.04221_jprb, 0.04224_jprb, 0.04228_jprb, 0.04230_jprb, 0.04233_jprb, 0.04236_jprb, 0.04239_jprb, 0.04242_jprb, &
    0.04245_jprb, 0.04248_jprb, 0.04251_jprb, 0.04254_jprb, 0.04256_jprb, 0.04259_jprb, 0.04261_jprb, 0.04263_jprb, 0.04266_jprb, &
    0.04268_jprb, 0.04270_jprb, 0.04273_jprb, 0.04275_jprb, 0.04277_jprb, 0.04279_jprb, 0.04281_jprb, 0.04282_jprb, 0.04284_jprb, &
    0.04286_jprb, 0.04288_jprb, 0.04290_jprb, 0.04292_jprb, 0.04293_jprb, 0.04294_jprb, 0.04295_jprb, 0.04296_jprb, 0.04297_jprb, &
    0.04298_jprb, 0.04299_jprb, 0.04300_jprb, 0.04301_jprb, 0.04302_jprb, 0.04303_jprb, 0.04303_jprb, 0.04304_jprb, 0.04304_jprb, &
    0.04305_jprb, 0.04305_jprb, 0.04306_jprb, 0.04307_jprb, 0.04307_jprb, 0.04309_jprb, 0.04310_jprb, 0.04311_jprb, 0.04312_jprb, &
    0.04313_jprb, 0.04315_jprb, 0.04316_jprb, 0.04317_jprb, 0.04318_jprb, 0.04320_jprb, 0.04322_jprb, 0.04324_jprb, 0.04326_jprb, &
    0.04328_jprb, 0.04330_jprb, 0.04331_jprb, 0.04333_jprb, 0.04335_jprb, 0.04337_jprb, 0.04338_jprb, 0.04340_jprb, 0.04341_jprb, &
    0.04342_jprb, 0.04344_jprb, 0.04345_jprb, 0.04346_jprb, 0.04348_jprb, 0.04349_jprb, 0.04349_jprb, 0.04350_jprb, 0.04350_jprb, &
    0.04350_jprb, 0.04350_jprb, 0.04350_jprb, 0.04350_jprb, 0.04351_jprb, 0.04351_jprb, 0.04351_jprb, 0.04352_jprb, 0.04353_jprb, &
    0.04354_jprb, 0.04355_jprb, 0.04357_jprb, 0.04358_jprb, 0.04359_jprb, 0.04360_jprb, 0.04361_jprb, 0.04363_jprb, 0.04365_jprb, &
    0.04368_jprb, 0.04371_jprb, 0.04373_jprb, 0.04376_jprb, 0.04379_jprb, 0.04381_jprb, 0.04384_jprb, 0.04387_jprb, 0.04389_jprb, &
    0.04392_jprb, 0.04395_jprb, 0.04397_jprb, 0.04400_jprb, 0.04403_jprb, 0.04406_jprb, 0.04409_jprb, 0.04412_jprb, 0.04415_jprb, &
    0.04418_jprb, 0.04421_jprb, 0.04424_jprb, 0.04427_jprb, 0.04430_jprb, 0.04433_jprb, 0.04436_jprb, 0.04439_jprb, 0.04442_jprb, &
    0.04445_jprb, 0.04449_jprb, 0.04453_jprb, 0.04456_jprb, 0.04460_jprb, 0.04463_jprb, 0.04467_jprb, 0.04470_jprb, 0.04474_jprb, &
    0.04477_jprb, 0.04481_jprb, 0.04484_jprb, 0.04488_jprb, 0.04491_jprb, 0.04495_jprb, 0.04499_jprb, 0.04503_jprb, 0.04506_jprb, &
    0.04510_jprb, 0.04514_jprb, 0.04518_jprb, 0.04521_jprb, 0.04525_jprb, 0.04529_jprb, 0.04533_jprb, 0.04537_jprb, 0.04540_jprb, &
    0.04544_jprb, 0.04548_jprb, 0.04552_jprb, 0.04555_jprb, 0.04558_jprb, 0.04561_jprb, 0.04565_jprb, 0.04568_jprb, 0.04571_jprb, &
    0.04574_jprb, 0.04578_jprb, 0.04581_jprb, 0.04584_jprb, 0.04587_jprb, 0.04591_jprb, 0.04594_jprb, 0.04597_jprb, 0.04600_jprb, &
    0.04604_jprb, 0.04607_jprb, 0.04610_jprb, 0.04614_jprb, 0.04617_jprb, 0.04620_jprb, 0.04624_jprb, 0.04627_jprb, 0.04630_jprb, &
    0.04634_jprb, 0.04637_jprb, 0.04640_jprb, 0.04644_jprb, 0.04647_jprb, 0.04650_jprb, 0.04654_jprb, 0.04657_jprb, 0.04661_jprb, &
    0.04664_jprb, 0.04667_jprb, 0.04671_jprb, 0.04674_jprb, 0.04678_jprb, 0.04681_jprb, 0.04684_jprb, 0.04688_jprb, 0.04691_jprb, &
    0.04695_jprb, 0.04698_jprb, 0.04701_jprb, 0.04705_jprb, 0.04708_jprb, 0.04711_jprb, 0.04713_jprb, 0.04715_jprb, 0.04718_jprb, &
    0.04720_jprb, 0.04722_jprb, 0.04725_jprb, 0.04727_jprb, 0.04729_jprb, 0.04732_jprb, 0.04734_jprb, 0.04736_jprb, 0.04739_jprb, &
    0.04741_jprb, 0.04743_jprb, 0.04745_jprb, 0.04747_jprb, 0.04749_jprb, 0.04751_jprb, 0.04752_jprb, 0.04754_jprb, 0.04755_jprb, &
    0.04757_jprb, 0.04759_jprb, 0.04760_jprb, 0.04762_jprb, 0.04763_jprb, 0.04765_jprb, 0.04767_jprb, 0.04768_jprb, 0.04770_jprb, &
    0.04771_jprb, 0.04773_jprb, 0.04775_jprb, 0.04776_jprb, 0.04778_jprb, 0.04780_jprb, 0.04781_jprb, 0.04783_jprb, 0.04784_jprb, &
    0.04786_jprb, 0.04788_jprb, 0.04789_jprb, 0.04791_jprb, 0.04793_jprb, 0.04794_jprb, 0.04796_jprb, 0.04798_jprb, 0.04799_jprb, &
    0.04801_jprb, 0.04803_jprb, 0.04804_jprb, 0.04806_jprb, 0.04807_jprb, 0.04809_jprb, 0.04811_jprb, 0.04812_jprb, 0.04814_jprb, &
    0.04815_jprb, 0.04817_jprb, 0.04819_jprb, 0.04820_jprb, 0.04822_jprb, 0.04824_jprb, 0.04825_jprb, 0.04827_jprb, 0.04828_jprb, &
    0.04830_jprb, 0.04832_jprb, 0.04833_jprb, 0.04835_jprb, 0.04836_jprb, 0.04838_jprb, 0.04839_jprb, 0.04841_jprb, 0.04843_jprb, &
    0.04844_jprb, 0.04846_jprb, 0.04847_jprb, 0.04849_jprb, 0.04850_jprb, 0.04852_jprb, 0.04854_jprb, 0.04855_jprb, 0.04856_jprb, &
    0.04856_jprb, 0.04857_jprb, 0.04857_jprb, 0.04858_jprb, 0.04858_jprb, 0.04859_jprb, 0.04859_jprb, 0.04859_jprb, 0.04860_jprb, &
    0.04860_jprb, 0.04861_jprb, 0.04861_jprb, 0.04862_jprb, 0.04862_jprb, 0.04863_jprb, 0.04863_jprb, 0.04864_jprb, 0.04864_jprb, &
    0.04864_jprb, 0.04864_jprb, 0.04864_jprb, 0.04864_jprb, 0.04864_jprb, 0.04864_jprb, 0.04864_jprb, 0.04864_jprb, 0.04864_jprb, &
    0.04865_jprb, 0.04866_jprb, 0.04866_jprb, 0.04867_jprb, 0.04867_jprb, 0.04868_jprb, 0.04869_jprb, 0.04869_jprb, 0.04870_jprb, &
    0.04870_jprb, 0.04871_jprb, 0.04871_jprb, 0.04872_jprb, 0.04872_jprb, 0.04872_jprb, 0.04873_jprb, 0.04873_jprb, 0.04874_jprb, &
    0.04874_jprb, 0.04874_jprb, 0.04875_jprb, 0.04875_jprb, 0.04876_jprb, 0.04875_jprb, 0.04874_jprb, 0.04874_jprb, 0.04873_jprb, &
    0.04872_jprb, 0.04871_jprb, 0.04870_jprb, 0.04869_jprb, 0.04869_jprb, 0.04868_jprb, 0.04867_jprb, 0.04866_jprb, 0.04865_jprb, &
    0.04863_jprb, 0.04862_jprb, 0.04860_jprb, 0.04859_jprb, 0.04857_jprb, 0.04856_jprb, 0.04854_jprb, 0.04853_jprb, 0.04851_jprb, &
    0.04850_jprb, 0.04848_jprb, 0.04847_jprb, 0.04844_jprb, 0.04842_jprb, 0.04840_jprb, 0.04838_jprb, 0.04836_jprb, 0.04834_jprb, &
    0.04832_jprb, 0.04830_jprb, 0.04828_jprb, 0.04826_jprb, 0.04824_jprb, 0.04822_jprb, 0.04820_jprb, 0.04818_jprb, 0.04816_jprb, &
    0.04814_jprb, 0.04812_jprb, 0.04810_jprb, 0.04808_jprb /

CONTAINS

  !> Linearly interpolates the reflectance data onto given wavenumbers
  !! @param[in]       wvn       Wavenumbers onto which reflectances are to be interpolated
  !! @param[in,out]   refl_ow   Interpolated ocean water reflectance
  !! @param[in,out]   refl_fw   Interpolated coastal water reflectance
  SUBROUTINE rttov_refl_water_interp(wvn, refl_ow, refl_fw)

    REAL(jprb), INTENT(IN)    :: wvn(:)
    REAL(jprb), INTENT(INOUT) :: refl_ow(SIZE(wvn))
    REAL(jprb), INTENT(INOUT) :: refl_fw(SIZE(wvn))

    REAL(jprb) :: wvn_waters_ref(numwave)

    INTEGER(jpim) :: nwvn, i, j
    REAL(jprb)    :: z, x1, x2

    !---------------------------------------------------------
    ! create wvn_waters_ref array (copied from mod_brdf_atlas)
    !---------------------------------------------------------
    !  numwave = 2101, hsr_wavenum(4000:25000:10)

    DO i = 1, numwave
      wvn_waters_ref(i) = 4000._jprb + (i - 1) * 10._jprb
    ENDDO


    nwvn = SIZE(wvn)
    DO i = 1, nwvn
      IF (wvn(i) <= wvn_waters_ref(1)) THEN
        refl_ow(i) = ocean_waters_ref(1)
        refl_fw(i) = coastal_waters_ref(1)
      ELSEIF (wvn(i) >= wvn_waters_ref(numwave)) THEN
        refl_ow(i) = ocean_waters_ref(numwave)
        refl_fw(i) = coastal_waters_ref(numwave)
      ELSE
        DO j = 2, numwave
          IF (wvn(i) < wvn_waters_ref(j)) EXIT
        ENDDO
        z = 1._jprb / (wvn_waters_ref(j) - wvn_waters_ref(j-1))
        x1 = (wvn_waters_ref(j) - wvn(i)) * z
        x2 = (wvn(i) - wvn_waters_ref(j-1)) * z
        refl_ow(i) = x1 * ocean_waters_ref(j-1) + x2 * ocean_waters_ref(j)
        refl_fw(i) = x1 * coastal_waters_ref(j-1) + x2 * coastal_waters_ref(j)
      ENDIF
    ENDDO

  END SUBROUTINE rttov_refl_water_interp

END MODULE rttov_solar_refl_mod
