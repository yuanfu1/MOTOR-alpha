MODULE mem_vol

  INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: Reflectivity, Reflectivity_HI, RadialVelocity, RadialVelocity_HI

  REAL, ALLOCATABLE, DIMENSION(:, :) :: elevationR, elevationR_HI, elevationV, elevationV_HI

  REAL, ALLOCATABLE, DIMENSION(:, :) :: azimuthR, azimuthR_HI, azimuthV, azimuthV_HI

  REAL, ALLOCATABLE, DIMENSION(:) :: distanceR, distanceR_HI, distanceV, distanceV_HI, nyquistVelocityV, nyquistVelocityV_HI

  INTEGER :: gateR, gateR_HI, gateV, gateV_HI, radialR, radialR_HI, &
             radialV, radialV_HI, scanR, scanR_HI, scanV, &
             scanV_HI, nf_fid, nf_vid, nf_status

  CHARACTER*8 :: c8_fname_format

END MODULE
