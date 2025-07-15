subroutine convert_hydrotable(input_file, cconfig, nconfig, ext_tab, ssa_tab, asm_tab, zef_tab)

! Description:
!
!    Creates bulk-scattering coefficient files for use by RTTOV_SCATT
!
!    IN: input_file - input configuration file
!        ext_tab    - Extinction               [km^-1]
!        ssa_tab    - Single scattering albedo [ ]
!        asm_tab    - Asymmetry paramter       [ ]
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date        Comment
! -------   ----        -------
!           06/07/2011  Flexible length config section. Add comments (Alan Geer)

use parkind1, only: jprb, jpim
use mod_scattering, only: lmax, nmaxconfig
!INTF_OFF
use mod_scattering, only: n_temp, n_lwc, &
 & wc_slope, wc_offset, temp_offset, ctype, is_frozen
use rttov_const, only: scattering_coef_filestem, version, len_platform, len_instrument
!INTF_ON
implicit none

character(len=*), intent(in) :: input_file
character (len=lmax), intent(in) :: cconfig(nmaxconfig)
integer (kind=jpim) :: nconfig
real (kind=jprb), dimension (:,:,:,:), intent(in) :: ext_tab, ssa_tab, asm_tab, zef_tab
!INTF_END


character (len=len_platform) :: csat_in  
character (len=len_instrument) :: csens_in 
character (len=lmax) :: cin, cfile  

integer, allocatable :: id_freq(:)
real   , allocatable :: freq_c(:), freq_o(:)
integer, allocatable :: mpol(:)
integer (kind=jpim) :: n_freq, n_sensor
integer (kind=jpim) :: m_type, m_freq, i
integer (kind=jpim), allocatable, dimension (:) :: id_type 
integer, dimension (8) :: d
logical :: l_reflectivity

character (len=lmax) :: cdir_table
character (len=lmax) :: sensor_in_fmt, sensor_out_fmt

integer :: ifreq, jfreq, itype, itemp, isensor, nside, ii, jj, i_plat, i_inst, iostat

!-------------------------------------------------------------------------------

open(11, file = input_file, status = 'old')

do ! Find header
  read (11,*) cin
  if (trim(cin) == 'HEADER') exit
enddo
read (11,*) cin
read (11,'(a500)') cdir_table
read (11,*) n_freq
read (11,*) n_sensor
read (11,*) m_type 
allocate(id_type(m_type), id_freq(10*n_freq), freq_c(n_freq), freq_o(n_freq), mpol(n_freq))
read (11,*) id_type 
read (11,*) cin

do ! Find frequencies and read config lines
  read (11,'(a500)') cin
  if (trim(cin) == 'FREQUENCIES') then
    exit
  endif
enddo
read (11,*) cin
read (11,*) cin
read (11,*) cin
read (11,*) cin
read (11,*) cin

do ifreq = 1, n_freq
   read (11,*) jfreq, freq_c (ifreq), nside, freq_o (ifreq), mpol(ifreq)
end do
do i = 1, 5
   read (11,*) cin
end do

! Allow sensor / platform format to evolve if these strings get made longer
write(sensor_in_fmt,'(a2,i0,a4,i0,a10)') '(a' , len_platform , ',x,a' , len_instrument , ',x,2i4,l4)'
write(sensor_out_fmt,'(a2,i0,a4,i0,a1)') '(a' , len_platform , ',x,a' , len_instrument , ')'

do isensor = 1, n_sensor

  read (11,sensor_in_fmt) csat_in, csens_in, i_plat, i_inst, l_reflectivity
  read (11,*)    m_freq, (id_freq (ii), ii = 1, m_freq)   

  cfile = trim(cdir_table) // '/' // scattering_coef_filestem // '_' // trim(csat_in) //  &
                      & '_' // trim(csens_in) // '.dat'
  open (20, file = cfile, iostat=iostat)
  if(iostat /= 0) then
    write(*,*)' Problem opening '//trim(cfile)
    exit
  endif

  !* Header
  write (20,'(a16,a8)') ' ! Hydro-tables   ', csens_in
  call write_divider(20)

  !* Identification section
  write (20,'(a)') 'IDENTIFICATION'
  write (20,'(a)') ' !'
  write (20,'(2i5,a)') i_plat, i_inst,'         ! platform instrument'
  write (20,sensor_out_fmt) csat_in, csens_in
  write (20,'(a)') ' mw                ! sensor type [ir,mw,hi]' 
  write (20,'(i5,a)') version, '              ! RTTOV compatibility version'      
  write (20,'(a)') '  2                ! version'
  call date_and_time(VALUES=d)
  write(20, '(i5,i3,i3,a)') d(1), d(2), d(3), '        ! creation date'
  call write_divider(20)

  !* Dimensions section
  write (20,'(a)') 'DIMENSIONS'
  write (20,'(a)') ' !'               
  write (20,'(i3,i2,i3,i4,l2,a)') m_freq, m_type, n_temp, n_lwc, l_reflectivity, '        ! frequencies hydrometeor-types temperatures water-contents reflectivity(T/F)'
  call write_divider(20)

  !* Frequencies
  write (20,'(a)') 'FREQUENCIES'
  write (20,'(30f9.4)') (freq_c (id_freq (jj)), jj = 1, m_freq)
  call write_divider(20)

  !* Polarisation
  write (20,'(a)') 'POLARISATION'
  write (20,'(30I9)') (mpol   (id_freq (jj)), jj = 1, m_freq)
  call write_divider(20)
  
  !* Hydrometeor types  
  write (20,'(a)') 'HYDROMETEOR'    
  write (20,'(a)') ' ! name, ID, is_frozen'   
  write (20,'(100a)') (trim(ctype(id_type(ii)))//' ', ii = 1, m_type)
  write (20,'(30I9)') (id_type(ii), ii = 1, m_type)
  write (20,'(30L9)') (is_frozen(id_type(ii)), ii = 1, m_type)
  call write_divider(20)  
  
  write (20,'(a)') 'CONVERSIONS'
  write (20,'(a)') ' ! offset          ! T [K] = index [1..70] + offset  '
  write (20,'(30I4)') (int(temp_offset(id_type(jj))), jj = 1, m_type)
  write (20,'(a)') ' ! slope offset    ! LWC [g/m^3] = 10.0 ** (slope * (index [1..nwc] - offset))'
  write (20,'(F4.2,I5)' ) wc_slope, wc_offset             
  call write_divider(20)  

  !* Note the configuration options used in the hydrotable calculation (text copied straight from channels.dat) 
  do i=1,nconfig
    write (20,'(a)') trim(cconfig(i))  
  enddo

  !* Write selected output    
  write (20,'(a)') 'EXTINCTION'
  do ifreq = 1, m_freq 
    do itype = 1, m_type
      do itemp = 1, n_temp
        write (20,1111) ext_tab (id_freq (ifreq),itype,itemp,:)
      end do
    end do
  end do

  write (20,'(a)') 'ALBEDO'
  do ifreq = 1, m_freq 
    do itype = 1, m_type
      do itemp = 1, n_temp
        write (20,1111) ssa_tab (id_freq (ifreq),itype,itemp,:)
      end do
    end do
  end do

  write (20,'(a)') 'ASYMMETRY'
  do ifreq = 1, m_freq 
    do itype = 1, m_type
      do itemp = 1, n_temp
        write (20,1111) asm_tab (id_freq (ifreq),itype,itemp,:)
      end do
    end do
  end do

  if(l_reflectivity) then
    write (20,'(a)') 'REFLECTIVITY'
    do ifreq = 1, m_freq
      do itype = 1, m_type
        do itemp = 1, n_temp
          write (20,1111) zef_tab (id_freq (ifreq),itype,itemp,:)
        end do
      end do
    end do
  endif

  write (*,*) 'Written: ', trim(cfile)

  close (20)
end do

1111    format (5(1x,e23.16))

return

contains

subroutine write_divider(unit)
integer, intent(in) :: unit
write(unit,'(a)')' ! ------------------------------------------------------'
end subroutine write_divider

end subroutine convert_hydrotable


