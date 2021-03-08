!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROSOXI oi para kato maps prokiptoun apo Xamiltonianes oi opies:
!
! A)
! einai anexartites apo ti longitutinal thesi l opote panta
! exo delta_f=delta_i
! se diaforetiki periptosi (exatrisi tis xamiltonianis apo to l) ta
! maps tha ine diaforetika
!
! B)
! einai gia somatidia opou o Lorentz fuctor beta (β) ine poli
! konta sti monada beta_relativistic --> 1 (β -->1) ; ultra-relativistic limit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ta para kato ine 6D ALLA
! me to periorismo delta_f=delta_i (des pio pano)

! mia diafora me to MAD-X ine oti xrisimopio to n_steps kai n_ste
! pou ine o sinolikos arithos to maps se kathe integrator eno sto
! MAD-X xrisimopiite o arithmos ton slices n_kick gia na orisis
! ton olokliroti pou tha xrisimopiis



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! module ACCURACY_CONSTANTS_STRUCTURES !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module ACCURACY_CONSTANTS_STRUCTURES

  implicit none

  !!! for integer values that are larger than selected_int_kind(6) (+-2147483647) you must declare it explicitly ex. 10000000000_iac
  !!! for real values that are more accurate than selected_int_kind(6) (8 desimals) you must declare it explicitly ex. 1._rac = 1.000000000000000

  integer, parameter :: iac = selected_int_kind(15)
  integer, parameter :: rac = selected_real_kind(15, 307)

  !!! SELECTED_INT_KIND(R) return the kind value of the smallest integer type that can represent all values ranging from -10^R (exclusive) to 10^R (exclusive)
  !integer, parameter :: k8 = selected_int_kind(8)

  !!! SELECTED_REAL_KIND(P,R) returns the kind value of a real data type with decimal precision of at least P digits, exponent range of at least R
  !!! precision of 32-, 64-, and 128-bit reals
  !integer, parameter :: sp = selected_real_kind(6, 37)
  !integer, parameter :: dp = selected_real_kind(15, 307)
  !integer, parameter :: qp = selected_real_kind(33, 4931)

  ! precision up to that of the machine-compiler-specifics,
  ! ensures that the double and quad types are actually twice and four times the precision of a single
  !!!integer, parameter ::                            &
  !!!  sp = kind(1.0),                                &
  !!!  dp = selected_real_kind(2*precision(1.0_sp)),  &
  !!!  qp = selected_real_kind(2*precision(1.0_dp))

  ! CONSTANDS
  real (kind=rac), parameter :: Pi = acos(-1._rac)
  real (kind=rac), parameter :: clight = 299792458._rac ![m/s]
  real (kind=rac), parameter :: proton_mass = 0.9382720813_rac ![GeV/c^2]

  type conjugate_variables

    real (kind=rac) :: X
    real (kind=rac) :: PX
    real (kind=rac) :: Y
    real (kind=rac) :: PY
    real (kind=rac) :: l
    real (kind=rac) :: delta

  end type conjugate_variables

end module ACCURACY_CONSTANTS_STRUCTURES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! module LATTICE ATTRIBUTES !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module LATTICE_ATTRIBUTES

  use ACCURACY_CONSTANTS_STRUCTURES

  implicit none

  type lattice_description

    character (len=4) :: element_name
    character (len=60) :: element_id
    real (kind=rac) :: element_length
    real (kind=rac) :: hkick_normalised_strenqth
    real (kind=rac) :: vkick_normalised_strenqth
    real (kind=rac) :: dipole_pole_fase_rotation_in
    real (kind=rac) :: dipole_pole_fase_rotation_out
    real (kind=rac) :: bending_angle
    real (kind=rac) :: quadrupole_normalised_strenqth
    real (kind=rac) :: sextupole_normalised_strenqth
    real (kind=rac) :: octupole_normalised_strenqth
    real (kind=rac) :: decapole_normalised_strenqth
    real (kind=rac) :: dodecapole_normalised_strenqth
    real (kind=rac) :: rf_normalised_frequency

  end type lattice_description

end module LATTICE_ATTRIBUTES


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! module LATTICE_PARTICLES_INITIATION !!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module LATTICE_PARTICLES_INITIATION

  use ACCURACY_CONSTANTS_STRUCTURES
  use LATTICE_ATTRIBUTES

  implicit none

  type (lattice_description), dimension(:), allocatable :: lattice
  integer (kind=iac) :: number_elements
  type (conjugate_variables), dimension(:), allocatable :: particle_in_dy_va
  integer (kind=iac) :: number_particles
  character (len=4), dimension(11) :: supported_elements = ["MARK", "DRIF", "HVCO", "MULT", "DFRF", &
															"CSDI", "CRDI", "QUAD", "SEXT", "OCTU", "RFCA"]
  character (len=3) :: comment_symbol='#$%'

  contains

    subroutine GetLines (unit, n_lines, n_active_lines)

      implicit none

      integer , intent(in) :: unit
      integer (kind=iac), intent(out) :: n_lines, n_active_lines
      integer :: io
      character (len=3) :: comment_check

      ! in order to use GetLines by its owne I mast add the file_name in the arguments and use the open

      !open(unit=unit,file=file_name,status='old', action='read')

      n_lines = 0
      n_active_lines = 0
      do

        read(unit=unit, fmt=*, iostat=io) comment_check

        if (io > 0) then

          print*, "ERROR - subroutine GetLines !!!"
          print*, "Some error occured while reading the file with unit: ", unit, " !!!"
          rewind(unit)
          return

        elseif (io == 0 .and. comment_check /= comment_symbol) then

          n_active_lines = n_active_lines + 1 ! all the none commented lines

          n_lines = n_lines + 1

        elseif (io == 0 .and. comment_check == comment_symbol) then

          n_lines = n_lines + 1

        else

          print*, "The file with unit: ",unit ," is read normally."
          print*, "It has ", n_lines, "lines, ", n_lines-n_active_lines, "of them are comented out and the rest", &
                  n_active_lines, "are active."
          exit

        end if

      end do

      rewind (unit)

      !close (unit)

    end subroutine GetLines

    subroutine GetLattice (unit, file_name)

      implicit none

      integer, intent(in) :: unit
      character (len=100), intent(in) :: file_name
      integer (kind=iac) :: number_lines, number_active_lines, ii, nn
      integer :: io
      character (len=1000) :: check_element

      open (unit=unit, file=file_name, status='old', action='read')

      call GetLines (unit, number_lines, number_active_lines)
      number_elements = number_active_lines

      allocate (lattice(number_elements))

      nn = 1
      do ii=1,number_lines

        read (unit=unit, fmt="(A)", iostat=io) check_element

        if (io > 0) then

          print*, "ERROR - subroutine GetLattice !!!"
          print*, "Some error occured while reading the file with unit: ", unit, " !!!"
          rewind(unit)
          return

        else if (io == 0) then

          if (ANY(supported_elements == check_element(1:4))) then

            read (check_element, fmt=*) lattice(nn)%element_name, &
										lattice(nn)%element_id, &
                                        lattice(nn)%element_length, &
                                        lattice(nn)%hkick_normalised_strenqth, &
                                        lattice(nn)%vkick_normalised_strenqth, &
                                        lattice(nn)%dipole_pole_fase_rotation_in, &
                                        lattice(nn)%dipole_pole_fase_rotation_out, &
                                        lattice(nn)%bending_angle, &
                                        lattice(nn)%quadrupole_normalised_strenqth, &
                                        lattice(nn)%sextupole_normalised_strenqth, &
                                        lattice(nn)%octupole_normalised_strenqth, &
                                        lattice(nn)%decapole_normalised_strenqth, &
                                        lattice(nn)%dodecapole_normalised_strenqth, &
                                        lattice(nn)%rf_normalised_frequency
            nn = nn+1

          else if (check_element(1:3) == comment_symbol) then

            ! print*, "The line: ", ii, " is skipt since it is a comment!!!"

          else

            print*, "ERROR - subroutine GetLattice !!!"
            print*, "The element ", check_element, " is not a supported element!!!"
            print*, "The element ", check_element, " is not added in the lattice for tracking!!!"

          endif

        elseif (io < 0) then

          print*, "The file with unit: ",unit ," is read normally."

        endif

      enddo

      rewind(unit)
      close (unit)

    end subroutine GetLattice

    subroutine GetParticlesInitialDistribution (unit, file_name)

      implicit none

      integer, intent(in) :: unit
      character (len=100), intent(in) :: file_name
      integer (kind=iac) :: number_lines, number_active_lines, ii, nn
      integer :: io
      character (len=1000) :: check_ini_conditions

      open (unit=unit, file=file_name, status='old', action='read')

      call GetLines (unit, number_lines, number_active_lines)
      number_particles = number_active_lines

      allocate (particle_in_dy_va(number_particles))

      nn = 1
      do ii=1,number_lines

        read (unit=unit, fmt="(A)", iostat=io) check_ini_conditions

        if (io > 0) then

          print*, "ERROR - subroutine GetParticlesInitialDistribution !!!"
          print*, "Some error occured while reading the file with unit: ", unit, " !!!"
          rewind(unit)
          return

        else if (io == 0) then

          if (check_ini_conditions(1:3) == comment_symbol) then

            ! print*, "The line: ", ii, " is skipt since it is a comment!!!"

          else

            read (check_ini_conditions, fmt=*) particle_in_dy_va(nn)%X, particle_in_dy_va(nn)%PX, &
                                               particle_in_dy_va(nn)%Y, particle_in_dy_va(nn)%PY, &
                                               particle_in_dy_va(nn)%l, particle_in_dy_va(nn)%delta
            nn = nn+1

          endif

        elseif (io < 0) then

          print*, "The file with unit: ",unit ," is read normally."

        endif

      enddo

      rewind(unit)
      close (unit)

    end subroutine GetParticlesInitialDistribution

    subroutine GenerateParticlesInitialPolarDistribution (min_angle_d, max_angle_d, angle_slices, &
                                                          min_sigma, max_sigma, sigma_slices, initial_delta)

      implicit none

      real (kind=rac), intent(in) :: min_angle_d, max_angle_d, min_sigma, max_sigma, initial_delta
      integer (kind=iac), intent(in) :: angle_slices, sigma_slices
      integer (kind=iac) :: ii, jj, nn
      real (kind=rac) :: min_angle, max_angle, radius, theta

      min_angle = min_angle_d*Pi/180._rac
      max_angle = max_angle_d*Pi/180._rac

      if (min_angle > max_angle .or. min_sigma > max_sigma .or. min_angle <= 0._rac .or. max_angle <= 0._rac &
          .or. min_angle >= Pi/2._rac .or. max_angle >= Pi/2._rac .or. min_sigma <= 0._rac .or. max_sigma <= 0._rac &
          .or. angle_slices < 0 .or. sigma_slices < 0) then

        print*, "ERROR - subroutine GenerateParticlesInitialPolarDistribution !!!"
        print*, "Some of the following are violated: "
        print*, "min_angle > max_angle, min_sigma > max_sigma, min_angle <= 0, &
                   max_angle <= 0, min_angle >= Pi/2, max_angle >= Pi/2, min_sigma <= 0, max_sigma <= 0, &
                   angle_slices < 0, sigma_slices < 0 "

      elseif (min_angle == max_angle .or. min_sigma == max_sigma) then

        if (min_angle == max_angle .and. angle_slices > 0) then

          print*, "ERROR - subroutine GenerateParticlesInitialPolarDistribution !!!"
          print*, "The angle_slices must be zero !!!"

          return

        elseif (min_sigma == max_sigma .and. sigma_slices > 0) then

          print*, "ERROR - subroutine GenerateParticlesInitialPolarDistribution !!!"
          print*, "The sigma_slices must be zero !!!"

          return

        elseif (min_angle == max_angle .and. min_sigma == max_sigma) then

          number_particles = 1
          allocate (particle_in_dy_va(number_particles))

          particle_in_dy_va(1)%X = min_sigma*cos(min_angle)
          particle_in_dy_va(1)%PX = 0.0_rac
          particle_in_dy_va(1)%Y = min_sigma*sin(min_angle)
          particle_in_dy_va(1)%PY = 0.0_rac
          particle_in_dy_va(1)%l = 0.0_rac
          particle_in_dy_va(1)%delta = initial_delta

          print*, 'The particles initial conditions are calculated correct!'

        elseif (min_angle == max_angle) then

          number_particles = sigma_slices + 2
          allocate (particle_in_dy_va(number_particles))

          nn=1
          do ii=0,(sigma_slices + 1)

            radius = min_sigma + real(ii,kind=rac)*(max_sigma-min_sigma)/(sigma_slices + 1._rac)

            particle_in_dy_va(nn)%X = radius*cos(min_angle)
            particle_in_dy_va(nn)%PX = 0.0_rac
            particle_in_dy_va(nn)%Y = radius*sin(min_angle)
            particle_in_dy_va(nn)%PY = 0.0_rac
            particle_in_dy_va(nn)%l = 0.0_rac
            particle_in_dy_va(nn)%delta = initial_delta

            nn=nn+1

          enddo

          print*, 'The particles initial conditions are calculated correct!'

        elseif (min_sigma == max_sigma) then

          if (min_sigma == 0._rac) then

            number_particles = 1
            allocate (particle_in_dy_va(number_particles))

            theta = min_angle

            particle_in_dy_va(1)%X = min_sigma*cos(theta)
            particle_in_dy_va(1)%PX = 0.0_rac
            particle_in_dy_va(1)%Y = min_sigma*sin(theta)
            particle_in_dy_va(1)%PY = 0.0_rac
            particle_in_dy_va(1)%l = 0.0_rac
            particle_in_dy_va(1)%delta = initial_delta

            print*, 'The particles initial conditions are calculated correct!'

          else

            number_particles = angle_slices + 2
            allocate (particle_in_dy_va(number_particles))

            nn=1
            do jj=0,(angle_slices + 1)

              theta = min_angle + real(jj,kind=rac)*(max_angle-min_angle)/(angle_slices + 1._rac)

              particle_in_dy_va(nn)%X = min_sigma*cos(theta)
              particle_in_dy_va(nn)%PX = 0.0_rac
              particle_in_dy_va(nn)%Y = min_sigma*sin(theta)
              particle_in_dy_va(nn)%PY = 0.0_rac
              particle_in_dy_va(nn)%l = 0.0_rac
              particle_in_dy_va(nn)%delta = initial_delta

              nn=nn+1

            enddo

            print*, 'The particles initial conditions are calculated correct!'

          endif

        endif

      else

        if (min_sigma == 0._rac) then

          number_particles = (angle_slices + 2)*(sigma_slices + 1) + 1
          allocate (particle_in_dy_va(number_particles))

        else

          number_particles = (angle_slices + 2)*(sigma_slices + 2)
          allocate (particle_in_dy_va(number_particles))

        endif

        nn=1
        do ii=0,(sigma_slices + 1)

          radius = min_sigma + real(ii,kind=rac)*(max_sigma-min_sigma)/(sigma_slices + 1._rac)

          do jj=0,(angle_slices + 1)

            theta = min_angle + real(jj,kind=rac)*(max_angle-min_angle)/(angle_slices + 1._rac)

            if (radius == 0._rac) then

              particle_in_dy_va(nn)%X = radius*cos(theta)
              particle_in_dy_va(nn)%PX = 0.0_rac
              particle_in_dy_va(nn)%Y = radius*sin(theta)
              particle_in_dy_va(nn)%PY = 0.0_rac
              particle_in_dy_va(nn)%l = 0.0_rac
              particle_in_dy_va(nn)%delta = initial_delta

              nn=nn+1

              exit

            else

              particle_in_dy_va(nn)%X = radius*cos(theta)
              particle_in_dy_va(nn)%PX = 0.0_rac
              particle_in_dy_va(nn)%Y = radius*sin(theta)
              particle_in_dy_va(nn)%PY = 0.0_rac
              particle_in_dy_va(nn)%l = 0.0_rac
              particle_in_dy_va(nn)%delta = initial_delta

              nn=nn+1

            endif

          enddo

        enddo

        print*, 'The particles initial conditions are calculated correct!'

      endif

    end subroutine GenerateParticlesInitialPolarDistribution

end module LATTICE_PARTICLES_INITIATION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! module THICK_ELEMENTS !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module THICK_ELEMENTS

  use ACCURACY_CONSTANTS_STRUCTURES

  implicit none

  contains

    function UnitaryTransformation (un_in_dy_va) result (un_fi_dy_va)

      implicit none

      type (conjugate_variables), intent(in) :: un_in_dy_va
      type (conjugate_variables) :: un_fi_dy_va

      un_fi_dy_va%X = un_in_dy_va%X
      un_fi_dy_va%PX = un_in_dy_va%PX
      un_fi_dy_va%Y = un_in_dy_va%Y
      un_fi_dy_va%PY = un_in_dy_va%PY
      un_fi_dy_va%l = un_in_dy_va%l
      un_fi_dy_va%delta = un_in_dy_va%delta

    end function UnitaryTransformation

    ! exact solution, without expanding the sqrt
    function Drift (dr_in_dy_va, drift_len) result (dr_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: drift_len
      type (conjugate_variables), intent(in) :: dr_in_dy_va
      type (conjugate_variables) :: dr_fi_dy_va

      if (drift_len <= 0.0_rac) then

        print*, "ERROR - Drift !!!"
        print*, "The length of the drift must be grater than zero !!!"

      end if

      dr_fi_dy_va%X = dr_in_dy_va%X &
                      + drift_len * (dr_in_dy_va%PX / sqrt((1._rac + dr_in_dy_va%delta)**2 - dr_in_dy_va%PX**2 - dr_in_dy_va%PY**2))
      dr_fi_dy_va%PX = dr_in_dy_va%PX
      dr_fi_dy_va%Y = dr_in_dy_va%Y &
                      + drift_len * (dr_in_dy_va%Py / sqrt((1._rac + dr_in_dy_va%delta)**2 - dr_in_dy_va%PX**2 - dr_in_dy_va%PY**2))
      dr_fi_dy_va%PY = dr_in_dy_va%PY
      dr_fi_dy_va%l = dr_in_dy_va%l + drift_len &
                      * (1._rac - (1._rac + dr_in_dy_va%delta)/sqrt((1._rac + dr_in_dy_va%delta)**2 &
                      - dr_in_dy_va%PX**2 - dr_in_dy_va%PY**2))
      dr_fi_dy_va%delta = dr_in_dy_va%delta

    end function Drift

    function CombinedSectorDipole (csd_in_dy_va, sec_dipole_len, angle, sec_E1, sec_E2, quad_norm_stren) result (csd_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: sec_dipole_len, angle, sec_E1, sec_E2, quad_norm_stren
      type (conjugate_variables), intent(in) :: csd_in_dy_va
      type (conjugate_variables) :: csd_st_dy_va, csd_sts_dy_va, csd_fi_dy_va
      real (kind=rac) :: ash, alf, gam, fie


      if (sec_dipole_len <= 0.0_rac) then ! .or. angle <=0.0_rac ?

        print*, 'ERROR - CombinedSectorDipole !!!'
        print*, "The length of the sector dipole must be a possitive number (length > 0) !!!"
        return

      endif

      ! bending_radius = sec_dipole_len/angle
      ash = angle/sec_dipole_len
      alf = (ash**2 + quad_norm_stren)/(csd_in_dy_va%delta + 1._rac)
      gam = ash*csd_in_dy_va%delta/(csd_in_dy_va%delta + 1._rac)
      fie = quad_norm_stren/(csd_in_dy_va%delta + 1._rac)
      
      csd_st_dy_va%X = csd_in_dy_va%X
      csd_st_dy_va%PX = csd_in_dy_va%PX + ash*tan(sec_E1)*csd_in_dy_va%X
      csd_st_dy_va%Y = csd_in_dy_va%Y
      csd_st_dy_va%PY = csd_in_dy_va%PY - ash*tan(sec_E1)*csd_in_dy_va%Y
      csd_st_dy_va%l = csd_in_dy_va%l
      csd_st_dy_va%delta = csd_in_dy_va%delta

      if ( alf >= 0.0_rac .and. fie >= 0.0_rac ) then

        csd_sts_dy_va%X = cos(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%X &
                         + sin(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%PX / (sqrt(alf)*(1._rac + csd_st_dy_va%delta)) &
                         + (1._rac - cos(sqrt(alf)*sec_dipole_len))*gam/alf
        csd_sts_dy_va%PX = cos(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%PX &
                          - (1._rac + csd_st_dy_va%delta)*sqrt(alf)*sin(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%X &
                          + sin(sqrt(alf)*sec_dipole_len)*gam*(1._rac + csd_st_dy_va%delta)/sqrt(alf)
        if (fie == 0.0_rac) then

          csd_sts_dy_va%Y = csd_st_dy_va%Y + sec_dipole_len * csd_st_dy_va%Py / (1._rac + csd_st_dy_va%delta)
          csd_sts_dy_va%PY = csd_st_dy_va%PY
          csd_sts_dy_va%l = csd_st_dy_va%l + (1._rac/4._rac) * ( 2._rac*gam*sec_dipole_len*csd_st_dy_va%X &
                         - alf*sec_dipole_len*csd_st_dy_va%X**2 &
                         - sec_dipole_len*(csd_st_dy_va%PX**2 + 2._rac*csd_st_dy_va%PY**2)/((1._rac + csd_st_dy_va%delta)**2) &
                         + csd_st_dy_va%PX*csd_st_dy_va%X/(1._rac + csd_st_dy_va%delta) &
                         - (gam+4._rac*ash)*(csd_st_dy_va%PX+gam*sec_dipole_len*(1._rac + csd_st_dy_va%delta)) &
                         / (alf*(1._rac + csd_st_dy_va%delta)) + ((4._rac*ash*(1._rac + csd_st_dy_va%delta) &
                         + (gam-alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta)*cos(sqrt(alf)*sec_dipole_len) &
                         - sqrt(alf)*csd_st_dy_va%PX*sin(sqrt(alf)*sec_dipole_len)) &
                         * (sqrt(alf)*csd_st_dy_va%PX*cos(sqrt(alf)*sec_dipole_len) &
                         + (gam-alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta)*sin(sqrt(alf)*sec_dipole_len) )) &
                         / (sqrt(alf**3)*(1._rac + csd_st_dy_va%delta)**2) )

        else

          csd_sts_dy_va%Y = cosh(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%Y &
                         + sinh(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%PY / (sqrt(fie)*(1._rac + csd_st_dy_va%delta))
          csd_sts_dy_va%PY = cosh(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%PY &
                          + (1._rac + csd_st_dy_va%delta)*sqrt(fie)*sinh(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%Y
          csd_sts_dy_va%l = csd_st_dy_va%l - (1._rac/(8._rac*alf*(1._rac+csd_st_dy_va%delta)**2)) &
                         * ( 2._rac*csd_st_dy_va%PX*(gam-alf*csd_st_dy_va%X)*(1._rac+csd_st_dy_va%delta) &
                         + 2._rac*sec_dipole_len*(alf*csd_st_dy_va%PX**2 &
                         + ((gam-alf*csd_st_dy_va%X)*(1._rac+csd_st_dy_va%delta))**2) &
                         - 2._rac*csd_st_dy_va%PX*(gam-alf*csd_st_dy_va%X) &
                         * (1._rac+csd_st_dy_va%delta)*cos(2._rac*sqrt(alf)*sec_dipole_len) &
                         + 8._rac*ash*(1._rac+csd_st_dy_va%delta)*(csd_st_dy_va%PX+gam*sec_dipole_len*(1._rac+csd_st_dy_va%delta) &
                         - csd_st_dy_va%PX*cos(sqrt(alf)*sec_dipole_len) &
                         + ((alf*csd_st_dy_va%X-gam)*(1._rac+csd_st_dy_va%delta)*sin(sqrt(alf)*sec_dipole_len))/(sqrt(alf))) &
                         - ((-alf*csd_st_dy_va%PX**2 + ((gam-alf*csd_st_dy_va%X)*(1._rac+csd_st_dy_va%delta))**2) &
                         * sin(2._rac*sqrt(alf)*sec_dipole_len))/(sqrt(alf)) &
                         + alf*(2._rac*sec_dipole_len*csd_st_dy_va%PY**2 &
                         - 2._rac*csd_st_dy_va%PY*csd_st_dy_va%Y*(1._rac+csd_st_dy_va%delta) &
                         - 2._rac*fie*sec_dipole_len*(csd_st_dy_va%Y*(1._rac+csd_st_dy_va%delta))**2 &
                         + 2._rac*csd_st_dy_va%PY*csd_st_dy_va%Y*(1._rac+csd_st_dy_va%delta)*cosh(2._rac*sqrt(fie)*sec_dipole_len) &
                         + ((csd_st_dy_va%PY**2 + fie*(csd_st_dy_va%Y*(1._rac+csd_st_dy_va%delta))**2) &
                         * sinh(2._rac*sqrt(fie)*sec_dipole_len) )/(sqrt(fie))) )

        endif

        csd_sts_dy_va%delta = csd_st_dy_va%delta


      elseif ( alf >= 0.0_rac .and. fie < 0.0_rac ) then

        fie = abs(fie)

        csd_sts_dy_va%X = cos(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%X &
                         + sin(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%PX / (sqrt(alf)*(1._rac + csd_st_dy_va%delta)) &
                         + (1._rac - cos(sqrt(alf)*sec_dipole_len))*gam/alf
        csd_sts_dy_va%PX =  cos(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%PX &
                           - (1._rac + csd_st_dy_va%delta)*sqrt(alf)*sin(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%X &
                           + sin(sqrt(alf)*sec_dipole_len)*gam*(1._rac + csd_st_dy_va%delta)/sqrt(alf)
        csd_sts_dy_va%Y = cos(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%Y &
                          + sin(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%PY / (sqrt(fie)*(1._rac + csd_st_dy_va%delta))
        csd_sts_dy_va%PY = cos(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%PY &
                          - (1._rac + csd_st_dy_va%delta)*sqrt(fie)*sin(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%Y
        csd_sts_dy_va%l = csd_st_dy_va%l + (1._rac/8._rac) * ( 4._rac*gam*sec_dipole_len*csd_st_dy_va%X &
                         - 2._rac*alf*sec_dipole_len*csd_st_dy_va%X**2 - 2._rac*fie*sec_dipole_len*csd_st_dy_va%Y**2 &
                         - (2._rac*(csd_st_dy_va%PX**2 + csd_st_dy_va%PY**2)*sec_dipole_len)/((1._rac + csd_st_dy_va%delta)**2) &
                         + (2._rac*(csd_st_dy_va%PX*csd_st_dy_va%X +csd_st_dy_va%PY*csd_st_dy_va%Y))/(1._rac + csd_st_dy_va%delta) &
                         - (2._rac*(gam+4._rac*ash)*(csd_st_dy_va%PX+gam*sec_dipole_len*(1._rac + csd_st_dy_va%delta))) &
                         / (alf*(1._rac + csd_st_dy_va%delta)) + (1._rac/(1._rac + csd_st_dy_va%delta)**2) &
                         * (-2._rac*csd_st_dy_va%PY*csd_st_dy_va%Y*(1._rac + csd_st_dy_va%delta) &
                         * cos(2._rac*sqrt(fie)*sec_dipole_len) + (2._rac/sqrt(alf**3)) &
                         * (4._rac*ash*(1._rac + csd_st_dy_va%delta) + (gam-alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta) &
                         * cos(sqrt(alf)*sec_dipole_len) - sqrt(alf)*csd_st_dy_va%PX*sin(sqrt(alf)*sec_dipole_len) ) &
                         * (sqrt(alf)*csd_st_dy_va%PX*cos(sqrt(alf)*sec_dipole_len) &
                         + (gam-alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta)*sin(sqrt(alf)*sec_dipole_len) ) &
                         + (1._rac/sqrt(fie))*(-csd_st_dy_va%PY**2 + fie*(csd_st_dy_va%Y*(1._rac + csd_st_dy_va%delta))**2 ) &
                         * sin(2._rac*sqrt(fie)*sec_dipole_len) ) )
        csd_sts_dy_va%delta = csd_st_dy_va%delta


      elseif ( alf < 0.0_rac .and. fie >= 0.0_rac ) then

        alf = abs(alf)

        csd_sts_dy_va%X = (1._rac/alf) * (-gam + (gam+alf*csd_st_dy_va%X)*cosh(sqrt(alf)*sec_dipole_len) &
                         + sqrt(alf)*sinh(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%PX / (1._rac + csd_st_dy_va%delta) )
        csd_sts_dy_va%PX = cosh(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%PX &
                          + (gam+alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta)*sinh(sqrt(alf)*sec_dipole_len)/sqrt(alf)

        if (fie == 0.0_rac) then

          csd_sts_dy_va%Y = csd_st_dy_va%Y + sec_dipole_len * csd_st_dy_va%Py / (1._rac + csd_st_dy_va%delta)
          csd_sts_dy_va%PY = csd_st_dy_va%PY
          csd_sts_dy_va%l = csd_st_dy_va%l + (1._rac/(4._rac*sqrt(alf**3)*(1._rac + csd_st_dy_va%delta)**2)) &
                         * ( sqrt(alf) * ( sec_dipole_len*(alf*csd_st_dy_va%X*(1._rac + csd_st_dy_va%delta))**2 &
                         + (gam+4._rac*ash)*(1._rac + csd_st_dy_va%delta)*(csd_st_dy_va%PX &
                         + gam*sec_dipole_len*(1._rac + csd_st_dy_va%delta)) + alf*(-sec_dipole_len*csd_st_dy_va%PX**2 &
                         + csd_st_dy_va%PX*csd_st_dy_va%X*(1._rac + csd_st_dy_va%delta) + 2._rac*sec_dipole_len &
                         * (-csd_st_dy_va%PY**2 +gam*csd_st_dy_va%X*(1._rac + csd_st_dy_va%delta)**2) ) ) &
                         - (4._rac*ash*(1._rac + csd_st_dy_va%delta)+(gam+alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta) &
                         * cosh(sqrt(alf)*sec_dipole_len) +sqrt(alf)*csd_st_dy_va%PX*sinh(sqrt(alf)*sec_dipole_len)) &
                         * ( sqrt(alf)*csd_st_dy_va%PX*cosh(sqrt(alf)*sec_dipole_len) &
                         + (gam+alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta)*sinh(sqrt(alf)*sec_dipole_len) ) )

        else

          csd_sts_dy_va%Y = cosh(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%Y &
                         + sinh(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%PY / (sqrt(fie)*(1._rac + csd_st_dy_va%delta))
          csd_sts_dy_va%PY = cosh(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%PY &
                          + (1._rac + csd_st_dy_va%delta)*sqrt(fie)*sinh(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%Y
          csd_sts_dy_va%l = csd_st_dy_va%l + (1._rac/4._rac) * ( gam*(gam+4._rac*ash)*sec_dipole_len/alf &
                         + sec_dipole_len*(2._rac*gam*csd_st_dy_va%X+alf*csd_st_dy_va%X**2 + fie*csd_st_dy_va%Y**2) &
                         - (csd_st_dy_va%PX**2 + csd_st_dy_va%PY**2)*sec_dipole_len/((1._rac + csd_st_dy_va%delta)**2) &
                         + (gam*csd_st_dy_va%PX+4._rac*ash*csd_st_dy_va%PX+alf*csd_st_dy_va%PX*csd_st_dy_va%X &
                         + alf*csd_st_dy_va%PY*csd_st_dy_va%Y)/(alf*(1._rac + csd_st_dy_va%delta)) ) &
                         - (1._rac/(8._rac*(1._rac + csd_st_dy_va%delta)**2)) &
                         * ( 2._rac*csd_st_dy_va%PY*csd_st_dy_va%Y*(1._rac + csd_st_dy_va%delta) &
                         * cosh(2._rac*sqrt(fie)*sec_dipole_len) + (2._rac/sqrt(alf**3))*(4._rac*ash*(1._rac + csd_st_dy_va%delta) &
                         + (gam+alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta)*cosh(sqrt(alf)*sec_dipole_len) &
                         + sqrt(alf)*csd_st_dy_va%PX*sinh(sqrt(alf)*sec_dipole_len)) &
                         * ( sqrt(alf)*csd_st_dy_va%PX*cosh(sqrt(alf)*sec_dipole_len) &
                         + (gam+alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta)*sinh(sqrt(alf)*sec_dipole_len) ) &
                         + (csd_st_dy_va%PY**2 + fie*(csd_st_dy_va%Y*(1._rac + csd_st_dy_va%delta))**2) &
                         * sinh(2._rac*sqrt(fie)*sec_dipole_len)/sqrt(fie) )

        endif

        csd_sts_dy_va%delta = csd_st_dy_va%delta


      elseif ( alf < 0.0_rac .and. fie < 0.0_rac ) then

          alf = abs(alf)
          fie = abs(fie)

          csd_sts_dy_va%X = (1._rac/alf) * (-gam + (gam+alf*csd_st_dy_va%X)*cosh(sqrt(alf)*sec_dipole_len) &
                         + sqrt(alf)*sinh(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%PX / (1._rac + csd_st_dy_va%delta) )
          csd_sts_dy_va%PX = cosh(sqrt(alf)*sec_dipole_len)*csd_st_dy_va%PX &
                          + (gam+alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta)*sinh(sqrt(alf)*sec_dipole_len)/sqrt(alf)
          csd_sts_dy_va%Y = cos(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%Y &
                          + sin(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%PY / (sqrt(fie)*(1._rac + csd_st_dy_va%delta))
          csd_sts_dy_va%PY = cos(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%PY &
                          - (1._rac + csd_st_dy_va%delta)*sqrt(fie)*sin(sqrt(fie)*sec_dipole_len)*csd_st_dy_va%Y
          csd_sts_dy_va%l = csd_st_dy_va%l + (1._rac/(8._rac*(1._rac + csd_st_dy_va%delta)**2)) &
                         * (2._rac*( -sec_dipole_len*csd_st_dy_va%PX**2 -sec_dipole_len*csd_st_dy_va%PY**2 &
                         + csd_st_dy_va%PX*csd_st_dy_va%X + 2._rac*gam*sec_dipole_len*csd_st_dy_va%X &
                         + csd_st_dy_va%PY*csd_st_dy_va%Y - fie*sec_dipole_len*csd_st_dy_va%Y**2 &
                         + alf*sec_dipole_len*(csd_st_dy_va%X*(1._rac + csd_st_dy_va%delta))**2 &
                         + (gam+4._rac*ash)*(1._rac + csd_st_dy_va%delta) &
                         * (csd_st_dy_va%PX+gam*sec_dipole_len*(1._rac + csd_st_dy_va%delta))/alf &
                         + csd_st_dy_va%delta*(csd_st_dy_va%PX*csd_st_dy_va%X + csd_st_dy_va%PY*csd_st_dy_va%Y &
                         + sec_dipole_len*(2._rac*gam*csd_st_dy_va%X-fie*csd_st_dy_va%Y**2)*(2._rac + csd_st_dy_va%delta)) &
                         - csd_st_dy_va%PY*csd_st_dy_va%Y*(1._rac + csd_st_dy_va%delta)*cos(2._rac*sqrt(fie)*sec_dipole_len) ) &
                         + (-csd_st_dy_va%PY**2 + fie*(csd_st_dy_va%Y*(1._rac + csd_st_dy_va%delta))**2) &
                         * sin(2._rac*sqrt(fie)*sec_dipole_len)/sqrt(fie) - (2._rac/sqrt(alf**3)) &
                         * (4._rac*ash*(1._rac + csd_st_dy_va%delta) + (gam+alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta) &
                         * cosh(sqrt(alf)*sec_dipole_len) + sqrt(alf)*csd_st_dy_va%PX*sinh(sqrt(alf)*sec_dipole_len) ) &
                         * (sqrt(alf)*csd_st_dy_va%PX*cosh(sqrt(alf)*sec_dipole_len) &
                         + (gam+alf*csd_st_dy_va%X)*(1._rac + csd_st_dy_va%delta)*sinh(sqrt(alf)*sec_dipole_len)) )
          csd_sts_dy_va%delta = csd_st_dy_va%delta

      endif
      
      csd_fi_dy_va%X = csd_sts_dy_va%X
      csd_fi_dy_va%PX = csd_sts_dy_va%PX + ash*tan(sec_E2)*csd_sts_dy_va%X
      csd_fi_dy_va%Y = csd_sts_dy_va%Y
      csd_fi_dy_va%PY = csd_sts_dy_va%PY - ash*tan(sec_E2)*csd_sts_dy_va%Y
      csd_fi_dy_va%l = csd_sts_dy_va%l
      csd_fi_dy_va%delta = csd_sts_dy_va%delta

    end function CombinedSectorDipole
    
    function CombinedRectangularDipole (crd_in_dy_va, rec_dipole_len, angle, rec_E1, rec_E2, quad_norm_stren) result (crd_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: rec_dipole_len, angle, rec_E1, rec_E2, quad_norm_stren
      type (conjugate_variables), intent(in) :: crd_in_dy_va
      type (conjugate_variables) :: crd_st_dy_va, crd_fi_dy_va
      real (kind=rac) :: half_angle, ash


      if (rec_dipole_len <= 0.0_rac) then ! .or. angle <=0.0_rac ?

        print*, 'ERROR - CombinedRectangularDipole !!!'
        print*, "The length of the rectangular dipole must be a possitive number (length > 0) !!!"
        return

      endif

      ! bending_radius = rec_dipole_len/angle
      half_angle = 0.5_rac*angle
      ash = angle/rec_dipole_len
      
      crd_st_dy_va%X = crd_in_dy_va%X
      crd_st_dy_va%PX = crd_in_dy_va%PX + ash*tan(rec_E1)*crd_in_dy_va%X
      crd_st_dy_va%Y = crd_in_dy_va%Y
      crd_st_dy_va%PY = crd_in_dy_va%PY - ash*tan(rec_E1)*crd_in_dy_va%Y
      crd_st_dy_va%l = crd_in_dy_va%l
      crd_st_dy_va%delta = crd_in_dy_va%delta
      
      crd_st_dy_va = CombinedSectorDipole (crd_st_dy_va, rec_dipole_len, angle, half_angle, half_angle, quad_norm_stren)
      
      crd_fi_dy_va%X = crd_st_dy_va%X
      crd_fi_dy_va%PX = crd_st_dy_va%PX + ash*tan(rec_E2)*crd_st_dy_va%X
      crd_fi_dy_va%Y = crd_st_dy_va%Y
      crd_fi_dy_va%PY = crd_st_dy_va%PY - ash*tan(rec_E2)*crd_st_dy_va%Y
      crd_fi_dy_va%l = crd_st_dy_va%l
      crd_fi_dy_va%delta = crd_st_dy_va%delta
      
    end function CombinedRectangularDipole
    
    function HorizontalCorrector (hc_in_dy_va, core_len, core_norm_stren) result (hc_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: core_len, core_norm_stren
      type (conjugate_variables), intent(in) :: hc_in_dy_va
      type (conjugate_variables) :: hc_fi_dy_va
      real (kind=rac) :: eff_delta, phi

	  if (core_norm_stren == 0.0_rac .or. core_len <= 0.0_rac) then

        print*, "ERROR - function HorizontalCorrector !!!"
        print*, "The stength of the horizontal corrector can not be equal to zero !!!"
        print*, "The length of the horizontal corrector must be greater than zero !!!"

      end if

	  eff_delta = sqrt(1._rac + hc_in_dy_va%delta)
	  phi = core_len*core_norm_stren/eff_delta
	  
	  hc_fi_dy_va%X = cos(phi)*hc_in_dy_va%X &
				    + sin(phi)*hc_in_dy_va%PX/(eff_delta*core_norm_stren) &
				    +(1._rac + cos(phi))*hc_in_dy_va%delta/core_norm_stren
	  hc_fi_dy_va%PX = cos(phi)*hc_in_dy_va%PX &
				     - sin(phi)*core_norm_stren*eff_delta*hc_in_dy_va%X &
				     + sin(phi)*eff_delta*hc_in_dy_va%delta
	  hc_fi_dy_va%Y = hc_in_dy_va%Y + core_len*hc_in_dy_va%PY/(1._rac + hc_in_dy_va%delta)
	  hc_fi_dy_va%PY = hc_in_dy_va%PY 
	  hc_fi_dy_va%l = hc_in_dy_va%l -core_norm_stren*hc_fi_dy_va%X &
					- (hc_fi_dy_va%PX**2 + hc_fi_dy_va%PY**2)/(2*(1._rac + hc_in_dy_va%delta)**2)
	  hc_fi_dy_va%delta = hc_in_dy_va%delta

    end function HorizontalCorrector
    
    function VerticalCorrector (vc_in_dy_va, core_len, core_norm_stren) result (vc_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: core_len, core_norm_stren
      type (conjugate_variables), intent(in) :: vc_in_dy_va
      type (conjugate_variables) :: vc_fi_dy_va
      real (kind=rac) :: eff_delta, phi

	  if (core_norm_stren == 0.0_rac .or. core_len <= 0.0_rac) then

        print*, "ERROR - function VerticalCorrector !!!"
        print*, "The stength of the vertical corrector can not be equal to zero !!!"
        print*, "The length of the vertical corrector must be greater than zero !!!"

      end if

	  eff_delta = sqrt(1._rac + vc_in_dy_va%delta)
	  phi = core_len*core_norm_stren/eff_delta
	  
	  vc_fi_dy_va%X = vc_in_dy_va%X + core_len*vc_in_dy_va%PX/(1._rac + vc_in_dy_va%delta)
	  vc_fi_dy_va%PX = vc_in_dy_va%PX 
	  vc_fi_dy_va%Y = cos(phi)*vc_in_dy_va%Y &
				    + sin(phi)*vc_in_dy_va%PY/(eff_delta*core_norm_stren) &
				    +(1._rac + cos(phi))*vc_in_dy_va%delta/core_norm_stren
	  vc_fi_dy_va%PY = cos(phi)*vc_in_dy_va%PY &
				     - sin(phi)*core_norm_stren*eff_delta*vc_in_dy_va%Y &
				     + sin(phi)*eff_delta*vc_in_dy_va%delta	  
	  vc_fi_dy_va%l = vc_in_dy_va%l -core_norm_stren*vc_fi_dy_va%Y &
					- (vc_fi_dy_va%PX**2 + vc_fi_dy_va%PY**2)/(2*(1._rac + vc_in_dy_va%delta)**2)
	  vc_fi_dy_va%delta = vc_in_dy_va%delta

    end function VerticalCorrector

    function Quadrupole (q_in_dy_va, quad_len, quad_norm_stren) result (q_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: quad_len, quad_norm_stren
      type (conjugate_variables), intent(in) :: q_in_dy_va
      type (conjugate_variables) :: q_fi_dy_va
      real (kind=rac) :: eff_qns

      ! gia thetiko strength exo foc ston x
      if (quad_norm_stren > 0.0_rac .and. quad_len > 0.0_rac) then

        eff_qns = sqrt(quad_norm_stren/(1._rac + q_in_dy_va%delta))
        q_fi_dy_va%X = cos(eff_qns*quad_len)*q_in_dy_va%X &
                       + sin(eff_qns*quad_len)*q_in_dy_va%PX/(eff_qns*(1._rac + q_in_dy_va%delta))
        q_fi_dy_va%PX = cos(eff_qns*quad_len)*q_in_dy_va%PX &
                       - sin(eff_qns*quad_len)*eff_qns*(1._rac + q_in_dy_va%delta)*q_in_dy_va%X
        q_fi_dy_va%Y = cosh(eff_qns*quad_len)*q_in_dy_va%Y &
                       + sinh(eff_qns*quad_len)*q_in_dy_va%PY/(eff_qns*(1._rac+q_in_dy_va%delta))
        q_fi_dy_va%PY = cosh(eff_qns*quad_len)*q_in_dy_va%PY &
                       + sinh(eff_qns*quad_len)*eff_qns*(1._rac +q_in_dy_va%delta)*q_in_dy_va%Y
        q_fi_dy_va%l = q_in_dy_va%l - (1._rac/(8._rac*(1._rac + q_in_dy_va%delta)**2)) * ( ((sin(eff_qns*2._rac*quad_len) &
                       * (q_in_dy_va%PX**2 - quad_norm_stren*(1._rac + q_in_dy_va%delta)*q_in_dy_va%X**2))/eff_qns) &
                       + ((sinh(eff_qns*2._rac*quad_len)*(q_in_dy_va%PY**2 &
                       + quad_norm_stren*(1._rac + q_in_dy_va%delta)*q_in_dy_va%Y**2))/eff_qns) &
                       + 2._rac*q_in_dy_va%PX*q_in_dy_va%X*(1._rac + q_in_dy_va%delta)*cos(eff_qns*2._rac*quad_len) &
                       + 2._rac*q_in_dy_va%PY*q_in_dy_va%Y*(1._rac + q_in_dy_va%delta)*cosh(eff_qns*2._rac*quad_len) &
                       + 2._rac*quad_len*(q_in_dy_va%PX**2 + quad_norm_stren*(1._rac + q_in_dy_va%delta)*q_in_dy_va%X**2) &
                       + 2._rac*quad_len*(q_in_dy_va%PY**2 - quad_norm_stren*(1._rac + q_in_dy_va%delta)*q_in_dy_va%Y**2) &
                       - 2._rac*q_in_dy_va%PX*q_in_dy_va%X*(1._rac + q_in_dy_va%delta) &
                       - 2._rac*q_in_dy_va%PY*q_in_dy_va%Y*(1._rac + q_in_dy_va%delta) )
        q_fi_dy_va%delta = q_in_dy_va%delta

      ! gia arnitiko strength exo foc ston y
      elseif (quad_norm_stren < 0.0_rac .and. quad_len > 0.0_rac) then

        eff_qns = sqrt(abs(quad_norm_stren)/(1._rac + q_in_dy_va%delta))
        q_fi_dy_va%X = cosh(eff_qns*quad_len)*q_in_dy_va%X &
                       + sinh(eff_qns*quad_len)*q_in_dy_va%PX/(eff_qns*(1._rac+q_in_dy_va%delta))
        q_fi_dy_va%PX = cosh(eff_qns*quad_len)*q_in_dy_va%PX &
                       + sinh(eff_qns*quad_len)*eff_qns*(1._rac +q_in_dy_va%delta)*q_in_dy_va%X
        q_fi_dy_va%Y = cos(eff_qns*quad_len)*q_in_dy_va%Y &
                       + sin(eff_qns*quad_len)*q_in_dy_va%PY/(eff_qns*(1._rac + q_in_dy_va%delta))
        q_fi_dy_va%PY = cos(eff_qns*quad_len)*q_in_dy_va%PY &
                       - sin(eff_qns*quad_len)*eff_qns*(1._rac + q_in_dy_va%delta)*q_in_dy_va%Y
        q_fi_dy_va%l = q_in_dy_va%l - (1._rac/(8._rac*(1._rac + q_in_dy_va%delta)**2)) * ( ((sin(eff_qns*2._rac*quad_len) &
                       * (q_in_dy_va%PY**2 - quad_norm_stren*(1._rac + q_in_dy_va%delta)*q_in_dy_va%Y**2))/eff_qns) &
                       + ((sinh(eff_qns*2._rac*quad_len)*(q_in_dy_va%PX**2 &
                       + quad_norm_stren*(1._rac + q_in_dy_va%delta)*q_in_dy_va%X**2))/eff_qns) &
                       + 2._rac*q_in_dy_va%PY*q_in_dy_va%Y*(1._rac + q_in_dy_va%delta)*cos(eff_qns*2._rac*quad_len) &
                       + 2._rac*q_in_dy_va%PX*q_in_dy_va%X*(1._rac + q_in_dy_va%delta)*cosh(eff_qns*2._rac*quad_len) &
                       + 2._rac*quad_len*(q_in_dy_va%PY**2 + quad_norm_stren*(1._rac + q_in_dy_va%delta)*q_in_dy_va%Y**2) &
                       + 2._rac*quad_len*(q_in_dy_va%PX**2 - quad_norm_stren*(1._rac + q_in_dy_va%delta)*q_in_dy_va%X**2) &
                       - 2._rac*q_in_dy_va%PY*q_in_dy_va%Y*(1._rac + q_in_dy_va%delta) &
                       - 2._rac*q_in_dy_va%PX*q_in_dy_va%X*(1._rac + q_in_dy_va%delta) )
        q_fi_dy_va%delta = q_in_dy_va%delta

      elseif (quad_norm_stren == 0.0_rac .or. quad_len <= 0.0_rac) then

        print*, "ERROR - function Quadrupole !!!"
        print*, "The stength of the quadrupoles can not be equal to zero !!!"
        print*, "The length of the quadrupoles must be greater than zero !!!"

      end if

    end function Quadrupole
    
    function RFcavite (rf_in_dy_va, rf_len, k_rf, phase_slip) result (rf_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: rf_len, k_rf, phase_slip
      type (conjugate_variables), intent(in) :: rf_in_dy_va
      type (conjugate_variables) :: rf_fi_dy_va

        rf_fi_dy_va%X = rf_in_dy_va%X + rf_len*rf_in_dy_va%PX
        rf_fi_dy_va%PX = rf_in_dy_va%PX
        rf_fi_dy_va%Y = rf_in_dy_va%Y + rf_len*rf_in_dy_va%PY
        rf_fi_dy_va%PY = rf_in_dy_va%PY
        rf_fi_dy_va%l = cos(k_rf*rf_len)*rf_in_dy_va%l - rf_in_dy_va%delta*sin(k_rf*rf_len)*phase_slip/k_rf
        rf_fi_dy_va%delta = cos(k_rf*rf_len)*rf_in_dy_va%delta - rf_in_dy_va%l*sin(k_rf*rf_len)*k_rf/phase_slip

      if (k_rf == 0.0_rac .or. rf_len <= 0.0_rac) then

        print*, "ERROR - function RFcavite !!!"
        print*, "The stength of the RF cavite can not be equal to zero !!!"
        print*, "The length of the RF cavite must be greater than zero !!!"

      end if

    end function RFcavite

end module THICK_ELEMENTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!! module THIN_ELEMENTS !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



module THIN_ELEMENTS

  use ACCURACY_CONSTANTS_STRUCTURES

  implicit none

  contains

    function CSDIDriftInt (csdidk_in_dy_va, csdi_len, angle, weight) result (csdidk_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: csdi_len, angle, weight
      type (conjugate_variables), intent(in) :: csdidk_in_dy_va
      type (conjugate_variables) :: csdidk_fi_dy_va
      real (kind=rac) :: ash, phi

      if (csdi_len <= 0.0_rac) then

        print*, "ERROR - function CSDIDriftInt !!!"
        print*, "The length of the CSDI must be grater than zero !!!"

      end if

      ! bending_radius = csdi_len/angle
      ash = angle/csdi_len
      phi = csdidk_in_dy_va%PY*ash*csdi_len*weight/(2._rac*(csdidk_in_dy_va%delta+1._rac))

      if ( csdidk_in_dy_va%PY == 0._rac ) then
        csdidk_fi_dy_va%X = (-1._rac + (1._rac+ash*csdidk_in_dy_va%X) &
                          * (1._rac+(ash*csdi_len*weight*csdidk_in_dy_va%PX/(2._rac*(csdidk_in_dy_va%delta+1._rac))))**2)/ash
      else
        csdidk_fi_dy_va%X = (1._rac/ash)*((1._rac+ash*csdidk_in_dy_va%X) &
                          * (cos(phi)+(csdidk_in_dy_va%PX/csdidk_in_dy_va%PY)*sin(phi))**2 - 1._rac)
      endif
      if ( csdidk_in_dy_va%PY == 0._rac ) then
        csdidk_fi_dy_va%PX = csdidk_in_dy_va%PX &
                           / (1._rac + (ash*csdi_len*weight/(2._rac*(csdidk_in_dy_va%delta+1._rac)))*csdidk_in_dy_va%PX)
      else
        csdidk_fi_dy_va%PX = csdidk_in_dy_va%PY * ((csdidk_in_dy_va%PX - csdidk_in_dy_va%PY*tan(phi)) &
                          / (csdidk_in_dy_va%PY + csdidk_in_dy_va%PX*tan(phi)))
      endif
      if ( csdidk_in_dy_va%PY == 0._rac ) then
        csdidk_fi_dy_va%Y = csdidk_in_dy_va%Y
      else
        csdidk_fi_dy_va%Y = csdidk_in_dy_va%Y + ((1._rac + ash*csdidk_in_dy_va%X)/ash) &
                          * ((phi/csdidk_in_dy_va%PY**2)*(csdidk_in_dy_va%PX**2 + csdidk_in_dy_va%PY**2) &
                          + (sin(2._rac*phi)/(2._rac*csdidk_in_dy_va%PY**2))*(csdidk_in_dy_va%PY**2 - csdidk_in_dy_va%PX**2) &
                          + 2._rac*sin(phi)*sin(phi)*(csdidk_in_dy_va%PX/csdidk_in_dy_va%PY))
      endif
      csdidk_fi_dy_va%PY = csdidk_in_dy_va%PY
      csdidk_fi_dy_va%l = csdidk_in_dy_va%l - phi*(((csdidk_in_dy_va%PX**2 + csdidk_in_dy_va%PY**2) &
                          * (1._rac + ash*csdidk_in_dy_va%X))/(ash*csdidk_in_dy_va%PY*(csdidk_in_dy_va%delta + 1._rac)))
      csdidk_fi_dy_va%l = csdidk_in_dy_va%l - csdi_len*weight*(((csdidk_in_dy_va%PX**2 + csdidk_in_dy_va%PY**2) &
                          * (1._rac + ash*csdidk_in_dy_va%X))/(2._rac*(csdidk_in_dy_va%delta + 1._rac)**2))
      csdidk_fi_dy_va%delta = csdidk_in_dy_va%delta

    end function CSDIDriftInt

    function CSDIKick (csdik_in_dy_va, csdi_len, angle, quad_norm_stren, weight) result (csdik_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: csdi_len, angle, quad_norm_stren, weight
      type (conjugate_variables), intent(in) :: csdik_in_dy_va
      type (conjugate_variables) :: csdik_fi_dy_va
      real (kind=rac) :: ash

      if (csdi_len <= 0.0_rac .or. angle == 0.0_rac) then

        print*, "ERROR - function CSDIKick !!!"
        print*, "The length of the CSDI must be greater than zero !!!"
        print*, "The angle of the CSDI can not be equal to zero !!!"

      end if

      ! bending_radius = csdi_len/angle
      ash = angle/csdi_len

      csdik_fi_dy_va%X = csdik_in_dy_va%X
      csdik_fi_dy_va%PX = csdik_in_dy_va%PX + (ash*csdik_in_dy_va%delta-(ash**2 + quad_norm_stren)*csdik_in_dy_va%X)*weight*csdi_len
      csdik_fi_dy_va%Y = csdik_in_dy_va%Y
      csdik_fi_dy_va%PY = csdik_in_dy_va%PY + (quad_norm_stren*csdi_len*weight) * csdik_in_dy_va%Y
      csdik_fi_dy_va%l = csdik_in_dy_va%l
      csdik_fi_dy_va%delta = csdik_in_dy_va%delta

    end function CSDIKick

    function CSDICorr (csdic_in_dy_va, csdi_len, angle, quad_norm_stren, weight) result (csdic_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: csdi_len, angle, quad_norm_stren, weight
      type (conjugate_variables), intent(in) :: csdic_in_dy_va
      type (conjugate_variables) :: csdic_fi_dy_va
      real (kind=rac) :: ash, ess

      if (csdi_len <= 0.0_rac .or. angle == 0.0_rac) then

        print*, "ERROR - function CSDICorr !!!"
        print*, "The length of the CSDI must be greater than zero !!!"
        print*, "The angle of the CSDI can not be equal to zero !!!"

      end if

      ! bending_radius = csdi_len/angle
      ash = angle/csdi_len
      ess = -(weight*csdi_len**3)/2._rac

      csdic_fi_dy_va%X = csdic_in_dy_va%X
      csdic_fi_dy_va%PX = csdic_in_dy_va%PX - (ash/(csdic_in_dy_va%delta+1._rac) * ((quad_norm_stren*csdic_in_dy_va%Y)**2 &
                          + (csdic_in_dy_va%X*(quad_norm_stren+ash**2)-ash*csdic_in_dy_va%delta)**2))*ess &
                          - (((1._rac+csdic_in_dy_va%X*ash)/(csdic_in_dy_va%delta+1._rac)) * 2._rac*(quad_norm_stren + ash**2) &
                          * (csdic_in_dy_va%X*(quad_norm_stren + ash**2)-ash*csdic_in_dy_va%delta))*ess
      csdic_fi_dy_va%Y = csdic_in_dy_va%Y
      csdic_fi_dy_va%PY = csdic_in_dy_va%PY - ((1._rac+csdic_in_dy_va%X*ash)/(csdic_in_dy_va%delta+1._rac)) &
                          * 2._rac*csdic_in_dy_va%Y*ess*quad_norm_stren**2
      csdic_fi_dy_va%l = csdic_in_dy_va%l - ess*((1._rac+csdic_in_dy_va%X*ash)/(csdic_in_dy_va%delta+1._rac)**2) &
                          * (csdic_in_dy_va%X*(quad_norm_stren + ash**2)*(csdic_in_dy_va%X*(quad_norm_stren + ash**2)+2._rac*ash) &
                          - (ash*csdic_in_dy_va%delta)**2 + (quad_norm_stren*csdic_in_dy_va%Y)**2 - 2*csdic_in_dy_va%delta*ash**2)
      csdic_fi_dy_va%delta = csdic_in_dy_va%delta

    end function CSDICorr
    
    function DFRFKick (dfrfk_in_dy_va, dipole_len, pole_face_angle, dipole_angle, weight) result (dfrfk_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: dipole_len, pole_face_angle, dipole_angle, weight
      type (conjugate_variables), intent(in) :: dfrfk_in_dy_va
      type (conjugate_variables) :: dfrfk_fi_dy_va
      real (kind=rac) :: ash

      if (dipole_angle == 0.0_rac) then

        print*, "ERROR - function DFRFKick !!!"
        print*, "The angle of the DFRFKick can not be equal to zero !!!"

      end if

      ! bending_radius = dipole_len/dipole_angle
      ash = dipole_angle/dipole_len
      
      dfrfk_fi_dy_va%X = dfrfk_in_dy_va%X
      dfrfk_fi_dy_va%PX = dfrfk_in_dy_va%PX + weight*ash*tan(pole_face_angle)*dfrfk_in_dy_va%X
      dfrfk_fi_dy_va%Y = dfrfk_in_dy_va%Y
      dfrfk_fi_dy_va%PY = dfrfk_in_dy_va%PY - weight*ash*tan(pole_face_angle)*dfrfk_in_dy_va%Y
      dfrfk_fi_dy_va%l = dfrfk_in_dy_va%l
      dfrfk_fi_dy_va%delta = dfrfk_in_dy_va%delta

    end function DFRFKick

    function DriftInt (dk_in_dy_va, drift_len, weight) result (dk_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: drift_len, weight
      type (conjugate_variables), intent(in) :: dk_in_dy_va
      type (conjugate_variables) :: dk_fi_dy_va

      if (drift_len <= 0.0_rac) then

        print*, "ERROR - function DriftInt !!!"
        print*, "The length of the drift must be grater than zero !!!"

      end if

      dk_fi_dy_va%X = dk_in_dy_va%X + (drift_len*weight) * (dk_in_dy_va%PX / (1._rac + dk_in_dy_va%delta))
      dk_fi_dy_va%PX = dk_in_dy_va%PX
      dk_fi_dy_va%Y = dk_in_dy_va%Y + (drift_len*weight) * (dk_in_dy_va%Py / (1._rac + dk_in_dy_va%delta))
      dk_fi_dy_va%PY = dk_in_dy_va%PY
      dk_fi_dy_va%l = dk_in_dy_va%l &
                      - (drift_len*weight) * ((dk_in_dy_va%PX**2 + dk_in_dy_va%PY**2)/(2._rac*(1._rac+dk_in_dy_va%delta)**2))
      dk_fi_dy_va%delta = dk_in_dy_va%delta

    end function DriftInt
    
    
    function HoVeKick (hvk_in_dy_va, hvkic_len, hkic_norm_stren, vkic_norm_stren, weight) result (hvk_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: hvkic_len, hkic_norm_stren, vkic_norm_stren, weight
      type (conjugate_variables), intent(in) :: hvk_in_dy_va
      real (kind=rac) :: effective_len
      type (conjugate_variables) :: hvk_fi_dy_va

      if (hkic_norm_stren == 0.0_rac .and. vkic_norm_stren == 0.0_rac) then

        print*, "ERROR - function HoVeKick !!!"
        print*, "The stength of the horizontal/vertical kick can not be equal to zero !!!"

      end if
      
      effective_len = hvkic_len
      
      if (hvkic_len < 0.0_rac) then
      
        print*, "ERROR - function HoVeKick !!!"
        print*, "The length of the multipole must be positive !!!"
      
      elseif (hvkic_len == 0.0_rac) then
      
		effective_len = 1.0_rac
		
	  endif

      hvk_fi_dy_va%X = hvk_in_dy_va%X
      hvk_fi_dy_va%PX = hvk_in_dy_va%PX + (hkic_norm_stren*effective_len*weight)
      hvk_fi_dy_va%Y = hvk_in_dy_va%Y
      hvk_fi_dy_va%PY = hvk_in_dy_va%PY + (vkic_norm_stren*effective_len*weight)
      hvk_fi_dy_va%l = hvk_in_dy_va%l
      hvk_fi_dy_va%delta = hvk_in_dy_va%delta

    end function HoVeKick
    

    function QuadKick (qk_in_dy_va, quad_len, quad_norm_stren, weight) result (qk_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: quad_len, quad_norm_stren, weight
      type (conjugate_variables), intent(in) :: qk_in_dy_va
      type (conjugate_variables) :: qk_fi_dy_va

      if (quad_norm_stren == 0.0_rac .or. quad_len <= 0.0_rac) then

        print*, "ERROR - function QuadKick !!!"
        print*, "The stength of the quadrupoles can not be equal to zero !!!"
        print*, "The length of the quadrupoles must be greater than zero !!!"

      end if

      qk_fi_dy_va%X = qk_in_dy_va%X
      qk_fi_dy_va%PX = qk_in_dy_va%PX - (quad_norm_stren*quad_len*weight) * qk_in_dy_va%X
      qk_fi_dy_va%Y = qk_in_dy_va%Y
      qk_fi_dy_va%PY = qk_in_dy_va%PY + (quad_norm_stren*quad_len*weight) * qk_in_dy_va%Y
      qk_fi_dy_va%l = qk_in_dy_va%l
      qk_fi_dy_va%delta = qk_in_dy_va%delta

    end function QuadKick


    function QuadCorr (qc_in_dy_va, quad_len, quad_norm_stren, weight) result (qc_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: quad_len, quad_norm_stren, weight
      type (conjugate_variables), intent(in) :: qc_in_dy_va
      type (conjugate_variables) :: qc_fi_dy_va

      if (quad_norm_stren == 0.0_rac .or. quad_len <= 0.0_rac) then

        print*, "ERROR - function QuadCorr !!!"
        print*, "The stength of the quadrupoles can not be equal to zero !!!"
        print*, "The length of the quadrupoles must be greater than zero !!!"

      end if

      qc_fi_dy_va%X = qc_in_dy_va%X
      qc_fi_dy_va%PX = qc_in_dy_va%PX + (quad_norm_stren**2 * quad_len**3 * weight) * (qc_in_dy_va%X / (1._rac + qc_in_dy_va%delta))
      qc_fi_dy_va%Y = qc_in_dy_va%Y
      qc_fi_dy_va%PY = qc_in_dy_va%PY + (quad_norm_stren**2 * quad_len**3 * weight) * (qc_in_dy_va%Y / (1._rac + qc_in_dy_va%delta))
      qc_fi_dy_va%l = qc_in_dy_va%l + quad_norm_stren**2 * quad_len**3 * weight &
                      * (qc_in_dy_va%X**2 + qc_in_dy_va%Y**2)/(2._rac*(1._rac + qc_in_dy_va%delta)**2)
      qc_fi_dy_va%delta = qc_in_dy_va%delta

    end function QuadCorr


    function SextKick (sk_in_dy_va, sex_len, sex_norm_stren, weight) result (sk_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: sex_len, sex_norm_stren, weight
      type (conjugate_variables), intent(in) :: sk_in_dy_va
      type (conjugate_variables) :: sk_fi_dy_va

      if (sex_norm_stren == 0.0_rac .or. sex_len <= 0.0_rac) then

        print*, "ERROR - function SextKick !!!"
        print*, "The stength of the sextupoles can not be equal to zero !!!"
        print*, "The length of the sextupoles must be greater than zero !!!"

      end if

      sk_fi_dy_va%X = sk_in_dy_va%X
      sk_fi_dy_va%PX = sk_in_dy_va%PX - (sex_norm_stren*sex_len*weight/2._rac) * (sk_in_dy_va%X**2 - sk_in_dy_va%Y**2)
      sk_fi_dy_va%Y = sk_in_dy_va%Y
      sk_fi_dy_va%PY = sk_in_dy_va%PY + (sex_norm_stren*sex_len*weight) * sk_in_dy_va%X*sk_in_dy_va%Y
      sk_fi_dy_va%l = sk_in_dy_va%l
      sk_fi_dy_va%delta = sk_in_dy_va%delta

    end function SextKick


    function SextCorr (sc_in_dy_va, sex_len, sex_norm_stren, weight) result (sc_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: sex_len, sex_norm_stren, weight
      type (conjugate_variables), intent(in) :: sc_in_dy_va
      type (conjugate_variables) :: sc_fi_dy_va

      if (sex_norm_stren == 0.0_rac .or. sex_len <= 0.0_rac) then

        print*, "ERROR - function SextCorr !!!"
        print*, "The stength of the sextupoles can not be equal to zero !!!"
        print*, "The length of the sextupoles must be greater than zero !!!"

      end if

      sc_fi_dy_va%X = sc_in_dy_va%X
      sc_fi_dy_va%PX = sc_in_dy_va%PX + (sex_norm_stren**2 * sex_len**3 * weight) &
                       * ((sc_in_dy_va%X**3 + sc_in_dy_va%X*sc_in_dy_va%Y**2) / (2._rac*(1._rac + sc_in_dy_va%delta)))
      sc_fi_dy_va%Y = sc_in_dy_va%Y
      sc_fi_dy_va%PY = sc_in_dy_va%PY + (sex_norm_stren**2 * sex_len**3 * weight) &
                       * ((sc_in_dy_va%Y**3 + sc_in_dy_va%Y*sc_in_dy_va%X**2) / (2._rac*(1._rac + sc_in_dy_va%delta)))
      sc_fi_dy_va%l = sc_in_dy_va%l + sex_norm_stren**2 * sex_len**3 * weight &
                      * ((sc_in_dy_va%X**2 + sc_in_dy_va%Y**2)**2)/(8._rac*(1._rac + sc_in_dy_va%delta)**2)
      sc_fi_dy_va%delta = sc_in_dy_va%delta

    end function SextCorr


    function OctuKick (ok_in_dy_va, oct_len, oct_norm_stren, weight) result (ok_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: oct_len, oct_norm_stren, weight
      type (conjugate_variables), intent(in) :: ok_in_dy_va
      type (conjugate_variables) :: ok_fi_dy_va

      if (oct_norm_stren == 0.0_rac .or. oct_len <= 0.0_rac) then

        print*, "ERROR - function OctuKick !!!"
        print*, "The stength of the octupoles can not be equal to zero !!!"
        print*, "The length of the octupoles must be greater than zero !!!"

      end if

      ok_fi_dy_va%X = ok_in_dy_va%X
      ok_fi_dy_va%PX = ok_in_dy_va%PX &
                       - (oct_norm_stren*oct_len*weight/6._rac) * (ok_in_dy_va%X**3 - 3._rac*ok_in_dy_va%X*ok_in_dy_va%Y**2)
      ok_fi_dy_va%Y = ok_in_dy_va%Y
      ok_fi_dy_va%PY = ok_in_dy_va%PY &
                       + (oct_norm_stren*oct_len*weight/6._rac) * (3._rac*ok_in_dy_va%Y*ok_in_dy_va%X**2 - ok_in_dy_va%Y**3)
      ok_fi_dy_va%l = ok_in_dy_va%l
      ok_fi_dy_va%delta = ok_in_dy_va%delta

    end function OctuKick


    function OctuCorr (oc_in_dy_va, oct_len, oct_norm_stren, weight) result (oc_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: oct_len, oct_norm_stren, weight
      type (conjugate_variables), intent(in) :: oc_in_dy_va
      type (conjugate_variables) :: oc_fi_dy_va

      if (oct_norm_stren == 0.0_rac .or. oct_len <= 0.0_rac) then

        print*, "ERROR - function OctuCorr !!!"
        print*, "The stength of the octupoles can not be equal to zero !!!"
        print*, "The length of the octupoles must be greater than zero !!!"

      end if

      oc_fi_dy_va%X = oc_in_dy_va%X
      oc_fi_dy_va%PX = oc_in_dy_va%PX + (oct_norm_stren**2 * oct_len**3 * weight) &
                       * ((oc_in_dy_va%X*(oc_in_dy_va%X**2 + oc_in_dy_va%Y**2)**2) / (12._rac*(1._rac + oc_in_dy_va%delta)))
      oc_fi_dy_va%Y = oc_in_dy_va%Y
      oc_fi_dy_va%PY = oc_in_dy_va%PY + (oct_norm_stren**2 * oct_len**3 * weight) &
                       * ((oc_in_dy_va%Y*(oc_in_dy_va%X**2 + oc_in_dy_va%Y**2)**2) / (12._rac*(1._rac + oc_in_dy_va%delta)))
      oc_fi_dy_va%l = oc_in_dy_va%l + oct_norm_stren**2 * oct_len**3 * weight &
                      * ((oc_in_dy_va%X**2 + oc_in_dy_va%Y**2)**3)/(72._rac*(1._rac + oc_in_dy_va%delta)**2)
      oc_fi_dy_va%delta = oc_in_dy_va%delta

    end function OctuCorr
    
    
    function DecaKick (dk_in_dy_va, dec_len, dec_norm_stren, weight) result (dk_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: dec_len, dec_norm_stren, weight
      type (conjugate_variables), intent(in) :: dk_in_dy_va
      type (conjugate_variables) :: dk_fi_dy_va

      if (dec_norm_stren == 0.0_rac .or. dec_len <= 0.0_rac) then

        print*, "ERROR - function DecaKick !!!"
        print*, "The stength of the decapoles can not be equal to zero !!!"
        print*, "The length of the decapoles must be greater than zero !!!"

      end if

      dk_fi_dy_va%X = dk_in_dy_va%X
      dk_fi_dy_va%PX = dk_in_dy_va%PX - (dec_norm_stren*dec_len*weight/24._rac) &
                       * (dk_in_dy_va%X**4 - 6._rac*(dk_in_dy_va%X*dk_in_dy_va%Y)**2 + dk_in_dy_va%Y**4)
      dk_fi_dy_va%Y = dk_in_dy_va%Y
      dk_fi_dy_va%PY = dk_in_dy_va%PY + (dec_norm_stren*dec_len*weight/6._rac) &
                       * (dk_in_dy_va%Y*dk_in_dy_va%X**3 - dk_in_dy_va%X*dk_in_dy_va%Y**3)
      dk_fi_dy_va%l = dk_in_dy_va%l
      dk_fi_dy_va%delta = dk_in_dy_va%delta

    end function DecaKick
    
    function DodeKick (dok_in_dy_va, dod_len, dod_norm_stren, weight) result (dok_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: dod_len, dod_norm_stren, weight
      type (conjugate_variables), intent(in) :: dok_in_dy_va
      type (conjugate_variables) :: dok_fi_dy_va

      if (dod_norm_stren == 0.0_rac .or. dod_len <= 0.0_rac) then

        print*, "ERROR - function DodeKick !!!"
        print*, "The stength of the dodecapoles can not be equal to zero !!!"
        print*, "The length of the dodecapoles must be greater than zero !!!"

      end if

      dok_fi_dy_va%X = dok_in_dy_va%X
      dok_fi_dy_va%PX = dok_in_dy_va%PX - (dod_norm_stren*dod_len*weight/120._rac) * (dok_in_dy_va%X**5 &
                        - 10._rac*dok_in_dy_va%X**3*dok_in_dy_va%Y**2 + 5._rac*dok_in_dy_va%X*dok_in_dy_va%Y**4)
      dok_fi_dy_va%Y = dok_in_dy_va%Y
      dok_fi_dy_va%PY = dok_in_dy_va%PY + (dod_norm_stren*dod_len*weight/120._rac) * (dok_in_dy_va%Y**5 &
                        - 10._rac*dok_in_dy_va%X**2*dok_in_dy_va%Y**3 + 5._rac*dok_in_dy_va%Y*dok_in_dy_va%X**4)
      dok_fi_dy_va%l = dok_in_dy_va%l
      dok_fi_dy_va%delta = dok_in_dy_va%delta

    end function DodeKick
    
    
    function MultKick (mk_in_dy_va, mul_len, angle, k1, k2, k3, k4, k5, weight) result (mk_fi_dy_va)

      implicit none

      real (kind=rac), intent(in) :: mul_len, angle, k1, k2, k3, k4, k5, weight
      type (conjugate_variables), intent(in) :: mk_in_dy_va
      real (kind=rac) :: effective_len
      type (conjugate_variables) :: mk_fi_dy_va, st_fi_dy_va

      if (angle == 0.0_rac .and. k1 == 0.0_rac .and. k2 == 0.0_rac .and. &
		   k3 == 0.0_rac .and. k4 == 0.0_rac .and. k5 == 0.0_rac) then

        print*, "ERROR - function MultKick !!!"
        print*, "The stength of all the multipole components can not be equal to zero !!!"

      end if
      
      effective_len = mul_len
      
      if (mul_len < 0.0_rac) then
      
        print*, "ERROR - function MultKick !!!"
        print*, "The length of the multipole must be positive !!!"
      
      elseif (mul_len == 0.0_rac) then
      
		effective_len = 1.0_rac
		
	  endif
      
      st_fi_dy_va = mk_in_dy_va

      if (angle /= 0.0_rac .or. k1 /= 0.0_rac) then
		st_fi_dy_va = CSDIKick (st_fi_dy_va, effective_len, angle, k1, weight)
	  endif
	  if (k2 /= 0.0_rac) then
		st_fi_dy_va = SextKick (st_fi_dy_va, effective_len, k2, weight)
	  endif
	  if (k3 /= 0.0_rac) then
		st_fi_dy_va = OctuKick (st_fi_dy_va, effective_len, k3, weight)
	  endif
	  if (k4 /= 0.0_rac) then
		st_fi_dy_va = DecaKick (st_fi_dy_va, effective_len, k4, weight)
	  endif
	  if (k5 /= 0.0_rac) then
		st_fi_dy_va = DodeKick (st_fi_dy_va, effective_len, k5, weight)
	  endif
	  
	  mk_fi_dy_va = st_fi_dy_va

    end function MultKick

end module THIN_ELEMENTS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! module EXACT_INTEGRATOR !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module EXACT_INTEGRATOR

  use THICK_ELEMENTS
  use THIN_ELEMENTS
  use LATTICE_ATTRIBUTES

  implicit none

  contains

    function ExactMap (lattice_traits, ini_dyn_var, max_pipe_radius) result (ex_fi_dy_va)

      implicit none

      type (lattice_description), intent(in) :: lattice_traits
      type (conjugate_variables), dimension (:), intent(in) :: ini_dyn_var
      real (kind=rac), intent(in) :: max_pipe_radius
      type (conjugate_variables), dimension (size(ini_dyn_var,1)) :: ex_fi_dy_va
      character (len=4) :: element
      real (kind=rac) :: length, hkick, vkick, E1, E2, angle, k1, k2, k3, k4, k5, krf
      !type (conjugate_variables) :: st_d_fi_dy_va, st_k_fi_dy_va
      integer (kind=iac) :: n_particles, jj

      n_particles = size(ini_dyn_var,1)

      element = lattice_traits%element_name
      length = lattice_traits%element_length
      hkick = lattice_traits%hkick_normalised_strenqth
      vkick = lattice_traits%vkick_normalised_strenqth
      E1 = lattice_traits%dipole_pole_fase_rotation_in
      E2 = lattice_traits%dipole_pole_fase_rotation_out
      angle = lattice_traits%bending_angle
      k1 = lattice_traits%quadrupole_normalised_strenqth
      k2 = lattice_traits%sextupole_normalised_strenqth
      k3 = lattice_traits%octupole_normalised_strenqth
      k4 = lattice_traits%decapole_normalised_strenqth
      k5 = lattice_traits%dodecapole_normalised_strenqth
      krf = lattice_traits%rf_normalised_frequency

      if (element=="MARK") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              ex_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          ex_fi_dy_va(jj) = UnitaryTransformation (ini_dyn_var(jj))

        end do

      elseif (element=="DRIF") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              ex_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          ex_fi_dy_va(jj) = Drift (ini_dyn_var(jj), length)

        end do
        
      elseif (element=="HVCO") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              ex_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif
          
          if (length==0._rac) then
          
            ex_fi_dy_va(jj) = HoVeKick (ini_dyn_var(jj), length, hkick, vkick, 1.0_rac)
            
          else 
          
			if (hkick /= 0._rac .and. vkick == 0._rac) then
			
			  ex_fi_dy_va(jj) = HorizontalCorrector (ini_dyn_var(jj), length, hkick)
			  
			elseif (vkick /= 0._rac .and. hkick == 0._rac) then
			
			  ex_fi_dy_va(jj) = VerticalCorrector (ini_dyn_var(jj), length, vkick)
			
			endif
			
		  endif

        end do

      elseif (element=="CSDI") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              ex_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          ex_fi_dy_va(jj) = CombinedSectorDipole (ini_dyn_var(jj), length, angle, E1, E2, k1)

        end do
        
      elseif (element=="CRDI") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              ex_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          ex_fi_dy_va(jj) = CombinedRectangularDipole (ini_dyn_var(jj), length, angle, E1, E2, k1)

        end do

      elseif (element=="QUAD") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              ex_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          ex_fi_dy_va(jj) = Quadrupole (ini_dyn_var(jj), length, k1)

        end do
        
      elseif (element=="RFCA") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              ex_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          ex_fi_dy_va(jj) = RFcavite (ini_dyn_var(jj), length, krf, 3.199791632068303_rac*1E-4)  

        end do
        
      elseif (element=="MULT") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              ex_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              ex_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          ex_fi_dy_va(jj) = MultKick (ini_dyn_var(jj), length, angle, k1, k2, k3, k4, k5, 1.0_rac)

        end do

      elseif (element=="SEXT" .or. element=="OCTU") then

        print *, "ERROR - function ExactMap !!!"
        print *, "There is no exact solution for sextupoles or octupoles !!!"

      else

        print*, "ERROR - function ExactMap !!!"
        print*, element, " element is not supported !!!"

      end if

    end function ExactMap

end module EXACT_INTEGRATOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! module SIMPLE_INTEGRATOR !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module SIMPLE_INTEGRATOR

  use THIN_ELEMENTS
  use THICK_ELEMENTS
  use LATTICE_ATTRIBUTES

  implicit none

  type simple_n_strength_length_weights

    real (kind=rac) :: kick_w
    real (kind=rac) :: drift_w

  end type simple_n_strength_length_weights

  type (simple_n_strength_length_weights), private :: simple_n

  private :: AssignWeightSimple

  contains

    subroutine AssignWeightSimple (n_ste,n_kic,simp)

      implicit none

      integer (kind=iac), intent(in) :: n_ste
      integer (kind=iac), intent(out) :: n_kic
      type (simple_n_strength_length_weights), intent(out) :: simp

      n_kic = (n_ste-1)/2
      simp%kick_w = 1._rac/n_kic
      simp%drift_w = 1._rac/(n_kic + 1._rac)

    end subroutine AssignWeightSimple

    function Simple (n_steps, lattice_traits, ini_dyn_var, max_pipe_radius) result (si_fi_dy_va)

      implicit none

      integer (kind=iac), intent(in) :: n_steps
      type (lattice_description), intent(in) :: lattice_traits
      type (conjugate_variables), dimension (:), intent(in) :: ini_dyn_var
      real (kind=rac), intent(in) :: max_pipe_radius
      character (len=4) :: element
      real (kind=rac) :: length, hkick, vkick, E1, E2, angle, h_angle, k1, k2, k3, k4, k5, krf
      type (conjugate_variables), dimension (size(ini_dyn_var,1)) :: si_fi_dy_va
      type (conjugate_variables) :: st_d_fi_dy_va, st_k_fi_dy_va
      integer (kind=iac) :: n_particles, n_kicks, ii, jj

      call AssignWeightSimple (n_steps,n_kicks,simple_n)

      n_particles = size(ini_dyn_var,1)

      element = lattice_traits%element_name
      length = lattice_traits%element_length
      hkick = lattice_traits%hkick_normalised_strenqth
      vkick = lattice_traits%vkick_normalised_strenqth
      E1 = lattice_traits%dipole_pole_fase_rotation_in
      E2 = lattice_traits%dipole_pole_fase_rotation_out
      angle = lattice_traits%bending_angle
      k1 = lattice_traits%quadrupole_normalised_strenqth
      k2 = lattice_traits%sextupole_normalised_strenqth
      k3 = lattice_traits%octupole_normalised_strenqth
      k4 = lattice_traits%decapole_normalised_strenqth
      k5 = lattice_traits%dodecapole_normalised_strenqth
      krf = lattice_traits%rf_normalised_frequency

      if (element=="CSDI" .or. element=="CRDI") then

        h_angle = 0.5_rac*angle

        if (n_steps<3 .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print *, "ERROR - function Simple !!!"
          print *, "The number of steps must be greater or equal to 3 !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,... !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_k_fi_dy_va = DFRFKick (ini_dyn_var(jj), length, E1, angle, 1.0_rac) 
            if (element=="CRDI") then
			  st_k_fi_dy_va = DFRFKick (st_k_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, 0.5_rac)
            st_k_fi_dy_va = CSDIKick (st_d_fi_dy_va, length, angle, k1, 1.0_rac)
            st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, 0.5_rac)
            if (element=="CRDI") then
			  st_d_fi_dy_va = DFRFKick (st_d_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            si_fi_dy_va(jj) = DFRFKick (st_d_fi_dy_va, length, E2, angle, 1.0_rac) 

          end do

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

			st_k_fi_dy_va = DFRFKick (ini_dyn_var(jj), length, E1, angle, 1.0_rac)
			if (element=="CRDI") then
			  st_k_fi_dy_va = DFRFKick (st_k_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, simple_n%drift_w)
            st_k_fi_dy_va = CSDIKick (st_d_fi_dy_va, length, angle, k1, simple_n%kick_w)

            do ii=2, n_kicks

              st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, simple_n%drift_w)
              st_k_fi_dy_va = CSDIKick (st_d_fi_dy_va, length, angle, k1, simple_n%kick_w)

            end do

            st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, simple_n%drift_w)
            if (element=="CRDI") then
			  st_d_fi_dy_va = DFRFKick (st_d_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            si_fi_dy_va(jj) = DFRFKick (st_d_fi_dy_va, length, E2, angle, 1.0_rac) 

          end do

        end if

      elseif (element=="QUAD") then

        if (n_steps<3 .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print *, "ERROR - function Simple !!!"
          print *, "The number of steps must be greater or equal to 3 !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,... !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, 0.5_rac)
            st_k_fi_dy_va = QuadKick (st_d_fi_dy_va, length, k1, 1.0_rac)
            si_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, 0.5_rac)

          end do

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, simple_n%drift_w)
            st_k_fi_dy_va = QuadKick (st_d_fi_dy_va, length, k1, simple_n%kick_w)

            do ii=2, n_kicks

              st_d_fi_dy_va = DriftInt (st_k_fi_dy_va, length, simple_n%drift_w)
              st_k_fi_dy_va = QuadKick (st_d_fi_dy_va, length, k1, simple_n%kick_w)

            end do

            si_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, simple_n%drift_w)

          end do

        end if

      elseif (element=="SEXT") then

        if (n_steps<3 .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print*, "ERROR - function Simple !!!"
          print *, "The number of steps must be greater or equal to 3 !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,... !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, 0.5_rac)
            st_k_fi_dy_va = SextKick (st_d_fi_dy_va, length, k2, 1.0_rac)
            si_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, 0.5_rac)

          end do

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, simple_n%drift_w)
            st_k_fi_dy_va = SextKick (st_d_fi_dy_va, length, k2, simple_n%kick_w)

            do ii=2, n_kicks

              st_d_fi_dy_va = DriftInt (st_k_fi_dy_va, length, simple_n%drift_w)
              st_k_fi_dy_va = SextKick (st_d_fi_dy_va, length, k2, simple_n%kick_w)

            end do

            si_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, simple_n%drift_w)

          end do

        end if

      elseif (element=="OCTU") then

        if (n_steps<3 .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print*, "ERROR - function Simple !!!"
          print *, "The number of steps must be greater or equal to 3 !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,... !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, 0.5_rac)
            st_k_fi_dy_va = OctuKick (st_d_fi_dy_va, length, k3, 1.0_rac)
            si_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, 0.5_rac)

          end do

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, simple_n%drift_w)
            st_k_fi_dy_va = OctuKick (st_d_fi_dy_va, length, k3, simple_n%kick_w)

            do ii=2, n_kicks

              st_d_fi_dy_va = DriftInt (st_k_fi_dy_va, length, simple_n%drift_w)
              st_k_fi_dy_va = OctuKick (st_d_fi_dy_va, length, k3, simple_n%kick_w)

            end do

            si_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, simple_n%drift_w)

          end do

        end if
        
      elseif (element=="MULT") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          si_fi_dy_va(jj) = MultKick (ini_dyn_var(jj), length, angle, k1, k2, k3, k4, k5, 1.0_rac)

        end do

      elseif (element=="MARK") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          si_fi_dy_va(jj) = UnitaryTransformation (ini_dyn_var(jj))

        end do

      elseif (element=="DRIF") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              si_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              si_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          si_fi_dy_va(jj) = Drift (ini_dyn_var(jj), length)

        end do

      else

        print*, "ERROR - function Simple !!!"
        print*, element, " element is not supported !!!"

      end if

    end function Simple

end module SIMPLE_INTEGRATOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! module TEAPOT_INTEGRATOR !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module TEAPOT_INTEGRATOR

  use THIN_ELEMENTS
  use THICK_ELEMENTS
  use LATTICE_ATTRIBUTES

  implicit none

  type teapot_n_strength_length_weights

    real (kind=rac) :: kick_w
    real (kind=rac) :: edge_drift_w
    real (kind=rac) :: center_drift_w

  end type teapot_n_strength_length_weights

  type (teapot_n_strength_length_weights), private :: teapot_n

  private :: AssignWeightTeapot

  contains

    subroutine AssignWeightTeapot (n_ste,n_kic,teap)

      implicit none

      integer (kind=iac), intent(in) :: n_ste
      integer (kind=iac), intent(out) :: n_kic
      type (teapot_n_strength_length_weights), intent(out) :: teap

      n_kic = (n_ste-1)/2
      teap%kick_w = 1._rac/n_kic
      teap%edge_drift_w = 1._rac/(2._rac + 2._rac*n_kic)
      teap%center_drift_w = n_kic/(-1._rac + n_kic**2)

    end subroutine AssignWeightTeapot

    function Teapot (n_steps, lattice_traits, ini_dyn_var, max_pipe_radius) result (te_fi_dy_va)

      implicit none

      integer (kind=iac), intent(in) :: n_steps
      type (lattice_description), intent(in) :: lattice_traits
      type (conjugate_variables), dimension(:), intent(in) :: ini_dyn_var
      real (kind=rac), intent(in) :: max_pipe_radius
      character (len=4) :: element
      real (kind=rac) :: length, hkick, vkick, E1, E2, angle, h_angle, k1, k2, k3, k4, k5, krf
      type (conjugate_variables), dimension(size(ini_dyn_var,1)) :: te_fi_dy_va
      type (conjugate_variables) :: st_d_fi_dy_va, st_k_fi_dy_va
      integer (kind=iac) :: n_particles, n_kicks, ii, jj

      call AssignWeightTeapot (n_steps,n_kicks,teapot_n)

      n_particles = size(ini_dyn_var,1)

      element = lattice_traits%element_name
      length = lattice_traits%element_length
      hkick = lattice_traits%hkick_normalised_strenqth
      vkick = lattice_traits%vkick_normalised_strenqth
      E1 = lattice_traits%dipole_pole_fase_rotation_in
      E2 = lattice_traits%dipole_pole_fase_rotation_out
      angle = lattice_traits%bending_angle
      k1 = lattice_traits%quadrupole_normalised_strenqth
      k2 = lattice_traits%sextupole_normalised_strenqth
      k3 = lattice_traits%octupole_normalised_strenqth
      k4 = lattice_traits%decapole_normalised_strenqth
      k5 = lattice_traits%dodecapole_normalised_strenqth
      krf = lattice_traits%rf_normalised_frequency

      if (element=="CSDI" .or. element=="CRDI") then

		h_angle = 0.5_rac*angle

        if (n_steps<3 .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print*, "ERROR - function Teapot !!!"
          print *, "The number of steps must be greater or equal to 3 !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,... !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

			st_k_fi_dy_va = DFRFKick (ini_dyn_var(jj), length, E1, angle, 1.0_rac) 
            if (element=="CRDI") then
			  st_k_fi_dy_va = DFRFKick (st_k_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, 0.5_rac)
            st_k_fi_dy_va = CSDIKick (st_d_fi_dy_va, length, angle, k1, 1.0_rac)
            st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, 0.5_rac)
            if (element=="CRDI") then
			  st_d_fi_dy_va = DFRFKick (st_d_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            te_fi_dy_va(jj) = DFRFKick (st_d_fi_dy_va, length, E2, angle, 1.0_rac) 

          end do

          print *, "There is no TEAPOT that consists of 3 steps. Leapfrog (DKD) is used instead !!!"

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

			st_k_fi_dy_va = DFRFKick (ini_dyn_var(jj), length, E1, angle, 1.0_rac) 
            if (element=="CRDI") then
			  st_k_fi_dy_va = DFRFKick (st_k_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, teapot_n%edge_drift_w)
            st_k_fi_dy_va = CSDIKick (st_d_fi_dy_va, length, angle, k1, teapot_n%kick_w)

            do ii=2, n_kicks

              st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, teapot_n%center_drift_w)
              st_k_fi_dy_va = CSDIKick (st_d_fi_dy_va, length, angle, k1, teapot_n%kick_w)

            end do

            st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, teapot_n%edge_drift_w)
            if (element=="CRDI") then
			  st_d_fi_dy_va = DFRFKick (st_d_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            te_fi_dy_va(jj) = DFRFKick (st_d_fi_dy_va, length, E2, angle, 1.0_rac) 

          end do

        end if

      elseif (element=="QUAD") then

        if (n_steps<3 .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print*, "ERROR - function Teapot !!!"
          print *, "The number of steps must be greater or equal to 3 !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,... !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, 0.5_rac)
            st_k_fi_dy_va = QuadKick (st_d_fi_dy_va, length, k1, 1.0_rac)
            te_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, 0.5_rac)

          end do

          print *, "There is no TEAPOT that consists of 3 steps. Leapfrog (DKD) is used instead !!!"

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, teapot_n%edge_drift_w)
            st_k_fi_dy_va = QuadKick (st_d_fi_dy_va, length, k1, teapot_n%kick_w)

            do ii=2, n_kicks

              st_d_fi_dy_va = DriftInt (st_k_fi_dy_va, length, teapot_n%center_drift_w)
              st_k_fi_dy_va = QuadKick (st_d_fi_dy_va, length, k1, teapot_n%kick_w)

            end do

            te_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, teapot_n%edge_drift_w)

          end do

        end if

      elseif (element=="SEXT") then

        if (n_steps<3 .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print*, "ERROR - function Teapot !!!"
          print *, "The number of steps must be greater or equal to 3 !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,... !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, 0.5_rac)
            st_k_fi_dy_va = SextKick (st_d_fi_dy_va, length, k2, 1.0_rac)
            te_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, 0.5_rac)

          end do

          print *, "There is no TEAPOT that consists of 3 steps. Leapfrog (DKD) is used instead !!!"

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, teapot_n%edge_drift_w)
            st_k_fi_dy_va = SextKick (st_d_fi_dy_va, length, k2, teapot_n%kick_w)

            do ii=2, n_kicks

              st_d_fi_dy_va = DriftInt (st_k_fi_dy_va, length, teapot_n%center_drift_w)
              st_k_fi_dy_va = SextKick (st_d_fi_dy_va, length, k2, teapot_n%kick_w)

            end do

            te_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, teapot_n%edge_drift_w)

          end do

        end if

      elseif (element=="OCTU") then

        if (n_steps<3 .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print*, "ERROR - function Teapot !!!"
          print *, "The number of steps must be greater or equal to 3 !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,... !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, 0.5_rac)
            st_k_fi_dy_va = OctuKick (st_d_fi_dy_va, length, k3, 1.0_rac)
            te_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, 0.5_rac)

          end do

          print *, "There is no TEAPOT that consists of 3 steps. Leapfrog (DKD) is used instead !!!"

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, teapot_n%edge_drift_w)
            st_k_fi_dy_va = OctuKick (st_d_fi_dy_va, length, k3, teapot_n%kick_w)

            do ii=2, n_kicks

              st_d_fi_dy_va = DriftInt (st_k_fi_dy_va, length, teapot_n%center_drift_w)
              st_k_fi_dy_va = OctuKick (st_d_fi_dy_va, length, k3, teapot_n%kick_w)

            end do

            te_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, teapot_n%edge_drift_w)

          end do

        end if
        
      elseif (element=="MULT") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          te_fi_dy_va(jj) = MultKick (ini_dyn_var(jj), length, angle, k1, k2, k3, k4, k5, 1.0_rac)

        end do

      elseif (element=="MARK") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          te_fi_dy_va(jj) = UnitaryTransformation (ini_dyn_var(jj))

        end do

      elseif (element=="DRIF") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              te_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              te_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          te_fi_dy_va(jj) = Drift (ini_dyn_var(jj), length)

        end do

      else

        print*, "ERROR - function Teapot !!!"
        print*, element, " element is not supported !!!"

      end if

    end function Teapot

end module TEAPOT_INTEGRATOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! module SABAC_INTEGRATOR !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module SABAC_INTEGRATOR

  use THIN_ELEMENTS
  use THICK_ELEMENTS
  use LATTICE_ATTRIBUTES

  implicit none

  integer (kind=iac), private :: sabac_max_steps=23

  type sabac_n_strength_length_weights

    ! up to the SABAC10 (12,23) ~~> (12 kicks --> 10 kicks + 2 correctors, 23 total maps)
    real (kind=rac) :: corrector_w
    real (kind=rac), dimension (10) :: kick_w
    real (kind=rac), dimension (11) :: drift_w

  end type sabac_n_strength_length_weights

  type (sabac_n_strength_length_weights), dimension (10), private :: sabac_n

  private :: AssignWeightSabac

  contains

    subroutine AssignWeightSabac (n_ste,n_kic,sabac)

      implicit none

      integer (kind=iac), intent(in) :: n_ste
      integer (kind=iac), intent(out) :: n_kic
      type (sabac_n_strength_length_weights), dimension (10), intent(out) :: sabac

      !(number of kicks without taking into account the correctors)
      n_kic = (n_ste - 3)/2

      ! SABAC1 (3,5)
      if (n_ste == 5) then

        sabac(1)%corrector_w = 1._rac/12._rac
        sabac(1)%kick_w(1) = 1._rac
        sabac(1)%kick_w(2:) = 0._rac
        sabac(1)%drift_w(1:2) = 1._rac/2._rac
        sabac(1)%drift_w(3:) = 0._rac

      ! SABAC2 (4,7)
      elseif (n_ste == 7) then

        sabac(2)%corrector_w = (2._rac - sqrt(3._rac))/24._rac
        sabac(2)%kick_w(1:2) = 1._rac/2._rac
        sabac(2)%kick_w(3:) = 0._rac
        sabac(2)%drift_w(1:3:2) = 1._rac/2._rac - sqrt(3._rac)/6._rac
        sabac(2)%drift_w(2) = sqrt(3._rac)/3._rac
        sabac(2)%drift_w(4:) = 0._rac

      ! SABAC3 (5,9)
      elseif (n_ste == 9) then

        sabac(3)%corrector_w = (54._rac - 13._rac*sqrt(15._rac))/648._rac
        sabac(3)%kick_w(1:3:2) = 5._rac/18._rac
        sabac(3)%kick_w(2) = 4._rac/9._rac
        sabac(3)%kick_w(4:) = 0._rac
        sabac(3)%drift_w(1:4:3) = 1._rac/2._rac - sqrt(15._rac)/10._rac
        sabac(3)%drift_w(2:3) = sqrt(15._rac)/10._rac
        sabac(3)%drift_w(5:) = 0._rac

      ! SABAC4 (6,11)
      elseif (n_ste == 11) then

        sabac(4)%corrector_w = 0.003396775048208601331532157783492144_rac
        sabac(4)%kick_w(1:4:3) = 1._rac/4._rac - sqrt(30._rac)/72._rac
        sabac(4)%kick_w(2:3) = 1._rac/4._rac + sqrt(30._rac)/72._rac
        sabac(4)%kick_w(5:) = 0._rac
        sabac(4)%drift_w(1:5:4) = 1._rac/2._rac - sqrt(525._rac + 70._rac*sqrt(30._rac))/70._rac
        sabac(4)%drift_w(2:4:2) = (sqrt(525._rac + 70._rac*sqrt(30._rac)) - sqrt(525._rac - 70._rac*sqrt(30._rac)))/70._rac
        sabac(4)%drift_w(3) = sqrt(525._rac - 70._rac*sqrt(30._rac))/35._rac
        sabac(4)%drift_w(6:) = 0._rac

      ! SABAC5 (7,13)
      elseif (n_ste == 13) then

        sabac(5)%corrector_w = 0.002270543121419264819434955050039130_rac
        sabac(5)%kick_w(1:5:4) = (322._rac - 13._rac*sqrt(70._rac))/1800._rac
        sabac(5)%kick_w(2:4:2) = (322._rac + 13._rac*sqrt(70._rac))/1800._rac
        sabac(5)%kick_w(3) = 64._rac/225._rac
        sabac(5)%kick_w(6:) = 0._rac
        sabac(5)%drift_w(1:6:5) = 1._rac/2._rac &
                                  - (sqrt(490._rac + 42._rac*sqrt(105._rac)) +sqrt(490._rac - 42._rac*sqrt(105._rac)))/84._rac
        sabac(5)%drift_w(2:5:3) = sqrt(490._rac - 42._rac*sqrt(105._rac))/42._rac
        sabac(5)%drift_w(3:4) = (sqrt(490._rac + 42._rac*sqrt(105._rac)) - sqrt(490._rac - 42._rac*sqrt(105._rac)))/84._rac
        sabac(5)%drift_w(7:) = 0._rac

      ! SABAC6 (8,15)
      elseif (n_ste == 15) then

        sabac(6)%corrector_w = 0.001624459841624282521452258512463608_rac
        sabac(6)%kick_w(1:6:5) = 0.085662246189585172520148071086366447_rac
        sabac(6)%kick_w(2:5:3) = 0.180380786524069303784916756918858056_rac
        sabac(6)%kick_w(3:4) = 0.233956967286345523694935171994775497_rac
        sabac(6)%kick_w(7:) = 0._rac
        sabac(6)%drift_w(1:7:6) = 0.033765242898423986093849222753002695_rac
        sabac(6)%drift_w(2:6:4) = 0.135630063868443757075450979737044631_rac
        sabac(6)%drift_w(3:5:2) = 0.211295100191533802515448936669596706_rac
        sabac(6)%drift_w(4) = 0.238619186083196908630501721680711935_rac
        sabac(6)%drift_w(8:) = 0._rac

      ! SABAC7 (9,17)
      elseif (n_ste == 17) then

        sabac(7)%corrector_w = 0.001219643912760418472579211822331645_rac
        sabac(7)%kick_w(1:7:6) = 0.064742483084434846635305716339541009_rac
        sabac(7)%kick_w(2:6:4) = 0.139852695744638333950733885711889791_rac
        sabac(7)%kick_w(3:5:2) = 0.190915025252559472475184887744487567_rac
        sabac(7)%kick_w(4) = 256._rac/1225._rac
        sabac(7)%kick_w(8:) = 0._rac
        sabac(7)%drift_w(1:8:7) = 0.025446043828620737736905157976074369_rac
        sabac(7)%drift_w(2:7:5) = 0.103788363371682042331162455383531428_rac
        sabac(7)%drift_w(3:6:3) = 0.167843017110998636478629180601913472_rac
        sabac(7)%drift_w(4:5) = 0.202922575688698583453303206038480732_rac
        sabac(7)%drift_w(9:) = 0._rac

      ! SABAC8 (10,19)
      elseif (n_ste == 19) then

        sabac(8)%corrector_w = 0.000949308177745602234792177503535054_rac
        sabac(8)%kick_w(1:8:7) = 0.050614268145188129576265677154981095_rac
        sabac(8)%kick_w(2:7:5) = 0.111190517226687235272177997213120442_rac
        sabac(8)%kick_w(3:6:3) = 0.156853322938943643668981100993300657_rac
        sabac(8)%kick_w(4:5) = 0.18134189168918099148257522463859781_rac
        sabac(8)%kick_w(9:) = 0._rac
        sabac(8)%drift_w(1:9:8) = 0.019855071751231884158219565715263505_rac
        sabac(8)%drift_w(2:8:6) = 0.081811689541954746046003466046821277_rac
        sabac(8)%drift_w(3:7:4) = 0.135567033748648876886907443643292044_rac
        sabac(8)%drift_w(4:6:2) = 0.171048883710339590439131453414531184_rac
        sabac(8)%drift_w(5) = 0.183434642495649804939476142360183981_rac
        sabac(8)%drift_w(10:) = 0._rac

      ! SABAC9 (11,21)
      elseif (n_ste == 21) then

        sabac(9)%corrector_w = 0.000759846022860436646358196674176815_rac
        sabac(9)%kick_w(1:9:8) = 0.040637194180787205985946079055261825_rac
        sabac(9)%kick_w(2:8:6) = 0.090324080347428702029236015621456405_rac
        sabac(9)%kick_w(3:7:4) = 0.130305348201467731159371434709316425_rac
        sabac(9)%kick_w(4:6:2) = 0.156173538520001420034315203292221833_rac
        sabac(9)%kick_w(5) = 16384._rac/99225._rac
        sabac(9)%kick_w(10:) = 0._rac
        sabac(9)%drift_w(1:10:9) = 0.015919880246186955082211898548163565_rac
        sabac(9)%drift_w(2:9:7) = 0.066064566090495147768073207416968997_rac
        sabac(9)%drift_w(3:8:5) = 0.111329837313022698495363874364130346_rac
        sabac(9)%drift_w(4:7:3) = 0.144559004648390734135082012349068788_rac
        sabac(9)%drift_w(5:6) = 0.162126711701904464519269007321668304_rac
        sabac(9)%drift_w(11:) = 0._rac

      ! SABAC10 (12,23)
      elseif (n_ste == 23) then

        sabac(10)%corrector_w = 0.000621934331486166426497049845358646_rac
        sabac(10)%kick_w(1:10:9) = 0.033335672154344068796784404946665896_rac
        sabac(10)%kick_w(2:9:7) = 0.074725674575290296572888169828848666_rac
        sabac(10)%kick_w(3:8:5) = 0.109543181257991021997767467114081596_rac
        sabac(10)%kick_w(4:7:3) = 0.134633359654998177545613460784734677_rac
        sabac(10)%kick_w(5:6) = 0.147762112357376435086946497325669165_rac
        sabac(10)%drift_w(1:11:10) = 0.013046735741414139961017993957773973_rac
        sabac(10)%drift_w(2:10:8) = 0.054421580914093604672933661830479502_rac
        sabac(10)%drift_w(3:9:6) = 0.092826899194980052248884661654309736_rac
        sabac(10)%drift_w(4:8:4) = 0.123007087084888607717530710974544707_rac
        sabac(10)%drift_w(5:7:2) = 0.142260527573807989957219971018032089_rac
        sabac(10)%drift_w(6) = 0.148874338981631210884826001129719985_rac

      endif

    end subroutine AssignWeightSabac

    function Sabac (n_steps, lattice_traits, ini_dyn_var, max_pipe_radius) result (sac_fi_dy_va)

      implicit none

      integer (kind=iac), intent(in) :: n_steps
      type (lattice_description), intent(in) :: lattice_traits
      type (conjugate_variables), dimension(:), intent(in) :: ini_dyn_var
      real (kind=rac), intent(in) :: max_pipe_radius
      character (len=4) :: element
      real (kind=rac) :: length, hkick, vkick, E1, E2, angle, h_angle, k1, k2, k3, k4, k5, krf
      type (conjugate_variables), dimension(size(ini_dyn_var,1)) :: sac_fi_dy_va
      type (conjugate_variables) :: st_d_fi_dy_va, st_k_fi_dy_va, st_c_fi_dy_va
      integer (kind=iac) :: n_particles, n_kicks, ii, jj

      call AssignWeightSabac (n_steps,n_kicks,sabac_n)

      n_particles = size(ini_dyn_var,1)

      element = lattice_traits%element_name
      length = lattice_traits%element_length
      hkick = lattice_traits%hkick_normalised_strenqth
      vkick = lattice_traits%vkick_normalised_strenqth
      E1 = lattice_traits%dipole_pole_fase_rotation_in
      E2 = lattice_traits%dipole_pole_fase_rotation_out
      angle = lattice_traits%bending_angle
      k1 = lattice_traits%quadrupole_normalised_strenqth
      k2 = lattice_traits%sextupole_normalised_strenqth
      k3 = lattice_traits%octupole_normalised_strenqth
      k4 = lattice_traits%decapole_normalised_strenqth
      k5 = lattice_traits%dodecapole_normalised_strenqth
      krf = lattice_traits%rf_normalised_frequency

      if (element=="CSDI" .or. element=="CRDI") then

		h_angle = 0.5_rac*angle

        if (n_steps<3 .or. n_steps>sabac_max_steps .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print*, "ERROR - function Sabac !!!"
          print *, "The number of steps must be greater or equal to 3 and smaller or equal to", sabac_max_steps," !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,...,",sabac_max_steps," !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

			st_k_fi_dy_va = DFRFKick (ini_dyn_var(jj), length, E1, angle, 1.0_rac) 
            if (element=="CRDI") then
			  st_k_fi_dy_va = DFRFKick (st_k_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, 0.5_rac)
            st_k_fi_dy_va = CSDIKick (st_d_fi_dy_va, length, angle, k1, 1.0_rac)
            st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, 0.5_rac)
            if (element=="CRDI") then
			  st_d_fi_dy_va = DFRFKick (st_d_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            sac_fi_dy_va(jj) = DFRFKick (st_d_fi_dy_va, length, E2, angle, 1.0_rac)

          end do

          print *, "There is no SABAC that consists of 3 steps. SABA1 which it is identical to Leapfrog (DKD) is used instead !!!"

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

			st_k_fi_dy_va = DFRFKick (ini_dyn_var(jj), length, E1, angle, 1.0_rac) 
            if (element=="CRDI") then
			  st_k_fi_dy_va = DFRFKick (st_k_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            st_c_fi_dy_va = CSDICorr (st_k_fi_dy_va, length, angle, k1, sabac_n(n_kicks)%corrector_w)
            st_d_fi_dy_va = CSDIDriftInt (st_c_fi_dy_va, length, angle, sabac_n(n_kicks)%drift_w(1))

            do ii=1, n_kicks

              st_k_fi_dy_va = CSDIKick (st_d_fi_dy_va, length, angle, k1, sabac_n(n_kicks)%kick_w(ii))
              st_d_fi_dy_va = CSDIDriftInt (st_k_fi_dy_va, length, angle, sabac_n(n_kicks)%drift_w(ii+1))

            end do

            st_c_fi_dy_va = CSDICorr (st_d_fi_dy_va, length, angle, k1, sabac_n(n_kicks)%corrector_w)
            if (element=="CRDI") then
			  st_c_fi_dy_va = DFRFKick (st_c_fi_dy_va, length, h_angle, angle, 1.0_rac) 
			endif
            sac_fi_dy_va(jj) = DFRFKick (st_c_fi_dy_va, length, E2, angle, 1.0_rac)

          end do

        end if

      elseif (element=="QUAD") then

        if (n_steps<3 .or. n_steps>sabac_max_steps .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print*, "ERROR - function Sabac !!!"
          print *, "The number of steps must be greater or equal to 3 and smaller or equal to", sabac_max_steps," !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,...,",sabac_max_steps," !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, 0.5_rac)
            st_k_fi_dy_va = QuadKick (st_d_fi_dy_va, length, k1, 1.0_rac)
            sac_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, 0.5_rac)

          end do

          print *, "There is no SABAC that consists of 3 steps. SABA1 which it is identical to Leapfrog (DKD) is used instead !!!"

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_c_fi_dy_va = QuadCorr (ini_dyn_var(jj), length, k1, sabac_n(n_kicks)%corrector_w)
            st_d_fi_dy_va = DriftInt (st_c_fi_dy_va, length, sabac_n(n_kicks)%drift_w(1))

            do ii=1, n_kicks

              st_k_fi_dy_va = QuadKick (st_d_fi_dy_va, length, k1, sabac_n(n_kicks)%kick_w(ii))
              st_d_fi_dy_va = DriftInt (st_k_fi_dy_va, length, sabac_n(n_kicks)%drift_w(ii+1))

            end do

            sac_fi_dy_va(jj) = QuadCorr (st_d_fi_dy_va, length, k1, sabac_n(n_kicks)%corrector_w)

          end do

        end if

      elseif (element=="SEXT") then

        if (n_steps<3 .or. n_steps>sabac_max_steps .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print*, "ERROR - function Sabac !!!"
          print *, "The number of steps must be greater or equal to 3 and smaller or equal to", sabac_max_steps," !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,...,", sabac_max_steps,"  !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, 0.5_rac)
            st_k_fi_dy_va = SextKick (st_d_fi_dy_va, length, k2, 1.0_rac)
            sac_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, 0.5_rac)

          end do

          print *, "There is no SABAC that consists of 3 steps. SABA1 which it is identical to Leapfrog (DKD) is used instead !!!"

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_c_fi_dy_va = SextCorr (ini_dyn_var(jj), length, k2, sabac_n(n_kicks)%corrector_w)
            st_d_fi_dy_va = DriftInt (st_c_fi_dy_va, length, sabac_n(n_kicks)%drift_w(1))

            do ii=1, n_kicks

              st_k_fi_dy_va = SextKick (st_d_fi_dy_va, length, k2, sabac_n(n_kicks)%kick_w(ii))
              st_d_fi_dy_va = DriftInt (st_k_fi_dy_va, length, sabac_n(n_kicks)%drift_w(ii+1))

            end do

            sac_fi_dy_va(jj) = SextCorr (st_d_fi_dy_va, length, k2, sabac_n(n_kicks)%corrector_w)

          end do

        end if

      elseif (element=="OCTU") then

        if (n_steps<3 .or. n_steps>sabac_max_steps .or. modulo(n_steps,2)/=1 .or. n_particles <= 0) then

          print*, "ERROR - function Sabac !!!"
          print *, "The number of steps must be greater or equal to 3 and smaller or equal to", sabac_max_steps," !!!"
          print *, "The number of steps must be an odd integer number equal to 3,5,7,...,", sabac_max_steps," !!!"
          print *, "The number of particles must be an integer number greater than or equal to 1 !!!"

        elseif (n_steps == 3) then

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_d_fi_dy_va = DriftInt (ini_dyn_var(jj), length, 0.5_rac)
            st_k_fi_dy_va = OctuKick (st_d_fi_dy_va, length, k3, 1.0_rac)
            sac_fi_dy_va(jj) = DriftInt (st_k_fi_dy_va, length, 0.5_rac)

          end do

          print *, "There is no SABAC that consists of 3 steps. SABA1 which it is identical to Leapfrog (DKD) is used instead !!!"

        else

          do jj=1, n_particles

            if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
                isnan(ini_dyn_var(jj)%X)) then

                sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
                sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

                cycle

            endif

            st_c_fi_dy_va = OctuCorr (ini_dyn_var(jj), length, k3, sabac_n(n_kicks)%corrector_w)
            st_d_fi_dy_va = DriftInt (st_c_fi_dy_va, length, sabac_n(n_kicks)%drift_w(1))

            do ii=1, n_kicks

              st_k_fi_dy_va = OctuKick (st_d_fi_dy_va, length, k3, sabac_n(n_kicks)%kick_w(ii))
              st_d_fi_dy_va = DriftInt (st_k_fi_dy_va, length, sabac_n(n_kicks)%drift_w(ii+1))

            end do

            sac_fi_dy_va(jj) = OctuCorr (st_d_fi_dy_va, length, k3, sabac_n(n_kicks)%corrector_w)

          end do

        end if
        
      elseif (element=="MULT") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          sac_fi_dy_va(jj) = MultKick (ini_dyn_var(jj), length, angle, k1, k2, k3, k4, k5, 1.0_rac)

        end do

      elseif (element=="MARK") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          sac_fi_dy_va(jj) = UnitaryTransformation (ini_dyn_var(jj))

        end do

      elseif (element=="DRIF") then

        do jj=1, n_particles

          if (abs(ini_dyn_var(jj)%X)>=max_pipe_radius .or. abs(ini_dyn_var(jj)%Y)>=max_pipe_radius .or. &
              isnan(ini_dyn_var(jj)%X)) then

              sac_fi_dy_va(jj)%X = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%PX = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%Y = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%PY = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%l = transfer((huge(1_rac)), 1.0_rac)
              sac_fi_dy_va(jj)%delta = transfer((huge(1_rac)), 1.0_rac)

              cycle

          endif

          sac_fi_dy_va(jj) = Drift (ini_dyn_var(jj), length)

        end do

      else

        print*, "ERROR - function Sabac !!!"
        print*, element, " element is not supported !!!"

      end if

    end function Sabac

end module SABAC_INTEGRATOR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!! module TRACKING_DATA !!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module TRACKING_DATA

  use ACCURACY_CONSTANTS_STRUCTURES
  use LATTICE_ATTRIBUTES

  implicit none

  contains

    subroutine AllocateInitialize (total_data_lines, particles_initial_dy_va, data_array)

      implicit none

      integer (kind=iac), intent(in) :: total_data_lines
      type (conjugate_variables), dimension(:), intent(in) :: particles_initial_dy_va
      type (conjugate_variables), dimension(:), allocatable, intent(out) :: data_array

      ! Allocate (set array size)
      allocate(data_array(((total_data_lines+1)*size(particles_initial_dy_va,1))))

      ! Iitialize with NaN
      data_array(:)%X = transfer((huge(1_rac)), 1.0_rac)
      data_array(:)%PX = transfer((huge(1_rac)), 1.0_rac)
      data_array(:)%Y = transfer((huge(1_rac)), 1.0_rac)
      data_array(:)%PY = transfer((huge(1_rac)), 1.0_rac)
      data_array(:)%l = transfer((huge(1_rac)), 1.0_rac)
      data_array(:)%delta = transfer((huge(1_rac)), 1.0_rac)

      ! Save the initial contition for each particle
      data_array(1::total_data_lines+1) = particles_initial_dy_va(1:)

      print*, 'The allocation and the inialization is complited normally.'

    end subroutine AllocateInitialize


    subroutine StoreData (element_id_interest, current_element, current_turn_serial_number, turn_start_store, turn_stop_store, &
                          turn_step_store, current_dy_var, current_data_array, save_data_in_file, binary_file_name)

      implicit none

      character (len=*), intent(in) :: element_id_interest, binary_file_name
      type (lattice_description), intent(in) :: current_element
      integer (kind=iac), intent(in) :: current_turn_serial_number, turn_start_store, turn_stop_store, turn_step_store
      type (conjugate_variables), dimension(:), intent(in) :: current_dy_var
      type (conjugate_variables), dimension(:), intent(out) :: current_data_array
      logical, intent(in) :: save_data_in_file
      integer (kind=iac) :: n_tracked_lines, line_id
      integer :: io


      if ( current_element%element_id(1:len(element_id_interest)) == element_id_interest .and. &
           current_turn_serial_number >= turn_start_store .and. current_turn_serial_number <= turn_stop_store ) then

        n_tracked_lines = (size(current_data_array,1)-size(current_dy_var,1))/size(current_dy_var,1)

        if (current_turn_serial_number == turn_start_store .or. modulo(current_turn_serial_number,turn_step_store) == 0_iac ) then

          if (current_turn_serial_number == turn_start_store) then

            line_id = 1_iac

          else

            line_id = (current_turn_serial_number - turn_start_store)/turn_step_store + 1_iac

          endif

          current_data_array(line_id+1::n_tracked_lines+1) = current_dy_var(1:)

          if (current_turn_serial_number == turn_stop_store .and. save_data_in_file) then

            open (unit=19, file=binary_file_name, form='unformatted', access='stream', status='replace')

            write (unit=19, iostat=io) current_data_array

            if (io > 0) then

              print*, "ERROR - subroutine StoreData !!!"
              print*, "Some error occured while writing the file ", binary_file_name, " !!!"

            else

              print*, "The file ", binary_file_name, " is written normally."

            endif

            close (unit=19)

          end if

        end if

      end if

    end subroutine StoreData


end module  TRACKING_DATA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program TRACKING

  use ACCURACY_CONSTANTS_STRUCTURES ! real (kind=rac) & type (conjugate_variables)
  use EXACT_INTEGRATOR ! ExactMap (element, norm_strength, length, ini_dyn_var)
  use SIMPLE_INTEGRATOR ! Simple (n_steps, element, norm_strength, length, ini_dyn_var)
  use TEAPOT_INTEGRATOR ! Teapot (n_steps, element, norm_strength, length, ini_dyn_var)
  use SABAC_INTEGRATOR ! Sabac (n_steps, element, norm_strength, length, ini_dyn_var)
  use LATTICE_PARTICLES_INITIATION ! [lattice & call GetLattice (unit, file_name) in the main progra always use character (len=60) :: lattice_file_name = "name" and deallocate (lattice)]
                                   ! [particle_in_dy_va & call GetParticlesInitialDistribution (unit, file_name) in the main progra always use character (len=60) :: paricles_file_name = "name" and deallocate (particle_in_dy_va)]
  use TRACKING_DATA

  implicit none

  real (kind=rac) :: energy, gamma_rel, beta_rel
  real (kind=rac) :: emittance_norm_x, emittance_geo_x, emittance_norm_y, emittance_geo_y, beta_x_start, beta_y_start
  real (kind=rac) :: min_angle_deg, max_angle_deg, min_sigma, max_sigma, sigma, sigma_x, sigma_y, delta_in
  integer (kind=iac) :: n_angle_slices, n_sigma_slices
  real (kind=rac) :: beam_pipe_radius

  character (len=100) :: lattice_file_name= "LHC_like_lattice_data.txt" !"LHC_like_lattice_data.txt" "lattice_data.txt"
  character (len=100) :: particles_initial_condition_file_name="particles_initial_conditions_data.txt"

  type (conjugate_variables), dimension(:), allocatable :: in_dy_va, st_dy_va
  type (conjugate_variables), dimension(:), allocatable :: fi_ex_dy_va, fi_sac_dy_va, fi_te_dy_va
  type (conjugate_variables), dimension(:), allocatable :: tbt_dy_va_ex, tbt_dy_va_sac, tbt_dy_va_te
  character (len=100), dimension(4), parameter :: tbt_data_file_name = [character(len=200) :: &
                                                  'tbt_data_ex_1_1001.bin', 'tbt_data_csa2_1_1001.bin', &
                                                  'tbt_data_te3_1_1001.bin', 'tbt_data_te7_1_1001.bin']
  logical :: create_data_file = .true.

  integer (kind=iac) :: jj, ii, integrator_number_steps
  integer (kind=iac) :: number_rotations, rotation_start_store, rotation_stop_store, step_stored_rotations, tot_stored_rotations
  character (len=60) :: element_id_dump
  character (len=4), dimension(2) :: exact_int_elements = ["HVCO", "RFCA"]
  character (len=4), dimension(2) :: nonlinear_thick_elements = ["SEXT", "OCTU"]

  character (len=precision(1._rac)+7) :: string_ex, string_int
  integer (kind=iac) :: accuracy_order, n_var, n_pre


  beam_pipe_radius = 1.0_rac !0.03_rac ![m]
  energy = 6500._rac ![GeV]
  gamma_rel = energy/proton_mass
  beta_rel = sqrt(1-(1/gamma_rel**2))
  emittance_geo_x = 3.608738638461539_rac*1E-10 ![m]
  emittance_geo_y = 3.608738638461539_rac*1E-10 ![m]
  emittance_norm_x = emittance_geo_x*gamma_rel*beta_rel ![m]
  emittance_norm_y = emittance_geo_x*gamma_rel*beta_rel ![m]
  beta_x_start = 4.450431806072019_rac*1E1 ![m] 
  beta_y_start = 1.064066102410145_rac*1E2 ![m] 
  sigma_x = sqrt(beta_x_start*emittance_geo_x)
  sigma_y = sqrt(beta_y_start*emittance_geo_y)
  sigma = sqrt(sigma_x**2 + sigma_y**2)
  print*, 'The particles energy is:', energy, '[GeV]'
  print*, 'The horizontal normalize emittance is', emittance_norm_x, '[m]'
  print*, 'The vertical normalize emittance is', emittance_norm_x, '[m]'

  call GetLattice (17, lattice_file_name)
  !number_elements=SIZE(lattice,1)=SHAPE(lattice)
  print*, 'The number of the elements in the lattice are: ', number_elements
  if (number_elements<=10)then
    print*, 'The lattice is:'
    print*, lattice%element_name
  endif


  ! Angles
  min_angle_deg = 0.5_rac
  max_angle_deg = 89.5_rac !80.5_rac
  n_angle_slices = 88_iac !8_iac
  ! Sigmas
  min_sigma = 0.1_rac*sigma!1._rac*sigma
  max_sigma = 6.1_rac*sigma!10._rac*sigma
  n_sigma_slices = 59_iac!98_iac
  ! Initial delta
  delta_in = 0.00027_rac!0._rac!0.00027_rac
  call GenerateParticlesInitialPolarDistribution (min_angle_deg, max_angle_deg, n_angle_slices, &
                                                  min_sigma, max_sigma, n_sigma_slices, delta_in)
  !call GetParticlesInitialDistribution (18, particles_initial_condition_file_name)
  !number_particles=SIZE(particle_in_dy_va,1)=SHAPE(particle_in_dy_va)
  allocate (in_dy_va(number_particles))
  allocate (st_dy_va(number_particles))
  allocate (fi_ex_dy_va(number_particles))
  allocate (fi_sac_dy_va(number_particles))
  allocate (fi_te_dy_va(number_particles))
  in_dy_va = particle_in_dy_va

  if (number_particles<=10)then
    print*, 'The particles initial conditions are:'
    print*, in_dy_va
  endif
  print*, 'The number of the tracked particles is: ', number_particles


  ! Number of revolutions
  number_rotations = 1001_iac
  ! Stored (dump) revolutions 
  ! Turn that the storation start (1 <= rotation_start_store <= number_rotations) stop (rotation_start_store <= rotation_start_store <= number_rotations) and step
  rotation_start_store = 1_iac
  rotation_stop_store = number_rotations
  step_stored_rotations = 1_iac
  tot_stored_rotations = int((rotation_stop_store-rotation_start_store)/step_stored_rotations,kind=iac) + 1_iac

  if ( rotation_start_store > rotation_stop_store .or. rotation_start_store < 0_iac .or. rotation_stop_store < 0_iac .or. &
        rotation_start_store > number_rotations .or. rotation_stop_store > number_rotations .or. &
        tot_stored_rotations > number_rotations ) then

    print*, "ERROR - main program !!!"
    print*, "Some error occured with the range of the stored rotations !!!"
    print*, "The following are set:"
    print*, "Number of stored rotations: ", tot_stored_rotations
    print*, "Start storing at rotations: ", rotation_start_store
    print*, "Stop storing at rotations: ", rotation_stop_store
          
  else

    print*, 'The number of rotations is:', number_rotations
    print*, 'The data storing will start at:', rotation_start_store, 'rotation'
    print*, 'The data will be stored every', step_stored_rotations, 'rotations'
    print*, 'The data storing will stop at:', rotation_stop_store, 'rotation'
    print*, 'The total number of the stored rotations is', tot_stored_rotations

  endif

  ! Element for data collection
  element_id_dump = lattice(ubound(lattice,1))%element_id
  print*, 'Element type for data collection:', lattice(ubound(lattice,1))%element_name, ' with ID:', element_id_dump

  
  integrator_number_steps= 621_iac !5001_iac !!!!! odd number >=3 !!!!!

  print*, "Exact Map"

  st_dy_va=in_dy_va

  call AllocateInitialize (tot_stored_rotations, in_dy_va, tbt_dy_va_ex)

  loop1_rev : do ii=1,number_rotations

    do jj=1,number_elements

       if ( ANY(lattice(jj)%element_name==nonlinear_thick_elements) ) then

         st_dy_va = Simple (integrator_number_steps, lattice(jj), st_dy_va, beam_pipe_radius)

       else

         st_dy_va = ExactMap (lattice(jj), st_dy_va, beam_pipe_radius)

       endif

       call StoreData (element_id_dump, lattice(jj), ii, rotation_start_store, rotation_stop_store, step_stored_rotations, &
                       st_dy_va, tbt_dy_va_ex, create_data_file, tbt_data_file_name(1))

       if (all(isnan(st_dy_va(:)%X))) then

         print*, 'All the particles get lost!'
         exit loop1_rev

       endif

    enddo

    if (modulo(ii,100_iac)==0_iac) then
      print*, ii,'/',number_rotations
    endif

  enddo loop1_rev

  fi_ex_dy_va=st_dy_va
  !print*, fi_ex_dy_va


  integrator_number_steps = 7_iac   ! !!!!! odd number >=3 !!!!!

  print*, "SABAC Integrator Map"
  print*, 'Number of steps',integrator_number_steps

  st_dy_va=in_dy_va

  call AllocateInitialize (tot_stored_rotations, in_dy_va, tbt_dy_va_sac)

  loop2_rev : do ii=1,number_rotations

    do jj=1,number_elements

      if ( ANY(lattice(jj)%element_name==exact_int_elements) ) then

        st_dy_va = ExactMap (lattice(jj), st_dy_va, beam_pipe_radius)

      else

        st_dy_va = Sabac (integrator_number_steps, lattice(jj), st_dy_va, beam_pipe_radius)

      endif

      call StoreData (element_id_dump, lattice(jj), ii, rotation_start_store, rotation_stop_store, step_stored_rotations, &
                       st_dy_va, tbt_dy_va_sac, create_data_file, tbt_data_file_name(2))

      if (all(isnan(st_dy_va(:)%X))) then

        print*, 'All the particles get lost!'
        exit loop2_rev

      endif

    enddo

    if (modulo(ii,100_iac)==0_iac) then
      print*, ii,'/',number_rotations
    endif

  enddo loop2_rev

  fi_sac_dy_va=st_dy_va
  !print*, fi_sac_dy_va


  integrator_number_steps = 7_iac   ! !!!!! odd number >=3 !!!!!

  print*, "TEAPOT Integrator Map"
  print*, 'Number of steps',integrator_number_steps

  st_dy_va=in_dy_va

  call AllocateInitialize (tot_stored_rotations, in_dy_va, tbt_dy_va_te)

  loop3_rev : do ii=1,number_rotations

    do jj=1,number_elements

      if ( ANY(lattice(jj)%element_name==exact_int_elements) ) then

        st_dy_va = ExactMap (lattice(jj), st_dy_va, beam_pipe_radius)

      else

        st_dy_va = Teapot (integrator_number_steps, lattice(jj), st_dy_va, beam_pipe_radius)

      endif

      call StoreData (element_id_dump, lattice(jj), ii, rotation_start_store, rotation_stop_store, step_stored_rotations, &
                       st_dy_va, tbt_dy_va_te, create_data_file, tbt_data_file_name(3))

      if (all(isnan(st_dy_va(:)%X))) then

        print*, 'All the particles get lost!'
        exit loop3_rev

      endif

    enddo

    if (modulo(ii,100_iac)==0_iac) then
      print*, ii,'/',number_rotations
    endif

  enddo loop3_rev

  fi_te_dy_va=st_dy_va
  !print*, fi_te_dy_va


  integrator_number_steps = 15_iac   ! !!!!! odd number >=3 !!!!!

  print*, "TEAPOT Integrator Map"
  print*, 'Number of steps',integrator_number_steps

  st_dy_va=in_dy_va

  call AllocateInitialize (tot_stored_rotations, in_dy_va, tbt_dy_va_te)

  loop4_rev : do ii=1,number_rotations

    do jj=1,number_elements

      if ( ANY(lattice(jj)%element_name==exact_int_elements) ) then

        st_dy_va = ExactMap (lattice(jj), st_dy_va, beam_pipe_radius)

      else

        st_dy_va = Teapot (integrator_number_steps, lattice(jj), st_dy_va, beam_pipe_radius)

      endif

      call StoreData (element_id_dump, lattice(jj), ii, rotation_start_store, rotation_stop_store, step_stored_rotations, &
                       st_dy_va, tbt_dy_va_te, create_data_file, tbt_data_file_name(4))

      if (all(isnan(st_dy_va(:)%X))) then

        print*, 'All the particles get lost!'
        exit loop4_rev

      endif

    enddo

    if (modulo(ii,100_iac)==0_iac) then
      print*, ii,'/',number_rotations
    endif

  enddo loop4_rev

  fi_te_dy_va=st_dy_va
  !print*, fi_te_dy_va


end program TRACKING
