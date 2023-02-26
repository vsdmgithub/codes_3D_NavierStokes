! --------------------------------------------------------------
! -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
! CODE BY:
! --------   |         |   ---------        /\        |\      |
! |          |         |  |                /  \       | \     |
! |          |         |  |               /    \      |  \    |
! --------   |         |  |   ------|    /------\     |   \   |
!         |  |         |  |         |   /        \    |    \  |
!         |  |         |  |         |  /          \   |     \ |
! ---------   ----------  ----------  /            \  |      \|
! --------------------------------------------------------------

! #########################
! MODULE: system_basicfunctions
! LAST MODIFIED: 20 FEBRAURY 2023
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! BASIC FUNCTIONS MODULE FOR 3D NSE ANALYSIS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_basicfunctions
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major system basic functions for the code to run.
! THese are standard functions, more advanced are done in advanced functions module
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_initialcondition
  USE system_basicoutput

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE normalized_initial_condition
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS TO :
  ! Get a normalized initial condition, first gets from a list of
  ! available ones, then checks for error and initializations.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  I  N  I  T  I  A  L        C  O  N  D  I  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL init_initcondn
    ! Calls the subroutine to get a initial condition
    ! REF-> <<< system_initialcondition >>>

    CALL compute_energy_spectral_data
    ! REF-> <<< system_initialcondition >>>

    CALL perform_debug
    ! Checks for any compressibility and NaN in data

    IF ( ( num_nan .EQ. 0 ) .AND. ( inc_err .EQ. 0) ) THEN

      pre_status = 1

      CALL compute_vorticity
      ! Calculates the vorticity (for the first time)

      CALL compute_spectral_data
      ! Gets the energy,enstrophy from spectral space

      CALL fft_c2r( V_x, V_y, V_z, N, Nh, U_x, U_y, U_z )
      ! FFT spectral to real velocity

    END IF

  END

  SUBROUTINE compute_spectral_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This calculates energy, spectral shell wise. It goes through each
  ! spectral mode and puts the energy in the corresponding shell.
  ! This gives the ENERGY SPECTRUM.  When the time is right, it saves
  ! them too.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::sh_no

    Eng_k     = zero
    ! Reset the array

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  P  E  C  T  R  U  M     C  A   L   C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Gets the energy in that particular mode (i_x,i_y,i_z)
    ! Keeps adding energy to that particular shell (|k| fixed), from all solid angles

    DO j_x =    1 , k_tru
    DO j_y = - k_tru, k_tru
    DO j_z = - k_tru, k_tru
    IF ( Lapla ( j_x, j_y, j_z ) .LT. k_tru_sqr ) THEN

      sh_no                       = Shell( j_x, j_y, j_z )
      eng_mod                 = CDABS( V_x( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( V_y( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( V_z( j_x, j_y, j_z ) ) ** two
      Eng_k( sh_no )    = Eng_k( sh_no )    + eng_mod

    END IF
    END DO
    END DO
    END DO

    j_x    =   0
    DO j_y = - k_tru, k_tru
    DO j_z = - k_tru, -1
    IF ( Lapla ( j_x, j_y, j_z ) .LT. k_tru_sqr ) THEN

      sh_no                       = Shell( j_x, j_y, j_z )
      eng_mod                 = CDABS( V_x( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( V_y( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( V_z( j_x, j_y, j_z ) ) ** two
      Eng_k( sh_no )    = Eng_k( sh_no )    + eng_mod

    END IF
    END DO
    END DO

    j_z    = 0
    DO j_y = 1, k_tru

      sh_no                       = Shell( j_x, j_y, j_z )
      eng_mod                 = CDABS( V_x( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( V_y( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( V_z( j_x, j_y, j_z ) ) ** two
      Eng_k( sh_no )    = Eng_k( sh_no )      + eng_mod

    END DO

    j_y = 0
    sh_no                       = Shell( j_x, j_y, j_z )
    eng_mod                 = CDABS( V_x( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( V_y( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( V_z( j_x, j_y, j_z ) ) ** two
    Eng_k( sh_no )    = Eng_k( sh_no )        + hf * eng_mod
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  S  H  E  L  L      A  V  E  R  A  G  I  N  G
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Eng_k_avg( 0 )         = Eng_k( 0 )
    Eng_k_avg( 1 )         = qtr * ( thr * Eng_k( 1 )    + Eng_k( 2 ) )
    DO k_ind                          = 2, k_max - 1
      Eng_k_avg( k_ind )    = qtr * ( Eng_k( k_ind - 1 )    + Eng_k( k_ind + 1 ) ) + &
                                        hf * ( Eng_k( k_ind ) )
    END DO

    ! eng = SUM( Eng_k )
    ! ens = SUM( Lapl * Eng_k )
    ! Computes the net energy, enstrophy,

  END

  SUBROUTINE compute_temporal_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This calculates energy, enstrophy and helicity
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  N  E  T     E , Z , H , D
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL compute_energy
    CALL compute_vorticity
    CALL compute_enstrophy
    CALL compute_helicity

    dis_eng = ( eng_pre - eng ) / dt  ! Estimates the dissipation rate of energy
    dis     = two * vis * ens         ! Estimates the viscous dissipation
    eng_pre = eng                     ! For next step

    CALL write_temporal_data
    ! REF-> <<< system_basicoutput >>>

  END

  SUBROUTINE compute_forcing_velocity_related
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! To assign values to the modes to force energy
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::ct
    DOUBLE PRECISION::frc_fac

    eng_frc = zero
    LOOP_FORCING_MODES_301: DO ct = 1, num_mod_frc

      j_x = F_kx_ind( ct )
      j_y = F_ky_ind( ct )
      j_z = F_kz_ind( ct )
      KX_EQ_ZERO_CHECK: IF ( j_x .EQ. 0 ) THEN
      eng_frc = eng_frc + hf *( CDABS( V_x( j_x, j_y, j_z ) ) ** two + &
                                                          CDABS( V_y( j_x, j_y, j_z ) ) ** two + &
                                                          CDABS( V_z( j_x, j_y, j_z ) ) ** two )
      ELSE
      eng_frc = eng_frc +  CDABS( V_x( j_x, j_y, j_z ) ) ** two + &
                                                     CDABS( V_y( j_x, j_y, j_z ) ) ** two + &
                                                     CDABS( V_z( j_x, j_y, j_z ) ) ** two
      END IF KX_EQ_ZERO_CHECK

    END DO LOOP_FORCING_MODES_301

		frc_fac  = dis  * ( ( one - ( eng - eng_0 ) ) ** two )
    ! Matches with the viscous dissipation, ! If still the energy is decreasing, then diss_rate would increase the forcing , and if energy is ! decreasing, it would decrease the forcing.

    IF ( eng_frc .GT. tol_double ) THEN
      frc_fac = frc_fac / ( two * eng_frc )
    ELSE
      frc_fac = zero ! To prevent NaN
    END IF

    LOOP_FORCING_MODES_302: DO ct = 1, num_mod_frc
      j_x = F_kx_ind( ct )
      j_y = F_ky_ind( ct )
      j_z = F_kz_ind( ct )
      I_fac( j_x, j_y, j_z ) = DEXP( - ( vis * Lapla( j_x, j_y, j_z ) - frc_fac ) * dt )
    END DO LOOP_FORCING_MODES_302

  END

  SUBROUTINE compute_random_forcing_adjusted
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! To assign values to the modes to force energy
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::ct
    DOUBLE PRECISION::phi
    DOUBLE COMPLEX::r_x,r_y,r_z
    DOUBLE PRECISION::r_nrm

		frc_fac = dis  * ( ( one - 10.0D0 * ( eng - eng_0 ) ) ** 8.0D0 )

    nrm_fac = DSQRT( frc_fac / ( num_mod_frc * dt ) )

    CALL init_random_seed
    ! REF-> <<< system_auxilaries >>>

    LOOP_FORCING_MODES_309: DO ct = 1, num_mod_frc

      j_x = F_kx_ind( ct )
      j_y = F_ky_ind( ct )
      j_z = F_kz_ind( ct )

      r_x = K_y( j_x, j_y, j_z ) * DCONJG( V_z( j_x, j_y, j_z ) ) &
          - K_z( j_x, j_y, j_z ) * DCONJG( V_y( j_x, j_y, j_z ) )
      r_y = K_z( j_x, j_y, j_z ) * DCONJG( V_x( j_x, j_y, j_z ) ) &
          - K_x( j_x, j_y, j_z ) * DCONJG( V_z( j_x, j_y, j_z ) )
      r_z = K_x( j_x, j_y, j_z ) * DCONJG( V_y( j_x, j_y, j_z ) ) &
          - K_y( j_x, j_y, j_z ) * DCONJG( V_x( j_x, j_y, j_z ) )
      r_nrm = DSQRT( CDABS( r_x ) ** two + CDABS( r_y ) ** two + CDABS( r_z ) ** two )
      r_x = r_x / r_nrm
      r_y = r_y / r_nrm
      r_z = r_z / r_nrm

      CALL RANDOM_NUMBER(phi)
      phi     = two_pi * phi

      F_x( j_x, j_y, j_z ) = nrm_fac * DCMPLX( DCOS( phi ), DSIN( phi ) ) * r_x
      F_y( j_x, j_y, j_z ) = nrm_fac * DCMPLX( DCOS( phi ), DSIN( phi ) ) * r_y
      F_z( j_x, j_y, j_z ) = nrm_fac * DCMPLX( DCOS( phi ), DSIN( phi ) ) * r_z

    END DO LOOP_FORCING_MODES_309

  END

  SUBROUTINE compute_random_forcing
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! To assign values to the modes to force energy
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::ct
    DOUBLE PRECISION::phi,the
    DOUBLE COMPLEX::r_x,r_y,r_z
    DOUBLE PRECISION,DIMENSION(3)::phs

		frc_fac = dis_fix

    nrm_fac = DSQRT( frc_fac / ( num_mod_frc * dt ) )

    CALL init_random_seed
    ! REF-> <<< system_auxilaries >>>

    LOOP_FORCING_MODES_389: DO ct = 1, num_mod_frc

      j_x = F_kx_ind( ct )
      j_y = F_ky_ind( ct )
      j_z = F_kz_ind( ct )

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(the)
      CALL RANDOM_NUMBER(phs)

      phi     = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      the   = DACOS( one - two * the )! Polar angle of \hat{u}_k vector
      phs      = two_pi * phs ! Phases of \hat{u}_k components

      r_x = nrm_fac * DSIN( the ) * DCOS( phi ) * DCMPLX( DCOS( phs( 1 ) ), DSIN( phs( 1 ) ) )
      r_y = nrm_fac * DSIN( the ) * DSIN( phi ) * DCMPLX( DCOS( phs( 2 ) ), DSIN( phs( 2 ) ) )
      r_z = nrm_fac * DCOS( the ) * DCMPLX( DCOS( phs( 3 ) ),DSIN( phs( 3 ) ) )

      F_x( j_x, j_y, j_z ) = Proj_xx( j_x, j_y, j_z ) * r_x + Proj_xy( j_x, j_y, j_z) * r_y + &
                                   Proj_zx( j_x, j_y, j_z ) * r_z
      F_y( j_x, j_y, j_z ) = Proj_xy( j_x, j_y, j_z ) * r_x + Proj_yy( j_x, j_y, j_z) * r_y + &
                                   Proj_yz( j_x, j_y, j_z ) * r_z
      F_z( j_x, j_y, j_z ) = Proj_zx( j_x, j_y, j_z ) * r_x + Proj_yz( j_x, j_y, j_z) * r_y + &
                                   Proj_zz( j_x, j_y, j_z ) * r_z

    END DO LOOP_FORCING_MODES_389

  END

  SUBROUTINE perturb_forcing
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to make a perturbation in the forcing vector
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    INTEGER(KIND=4)::ct
    DOUBLE PRECISION::fac

    ct = num_mod_frc
    j_x = F_kx_ind( ct )
    j_y = F_ky_ind( ct )
    j_z = F_kz_ind( ct )

    fac = 0.98D0
    F_x( j_x, j_y, j_z ) = fac * F_x( j_x, j_y, j_z )
    F_y( j_x, j_y, j_z ) = fac * F_y( j_x, j_y, j_z )
    F_z( j_x, j_y, j_z ) = fac * F_z( j_x, j_y, j_z )

  END

  SUBROUTINE compute_energy
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the Kinetic energy in real space
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    eng = hf * SUM( U_x ** two + U_y ** two + U_z ** two ) / N3

  END

  SUBROUTINE compute_vorticity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! -------------
  ! CALL this to get vorticity field
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    W_kx = i * ( K_y * V_z - K_z * V_y )
    W_ky = i * ( K_z * V_x - K_x * V_z )
    W_kz = i * ( K_x * V_y - K_y * V_x )
    ! Spectral Vorticity

    CALL fft_c2r( W_kx, W_ky, W_kz, N, Nh, W_x, W_y, W_z )
    ! Real Vorticity

  END

  SUBROUTINE compute_helicity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get helicity
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    hel = hf * SUM( W_x * U_x + W_y * U_y + W_z * U_z ) / N3

  END

  SUBROUTINE compute_enstrophy
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the enstrophy in real space
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ens = hf * SUM( W_x ** two + W_y ** two + W_z ** two ) / N3

  END

  SUBROUTINE perform_debug
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check incompressibility criterion and Nan in data
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    CALL check_inc
    CALL check_nan
    CALL check_cfl

  END

  SUBROUTINE check_inc
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check incompressibility condition. Sums over all residues
  ! of incompressibility in L2 norm and prints it. Of order 10^(-12).
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    inc             = zero

    inc_err       = 0

    DO j_x          = 0, Nh
    DO j_y          = -Nh, Nh - 1
    DO j_z          = -Nh, Nh - 1
                inc = inc + CDABS( j_x * V_x( j_x, j_y, j_z ) + &
                                   j_y * V_y( j_x, j_y, j_z ) + &
                                   j_z * V_z( j_x, j_y, j_z ) ) ** two
    END DO
    END DO
    END DO

    j_x             = 0
    DO j_y          = - Nh, Nh - 1
    DO j_z          = - Nh, Nh - 1
      inc           = inc + hf* CDABS( j_x * V_x( j_x, j_y, j_z ) + &
                                       j_y * V_y( j_x, j_y, j_z ) + &
                                       j_z * V_z( j_x, j_y, j_z ) ) ** two
    END DO
    END DO

    j_x             = Nh
    DO j_y          = - Nh, Nh - 1
    DO j_z          = - Nh, Nh - 1
      inc           = inc + hf* CDABS( j_x * V_x( j_x, j_y, j_z ) + &
                                       j_y * V_y( j_x, j_y, j_z ) + &
                                       j_z * V_z( j_x, j_y, j_z ) ) ** two
    END DO
    END DO

    COMPRESSIBILITY_CHECK_101: IF (inc .GT. tol_float) THEN

      inc_err = 1

      deb_err = 1
      ! This will jump out of evolution loop, if caught during that.

      CALL print_error_inc
      ! This will prompt error when checked
      ! REF-> <<< system_output >>>

    END IF COMPRESSIBILITY_CHECK_101

  END

  SUBROUTINE check_nan
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check if there is any NaN in the data.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    num_nan  = 0

    DO j_x         = 0, Nh
    DO j_y         = -Nh, Nh - 1
    DO j_z         = -Nh, Nh - 1
      IF ( V_x( j_x, j_y, j_z ) .NE. V_x( j_x, j_y, j_z ) ) THEN
        num_nan = num_nan + 1
      END IF
    END DO
    END DO
    END DO

    j_x            = 0
    DO j_y         = - Nh, Nh - 1
    DO j_z         = - Nh, Nh - 1
      IF ( V_x( j_x, j_y, j_z ) .NE. V_x( j_x, j_y, j_z ) ) THEN
        num_nan = num_nan + 1
      END IF
    END DO
    END DO

    j_x            = Nh
    DO j_y         = - Nh, Nh - 1
    DO j_z         = - Nh, Nh - 1
      IF ( V_x( j_x, j_y, j_z ) .NE. V_x( j_x, j_y, j_z ) ) THEN
        num_nan = num_nan + 1
      END IF
    END DO
    END DO

    NAN_CHECK_101: IF (num_nan .NE. 0) THEN

      nan_err   = 1

      deb_err = 1
      ! This will jump out of evolution loop, if caught during that.

      CALL print_error_nan
      ! This will prompt error when checked
      ! REF-> <<< system_output >>>

    END IF NAN_CHECK_101

  END

  SUBROUTINE check_cfl
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check if there is any NaN in the data.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    DOUBLE PRECISION,DIMENSION(3)::cfl_min
    DOUBLE PRECISION::cfl_new

    cfl_min( 1 )= l_grd / ( MAXVAL( DABS( U_x ) ) * dt )
    cfl_min( 2 )= l_grd / ( MAXVAL( DABS( U_y ) ) * dt )
    cfl_min( 3 )= l_grd / ( MAXVAL( DABS( U_z ) ) * dt )

    cfl_new = MAXVAL( cfl_min )

    IF ( cfl_new .LT. cfl ) THEN
      dt_max = dt * cfl_new / cfl
      CALL find_timestep( t_min, dt_max, dt )
      ! REF-> <<< system_auxilaries >>>
    END IF

  END

END MODULE system_basicfunctions
