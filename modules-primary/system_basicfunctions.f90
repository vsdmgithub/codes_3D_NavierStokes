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
! LAST MODIFIED: 21 JUNE 2021
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

    IF ( ( nan_count .EQ. 0 ) .AND. ( incomp_error .EQ. 0) ) THEN

      precheck_status = 1

      CALL compute_vorticity
      ! Calculates the vorticity (for the first time)

      CALL compute_spectral_data
      ! Gets the energy,enstrophy from spectral space

      CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
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

    spectral_energy     = zero
    spectral_enstrophy  = zero
    spectral_helicity   = zero
    ! Reset the array

    energy              = zero
    enstrophy           = zero
    helicity            = zero
    ! Reset the variables

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  P  E  C  T  R  U  M     C  A   L   C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Gets the energy, enstrophy, helicity in that particular mode (i_x,i_y,i_z)
    ! Keeps adding energy to that particular shell (|k| fixed), from all solid angles

    CALL compute_vorticity
    ! Computes vorticity in spectral space

    DO j_x =    1 , k_G
    DO j_y = - k_G, k_G
    DO j_z = - k_G, k_G
    IF ( k_2 ( j_x, j_y, j_z ) .LT. k_G_2 ) THEN

      sh_no                       = shell_no( j_x, j_y, j_z )
      energy_mode                 = CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( v_z( j_x, j_y, j_z ) ) ** two
      enstrophy_mode              = CDABS( w_vx( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( w_vy( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( w_vz( j_x, j_y, j_z ) ) ** two
      helicity_mode_complex       = v_x( j_x, j_y, j_z ) * DCONJG( w_vx( j_x, j_y, j_z ) ) + &
                                    v_y( j_x, j_y, j_z ) * DCONJG( w_vy( j_x, j_y, j_z ) ) + &
                                    v_z( j_x, j_y, j_z ) * DCONJG( w_vz( j_x, j_y, j_z ) )
      helicity_mode               = DREAL( helicity_mode_complex )
      spectral_energy( sh_no )    = spectral_energy( sh_no )    + energy_mode
      spectral_enstrophy( sh_no ) = spectral_enstrophy( sh_no ) + enstrophy_mode
      spectral_helicity( sh_no )  = spectral_helicity( sh_no )  + helicity_mode

    END IF
    END DO
    END DO
    END DO

    j_x    =   0
    DO j_y = - k_G, k_G
    DO j_z = - k_G, -1
    IF ( k_2 ( j_x, j_y, j_z ) .LT. k_G_2 ) THEN

      sh_no                       = shell_no( j_x, j_y, j_z )
      energy_mode                 = CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( v_z( j_x, j_y, j_z ) ) ** two
      enstrophy_mode              = CDABS( w_vx( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( w_vy( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( w_vz( j_x, j_y, j_z ) ) ** two
      helicity_mode_complex       = v_x( j_x, j_y, j_z ) * DCONJG( w_vx( j_x, j_y, j_z ) ) + &
                                    v_y( j_x, j_y, j_z ) * DCONJG( w_vy( j_x, j_y, j_z ) ) + &
                                    v_z( j_x, j_y, j_z ) * DCONJG( w_vz( j_x, j_y, j_z ) )
      helicity_mode               = DREAL( helicity_mode_complex )
      spectral_energy( sh_no )    = spectral_energy( sh_no )    + energy_mode
      spectral_enstrophy( sh_no ) = spectral_enstrophy( sh_no ) + enstrophy_mode
      spectral_helicity( sh_no )  = spectral_helicity( sh_no )  + helicity_mode

    END IF
    END DO
    END DO

    j_z    = 0
    DO j_y = 1, k_G

      sh_no                       = shell_no( j_x, j_y, j_z )
      energy_mode                 = CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( v_z( j_x, j_y, j_z ) ) ** two
      enstrophy_mode              = CDABS( w_vx( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( w_vy( j_x, j_y, j_z ) ) ** two + &
                                    CDABS( w_vz( j_x, j_y, j_z ) ) ** two
      helicity_mode_complex       = v_x( j_x, j_y, j_z ) * DCONJG( w_vx( j_x, j_y, j_z ) ) + &
                                    v_y( j_x, j_y, j_z ) * DCONJG( w_vy( j_x, j_y, j_z ) ) + &
                                    v_z( j_x, j_y, j_z ) * DCONJG( w_vz( j_x, j_y, j_z ) )
      helicity_mode               = DREAL( helicity_mode_complex )
      spectral_energy( sh_no )    = spectral_energy( sh_no )      + energy_mode
      spectral_enstrophy( sh_no ) = spectral_enstrophy( sh_no )   + enstrophy_mode
      spectral_helicity( sh_no )  = spectral_helicity( sh_no )    + helicity_mode

    END DO

    j_y = 0
    sh_no                       = shell_no( j_x, j_y, j_z )
    energy_mode                 = CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( v_z( j_x, j_y, j_z ) ) ** two
    enstrophy_mode              = CDABS( w_vx( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( w_vy( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( w_vz( j_x, j_y, j_z ) ) ** two
    helicity_mode_complex       = v_x( j_x, j_y, j_z ) * DCONJG( w_vx( j_x, j_y, j_z ) ) + &
                                  v_y( j_x, j_y, j_z ) * DCONJG( w_vy( j_x, j_y, j_z ) ) + &
                                  v_z( j_x, j_y, j_z ) * DCONJG( w_vz( j_x, j_y, j_z ) )
    helicity_mode               = DREAL( helicity_mode_complex )
    spectral_energy( sh_no )    = spectral_energy( sh_no )        + hf * energy_mode
    spectral_enstrophy( sh_no ) = spectral_enstrophy( sh_no )     + hf * enstrophy_mode
    spectral_helicity( sh_no )  = spectral_helicity( sh_no )      + hf * helicity_mode

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  S  H  E  L  L      A  V  E  R  A  G  I  N  G
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    spectral_energy_avg( 0 )         = spectral_energy( 0 )
    spectral_enstrophy_avg( 0 )      = spectral_enstrophy( 0 )
    spectral_helicity_avg( 0 )       = spectral_helicity( 0 )
    spectral_energy_avg( 1 )         = qtr * ( thr * spectral_energy( 1 )    + spectral_energy( 2 ) )
    spectral_enstrophy_avg( 1 )      = qtr * ( thr * spectral_enstrophy( 1 ) + spectral_enstrophy( 2 ) )
    spectral_helicity_avg( 1 )       = qtr * ( thr * spectral_helicity( 1 )  + spectral_helicity( 2 ) )

    DO k_no                          = 2, k_max - 1

      spectral_energy_avg( k_no )    = qtr * ( spectral_energy( k_no - 1 )    + spectral_energy( k_no + 1 ) ) + &
                                        hf * ( spectral_energy( k_no ) )
      spectral_enstrophy_avg( k_no ) = qtr * ( spectral_enstrophy( k_no - 1 ) + spectral_enstrophy( k_no + 1 ) ) + &
                                        hf * ( spectral_enstrophy( k_no ) )
      spectral_helicity_avg( k_no )  = qtr * ( spectral_helicity( k_no - 1 )  + spectral_helicity( k_no + 1 ) ) + &
                                        hf * ( spectral_helicity( k_no ) )

    END DO

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  N  E  T     E , Z , H , D
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    energy    = SUM( spectral_energy(    : ) )
    enstrophy = SUM( spectral_enstrophy( : ) )
    helicity  = SUM( spectral_helicity(  : ) )
    ! Computes the net energy, enstrophy, helicity

    diss_rate         = ( energy_initial - energy ) / dt
    energy_old        = energy
    ! Estimates the dissipation rate of energy

    diss_rate_viscous = two * viscosity * enstrophy
    ! Estimates the viscous dissipation from enstrophy

  END

  SUBROUTINE compute_forcing_in_modes
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
    DOUBLE PRECISION::rd1,rd2,std,avg,noise
    DOUBLE PRECISION::pre_factor_forcing

    energy_forcing_modes = zero
    LOOP_FORCING_MODES_301: DO ct = 1, tot_forced_modes

      j_x = fkx( ct )
      j_y = fky( ct )
      j_z = fkz( ct )
      KX_EQ_ZERO_CHECK: IF ( j_x .EQ. 0 ) THEN
      energy_forcing_modes = energy_forcing_modes + hf *( CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                                                          CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                                                          CDABS( v_z( j_x, j_y, j_z ) ) ** two )
      energy_forcing_modes = energy_forcing_modes + hf *( CDABS( v2_x( j_x, j_y, j_z ) ) ** two + &
                                                          CDABS( v2_y( j_x, j_y, j_z ) ) ** two + &
                                                          CDABS( v2_z( j_x, j_y, j_z ) ) ** two )
      ELSE
      energy_forcing_modes = energy_forcing_modes +  CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                                                     CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                                                     CDABS( v_z( j_x, j_y, j_z ) ) ** two
      energy_forcing_modes = energy_forcing_modes +  CDABS( v2_x( j_x, j_y, j_z ) ) ** two + &
                                                     CDABS( v2_y( j_x, j_y, j_z ) ) ** two + &
                                                     CDABS( v2_z( j_x, j_y, j_z ) ) ** two
      END IF KX_EQ_ZERO_CHECK

    END DO LOOP_FORCING_MODES_301

    energy_forcing_modes = energy_forcing_modes / two

		pre_factor_forcing   = diss_rate_viscous + diss_rate
    ! Matches with the visous dissipation
    ! If still the energy is decreasing, then diss_rate would increase the forcing , and if energy is
    ! decreasing, it would decrease the forcing.

    pre_factor_forcing = pre_factor_forcing / ( two * energy_forcing_modes )

    ! pre_factor_forcing = 0.710D0
    LOOP_FORCING_MODES_302: DO ct = 1, tot_forced_modes
      j_x = fkx( ct )
      j_y = fky( ct )
      j_z = fkz( ct )
      integrating_factor( j_x, j_y, j_z ) = DEXP( - ( viscosity * k_2( j_x, j_y, j_z ) - pre_factor_forcing ) * dt )
    END DO LOOP_FORCING_MODES_302

  END

  SUBROUTINE compute_forcing_in_modes_with_perturbation
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
    DOUBLE PRECISION::rd1,rd2,std,avg,noise
    DOUBLE PRECISION::pre_factor_forcing

    energy_forcing_modes = zero
    LOOP_FORCING_MODES_301: DO ct = 1, tot_forced_modes

      j_x = fkx( ct )
      j_y = fky( ct )
      j_z = fkz( ct )
      KX_EQ_ZERO_CHECK: IF ( j_x .EQ. 0 ) THEN
      energy_forcing_modes = energy_forcing_modes + hf *( CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                                                          CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                                                          CDABS( v_z( j_x, j_y, j_z ) ) ** two )
      ELSE
      energy_forcing_modes = energy_forcing_modes +  CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                                                     CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                                                     CDABS( v_z( j_x, j_y, j_z ) ) ** two
      END IF KX_EQ_ZERO_CHECK

    END DO LOOP_FORCING_MODES_301

		pre_factor_forcing = diss_rate_viscous + diss_rate
    ! Matches with the visous dissipation
    ! If still the energy is decreasing, then diss_rate would increase the forcing , and if energy is
    ! decreasing, it would decrease the forcing.
    pre_factor_forcing = pre_factor_forcing / ( two * energy_forcing_modes )

    ! pre_factor_forcing = 0.710D0
    LOOP_FORCING_MODES_302: DO ct = 1, tot_forced_modes
      j_x = fkx( ct )
      j_y = fky( ct )
      j_z = fkz( ct )
      integrating_factor( j_x, j_y, j_z ) = DEXP( - ( viscosity * k_2( j_x, j_y, j_z ) - 0.90D0 * pre_factor_forcing ) * dt )
    END DO LOOP_FORCING_MODES_302

  END

  SUBROUTINE compute_energy
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the Kinetic energy in real space
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    energy = hf * SUM( u_x ** two + u_y ** two + u_z ** two ) / N3

  END

  SUBROUTINE compute_vorticity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get vorticity field
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    w_vx = i * ( k_y * v_z - k_z * v_y )
    w_vy = i * ( k_z * v_x - k_x * v_z )
    w_vz = i * ( k_x * v_y - k_y * v_x )
    ! Spectral Vorticity

    CALL fft_c2r( w_vx, w_vy, w_vz, N, Nh, w_ux, w_uy, w_uz )
    ! Real Vorticity

  END

  SUBROUTINE compute_helicity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get helicity
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    helicity = hf * SUM( w_ux * u_x + w_uy * u_y + w_uz * u_z ) / N3

  END

  SUBROUTINE compute_enstrophy
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to get the enstrophy in real space
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    enstrophy = hf * SUM( w_ux ** two + w_uy ** two + w_uz ** two ) / N3

  END

  SUBROUTINE perform_debug
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check incompressibility criterion and Nan in data
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    CALL check_nan

    CALL compute_compressibility

  END

  SUBROUTINE compute_compressibility
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check incompressibility condition. Sums over all residues
  ! of incompressibility and prints it. Of order 10^(-12).
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    k_dot_v_norm   = zero

    incomp_error   = zero

    DO j_x         = 0, Nh
    DO j_y         = -Nh, Nh - 1
    DO j_z         = -Nh, Nh - 1
      k_dot_v_norm = k_dot_v_norm + CDABS( j_x * v_x( j_x, j_y, j_z ) + &
                                           j_y * v_y( j_x, j_y, j_z ) + &
                                           j_z * v_z( j_x, j_y, j_z ) ) ** two
    END DO
    END DO
    END DO

    j_x            = 0
    DO j_y         = - Nh, Nh - 1
    DO j_z         = - Nh, Nh - 1
      k_dot_v_norm = k_dot_v_norm + hf* CDABS( j_x * v_x( j_x, j_y, j_z ) + &
                                               j_y * v_y( j_x, j_y, j_z ) + &
                                               j_z * v_z( j_x, j_y, j_z ) ) ** two
    END DO
    END DO

    j_x            = Nh
    DO j_y         = - Nh, Nh - 1
    DO j_z         = - Nh, Nh - 1
      k_dot_v_norm = k_dot_v_norm + hf* CDABS( j_x * v_x( j_x, j_y, j_z ) + &
                                               j_y * v_y( j_x, j_y, j_z ) + &
                                               j_z * v_z( j_x, j_y, j_z ) ) ** two
    END DO
    END DO

    k_dot_v_norm = DSQRT( k_dot_v_norm )

    COMPRESSIBILITY_CHECK_101: IF (k_dot_v_norm .GT. tol_float ) THEN

      incomp_error = 1

      debug_error = 1
      ! This will jump out of evolution loop, if caught during that.

      CALL print_error_incomp
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

    nan_count  = 0

    DO j_x         = 0, Nh
    DO j_y         = -Nh, Nh - 1
    DO j_z         = -Nh, Nh - 1
      IF ( v_x( j_x, j_y, j_z ) .NE. v_x( j_x, j_y, j_z ) ) THEN
        nan_count = nan_count + 1
      END IF
    END DO
    END DO
    END DO

    j_x            = 0
    DO j_y         = - Nh, Nh - 1
    DO j_z         = - Nh, Nh - 1
      IF ( v_x( j_x, j_y, j_z ) .NE. v_x( j_x, j_y, j_z ) ) THEN
        nan_count = nan_count + 1
      END IF
    END DO
    END DO

    j_x            = Nh
    DO j_y         = - Nh, Nh - 1
    DO j_z         = - Nh, Nh - 1
      IF ( v_x( j_x, j_y, j_z ) .NE. v_x( j_x, j_y, j_z ) ) THEN
        nan_count = nan_count + 1
      END IF
    END DO
    END DO

    NAN_CHECK_101: IF (nan_count .NE. 0) THEN

      nan_error   = 1

      debug_error = 1
      ! This will jump out of evolution loop, if caught during that.

      CALL print_error_nan
      ! This will prompt error when checked
      ! REF-> <<< system_output >>>

    END IF NAN_CHECK_101

  END

END MODULE system_basicfunctions
