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
  ! Get a normalized initial condition
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

    IF ( ( NaN_count .EQ. 0 ) .AND. ( k_dot_v_error .EQ. 0) ) THEN

      check_status = 1

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

    DO i_x =    1 , k_G
    DO i_y = - k_G, k_G
    DO i_z = - k_G, k_G
    IF ( k_2 ( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

      sh_no                       = shell_no( i_x, i_y, i_z )
      energy_mode                 = CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( v_z( i_x, i_y, i_z ) ) ** two
      enstrophy_mode              = CDABS( w_vx( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( w_vy( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( w_vz( i_x, i_y, i_z ) ) ** two
      helicity_mode_complex       = v_x( i_x, i_y, i_z ) * DCONJG( w_vx( i_x, i_y, i_z ) ) + &
                                    v_y( i_x, i_y, i_z ) * DCONJG( w_vy( i_x, i_y, i_z ) ) + &
                                    v_z( i_x, i_y, i_z ) * DCONJG( w_vz( i_x, i_y, i_z ) )
      helicity_mode               = DREAL( helicity_mode_complex )
      spectral_energy( sh_no )    = spectral_energy( sh_no )    + energy_mode
      spectral_enstrophy( sh_no ) = spectral_enstrophy( sh_no ) + enstrophy_mode
      spectral_helicity( sh_no )  = spectral_helicity( sh_no )  + helicity_mode

    END IF
    END DO
    END DO
    END DO

    i_x    =   0
    DO i_y = - k_G, k_G
    DO i_z = - k_G, -1
    IF ( k_2 ( i_x, i_y, i_z ) .LT. k_G_2 ) THEN

      sh_no                       = shell_no( i_x, i_y, i_z )
      energy_mode                 = CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( v_z( i_x, i_y, i_z ) ) ** two
      enstrophy_mode              = CDABS( w_vx( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( w_vy( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( w_vz( i_x, i_y, i_z ) ) ** two
      helicity_mode_complex       = v_x( i_x, i_y, i_z ) * DCONJG( w_vx( i_x, i_y, i_z ) ) + &
                                    v_y( i_x, i_y, i_z ) * DCONJG( w_vy( i_x, i_y, i_z ) ) + &
                                    v_z( i_x, i_y, i_z ) * DCONJG( w_vz( i_x, i_y, i_z ) )
      helicity_mode               = DREAL( helicity_mode_complex )
      spectral_energy( sh_no )    = spectral_energy( sh_no )    + energy_mode
      spectral_enstrophy( sh_no ) = spectral_enstrophy( sh_no ) + enstrophy_mode
      spectral_helicity( sh_no )  = spectral_helicity( sh_no )  + helicity_mode

    END IF
    END DO
    END DO

    i_z    = 0
    DO i_y = 1, k_G

      sh_no                       = shell_no( i_x, i_y, i_z )
      energy_mode                 = CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( v_z( i_x, i_y, i_z ) ) ** two
      enstrophy_mode              = CDABS( w_vx( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( w_vy( i_x, i_y, i_z ) ) ** two + &
                                    CDABS( w_vz( i_x, i_y, i_z ) ) ** two
      helicity_mode_complex       = v_x( i_x, i_y, i_z ) * DCONJG( w_vx( i_x, i_y, i_z ) ) + &
                                    v_y( i_x, i_y, i_z ) * DCONJG( w_vy( i_x, i_y, i_z ) ) + &
                                    v_z( i_x, i_y, i_z ) * DCONJG( w_vz( i_x, i_y, i_z ) )
      helicity_mode               = DREAL( helicity_mode_complex )
      spectral_energy( sh_no )    = spectral_energy( sh_no )      + energy_mode
      spectral_enstrophy( sh_no ) = spectral_enstrophy( sh_no )   + enstrophy_mode
      spectral_helicity( sh_no )  = spectral_helicity( sh_no )    + helicity_mode

    END DO

    i_y = 0
    sh_no                       = shell_no( i_x, i_y, i_z )
    energy_mode                 = CDABS( v_x( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_y( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( v_z( i_x, i_y, i_z ) ) ** two
    enstrophy_mode              = CDABS( w_vx( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( w_vy( i_x, i_y, i_z ) ) ** two + &
                                  CDABS( w_vz( i_x, i_y, i_z ) ) ** two
    helicity_mode_complex       = v_x( i_x, i_y, i_z ) * DCONJG( w_vx( i_x, i_y, i_z ) ) + &
                                  v_y( i_x, i_y, i_z ) * DCONJG( w_vy( i_x, i_y, i_z ) ) + &
                                  v_z( i_x, i_y, i_z ) * DCONJG( w_vz( i_x, i_y, i_z ) )
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

    diss_rate         = ( energy_old - energy ) / dt
    energy_old        = energy
    ! Estimates the dissipation rate of energy

    diss_rate_viscous = two * viscosity * enstrophy
    ! Estimates the viscous dissipation from enstrophy

  END

  SUBROUTINE forcing_spectrum
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Create a template for forcing, near integral wavenumbers.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION             ::phi,theta
    DOUBLE PRECISION             ::toss, threshold, energy_diff
    DOUBLE PRECISION,DIMENSION(3)::ph
    DOUBLE PRECISION,DIMENSION(1)::g_toss
    DOUBLE PRECISION,DIMENSION(3)::F_toss
    DOUBLE PRECISION             ::F_k_mod,k_ratio
    DOUBLE COMPLEX,DIMENSION(3)  ::F_k
    INTEGER(KIND=4)::fmodes

    CALL RANDOM_NUMBER( toss )

    energy_diff    = ( energy - energy_initial ) / ( 0.5 * energy_initial )
    threshold      = DEXP( energy_diff )
    forcing_status = 0
    fmodes         = 0

    IF ( ( toss .LT. threshold ) .AND. ( energy_diff .LT. zero ) ) THEN

      forcing_status = 1

      CALL init_random_seed
      ! Randomizes seed for random numbers (in 'auxilary_functions' module )
      ! REF-> <<< system_auxilaries >>>

      forced_energy  = zero
      ! Energy rate injected by forcing

      DO i_x =    1, k_F
      DO i_y = -k_F, k_F
      DO i_z = -k_F, k_F
      IF ( k_2( i_x, i_y, i_z ) .LT. k_F_2 ) THEN

        CALL RANDOM_NUMBER(ph)

        ph      = two_pi * ph ! Phases of \hat{u}_k components

        F_toss  = normal_dist( 3, zero, one )

        F_k(1)  = F_toss(1) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
        F_k(2)  = F_toss(2) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
        F_k(3)  = F_toss(3) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
        ! 3 COMPLEX values for spectral velocity

        F_k_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * F_k(1) + proj_xy( i_x, i_y, i_z) * F_k(2) + &
                                     proj_zx( i_x, i_y, i_z ) * F_k(3)
        F_k_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * F_k(1) + proj_yy( i_x, i_y, i_z) * F_k(2) + &
                                     proj_yz( i_x, i_y, i_z ) * F_k(3)
        F_k_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * F_k(1) + proj_yz( i_x, i_y, i_z) * F_k(2) + &
                                     proj_zz( i_x, i_y, i_z ) * F_k(3)

        forced_energy = forced_energy + DREAL( F_k_x( i_x, i_y, i_z ) * DCONJG( v_x( i_x, i_y, i_z) ) + &
                                               F_k_y( i_x, i_y, i_z ) * DCONJG( v_y( i_x, i_y, i_z) ) + &
                                               F_k_z( i_x, i_y, i_z ) * DCONJG( v_z( i_x, i_y, i_z) ) )

        ! ----------------------------------------------------------------------------------------
      END IF
      END DO
      END DO
      END DO

      i_x    = 0
      DO i_y = -k_F, k_F
      DO i_z = -k_F, -1
      IF ( k_2( i_x, i_y, i_z ) .LT. k_F_2 ) THEN

        CALL RANDOM_NUMBER(ph)

        ph      = two_pi * ph ! Phases of \hat{u}_k components

        F_toss  = normal_dist( 3, zero, one )

        F_k(1)  = F_toss(1) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
        F_k(2)  = F_toss(2) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
        F_k(3)  = F_toss(3) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
        ! 3 COMPLEX values for spectral velocity

        F_k_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * F_k(1) + proj_xy( i_x, i_y, i_z) * F_k(2) + &
                                    proj_zx( i_x, i_y, i_z ) * F_k(3)
        F_k_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * F_k(1) + proj_yy( i_x, i_y, i_z) * F_k(2) + &
                                     proj_yz( i_x, i_y, i_z ) * F_k(3)
        F_k_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * F_k(1) + proj_yz( i_x, i_y, i_z) * F_k(2) + &
                                     proj_zz( i_x, i_y, i_z ) * F_k(3)

        F_k_x( i_x, - i_y, - i_z ) = DCONJG( F_k_x( i_x, i_y, i_z ) )
        F_k_y( i_x, - i_y, - i_z ) = DCONJG( F_k_y( i_x, i_y, i_z ) )
        F_k_z( i_x, - i_y, - i_z ) = DCONJG( F_k_z( i_x, i_y, i_z ) )

        forced_energy = forced_energy + DREAL( F_k_x( i_x, i_y, i_z ) * DCONJG( v_x( i_x, i_y, i_z) ) + &
                                               F_k_y( i_x, i_y, i_z ) * DCONJG( v_y( i_x, i_y, i_z) ) + &
                                               F_k_z( i_x, i_y, i_z ) * DCONJG( v_z( i_x, i_y, i_z) ) )

        ! ----------------------------------------------------------------------------------------

      END IF
      END DO
      END DO

      i_z    = 0
      DO i_y = 1, k_F

        CALL RANDOM_NUMBER(ph)

        ph      = two_pi * ph ! Phases of \hat{u}_k components

        F_toss  = normal_dist( 3, zero, one )

        F_k(1)  = F_toss(1) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
        F_k(2)  = F_toss(2) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
        F_k(3)  = F_toss(3) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
        ! 3 COMPLEX values for spectral velocity

        F_k_x( i_x, i_y, i_z )     = proj_xx( i_x, i_y, i_z ) * F_k(1) + proj_xy( i_x, i_y, i_z) * F_k(2) + &
                                     proj_zx( i_x, i_y, i_z ) * F_k(3)
        F_k_y( i_x, i_y, i_z )     = proj_xy( i_x, i_y, i_z ) * F_k(1) + proj_yy( i_x, i_y, i_z) * F_k(2) + &
                                     proj_yz( i_x, i_y, i_z ) * F_k(3)
        F_k_z( i_x, i_y, i_z )     = proj_zx( i_x, i_y, i_z ) * F_k(1) + proj_yz( i_x, i_y, i_z) * F_k(2) + &
                                     proj_zz( i_x, i_y, i_z ) * F_k(3)

        F_k_x( i_x, - i_y, - i_z ) = DCONJG( F_k_x( i_x, i_y, i_z ) )
        F_k_y( i_x, - i_y, - i_z ) = DCONJG( F_k_y( i_x, i_y, i_z ) )
        F_k_z( i_x, - i_y, - i_z ) = DCONJG( F_k_z( i_x, i_y, i_z ) )

        forced_energy = forced_energy + DREAL( F_k_x( i_x, i_y, i_z ) * DCONJG( v_x( i_x, i_y, i_z) ) + &
                                               F_k_y( i_x, i_y, i_z ) * DCONJG( v_y( i_x, i_y, i_z) ) + &
                                               F_k_z( i_x, i_y, i_z ) * DCONJG( v_z( i_x, i_y, i_z) ) )

        ! ----------------------------------------------------------------------------------------
      END DO

      ! Making sure, that the zero mode forcing is zero.
      F_k_x( 0, 0, 0 )    =     zero
      F_k_y( 0, 0, 0 )    =     zero
      F_k_z( 0, 0, 0 )    =     zero

      forced_energy = forced_energy * two
      ! For the conjugate modes

      g_toss       = normal_dist( 1, 1.0D0, 0.5D0 ) * ( one - threshold )**two
      force_factor = g_toss(1) * diss_rate_viscous / forced_energy
      ! Recalibration of forcing according to dissipation with a noise

! print*,'forced',energy,g_toss(1)
      F_k_x         = force_factor * dt * F_k_x * diss_Ifactor
      F_k_y         = force_factor * dt * F_k_y * diss_Ifactor
      F_k_z         = force_factor * dt * F_k_z * diss_Ifactor
      ! Normalizing the forcing so that it equals the dissipation rate

    END IF

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

    DO i_x         = 0, Nh
    DO i_y         = -Nh, Nh - 1
    DO i_z         = -Nh, Nh - 1
      k_dot_v_norm = k_dot_v_norm + CDABS( i_x * v_x( i_x, i_y, i_z ) + &
                                           i_y * v_y( i_x, i_y, i_z ) + &
                                           i_z * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO
    END DO

    i_x            = 0
    DO i_y         = - Nh, Nh - 1
    DO i_z         = - Nh, Nh - 1
      k_dot_v_norm = k_dot_v_norm + hf* CDABS( i_x * v_x( i_x, i_y, i_z ) + &
                                               i_y * v_y( i_x, i_y, i_z ) + &
                                               i_z * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO

    i_x            = Nh
    DO i_y         = - Nh, Nh - 1
    DO i_z         = - Nh, Nh - 1
      k_dot_v_norm = k_dot_v_norm + hf* CDABS( i_x * v_x( i_x, i_y, i_z ) + &
                                               i_y * v_y( i_x, i_y, i_z ) + &
                                               i_z * v_z( i_x, i_y, i_z ) ) ** two
    END DO
    END DO

    k_dot_v_norm = DSQRT( k_dot_v_norm )

    IF (k_dot_v_norm .GT. tol_float ) THEN

      k_dot_v_error = 1

      debug_error = 1
      ! This will jump out of evolution loop, if caught during that.

      CALL print_error_incomp
      ! This will prompt error when checked
      ! REF-> <<< system_output >>>

    END IF
  END

  SUBROUTINE check_nan
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check if there is any NaN in the data.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    NaN_count  = 0

    DO i_x         = 0, Nh
    DO i_y         = -Nh, Nh - 1
    DO i_z         = -Nh, Nh - 1
      IF ( v_x( i_x, i_y, i_z ) .NE. v_x( i_x, i_y, i_z ) ) THEN
        NaN_count = NaN_count + 1
      END IF
    END DO
    END DO
    END DO

    i_x            = 0
    DO i_y         = - Nh, Nh - 1
    DO i_z         = - Nh, Nh - 1
      IF ( v_x( i_x, i_y, i_z ) .NE. v_x( i_x, i_y, i_z ) ) THEN
        NaN_count = NaN_count + 1
      END IF
    END DO
    END DO

    i_x            = Nh
    DO i_y         = - Nh, Nh - 1
    DO i_z         = - Nh, Nh - 1
      IF ( v_x( i_x, i_y, i_z ) .NE. v_x( i_x, i_y, i_z ) ) THEN
        NaN_count = NaN_count + 1
      END IF
    END DO
    END DO

    IF (NaN_count .NE. 0) THEN

      debug_error = 1
      ! This will jump out of evolution loop, if caught during that.

      CALL print_error_nan
      ! This will prompt error when checked
      ! REF-> <<< system_output >>>

    END IF

  END

END MODULE system_basicfunctions
