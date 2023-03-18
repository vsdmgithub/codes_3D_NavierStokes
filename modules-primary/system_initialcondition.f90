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
! MODULE: system_initialcondition
! LAST MODIFIED: 20 FEBRAURY 2023
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! INITIAL CONDITION FOR 3D NSE EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_initialcondition
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Different initial conditions are coded here
! 1. Exponentially decaying spectrum
! 2. Thermalized spectrum
! 3. Kolmogorov like spectrum
! 4. Taylor-Green Initial Condition
! 5. Kida Peltz Initial Condition
! 6. ABC Flow initial condition
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicdeclaration
  USE system_fftw
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  IMPLICIT  NONE

  CONTAINS

  SUBROUTINE init_initcondn
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Initialize initial condition for velocity, Choose one.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! Initializing the initial velocity (spectral) and projecting it so that the flow is incompressible.

    CALL IC_exp_decaying_spectrum(eng_0)
    ! Generic randomized initial condition, with energy mainly in integral scale (spectrally)

    ! CALL IC_TG(eng_0)
    ! TAYLOR_GREEN Initial condition - lots of symmetries (although solver is not using them)

    ! CALL IC_from_file
    ! Read from file (binary).
    ! *****Check whether file is available already.

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! UNCOMMENT TO TRUNCATE IF NEEDED - (most of I.C are already truncated)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    V_x = Trunc * V_x
    V_y = Trunc * V_y
    V_z = Trunc * V_z

  END

  SUBROUTINE IC_exp_decaying_spectrum(eng_inp)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! A typical initial condition IN which only few large modes are excited.
  ! CALL this to give 3 components of COMPLEX spectral velocity at spectral grid 'i_x,i_y,i_z'
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)  ::eng_inp
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION             ::phi,the
    DOUBLE PRECISION,DIMENSION(3)::phs
    DOUBLE PRECISION             ::V_k_mod,k_rat
    DOUBLE COMPLEX,DIMENSION(3)  ::V_k
    INTEGER(KIND=4)              ::ple_exp

    icn_type = 'EXP-DEC'

    CALL init_random_seed
    ! Randomizes seed for random numbers (in 'auxilary_functions' module )
    ! REF-> <<< system_auxilaries >>>

    ! k_int        = 2
    ! Integral scale wavenumber

    ple_exp = 2
    ! The power in the spectrum E(k) for k<k_int
    ! Generally either 2 or 4 or 6. But has to be a even number

    CALL normalization_exponential_spectrum( ple_exp, k_int, nrm_fac)
    ! Returns the 'nrm_fac', so that theoretically net energy is O(1).
    ! Additionally normalization to any energy can be done with nrm_fac at end
    ! REF-> <<< system_auxilaries >>>

    LOOP_FX: DO j_x = 0, Nh
    LOOP_FY: DO j_y = -Nh, Nh - 1
    LOOP_FZ: DO j_z = -Nh, Nh - 1
    IF ( Lapla( j_x, j_y, j_z ) .LT. k_tru_sqr ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(the)
      CALL RANDOM_NUMBER(phs)

      phi     = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      the     = DACOS( one - two * the )! Polar angle of \hat{u}_k vector
      phs     = two_pi * phs ! Phases of \hat{u}_k components

      k_rat   = DSQRT( Lapla( j_x, j_y, j_z) ) / DBLE(k_int)

      V_k_mod = nrm_fac * k_rat**( hf * ple_exp - 1 ) &
                * DEXP( - qtr * ple_exp * ( k_rat ** 4.0D0 ) )

      V_k(1)  = V_k_mod * DSIN( the ) * DCOS( phi ) * DCMPLX( DCOS( phs( 1 ) ), DSIN( phs( 1 ) ) )
      V_k(2)  = V_k_mod * DSIN( the ) * DSIN( phi ) * DCMPLX( DCOS( phs( 2 ) ), DSIN( phs( 2 ) ) )
      V_k(3)  = V_k_mod * DCOS( the ) * DCMPLX( DCOS( phs( 3 ) ),DSIN( phs( 3 ) ) )
      ! 3 COMPLEX values for spectral velocity

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  V  E  L  O  C  I  T  Y          P  R  O  J  E  C  T  O  R
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Assigning the velocity along with projection so that it is Incompressible -- u(k).k=0
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  N  O  T  E  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! The following steps ensure, that except for j_x=0 or Nh plane, all the other planes are given Initial velocity
      ! But those planes require special attention, becoz, their Inversion lies IN the same plane,
      ! so conjugates must be placed accordingly.

      IF (((j_x .NE. 0) .AND. (j_x .NE. Nh)) .OR. (j_z .LT. 0)) THEN
        V_x( j_x, j_y, j_z )     = Proj_xx( j_x, j_y, j_z ) * V_k(1) + Proj_xy( j_x, j_y, j_z) * V_k(2) + &
                                   Proj_zx( j_x, j_y, j_z ) * V_k(3)
        V_y( j_x, j_y, j_z )     = Proj_xy( j_x, j_y, j_z ) * V_k(1) + Proj_yy( j_x, j_y, j_z) * V_k(2) + &
                                   Proj_yz( j_x, j_y, j_z ) * V_k(3)
        V_z( j_x, j_y, j_z )     = Proj_zx( j_x, j_y, j_z ) * V_k(1) + Proj_yz( j_x, j_y, j_z) * V_k(2) + &
                                   Proj_zz( j_x, j_y, j_z ) * V_k(3)

      IF (((j_x .EQ. 0) .OR. (j_x .EQ. Nh)) .AND. ((j_y .NE. -Nh) .AND. (j_z .NE. -Nh))) THEN
        V_x( j_x, - j_y, - j_z ) = DCONJG( V_x( j_x, j_y, j_z ) )
        V_y( j_x, - j_y, - j_z ) = DCONJG( V_y( j_x, j_y, j_z ) )
        V_z( j_x, - j_y, - j_z ) = DCONJG( V_z( j_x, j_y, j_z ) )
      END IF

      ELSE IF ((j_z .EQ. 0) .AND. (j_y .LE. 0)) THEN
        V_x( j_x, j_y, j_z )     = Proj_xx( j_x, j_y, j_z ) * V_k(1) + Proj_xy( j_x, j_y, j_z ) * V_k(2) + &
                                   Proj_zx( j_x, j_y, j_z) * V_k(3)
        V_y( j_x, j_y, j_z )     = Proj_xy( j_x, j_y, j_z ) * V_k(1) + Proj_yy( j_x, j_y, j_z ) * V_k(2) + &
                                   Proj_yz( j_x, j_y, j_z) * V_k(3)
        V_z( j_x, j_y, j_z )     = Proj_zx( j_x, j_y, j_z ) * V_k(1) + Proj_yz( j_x, j_y, j_z ) * V_k(2) + &
                                   Proj_zz( j_x, j_y, j_z) * V_k(3)

        V_x( j_x, - j_y, j_z )   = DCONJG( V_x ( j_x, j_y, j_z ) )
        V_y( j_x, - j_y, j_z )   = DCONJG( V_y ( j_x, j_y, j_z ) )
        V_z( j_x, - j_y, j_z )   = DCONJG( V_z ( j_x, j_y, j_z ) )
      ELSE IF ((j_y .EQ. -Nh) .AND. (j_z .GT. 0)) THEN
        V_x(j_x,j_y,j_z)         = Proj_xx( j_x, j_y, j_z ) * V_k(1) + Proj_xy( j_x, j_y, j_z ) * V_k(2) + &
                                   Proj_zx( j_x, j_y, j_z) * V_k(3)
        V_y(j_x,j_y,j_z)         = Proj_xy( j_x, j_y, j_z ) * V_k(1) + Proj_yy( j_x, j_y, j_z ) * V_k(2) + &
                                   Proj_yz( j_x, j_y, j_z) * V_k(3)
        V_z(j_x,j_y,j_z)         = Proj_zx( j_x, j_y, j_z ) * V_k(1) + Proj_yz( j_x, j_y, j_z ) * V_k(2) + &
                                   Proj_zz( j_x, j_y, j_z) * V_k(3)
      END IF
      ! ----------------------------------------------------------------------------------------

    END IF
    END DO LOOP_FZ
    END DO LOOP_FY
    END DO LOOP_FX

    ! Making sure, that the average velocity is zero.
    V_x( 0, 0, 0 )    =     zero
    V_y( 0, 0, 0 )    =     zero
    V_z( 0, 0, 0 )    =     zero

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    nrm_fac = DSQRT( eng_inp / eng )
    ! Normalizing the nrm_fac, so that we get energy='eng_inp'

    V_x = V_x * nrm_fac
    V_y = V_y * nrm_fac
    V_z = V_z * nrm_fac

  END

  SUBROUTINE IC_TG(eng_inp)
  ! This is TAYLOR_GREEN Vortex initial condition which has a lot of symmetries
    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)::eng_inp
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX::v0

    icn_type = 'TAY-GRE'

    v0            = i / 8.0D0
    !-------- 'x' velocity--------------
    V_x(+1,+1,+1) = - v0
    V_x(+1,+1,-1) = - v0
    V_x(+1,-1,+1) = - v0
    V_x(+1,-1,-1) = - v0
    !-------- 'y' velocity--------------
    V_y(+1,+1,+1) = + v0
    V_y(+1,+1,-1) = + v0
    V_y(+1,-1,+1) = - v0
    V_y(+1,-1,-1) = - v0

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    nrm_fac = DSQRT( eng_inp / eng )
    ! Normalizing the nrm_fac, so that we get energy='eng_inp'

    v0           = nrm_fac * v0
    !-------- 'x' velocity--------------
    V_x(+1,+1,+1) = - v0
    V_x(+1,+1,-1) = - v0
    V_x(+1,-1,+1) = - v0
    V_x(+1,-1,-1) = - v0
    !-------- 'y' velocity--------------
    V_y(+1,+1,+1) = + v0
    V_y(+1,+1,-1) = + v0
    V_y(+1,-1,+1) = - v0
    V_y(+1,-1,-1) = - v0

  END

  SUBROUTINE IC_from_file
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Read a file for the velocity data in real space.
  ! Unformatted file
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN=80)::IC_file_name

    icn_type = 'INP-REA'

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! V   E  L  O  C  I  T  Y       I  N  P  U  T     F  I  L  E
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IC_file_name  = '../NSE_data/vel/velocity_' // TRIM( ADJUSTL( res_char ) ) // '.in'
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! OPEN( UNIT = 48, FILE = IC_file_name, FORM='unformatted', STATUS='old')
    ! '+++++++++++++++++++++'
    ! 'I.C FROM FILE'
    ! '+++++++++++++++++++++'
    ! READ(48) ((( U_x( i_x, i_y, i_z ), U_y( i_x, i_y, i_z ), U_z( i_x, i_y, i_z ), &
                 ! i_z = 0 , N - 1 ) i_y = 0 , N - 1 ) i_x = 0 , N - 1 )

    ! CLOSE(48)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    CALL fft_r2c( U_x, U_y, U_z, N, Nh, V_x, V_y, V_z )
    ! FFT spectral to real velocity

    !  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! UNCOMMENT TO USE WITHOUT NORMALIZATION
    !  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    nrm_fac = DSQRT( eng_0 / eng )
    ! Normalizing the nrm_fac, so that we get energy='eng_inp'

    V_x = V_x * nrm_fac
    V_y = V_y * nrm_fac
    V_z = V_z * nrm_fac
    !  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  END

  SUBROUTINE compute_energy_spectral_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to check the presence of NaN in your spectral velocity data (v(k)),
  ! and also the L2 norm or the Kinetic energy.
  ! NOTE: Count certain modes once, certain modes half (owing to 1/2 factor)
  ! in the first loop i_x=0 plane is left. later it is considered
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    eng      = zero
    hel   = zero

    DO j_x      =   1, Nh - 1
    DO j_y      = -Nh, Nh - 1
    DO j_z      = -Nh, Nh - 1
      eng    = eng + CDABS( V_x( j_x, j_y, j_z ) ) ** two + &
                           CDABS( V_y( j_x, j_y, j_z ) ) ** two + &
                           CDABS( V_z( j_x, j_y, j_z ) ) ** two
      hel = hel- two * &
      ( DIMAG( DCONJG( V_x( j_x, j_y, j_z ) ) * V_z( j_x, j_y, j_z ) ) * K_y( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( V_y( j_x, j_y, j_z ) ) * V_x( j_x, j_y, j_z ) ) * K_z( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( V_z( j_x, j_y, j_z ) ) * V_y( j_x, j_y, j_z ) ) * K_x( j_x, j_y, j_z ) )
    END DO
    END DO
    END DO

    j_x         =   0
    DO j_y      = - Nh, Nh - 1
    DO j_z      = - Nh, Nh - 1
      eng    = eng + hf * ( CDABS( V_x( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( V_y( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( V_z( j_x, j_y, j_z ) ) ** two )
      hel = hel- &
      ( DIMAG( DCONJG( V_x( j_x, j_y, j_z ) ) * V_z( j_x, j_y, j_z ) ) * K_y( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( V_y( j_x, j_y, j_z ) ) * V_x( j_x, j_y, j_z ) ) * K_z( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( V_z( j_x, j_y, j_z ) ) * V_y( j_x, j_y, j_z ) ) * K_x( j_x, j_y, j_z ) )
    END DO
    END DO

    j_x         =   Nh
    DO j_y      = - Nh, Nh - 1
    DO j_z      = - Nh, Nh - 1
      eng    = eng + hf * ( CDABS( V_x( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( V_y( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( V_z( j_x, j_y, j_z ) ) ** two )
      hel = hel- &
      ( DIMAG( DCONJG( V_x( j_x, j_y, j_z ) ) * V_z( j_x, j_y, j_z ) ) * K_y( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( V_y( j_x, j_y, j_z ) ) * V_x( j_x, j_y, j_z ) ) * K_z( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( V_z( j_x, j_y, j_z ) ) * V_y( j_x, j_y, j_z ) ) * K_x( j_x, j_y, j_z ) )
    END DO
    END DO

  END

  SUBROUTINE compute_projected_velocity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to project the spectral velocity, to make it incompressible
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   ::v_P_x,v_P_y,v_P_z
    ALLOCATE(v_P_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_P_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_P_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

    v_P_x = V_x
    v_P_y = V_y
    v_P_z = V_z

    V_x   = Proj_xx * v_P_x + Proj_xy * v_P_y + Proj_zx * v_P_z
    V_y   = Proj_xy * v_P_x + Proj_yy * v_P_y + Proj_yz * v_P_z
    V_z   = Proj_zx * v_P_x + Proj_yz * v_P_y + Proj_zz * v_P_z

    DEALLOCATE(v_P_x,v_P_y,v_P_z)

  END

END MODULE system_initialcondition
