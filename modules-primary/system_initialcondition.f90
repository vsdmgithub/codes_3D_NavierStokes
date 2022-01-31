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
! LAST MODIFIED: 21 June 2021
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

    ! CALL IC_exp_decaying_spectrum(energy_initial)
    ! Generic randomized initial condition, with energy mainly in integral scale (spectrally)

    ! CALL IC_Kolmogorov_spectrum(energy_initial)
    ! Creating Kolmogorov spectrum, inbuilt k^-5/3 spectrum

    ! CALL IC_from_file_spectral
    ! Read from file.
    ! *****Check whether file is available already.

    CALL IC_from_file_real
    ! Read from file.
    ! *****Check whether file is available already.

    ! CALL IC_from_file_real_unformatted
    ! Read from file (binary).
    ! *****Check whether file is available already.

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! UNCOMMENT TO TRUNCATE IF NEEDED - (most of I.C are already truncated)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    v_x = truncator * v_x
    v_y = truncator * v_y
    v_z = truncator * v_z

  END

  SUBROUTINE IC_exp_decaying_spectrum(energy_input)
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
    DOUBLE PRECISION,INTENT(IN)  ::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION             ::phi,theta
    DOUBLE PRECISION,DIMENSION(3)::ph
    DOUBLE PRECISION             ::V_k_mod,norm_const,k_ratio
    DOUBLE COMPLEX,DIMENSION(3)  ::V_k
    INTEGER(KIND=4)              ::ple_exp

    IC_type = 'EXP-DEC'

    CALL init_random_seed
    ! Randomizes seed for random numbers (in 'auxilary_functions' module )
    ! REF-> <<< system_auxilaries >>>

    k_int        = 2
    ! Integral scale wavenumber

    ple_exp = 2
    ! The power in the spectrum E(k) for k<k_int
    ! Generally either 2 or 4 or 6. But has to be a even number

    CALL normalization_exponential_spectrum( ple_exp, k_int, norm_const)
    ! Returns the 'norm_const', so that theoretically net energy is O(1).
    ! Additionally normalization to any energy can be done with norm_factor at end
    ! REF-> <<< system_auxilaries >>>

    LOOP_FX: DO j_x = 0, Nh
    LOOP_FY: DO j_y = -Nh, Nh - 1
    LOOP_FZ: DO j_z = -Nh, Nh - 1
    IF ( k_2( j_x, j_y, j_z ) .LT. k_G_2 ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(theta)
      CALL RANDOM_NUMBER(ph)

      phi     = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      theta   = DACOS( one - two * theta )! Polar angle of \hat{u}_k vector
      ph      = two_pi * ph ! Phases of \hat{u}_k components

      k_ratio = DSQRT( k_2( j_x, j_y, j_z) ) / DBLE(k_int)

      V_k_mod = norm_const * k_ratio**( hf * ple_exp - 1 ) &
                * DEXP( - qtr * ple_exp * ( k_ratio ** two ) )

      V_k(1)  = V_k_mod * DSIN( theta ) * DCOS( phi ) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
      V_k(2)  = V_k_mod * DSIN( theta ) * DSIN( phi ) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
      V_k(3)  = V_k_mod * DCOS( theta ) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
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
        v_x( j_x, j_y, j_z )     = proj_xx( j_x, j_y, j_z ) * V_k(1) + proj_xy( j_x, j_y, j_z) * V_k(2) + &
                                   proj_zx( j_x, j_y, j_z ) * V_k(3)
        v_y( j_x, j_y, j_z )     = proj_xy( j_x, j_y, j_z ) * V_k(1) + proj_yy( j_x, j_y, j_z) * V_k(2) + &
                                   proj_yz( j_x, j_y, j_z ) * V_k(3)
        v_z( j_x, j_y, j_z )     = proj_zx( j_x, j_y, j_z ) * V_k(1) + proj_yz( j_x, j_y, j_z) * V_k(2) + &
                                   proj_zz( j_x, j_y, j_z ) * V_k(3)

      IF (((j_x .EQ. 0) .OR. (j_x .EQ. Nh)) .AND. ((j_y .NE. -Nh) .AND. (j_z .NE. -Nh))) THEN
        v_x( j_x, - j_y, - j_z ) = DCONJG( v_x( j_x, j_y, j_z ) )
        v_y( j_x, - j_y, - j_z ) = DCONJG( v_y( j_x, j_y, j_z ) )
        v_z( j_x, - j_y, - j_z ) = DCONJG( v_z( j_x, j_y, j_z ) )
      END IF

      ELSE IF ((j_z .EQ. 0) .AND. (j_y .LE. 0)) THEN
        v_x( j_x, j_y, j_z )     = proj_xx( j_x, j_y, j_z ) * V_k(1) + proj_xy( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_zx( j_x, j_y, j_z) * V_k(3)
        v_y( j_x, j_y, j_z )     = proj_xy( j_x, j_y, j_z ) * V_k(1) + proj_yy( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_yz( j_x, j_y, j_z) * V_k(3)
        v_z( j_x, j_y, j_z )     = proj_zx( j_x, j_y, j_z ) * V_k(1) + proj_yz( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_zz( j_x, j_y, j_z) * V_k(3)

        v_x( j_x, - j_y, j_z )   = DCONJG( v_x ( j_x, j_y, j_z ) )
        v_y( j_x, - j_y, j_z )   = DCONJG( v_y ( j_x, j_y, j_z ) )
        v_z( j_x, - j_y, j_z )   = DCONJG( v_z ( j_x, j_y, j_z ) )
      ELSE IF ((j_y .EQ. -Nh) .AND. (j_z .GT. 0)) THEN
        v_x(j_x,j_y,j_z)         = proj_xx( j_x, j_y, j_z ) * V_k(1) + proj_xy( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_zx( j_x, j_y, j_z) * V_k(3)
        v_y(j_x,j_y,j_z)         = proj_xy( j_x, j_y, j_z ) * V_k(1) + proj_yy( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_yz( j_x, j_y, j_z) * V_k(3)
        v_z(j_x,j_y,j_z)         = proj_zx( j_x, j_y, j_z ) * V_k(1) + proj_yz( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_zz( j_x, j_y, j_z) * V_k(3)
      END IF
      ! ----------------------------------------------------------------------------------------

    END IF
    END DO LOOP_FZ
    END DO LOOP_FY
    END DO LOOP_FX

    ! Making sure, that the average velocity is zero.
    v_x( 0, 0, 0 )    =     zero
    v_y( 0, 0, 0 )    =     zero
    v_z( 0, 0, 0 )    =     zero

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

  END

  SUBROUTINE IC_Kolmogorov_spectrum(energy_input)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! A typical initial condition with kolmogorov spectrum model, referred from Pope's Turbulence.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)  ::energy_input
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION             ::phi,theta
    DOUBLE PRECISION,DIMENSION(3)::ph
    DOUBLE PRECISION             ::V_k_mod,norm_const,k_ratio
    DOUBLE PRECISION             ::factor_integral,factor_dissipation
    DOUBLE PRECISION             ::eleven_by_twelve
    DOUBLE COMPLEX,DIMENSION(3)  ::V_k
    INTEGER(KIND=4)              ::ple_exp

    IC_type = 'KOL-INE'

    CALL init_random_seed
    ! Randomizes seed for random numbers (in 'auxilary_functions' module )
    ! REF-> <<< system_auxilaries >>>

    ! k_int      = FLOOR( DLOG( DBLE( N ) ) / DLOG( 4.0D0 ) ) - 1
    ! Integral scale wavenumber

    ! k_kolmo         = FLOOR ( two * DBLE( N ) / 9.0D0 )
    ! End of kolmogorov spectrum

    ple_exp = 2
    ! The power in the spectrum E(k) for k<k_int
    ! Generally either 2 or 4 or 6. But has to be a even number

    eleven_by_twelve  = 11.0D0 / 12.0D0

    CALL normalization_kolmogorov_spectrum( k_int, k_kol, norm_const)
    ! Returns the 'norm_const', so that theoretically net energy is O(1).
    ! Additionally normalization to any energy can be done with norm_factor
    ! REF-> <<< system_auxilaries >>>
    LOOP_FX: DO j_x = 0, Nh
    LOOP_FY: DO j_y = -Nh, Nh - 1
    LOOP_FZ: DO j_z = -Nh, Nh - 1

    IF ( k_2( j_x, j_y, j_z ) .LT. k_G_2 ) THEN

      CALL RANDOM_NUMBER(phi)
      CALL RANDOM_NUMBER(theta)
      CALL RANDOM_NUMBER(ph)

      phi                = two_pi * phi ! Azimuthal angle of \hat{u}_k vector
      theta              = DACOS( one - two * theta )! Polar angle of \hat{u}_k vector
      ph                 = two_pi * ph ! Phases of \hat{u}_k components

      k_ratio            = DSQRT( k_2( j_x, j_y, j_z) ) / DBLE(k_int)
      factor_integral    = one
      factor_dissipation = one

      IF ( k_ratio .LT. 1 ) THEN
        CALL kolmogorov_spectrum_integralscale_subpart(k_ratio,ple_exp,factor_integral)
        ! REF-> <<< system_auxilaries >>>
      END IF

      k_ratio            = DSQRT( k_2( j_x, j_y, j_z) ) / DBLE(k_kol)

      IF ( k_ratio .GT. hf ) THEN
        CALL kolmogorov_spectrum_dissipationscale_subpart(k_ratio,factor_dissipation)
        ! REF-> <<< system_auxilaries >>>
      END IF

      V_k_mod            = norm_const * ( k_2( j_x, j_y, j_z ) ** ( - eleven_by_twelve ) ) &
                          * factor_integral * factor_dissipation

      V_k(1)             = V_k_mod * DSIN( theta ) * DCOS( phi ) * DCMPLX( DCOS( ph( 1 ) ), DSIN( ph( 1 ) ) )
      V_k(2)             = V_k_mod * DSIN( theta ) * DSIN( phi ) * DCMPLX( DCOS( ph( 2 ) ), DSIN( ph( 2 ) ) )
      V_k(3)             = V_k_mod * DCOS( theta ) * DCMPLX( DCOS( ph( 3 ) ),DSIN( ph( 3 ) ) )
      ! 3 COMPLEX values for spectral velocity

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  V  E  L  O  C  I  T  Y          P  R  O  J  E  C  T  O  R
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Assigning the velocity along with projection so that it is INcompressible -- u(k).k=0
      ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  N  O  T  E  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      ! The following steps ensure, that except for j_x=0 or Nh plane, all the other planes are given INitial velocity
      ! But those planes require special attention, becoz, their INversion lies IN the same plane,
      ! so conjugates must be placed accordINgly.

      IF (((j_x .NE. 0) .AND. (j_x .NE. Nh)) .OR. (j_z .LT. 0)) THEN
        v_x( j_x, j_y, j_z )     = proj_xx( j_x, j_y, j_z ) * V_k(1) + proj_xy( j_x, j_y, j_z) * V_k(2) + &
                                   proj_zx( j_x, j_y, j_z ) * V_k(3)
        v_y( j_x, j_y, j_z )     = proj_xy( j_x, j_y, j_z ) * V_k(1) + proj_yy( j_x, j_y, j_z) * V_k(2) + &
                                   proj_yz( j_x, j_y, j_z ) * V_k(3)
        v_z( j_x, j_y, j_z )     = proj_zx( j_x, j_y, j_z ) * V_k(1) + proj_yz( j_x, j_y, j_z) * V_k(2) + &
                                   proj_zz( j_x, j_y, j_z ) * V_k(3)

      IF (((j_x .EQ. 0) .OR. (j_x .EQ. Nh)) .AND. ((j_y .NE. -Nh) .AND. (j_z .NE. -Nh))) THEN
        v_x( j_x, - j_y, - j_z ) = DCONJG( v_x( j_x, j_y, j_z ) )
        v_y( j_x, - j_y, - j_z ) = DCONJG( v_y( j_x, j_y, j_z ) )
        v_z( j_x, - j_y, - j_z ) = DCONJG( v_z( j_x, j_y, j_z ) )
      END IF

      ELSE IF ((j_z .EQ. 0) .AND. (j_y .LE. 0)) THEN
        v_x( j_x, j_y, j_z )     = proj_xx( j_x, j_y, j_z ) * V_k(1) + proj_xy( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_zx( j_x, j_y, j_z) * V_k(3)
        v_y( j_x, j_y, j_z )     = proj_xy( j_x, j_y, j_z ) * V_k(1) + proj_yy( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_yz( j_x, j_y, j_z) * V_k(3)
        v_z( j_x, j_y, j_z )     = proj_zx( j_x, j_y, j_z ) * V_k(1) + proj_yz( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_zz( j_x, j_y, j_z) * V_k(3)

        v_x( j_x, - j_y, j_z )   = DCONJG( v_x ( j_x, j_y, j_z ) )
        v_y( j_x, - j_y, j_z )   = DCONJG( v_y ( j_x, j_y, j_z ) )
        v_z( j_x, - j_y, j_z )   = DCONJG( v_z ( j_x, j_y, j_z ) )
      ELSE IF ((j_y .EQ. -Nh) .AND. (j_z .GT. 0)) THEN
        v_x(j_x,j_y,j_z)         = proj_xx( j_x, j_y, j_z ) * V_k(1) + proj_xy( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_zx( j_x, j_y, j_z) * V_k(3)
        v_y(j_x,j_y,j_z)         = proj_xy( j_x, j_y, j_z ) * V_k(1) + proj_yy( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_yz( j_x, j_y, j_z) * V_k(3)
        v_z(j_x,j_y,j_z)         = proj_zx( j_x, j_y, j_z ) * V_k(1) + proj_yz( j_x, j_y, j_z ) * V_k(2) + &
                                   proj_zz( j_x, j_y, j_z) * V_k(3)
      END IF
      ! ----------------------------------------------------------------------------------------

    END IF
    END DO LOOP_FZ
    END DO LOOP_FY
    END DO LOOP_FX

    ! Making sure, that the average velocity is zero.
    v_x( 0, 0, 0 )    =     zero
    v_y( 0, 0, 0 )    =     zero
    v_z( 0, 0, 0 )    =     zero

    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_input / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor

  END

  SUBROUTINE helical_decomposition(A_pos,A_neg)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Decompose the spectral field into helical modes
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! TRANSFER` VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,INTENT(IN)  ::A_pos,A_neg
    ! _________________________
    ! LOCAL ARRAYS
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE:: v_pos,v_neg

    ALLOCATE(v_pos(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(v_neg(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

    ! v_pos = v_x * h_neg_x + v_y * h_neg_y + v_z * h_neg_z
    ! v_neg = v_x * h_pos_x + v_y * h_pos_y + v_z * h_pos_z

    ! v_x = A_pos * v_pos * h_pos_x + A_neg * v_neg * h_neg_x
    ! v_y = A_pos * v_pos * h_pos_y + A_neg * v_neg * h_neg_y
    ! v_z = A_pos * v_pos * h_pos_z + A_neg * v_neg * h_neg_z

    DEALLOCATE(v_pos,v_neg)

  END

  SUBROUTINE IC_from_file_spectral
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Read a file for the velocity data in spectral space.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::real_part,imag_part
    CHARACTER(LEN=200)::IC_file_name

    IC_type = 'INP-SPE'

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! V   E  L  O  C  I  T  Y       I  N  P  U  T     F  I  L  E
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IC_file_name  = '../NSE_data/vel/spectral_velocity_' // TRIM( ADJUSTL( N_char ) ) // '_input.dat'
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 43, FILE = IC_file_name )

    ! '+++++++++++++++++++++'
    ! 'I.C FROM FILE'
    ! '+++++++++++++++++++++'

    LOOP_FX: DO j_x = 0, Nh
    LOOP_FY: DO j_y = -Nh, Nh - 1
    LOOP_FZ: DO j_z = -Nh, Nh - 1

      READ(43,f_c32p17,ADVANCE='NO') real_part, imag_part
      v_x( j_x, j_y, j_z ) = DCMPLX( real_part, imag_part )
      READ(43,f_c32p17,ADVANCE='NO') real_part, imag_part
      v_y( j_x, j_y, j_z ) = DCMPLX( real_part, imag_part )
      READ(43,f_c32p17,ADVANCE='YES') real_part, imag_part
      v_z( j_x, j_y, j_z ) = DCMPLX( real_part, imag_part )

    END DO LOOP_FZ
    END DO LOOP_FY
    END DO LOOP_FX

    CLOSE(43)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE IC_from_file_real
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Read a file for the velocity data in real space.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN=200)::IC_file_name

    IC_type = 'INP-REA'

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! V   E  L  O  C  I  T  Y       I  N  P  U  T     F  I  L  E
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IC_file_name  = '../vel/velocity_N' // TRIM( ADJUSTL( N_char ) ) // '_V8.dat'
    ! IC_file_name  = '../data_NSE_steadystate/N256/run_EC/3D_data/velocity_256_t_10.0000.dat'
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 44, FILE = IC_file_name )

    ! '+++++++++++++++++++++'
    ! 'I.C FROM FILE'
    ! '+++++++++++++++++++++'

    DO i_x = 0 , N - 1
    DO i_y = 0 , N - 1
    DO i_z = 0 , N - 1

      READ(44,f_d32p17,ADVANCE = 'NO')  u_x( i_x, i_y, i_z)
      READ(44,f_d32p17,ADVANCE = 'NO')  u_y( i_x, i_y, i_z)
      READ(44,f_d32p17,ADVANCE = 'YES') u_z( i_x, i_y, i_z)

    END DO
    END DO
    END DO

    CLOSE(44)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    CALL fft_r2c( u_x, u_y, u_z, N, Nh, v_x, v_y, v_z )
    ! FFT spectral to real velocity

    !  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! UNCOMMENT TO USE WITHOUT NORMALIZATION
    !  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_initial / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor
    !  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  END

  SUBROUTINE IC_from_file_real_unformatted
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

    IC_type = 'INP-REA'

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! V   E  L  O  C  I  T  Y       I  N  P  U  T     F  I  L  E
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    IC_file_name  = '../NSE_data/vel/velocity_' // TRIM( ADJUSTL( N_char ) ) // '.in'
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! OPEN( UNIT = 48, FILE = IC_file_name, FORM='unformatted', STATUS='old')
    ! '+++++++++++++++++++++'
    ! 'I.C FROM FILE'
    ! '+++++++++++++++++++++'
    ! READ(48) ((( u_x( i_x, i_y, i_z ), u_y( i_x, i_y, i_z ), u_z( i_x, i_y, i_z ), &
                 ! i_z = 0 , N - 1 ) i_y = 0 , N - 1 ) i_x = 0 , N - 1 )

    ! CLOSE(48)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    CALL fft_r2c( u_x, u_y, u_z, N, Nh, v_x, v_y, v_z )
    ! FFT spectral to real velocity

    !  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! UNCOMMENT TO USE WITHOUT NORMALIZATION
    !  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    CALL compute_energy_spectral_data
    ! Gets the energy from spectral space

    norm_factor = DSQRT( energy_initial / energy )
    ! Normalizing the norm_factor, so that we get energy='energy_input'

    v_x = v_x * norm_factor
    v_y = v_y * norm_factor
    v_z = v_z * norm_factor
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

    energy      = zero
    helicity    = zero

    DO j_x      =   1, Nh - 1
    DO j_y      = -Nh, Nh - 1
    DO j_z      = -Nh, Nh - 1
      energy    = energy + CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                           CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                           CDABS( v_z( j_x, j_y, j_z ) ) ** two
      helicity  = helicity - two * &
      ( DIMAG( DCONJG( v_x( j_x, j_y, j_z ) ) * v_z( j_x, j_y, j_z ) ) * k_y( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( v_y( j_x, j_y, j_z ) ) * v_x( j_x, j_y, j_z ) ) * k_z( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( v_z( j_x, j_y, j_z ) ) * v_y( j_x, j_y, j_z ) ) * k_x( j_x, j_y, j_z ) )
    END DO
    END DO
    END DO

    j_x         =   0
    DO j_y      = - Nh, Nh - 1
    DO j_z      = - Nh, Nh - 1
      energy    = energy + hf * ( CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( v_z( j_x, j_y, j_z ) ) ** two )
      helicity  = helicity - &
      ( DIMAG( DCONJG( v_x( j_x, j_y, j_z ) ) * v_z( j_x, j_y, j_z ) ) * k_y( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( v_y( j_x, j_y, j_z ) ) * v_x( j_x, j_y, j_z ) ) * k_z( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( v_z( j_x, j_y, j_z ) ) * v_y( j_x, j_y, j_z ) ) * k_x( j_x, j_y, j_z ) )
    END DO
    END DO

    j_x         =   Nh
    DO j_y      = - Nh, Nh - 1
    DO j_z      = - Nh, Nh - 1
      energy    = energy + hf * ( CDABS( v_x( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( v_y( j_x, j_y, j_z ) ) ** two + &
                                  CDABS( v_z( j_x, j_y, j_z ) ) ** two )
      helicity  = helicity - &
      ( DIMAG( DCONJG( v_x( j_x, j_y, j_z ) ) * v_z( j_x, j_y, j_z ) ) * k_y( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( v_y( j_x, j_y, j_z ) ) * v_x( j_x, j_y, j_z ) ) * k_z( j_x, j_y, j_z ) + &
        DIMAG( DCONJG( v_z( j_x, j_y, j_z ) ) * v_y( j_x, j_y, j_z ) ) * k_x( j_x, j_y, j_z ) )
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

    v_P_x = v_x
    v_P_y = v_y
    v_P_z = v_z

    v_x   = proj_xx * v_P_x + proj_xy * v_P_y + proj_zx * v_P_z
    v_y   = proj_xy * v_P_x + proj_yy * v_P_y + proj_yz * v_P_z
    v_z   = proj_zx * v_P_x + proj_yz * v_P_y + proj_zz * v_P_z

    DEALLOCATE(v_P_x,v_P_y,v_P_z)

  END

END MODULE system_initialcondition
