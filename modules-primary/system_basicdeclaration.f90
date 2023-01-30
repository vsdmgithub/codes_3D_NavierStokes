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
! MODULE: system_basicdeclaration
! LAST MODIFIED: 2 June 2021
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! BASIC VARIABLES AND ARRAYS FOR 3D NSE EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_basicdeclaration
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the global variables and arrays for the simulation space given values
! here, wheras temporary variables (if necessary) are declared within the SUBROUTINEs
! Further, each variable is classified based on where its purpose suits apt.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_auxilaries
  USE system_basicvariables
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  IMPLICIT NONE
  !TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

  CONTAINS

  SUBROUTINE read_input
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Read simulation parameters from a file 'input_file'
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER( LEN = 60 )::input_file

    input_file  = 'parameters.dat'
    ! This file contains all major Input parameters to be fed from outside file

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! R  E  A  D  I  N  G       I  N  P  U  T       F  I  L  E
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    OPEN( UNIT = 1001, FILE = TRIM( ADJUSTL(input_file) ) )

    READ( 1001, f_d8p4, ADVANCE ='yes')
    READ( 1001, f_i6,   ADVANCE ='yes')  N
  	! No of collocation points in physical space in one Dimension

    READ( 1001, f_d8p4,  ADVANCE ='yes')
    READ( 1001, f_d8p4,  ADVANCE ='yes')  time_total
    ! Total time to simulate

    READ( 1001, f_d8p4, ADVANCE ='yes')
    READ( 1001, f_i4,   ADVANCE ='yes')  no_of_saves
    ! No of saves

    READ( 1001, f_d8p4, ADVANCE ='yes')
    READ( 1001, f_i4,   ADVANCE ='yes')  no_of_3D_saves
    ! No of 3D saves

    CLOSE(1001)
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

  END

  SUBROUTINE init_global_variables
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !       Initialize all the variables used in the code and
  ! are explained too. Only variables are initialized here.
  ! Arrays are initialized in  another SUBROUTINE.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! S P A C E    A N D     T I M E
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    l_sys             = two_pi
    ! Length of periodic cube

    Nh                = INT( N / 2 )
    ! Maximum wavemuber in the First Brillouin zone.

    WRITE( N_char, f_i8 ) N
    ! Converting resolution value to character

    dx                = l_sys / DBLE( N )
    dy                = dx
    dz                = dx
    ! Grid distance

    dxdydz            = dx * dy * dz
    ! Grid volume

    N3                = DBLE( N * N * N )
    ! No of points in real space

    vol               = l_sys ** thr
    ! Volume of domain

    k_G               = FLOOR( DBLE( N ) / thr ) - 1
    k_G_2             = DBLE( k_G * k_G ) - one
    ! Truncation wavenumber shell radius.
    ! For k_i>k_G, modes are truncated to remove dealiasing error.

    k_max             = CEILING( DSQRT( thr ) * DBLE( Nh ) ) + 1
    ! Maximum shell no a mode in first Brilloun zone can get.

    tot_modes         = N**3

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! F L U I D D A T A
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    norm_factor       = one
    ! Normalization factor for energy - later changed so that initial energy is obtained.

    energy_initial    = one
    energy            = energy_initial
    energy_old        = energy
    ! Initial energy

    ! XXXXXXXXXXXXXXXXXXXXXX
    ! VISCOSITY_SELECTION:
    ! XXXXXXXXXXXXXXXXXXXXXX
    ! viscosity         = 4.0D0 * 1.0E-3 ! For N=256 ! THIS IS THE REFERENCE VISCOSITY
    viscosity         = ( 128.0D0 / DBLE( N ) ) * viscosity * 1.0E-3
    ! Viscosity of the system

    k_int             = 2
    ! Integral scale wavnumber

    k_for             = 4
    ! Wavenumbers below which forcing is implemented

    CALL find_diss_rate_ref( energy_initial, k_int, diss_rate_ref )
    ! Gives a estimate of dissipation rate fitting the Kolmogorov spectrum approximately
    diss_rate_ref     = 0.50D0 ! For N=256

    diss_rate_viscous = diss_rate_ref
    ! Assuming viscous dissipation as reference for start

    cfl_min           = 10
    ! - Courant-Friedrichs-Lewy (CFL) condition - CFL no is inverse of the above ratio
    ! No of steps (minimum) that should take to cross a grid

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! A U X I L A R Y
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    no_of_debug       = 5
    ! No of times that the program looks for any 'NaN' while marching forward in time.

    simulation_status = 0
    ! Meaning it is initializing , '1' means final, will be changed after the time_marching is DOne.

    forcing_status    = 1
    ! Tick that controls the forcing, 1 - YES, 0 - NO

    helicity_comp     = 0
    ! Tick that controls the computation of helicity, 1 - YES, 0 - NO

    velocitygrad_comp = 0
    ! Tick that controls the computation of gradient of velocity, 1 - YES, 0 - NO

    precheck_status   = 0
    ! Status that is checked before the time evolution

    nan_count         = 0
    ! No of NaN in v_x (if so)

    incomp_error      = 0
    ! Incompressibility error ( Sum(k.v(k)) exceeding tolerance )

    nan_error         = 0
    ! NaN found in the data - error

    viscosity_error   = 0
    ! Viscosity too small for this resolution

    debug_error       = 0
    ! some error found during debug (could be either nan or incomp)

	END

  SUBROUTINE compute_system_details
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !   This calculates the state of the system. Done at start and end.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION :: time_min

    u_rms             = DSQRT( two * energy / thr )
    ! 1D RMS Velocity

    time_grid         = dx / u_rms
    ! Time for particle to cross a grid

    k_kol             = FLOOR( ( diss_rate_viscous / ( viscosity ** thr) ) ** qtr )
    ! Kolmogorov wavenumber

    u_kol             = ( diss_rate_viscous * viscosity ) ** qtr
    ! Kolmogorov velocity

    time_kol          = DSQRT( viscosity / diss_rate_viscous )
    ! Kolmogorov timescale

    time_tur          = DSQRT( energy / diss_rate_viscous )
    ! Turbulent time scale

    resolving_power   = DBLE( k_G ) / DBLE( k_kol )
    ! Extent beyond which the smallest scale is resolved

    rey_no            = FLOOR( u_rms * l_sys / ( k_int * viscosity ) )
    tay_rey_no        = FLOOR( DSQRT( 20.0D0 * DBLE( rey_no ) / 3.0D0 ) )
    ! Reynold's number and Taylor scale Reynold's number.

    FIND_MIN_TIME: IF ( time_grid .LT. time_kol ) THEN
      time_min        = time_grid
    ELSE
      time_min        = time_kol
    END IF FIND_MIN_TIME

    dt_max            = time_min / DBLE( cfl_min )
    ! Maximum value of time step

    CALL find_CFL_timestep( time_min, dt_max, dt )
    ! Finds a time smaller than 'dt_max' in terms of '0.0..0p' , where p being '1' or '5'
    ! REF-> <<< system_auxilaries >>>

    cfl_system        = FLOOR( time_min / dt )
    ! CFL number of the system

    CALL time_to_step_convert(time_total,t_step_total,dt)
    ! returns the no of time_steps (\delta t) in a given time
    ! REF-> <<< system_auxilaries >>>

    t_step_save       = t_step_total / no_of_saves
    ! Determines how many time steps after the save has to be made.

    t_step_3D_save    = t_step_total / no_of_3D_saves
    ! Determines how many time steps after the PVD save has to be made.

    t_step_debug      = t_step_total / no_of_debug
    ! Determines how many time steps after the checking has to be made

    CALL step_to_time_convert(t_step_save,time_save,dt)
    ! Determines the saving time intervals
    ! REF-> <<< system_auxilaries >>>

    CALL time_to_step_convert( time_kol,t_step_kol,dt)
    ! REF-> <<< system_auxilaries >>>

	END
  SUBROUTINE init_global_arrays
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !   This defines few constant arrays that do not evolve. k^2 matrix,k matrix, projection matrix,
  ! shell no matrix, density_modes matrix.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::diff_ceiling
    DOUBLE PRECISION::kx,ky,kz
    DOUBLE PRECISION::k_per_2,k_mod,k_mod_2
    DOUBLE PRECISION::h_real_mod,h_imag_mod
    INTEGER(KIND=4)::ct


    CALL allocate_operators
    ! Allocates the arrays declared here.

    shell_no           = 0
    density_modes      = 0
    tot_forced_modes   = 0
    tot_active_modes   = 0
    integrating_factor = 0

    !  +++++++++++++++++++++++++++++++++
    !  A  X  I  S      G   R  I   D   S
    !  +++++++++++++++++++++++++++++++++
    LOOP_RX_000: DO i_x = 0, N-1
      axis(i_x) = DBLE( i_x ) * dx
      ! Location of grid points along principal axis
    END DO LOOP_RX_000

    LOOP_FX_101: DO j_x = 0, Nh
    LOOP_FY_101: DO j_y = -Nh, Nh - 1
    LOOP_FZ_101: DO j_z = -Nh, Nh - 1

      kx                    = DBLE( j_x )
      ky                    = DBLE( j_y )
      kz                    = DBLE( j_z )
      k_x( j_x, j_y, j_z )  = kx
      k_y( j_x, j_y, j_z )  = ky
      k_z( j_x, j_y, j_z )  = kz
      ! Just the k component matrix storing its grid points.

      k_mod_2               = kx**two + ky**two + kz**two
      k_per_2               = kx**two + ky**two
      k_2( j_x, j_y, j_z )  = k_mod_2
      k_mod                 = DSQRT( k_mod_2 )
      ! Square of distance to origin

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  P  R  O  J  E  C  T  I  O  N             M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Projection matrix \mathbb{P}_{ij}=\delta_{ij}-\frac{k_ik_j}{k^2}}

      PROJECTION_TENSOR: IF ( k_2 ( j_x, j_y, j_z ) .GT. tol ) THEN
        ! Checking IF k^2 is not too low, to cause NaN (this will happen only for (0,0,0))
        proj_xx( j_x, j_y, j_z ) = one - ( kx * kx ) / k_mod_2
        proj_yy( j_x, j_y, j_z ) = one - ( ky * ky ) / k_mod_2
        proj_zz( j_x, j_y, j_z ) = one - ( kz * kz ) / k_mod_2
        proj_xy( j_x, j_y, j_z ) = - ( kx * ky ) / k_mod_2
        proj_yz( j_x, j_y, j_z ) = - ( ky * kz ) / k_mod_2
        proj_zx( j_x, j_y, j_z ) = - ( kz * kx ) / k_mod_2
      ELSE
        proj_xx( j_x, j_y, j_z ) = + twothird
        proj_yy( j_x, j_y, j_z ) = + twothird
        proj_zz( j_x, j_y, j_z ) = + twothird
        proj_xy( j_x, j_y, j_z ) = - onethird
        proj_yz( j_x, j_y, j_z ) = - onethird
        proj_zx( j_x, j_y, j_z ) = - onethird
      END IF PROJECTION_TENSOR

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  S  H  E  L  L     N  O      M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      diff_ceiling                                   =  DBLE( CEILING( k_mod ) ) - k_mod
      ! figuring the decimal part of the |k|, so that it will be alloted to shell k or k+1

      SHELL_ALLOCATION: IF (diff_ceiling .GE. hf) THEN
        shell_no( j_x, j_y, j_z )                    =  FLOOR( k_mod )
      ELSE
        shell_no( j_x, j_y, j_z )                    =  CEILING( k_mod )
      END IF SHELL_ALLOCATION

      density_modes( shell_no( j_x, j_y, j_z ) ) =  density_modes( shell_no( j_x, j_y, j_z ) ) + 1
      ! counts no of grid points that belong to a particular shell, it should go as ~s^2

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  T  R  U  N  C  A  T  I  O  N  ,   S  H  E  L  L           M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Truncation mask matrix (multiply this with any spectral matrix to DO the truncation)
      TRUNCATION_MASK: IF ( k_mod_2 .LT. k_G_2 ) THEN

        truncator( j_x, j_y, j_z ) = one
        ! Spherical truncation filter matrix

        tot_active_modes           = tot_active_modes + 1
        ! Total no of active modes inside the truncation sphere.

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  D I S S I P A T I O N    I N T . F A C T O R    M A T R I X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        integrating_factor( j_x, j_y, j_z ) = DEXP( - viscosity * k_mod_2 * dt )

      ELSE

        truncator( j_x, j_y, j_z )                     =  zero
        ! Outised the truncation sphere.

      END IF TRUNCATION_MASK

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  F O R C I N G     S H E L L    C O U N T
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      FORCING_CHECK_101: IF ( forcing_status .EQ. 1 ) THEN
      FORCING_SHELL_COUNT_101:  IF ( shell_no( j_x, j_y, j_z ) .EQ. k_for )  THEN

        tot_forced_modes = tot_forced_modes + 1

      END IF FORCING_SHELL_COUNT_101
    END IF FORCING_CHECK_101


      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  H E L I C A L   B A S I S    M A T R I X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ! HELICAL_BASIS: IF ( k_per_2 .GT. tol ) THEN
      !
      !   h_real_mod               = root_2 * DSQRT( k_per_2 )
      !   h_imag_mod               = root_2 * DSQRT( k_per_2 * k_mod_2 )
      !   h_pos_x( j_x, j_y, j_z ) = ( - ky / h_real_mod ) + i * ( - kx * kz / h_imag_mod )
      !   h_pos_y( j_x, j_y, j_z ) = ( + kx / h_real_mod ) + i * ( - ky * kz / h_imag_mod )
      !   h_pos_z( j_x, j_y, j_z ) =                       + i * ( + k_per_2 / h_imag_mod )
      !   h_neg_x( j_x, j_y, j_z ) = ( - ky / h_real_mod ) - i * ( - kx * kz / h_imag_mod )
      !   h_neg_y( j_x, j_y, j_z ) = ( + kx / h_real_mod ) - i * ( - ky * kz / h_imag_mod )
      !   h_neg_z( j_x, j_y, j_z ) =                       - i * ( + k_per_2 / h_imag_mod )
      !
      ! ELSE
      !
      !   h_pos_x( j_x, j_y, j_z ) = +     one / root_2
      !   h_pos_y( j_x, j_y, j_z ) = + i * one / root_2
      !   h_pos_z( j_x, j_y, j_z ) =       zero
      !   h_neg_x( j_x, j_y, j_z ) = +     one / root_2
      !   h_neg_y( j_x, j_y, j_z ) = - i * one / root_2
      !   h_neg_z( j_x, j_y, j_z ) =       zero
      !
      ! END IF HELICAL_BASIS

    END DO LOOP_FZ_101
    END DO LOOP_FY_101
    END DO LOOP_FX_101

    k_2( 0, 0, 0 ) = one
    ! Just to make sure , when something is divided by k_2 , to avoid NaN. Numerator would be zero anyways,
    ! in most of such cases. So nothing wrong here.

    FORCING_CHECK_201: IF ( forcing_status .EQ. 1 ) THEN

      ALLOCATE( fkx( tot_forced_modes ) )
      ALLOCATE( fky( tot_forced_modes ) )
      ALLOCATE( fkz( tot_forced_modes ) )

      ct = 0
      LOOP_FX_201: DO j_x = 0, Nh
      LOOP_FY_201: DO j_y = -Nh, Nh - 1
      LOOP_FZ_201: DO j_z = -Nh, Nh - 1

        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  F O R C I N G     S H E L L    D A T A
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        FORCING_SHELL_COUNT_201:  IF ( shell_no( j_x, j_y, j_z ) .EQ. k_for )  THEN
          ct        = ct + 1
          fkx( ct ) = j_x
          fky( ct ) = j_y
          fkz( ct ) = j_z
        END IF FORCING_SHELL_COUNT_201

      END DO LOOP_FZ_201
      END DO LOOP_FY_201
      END DO LOOP_FX_201

    END IF FORCING_CHECK_201

  END

  SUBROUTINE allocate_operators
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate arrays which are constants basically, k2,truncator, shell no etc.,
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y         A  L  L  O  C  A  T  I  O  N .
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(axis(0:N-1))
    ALLOCATE(k_2(               0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(truncator(         0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(shell_no(          0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(integrating_factor(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(k_x(    0:Nh,-Nh:Nh-1,-Nh:Nh-1),k_y(    0:Nh,-Nh:Nh-1,-Nh:Nh-1),k_z(    0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(proj_xx(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_yy(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_zz(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(proj_xy(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_yz(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_zx(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(density_modes(0:k_max))
    ! ALLOCATE(h_pos_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),h_pos_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),h_pos_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ! ALLOCATE(h_neg_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),h_neg_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),h_neg_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

  END

  SUBROUTINE allocate_velocity
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays related to velocity(real and spectral) and it spectrum
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N    -   V  E  L  O  C  I  T  Y
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(u_x(0:N-1,0:N-1,0:N-1),u_y(0:N-1,0:N-1,0:N-1),u_z(0:N-1,0:N-1,0:N-1))
    ! ALLOCATE(u2_x(0:N-1,0:N-1,0:N-1),u2_y(0:N-1,0:N-1,0:N-1),u2_z(0:N-1,0:N-1,0:N-1))
    ALLOCATE(v_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ! ALLOCATE(v2_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v2_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v2_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(spectral_energy(0       :k_max))
    ALLOCATE(spectral_energy_avg(0   :k_max))
    ! ALLOCATE(spectral_enstrophy(0    :k_max))
    ! ALLOCATE(spectral_enstrophy_avg(0:k_max))
    ! ALLOCATE(spectral_helicity(0     :k_max))
    ! ALLOCATE(spectral_helicity_avg(0 :k_max))

	END

  SUBROUTINE allocate_vorticity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate arrays related to vorticity(real and spectral)
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N    -   V  O  R  T  I  C  I  T  Y
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(w_ux(0:N-1,0:N-1,0:N-1),w_uy(0:N-1,0:N-1,0:N-1),w_uz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(w_vx(0:Nh,-Nh:Nh-1,-Nh:Nh-1),w_vy(0:Nh,-Nh:Nh-1,-Nh:Nh-1),w_vz(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

  END

	SUBROUTINE deallocate_velocity
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(v_x,v_y,v_z)
		! DEALLOCATE(v2_x,v2_y,v2_z)
		DEALLOCATE(u_x,u_y,u_z)
		! DEALLOCATE(u2_x,u2_y,u2_z)
		DEALLOCATE(spectral_energy)
    DEALLOCATE(spectral_energy_avg)
		! DEALLOCATE(spectral_enstrophy)
		! DEALLOCATE(spectral_enstrophy_avg)
		! DEALLOCATE(spectral_helicity_avg)
		! DEALLOCATE(spectral_helicity)

	END

	SUBROUTINE deallocate_vorticity
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(w_ux,w_uy,w_uz)
		DEALLOCATE(w_vx,w_vy,w_vz)

	END

  SUBROUTINE deallocate_operators
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE(axis)
    DEALLOCATE(k_2,k_x,k_y,k_z)
    DEALLOCATE(truncator)
    DEALLOCATE(integrating_factor)
    DEALLOCATE(proj_xx,proj_yy,proj_zz)
    DEALLOCATE(proj_xy,proj_yz,proj_zx)
    DEALLOCATE(shell_no,density_modes)
    ! DEALLOCATE(h_pos_x,h_pos_y,h_pos_z)
    ! DEALLOCATE(h_neg_x,h_neg_y,h_neg_z)

    FORCING_CHECK_202: IF ( forcing_status .EQ. 1 ) THEN
      DEALLOCATE(fkx,fky,fkz)
    END IF FORCING_CHECK_202

	END

END MODULE system_basicdeclaration
