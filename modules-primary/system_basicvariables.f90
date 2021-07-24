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
! MODULE: system_basicvariables
! LAST MODIFIED: 2 June 2021
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! BASIC VARIABLES AND ARRAYS FOR 3D NSE EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_basicvariables
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the global variables and arrays for the simulation space are declared and given values
! here, wheras temporary variables (IF necessary) are declared within the SUBROUTINEs
! Further, each variable is classified based on where its purpose suits apt.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_constants
  USE system_auxilaries
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  IMPLICIT NONE
  ! _________________________
  ! REAL SPACE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION ::length,dx,dy,dz
  DOUBLE PRECISION ::N3,vol,dxdydz
  DOUBLE PRECISION ::k_G_2
  INTEGER(KIND=4)  ::N,Nh,k_G,k_max
  INTEGER(KIND=4)  ::i_x,i_y,i_z
  INTEGER(KIND=4)  ::j_x,j_y,j_z
  INTEGER(KIND=4)  ::k_no
  INTEGER(KIND=4)  ::tot_active_modes
  ! _________________________
  ! TIME VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION ::dt,dt_max
  DOUBLE PRECISION ::time_total,time_now,time_grid
  DOUBLE PRECISION ::time_save
  INTEGER (KIND=4) ::t_step,t_step_total,t_step_save,no_of_saves
  INTEGER (KIND=4) ::t_step_PVD_save,no_of_PVD_saves
  INTEGER (KIND=4) ::CFL_min,CFL_system
  ! _________________________
  ! FLUID VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION ::initial_circ
  DOUBLE PRECISION ::resolving_power
  DOUBLE PRECISION ::viscosity
  DOUBLE PRECISION ::diss_rate_ref, diss_rate, diss_rate_viscous
  DOUBLE PRECISION ::U_kolmo, time_kolmo
  DOUBLE PRECISION ::time_turb
  DOUBLE PRECISION ::k_F_2, forced_energy, force_factor
  INTEGER(KIND=4)  ::k_F
  INTEGER(KIND=4)  ::forcing_status
  INTEGER(KIND=4)  ::t_step_kolmo
  INTEGER(KIND=4)  ::k_integral,k_kolmo
  INTEGER(KIND=4)  ::Rey_No,Tay_Rey_No
  ! _________________________
  ! FUNCTION VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION ::energy
  DOUBLE PRECISION ::energy_old
  DOUBLE PRECISION ::helicity
  DOUBLE PRECISION ::energy_initial,U_rms
  DOUBLE PRECISION ::k_dot_v_norm
  DOUBLE PRECISION ::enstrophy
  DOUBLE PRECISION ::norm_factor
  DOUBLE PRECISION ::energy_mode
  DOUBLE PRECISION ::enstrophy_mode
  DOUBLE PRECISION ::helicity_mode
  DOUBLE COMPLEX   ::helicity_mode_complex
  ! _________________________
  ! CHARACTERS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  CHARACTER( LEN = 10) :: N_char
  CHARACTER( LEN = 20) :: IC_type
  CHARACTER( LEN = 3 ) :: run_code
  CHARACTER( LEN = 3 ) :: test_code
  CHARACTER( LEN = 3 ) :: solver_alg
  ! _________________________
  ! DEBUGGING VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  ::state_sim,check_status,debug_error
  INTEGER(KIND=4)  ::NaN_count,k_dot_v_error
  INTEGER(KIND=4)  ::no_of_debug,t_step_debug
  ! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::axis
  ! Arrays that store the grid
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::u_x,u_y,u_z
  ! Real velocity matrix  (updated after every time step)
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::w_ux,w_uy,w_uz
  ! Real vorticity
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::proj_xx,proj_yy,proj_zz
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::proj_xy,proj_yz,proj_zx
  ! projection operators matrix \mathbb{P}_{ij}=\delta_{ij}-\frac{k_ik_j}{k^2}}
  ! _________________________________________
  ! FOURIER SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4),DIMENSION(:,:,:),ALLOCATABLE  ::shell_no
  ! Every grid has its modulus, |k|
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::k_x,k_y,k_z,k_2
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::truncator
  ! wavenumber,truncator matrix
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::diss_Ifactor
  ! Integration factor for every time step in the dissipation
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   :: F_k_x, F_k_y, F_k_z
  ! Forcing in spectral space to match the dissipation rate
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   ::h_pos_x,h_pos_y,h_pos_z
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   ::h_neg_x,h_neg_y,h_neg_z
  ! Helical basis
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   :: v_x,v_y,v_z
  ! Spectral velocity matrix (will be updated after every time step)
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE   :: w_vx,w_vy,w_vz
  ! Spectral vorticity
  ! _________________________________________
  ! FOURIER (SHELL AVG) - SPECTRAL ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::spectral_energy,spectral_energy_avg
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::spectral_enstrophy,spectral_enstrophy_avg
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::spectral_helicity,spectral_helicity_avg
  ! Spectral data Integrated over spherical shell of radius k_mod
  INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE      ::count_modes_shell
  ! This counts no of modes that have (k-1)<|\vec{k}|<k

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
    READ( 1001, f_i4,   ADVANCE ='yes')  no_of_PVD_saves
    ! No of PVD saves

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

    length          = two_pi
    ! Length of periodic cube

    Nh              = INT( N / 2 )
    ! Maximum wavemuber in the First Brillouin zone.

    WRITE( N_char, f_i8 ) N
    ! converting resolution value to CHARACTER

    dx              = length / DBLE( N )
    dy              = dx
    dz              = dx
    ! Grid distance

    dxdydz          = dx * dy * dz
    ! Grid volume

    N3              = DBLE( N * N * N )
    ! No of poINts IN real space

    vol             = length ** thr
    ! Volume of domain

    k_G             = FLOOR( DBLE( N ) / thr ) - 1
    k_G_2           = DBLE( k_G * k_G ) - one
    ! Truncation wavenumber shell radius. For k_i>k_G, modes are truncated to remove dealiasing error.

    k_max           = CEILING( DSQRT( thr ) * DBLE( Nh ) ) + 1
    ! Maximum shell no it can get.

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! F L U I D D A T A
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    norm_factor     = one
    ! Normalization factor for energy - later changed so that initial energy is obtained.

    energy_initial  = one
    energy_old      = energy_initial
    ! Initial energy of the system

    U_rms           = DSQRT( two * energy_initial / thr )
    ! RMS Velocity

    time_grid       = dx / U_rms
    ! Time scale for particle to cross a grid

    CFL_min         = 10
    ! - Courant-Friedrichs-Lewy (CFL) condition - CFL no is inverse of the above ratio
    ! No of steps (minimum) that should take to cross a grid

    viscosity       = 20.0D0 * 10.0D0 ** ( - 3.0D0)
    ! Viscosity of the system

    k_integral      = 2
    ! Integral scale

    k_F             = 3
    k_F_2           = DBLE( k_F * k_F )
    ! Wavenumbers below which forcing is implemented

    CALL find_diss_rate_ref( energy_initial, k_integral, diss_rate_ref )
    ! Gives a estimate of dissipation rate fitting the Kolmogorov spectrum approximately

    k_kolmo         = FLOOR( ( diss_rate_ref / ( viscosity ** thr) ) ** qtr )
    ! Kolmogorov wavenumber

    U_kolmo         = ( diss_rate_ref * viscosity ) ** qtr
    ! Kolmogorov velocity

    time_kolmo      = DSQRT( viscosity / diss_rate_ref )
    ! Kolmogorov timescale

    time_turb       = DSQRT( energy_initial / diss_rate_ref )
    ! Turbulent time scale

    resolving_power = DBLE( k_G ) / DBLE( k_kolmo )
    ! Extent beyond which the smallest scale is resolved

    Rey_No          = FLOOR( DSQRT( energy_initial ) * length / viscosity )
    Tay_Rey_No      = FLOOR( 20.0D0 * DBLE( Rey_no ) / 3.0D0 )
    ! Reynold's no

    IF ( time_grid .LT. time_kolmo ) THEN

      dt_max        = time_grid / DBLE( CFL_min )
      ! Maximum value of time step

      CALL find_CFL_timestep( time_grid, dt_max, dt )
      ! Finds a time smaller than 'dt_max' in terms of '0.0..0p' , where p being '1' or '5'
      ! REF-> <<< system_auxilaries >>>

      CFL_system    = FLOOR( time_grid / dt )

    ELSE

      dt_max        = time_kolmo / DBLE( CFL_min )
      ! Maximum value of time step

      CALL find_CFL_timestep( time_kolmo, dt_max, dt )
      ! Finds a time smaller than 'dt_max' in terms of '0.0..0p' , where p being '1' or '5'
      ! REF-> <<< system_auxilaries >>>

      CFL_system    = FLOOR( time_kolmo / dt )

    END IF

    CALL time_to_step_convert(time_total,t_step_total,dt)
    ! returns the no of time_steps (\delta t) in a given time
    ! REF-> <<< system_auxilaries >>>

    t_step_save     = t_step_total / no_of_saves
    ! Determines how many time steps after the save has to be made.

    t_step_PVD_save = t_step_total / no_of_PVD_saves
    ! Determines how many time steps after the PVD save has to be made.

    CALL step_to_time_convert(t_step_save,time_save,dt)
    ! Determines the saving time intervals
    ! REF-> <<< system_auxilaries >>>

    CALL time_to_step_convert( time_kolmo, t_step_kolmo, dt )
    ! REF-> <<< system_auxilaries >>>

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! A U X I L A R Y
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    no_of_debug     = 10
    ! No of times that the program looks for any 'NaN' while marching forward in time.

    t_step_debug    = t_step_total / no_of_debug
    ! Determines how many time steps after the checking has to be made

    state_sim       = 0
    ! Meaning it is initializing , '1' means final, will be changed after the time_marching is DOne.

    forcing_status  = 0
    ! Tick that controls the forcing

    check_status    = 0
    ! Status that is checked before the time evolution

    NaN_count       = 0
    ! No of NaN IN v_x

    k_dot_v_error   = 0
    ! Incompressibility error

    debug_error     = 0
    ! Error found in debug, either NaN_count or k_dot_v_error

	END

  SUBROUTINE init_global_arrays
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !   This defines few arrays that do not evolve. k^2 matrix,k matrix, projection matrix,
  ! shell no matrix, count_modes_shell matrix.
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

    CALL allocate_operators
    ! Allocates the arrays declared here.

    shell_no          =    0
    count_modes_shell =    0

    !  +++++++++++++++++++++++++++++++++
    !  A  X  I  S      G   R  I   D   S
    !  +++++++++++++++++++++++++++++++++
    DO i_y = 0, N-1
      axis(i_y) = DBLE( i_y ) * dy
      ! Location of grid points along principal axis
    END DO

    DO i_x = 0, Nh
    DO i_y = -Nh, Nh - 1
    DO i_z = -Nh, Nh - 1

      kx                    = DBLE( i_x )
      ky                    = DBLE( i_y )
      kz                    = DBLE( i_z )
      k_x( i_x, i_y, i_z )  = kx
      k_y( i_x, i_y, i_z )  = ky
      k_z( i_x, i_y, i_z )  = kz
      ! Just the k component matrix storing its grid points.

      k_mod_2               = kx**two + ky**two + kz**two
      k_per_2               = kx**two + ky**two
      k_2( i_x, i_y, i_z )  = k_mod_2
      k_mod                 = DSQRT( k_mod_2 )
      ! Square of distance to origin

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  P  R  O  J  E  C  T  I  O  N             M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Projection matrix \mathbb{P}_{ij}=\delta_{ij}-\frac{k_ik_j}{k^2}}

      IF ( k_2 ( i_x, i_y, i_z ) .GT. tol ) THEN
        ! Checking IF k^2 is not too low, to cause NaN (this will happen only for (0,0,0))
        proj_xx( i_x, i_y, i_z ) = one - ( kx * kx ) / k_mod_2
        proj_yy( i_x, i_y, i_z ) = one - ( ky * ky ) / k_mod_2
        proj_zz( i_x, i_y, i_z ) = one - ( kz * kz ) / k_mod_2
        proj_xy( i_x, i_y, i_z ) = - ( kx * ky ) / k_mod_2
        proj_yz( i_x, i_y, i_z ) = - ( ky * kz ) / k_mod_2
        proj_zx( i_x, i_y, i_z ) = - ( kz * kx ) / k_mod_2
      ELSE
        proj_xx( i_x, i_y, i_z ) = + twothird
        proj_yy( i_x, i_y, i_z ) = + twothird
        proj_zz( i_x, i_y, i_z ) = + twothird
        proj_xy( i_x, i_y, i_z ) = - onethird
        proj_yz( i_x, i_y, i_z ) = - onethird
        proj_zx( i_x, i_y, i_z ) = - onethird
      END IF

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  S  H  E  L  L     N  O      M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      diff_ceiling                                   =  DBLE( CEILING( k_mod ) ) - k_mod
      ! figuring the decimal part of the |k|, so that it will be alloted to shell k or k+1

      IF (diff_ceiling .GE. hf) THEN
        shell_no( i_x, i_y, i_z )                    =  FLOOR( k_mod )
      ELSE
        shell_no( i_x, i_y, i_z )                    =  CEILING( k_mod )
      END IF

      count_modes_shell( shell_no( i_x, i_y, i_z ) ) =  count_modes_shell( shell_no( i_x, i_y, i_z ) ) + 1
      ! counts no of grid points that belong to a particular shell, it should go as ~s^2

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  T  R  U  N  C  A  T  I  O  N  ,   S  H  E  L  L           M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Truncation mask matrix (multiply this with any spectral matrix to DO the truncation)
      IF ( k_mod_2 .LT. k_G_2 ) THEN

        truncator( i_x, i_y, i_z ) = one
        ! Spherical truncation filter matrix

        tot_active_modes           = tot_active_modes + 1
        ! Total no of active modes inside the truncation sphere.

      ELSE

        truncator( i_x, i_y, i_z )                     =  zero
        ! Outised the truncation sphere.

      END IF

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  H E L I C A L   B A S I S    M A T R I X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      ! IF ( k_per_2 .GT. tol ) THEN
      !
      !   h_real_mod               = root_2 * DSQRT( k_per_2 )
      !   h_imag_mod               = root_2 * DSQRT( k_per_2 * k_mod_2 )
      !   h_pos_x( i_x, i_y, i_z ) = ( - ky / h_real_mod ) + i * ( - kx * kz / h_imag_mod )
      !   h_pos_y( i_x, i_y, i_z ) = ( + kx / h_real_mod ) + i * ( - ky * kz / h_imag_mod )
      !   h_pos_z( i_x, i_y, i_z ) =                       + i * ( + k_per_2 / h_imag_mod )
      !   h_neg_x( i_x, i_y, i_z ) = ( - ky / h_real_mod ) - i * ( - kx * kz / h_imag_mod )
      !   h_neg_y( i_x, i_y, i_z ) = ( + kx / h_real_mod ) - i * ( - ky * kz / h_imag_mod )
      !   h_neg_z( i_x, i_y, i_z ) =                       - i * ( + k_per_2 / h_imag_mod )
      !
      ! ELSE
      !
      !   h_pos_x( i_x, i_y, i_z ) = +     one / root_2
      !   h_pos_y( i_x, i_y, i_z ) = + i * one / root_2
      !   h_pos_z( i_x, i_y, i_z ) =       zero
      !   h_neg_x( i_x, i_y, i_z ) = +     one / root_2
      !   h_neg_y( i_x, i_y, i_z ) = - i * one / root_2
      !   h_neg_z( i_x, i_y, i_z ) =       zero
      !
      ! END IF

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  D I S S I P A T I O N    I N T . F A C T O R    M A T R I X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      diss_Ifactor( i_x, i_y, i_z ) = DEXP( - viscosity * k_mod_2 * dt )

    END DO
    END DO
    END DO

    k_2( 0, 0, 0 ) = one
    ! Just to make sure , when something is divided by k_2 , to avoid NaN. Numerator would be zero anyways,
    ! in most of such cases. So nothing wrong here.

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
    ALLOCATE(k_2(0         :Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(truncator(0   :Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(diss_Ifactor(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(shell_no(0    :Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(k_x(0    :Nh,-Nh:Nh-1,-Nh:Nh-1),k_y(0    :Nh,-Nh:Nh-1,-Nh:Nh-1),k_z(0    :Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(proj_xx(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_yy(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_zz(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(proj_xy(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_yz(0:Nh,-Nh:Nh-1,-Nh:Nh-1),proj_zx(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(F_k_x(0  :Nh,-Nh:Nh-1,-Nh:Nh-1),F_k_y(0  :Nh,-Nh:Nh-1,-Nh:Nh-1),F_k_z(0  :Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(count_modes_shell(0:k_max))
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
    ALLOCATE(v_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

    ALLOCATE(spectral_energy(0       :k_max))
    ALLOCATE(spectral_enstrophy(0    :k_max))
    ALLOCATE(spectral_enstrophy_avg(0:k_max))
    ALLOCATE(spectral_helicity(0     :k_max))
    ALLOCATE(spectral_helicity_avg(0 :k_max))
    ALLOCATE(spectral_energy_avg(0   :k_max))

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
		DEALLOCATE(u_x,u_y,u_z)
		DEALLOCATE(spectral_energy)
    DEALLOCATE(spectral_energy_avg)
		DEALLOCATE(spectral_enstrophy)
		DEALLOCATE(spectral_enstrophy_avg)
		DEALLOCATE(spectral_helicity_avg)
		DEALLOCATE(spectral_helicity)

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
    DEALLOCATE(diss_Ifactor)
    DEALLOCATE(proj_xx,proj_yy,proj_zz)
    DEALLOCATE(proj_xy,proj_yz,proj_zx)
    DEALLOCATE(F_k_x,F_k_y,F_k_z)
    DEALLOCATE(shell_no,count_modes_shell)
    ! DEALLOCATE(h_pos_x,h_pos_y,h_pos_z)
    ! DEALLOCATE(h_neg_x,h_neg_y,h_neg_z)

	END

END MODULE system_basicvariables
