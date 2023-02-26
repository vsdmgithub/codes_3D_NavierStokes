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
! LAST MODIFIED: 20 FEBRAURY 2023
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
    CHARACTER( LEN = 60 )::inp_file

    inp_file  = 'parameters.dat'
    ! This file contains all major Input parameters to be fed from outside file

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! R  E  A  D  I  N  G       I  N  P  U  T       F  I  L  E
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    OPEN( UNIT = 1001, FILE = TRIM( ADJUSTL(inp_file) ) )

    READ( 1001, f_d8p4, ADVANCE ='yes')
    READ( 1001, f_i6,   ADVANCE ='yes')  N
  	! No of collocation points in physical space in one Dimension

    READ( 1001, f_d8p4,  ADVANCE ='yes')
    READ( 1001, f_d8p4,  ADVANCE ='yes')  t_tot
    ! Total time to simulate

    READ( 1001, f_d8p4, ADVANCE ='yes')
    READ( 1001, f_i4,   ADVANCE ='yes')  num_sav
    ! No of saves

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

    IMPLICIT NONE

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! S P A C E A N D T I M E
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    l_sys          = two_pi
    ! Length of periodic cube

    Nh             = INT( N / 2 )
    ! Maximum wavemuber in the First Brillouin zone.

    l_grd          = l_sys / DBLE( N )
    ! Grid distance

    N3             = DBLE( N * N * N )
    ! No of points in real space

    k_tru          = FLOOR( DBLE( N ) / thr ) - 1
    k_tru_sqr      = DBLE( k_tru * k_tru ) + tol_double
    ! Truncation wavenumber shell radius.
    ! For k_i>k_tru, modes are truncated to remove dealiasing error.

    k_max          = CEILING( DSQRT( thr ) * DBLE( Nh ) )
    ! Maximum shell no a mode in first Brilloun zone can get.

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! F L U I D D A T A
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    nrm_fac        = one
    ! Normalization factor for energy - later changed so that initial energy is obtained.

    eng_0          = one
    eng            = eng_0
    eng_pre        = eng
    ! Initial energy

    ! XXXXXXXXXXXXXXXXXXXXXX
    ! VISCOSITY_SELECTION:
    ! XXXXXXXXXXXXXXXXXXXXXX
    ! vis          = 4.0D0 * 1.0E-3 ! For N=256 ! THIS IS THE REFERENCE VISCOSITY
    vis            = ( 128.0D0 / DBLE( N ) ) * vis * 2.0E-3
    ! Viscosity of the system

    k_int          = 3
    ! Integral scale wavnumber - ASSUMPTION

    CALL find_diss_rate_ref( eng_0, k_int, dis_ref )
    ! REF-> <<< system_auxilaries >>>
    ! dis_ref = 0.50D0 ! CUSTOM For N=256

    dis            = dis_ref
    ! Assuming viscous dissipation as reference for start

    cfl_min        = 10
    ! - Courant-Friedrichs-Lewy (CFL) condition - CFL
    ! No of steps (minimum) that should take to cross a grid

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! A U X I L A R Y
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    num_deb        = 10
    ! No of times that the program looks for any 'NaN' while marching forward in time.

    sim_status     = 0
    ! Meaning it is initializing , will be set to '1' at the end

		tur_status     = 0
		! Preparing the state of turbulent flow with stationarity

		lyp_status     = 0
		! Once the perturbation is created, then set to '1'

    frc_status     = 1
    ! Tick that controls the forcing, 1 - YES, 0 - NO

    pre_status     = 0
    ! Status that is checked before the time evolution begins

    num_nan        = 0
    ! No of NaN in V_x (if so)

    inc_err      = 0
    ! Incompressibility error ( Sum(|k.v(k)|) exceeding tolerance )

    nan_err      = 0
    ! NaN found in the data - error

    vis_err      = 0
    ! Viscosity too small for this resolution

    deb_err      = 0
    ! some error found during debug (could be either nan or incomp)

    WRITE( res_char, f_i8 ) N
    ! Converting resolution value to character

	  CALL compute_system_details
		! REF-> <<< system_basicdeclaration >>>

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

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! KOLMOGOROV SCALES -
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		l_kol      = ( ( vis ** thr) / dis ) ** qtr
    ! Kolmogorov length scale

    u_kol      = ( dis * vis ) ** qtr
    ! Kolmogorov velocity

    t_kol      = DSQRT( vis / dis )
    ! Kolmogorov timescale

		k_kol      = FLOOR( l_sys / l_kol )
		! Kolmogorov wavenumber

    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    ! INTEGRAL SCALE
    ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		l_int      = ( eng ** 1.5D0 ) / dis
		! Integral length scale

    u_int      = DSQRT( two / thr ) * DSQRT( eng )
    ! 1D RMS Velocity

    t_int      = DSQRT( two / thr ) * eng / dis
    ! Turbulent time scale ( Large eddy turn-over timescale )

		k_int      = CEILING( l_sys / l_int )
		! Integral scale wavenumber

    res_pow    = DBLE( k_tru ) * l_kol
    ! Extent beyond which the smallest scale is resolved - atleast 1

    t_grd      = l_grd / u_int
    ! Time for particle to cross a grid

		l_tay      = DSQRT( two * fiv * vis * eng / dis )
		! Taylor scale
		k_tay      = CEILING( l_sys / l_tay )
		! Integral scale wavenumber

    rey_int    = FLOOR( DSQRT( eng ) * l_int / vis )
    rey_tay    = FLOOR( DSQRT( two / thr ) * DSQRT( eng ) * l_tay / vis )
    ! Reynold's number and Taylor scale Reynold's number.

    FIND_MIN_TIME: IF ( t_grd .LT. t_kol ) THEN
      t_min    = t_grd
    ELSE
      t_min    = t_kol
    END IF FIND_MIN_TIME

    dt_max     = t_min / DBLE( cfl_min )
    ! Maximum value of time step

    CALL find_timestep( t_min, dt_max, dt )
    ! Finds a time smaller than 'dt_max' in terms of '0.0..0p' , where p being '1' or '5'
    ! REF-> <<< system_auxilaries >>>

    cfl        = FLOOR( t_min / dt )
    ! CFL number of the system

    CALL time_to_step_convert(t_tot,t_step_tot,dt)
    ! returns the no of time_steps (\delta t) in a given time
    ! REF-> <<< system_auxilaries >>>

    t_step_sav = t_step_tot / num_sav
    ! Determines how many time steps after the save has to be made.

    t_step_deb = t_step_tot / num_deb
    ! Determines how many time steps after the checking has to be made

    CALL step_to_time_convert(t_step_sav,t_sav,dt)
    ! REF-> <<< system_auxilaries >>>

    CALL time_to_step_convert( t_kol,t_step_kol,dt)
    ! REF-> <<< system_auxilaries >>>

	END
  SUBROUTINE init_global_arrays
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !   This defines few constant arrays that do not evolve. k^2 matrix,k matrix, projection matrix,
  ! shell no matrix, Den_k matrix.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::dif_ceil
    DOUBLE PRECISION::kx,ky,kz
    DOUBLE PRECISION::k_mod,k_mod_sqr
    INTEGER(KIND=4)::ct


    CALL allocate_operators
    ! Allocates the arrays declared here.

    Shell       = 0
    Den_k       = 0
    I_fac       = 0
    num_mod_frc = 0
    num_mod     = 0

    k_frc       = 3
    ! Wavenumbers below which forcing is implemented

    ! +++++++++++++++++++++++++++++++++
    !  A  X  I  S      G   R  I   D   S
    !  +++++++++++++++++++++++++++++++++
    LOOP_RX_000: DO i_x = 0, N-1
      Axis(i_x) = DBLE( i_x ) * l_grd
      ! Location of grid points along principal Axis
    END DO LOOP_RX_000

    LOOP_FX_101: DO j_x = 0, Nh
    LOOP_FY_101: DO j_y = -Nh, Nh - 1
    LOOP_FZ_101: DO j_z = -Nh, Nh - 1

      kx                    = DBLE( j_x )
      ky                    = DBLE( j_y )
      kz                    = DBLE( j_z )
      K_x( j_x, j_y, j_z )  = kx
      K_y( j_x, j_y, j_z )  = ky
      K_z( j_x, j_y, j_z )  = kz
      ! Just the k component matrix storing its grid points.

      k_mod_sqr              = kx**two + ky**two + kz**two
      Lapla( j_x, j_y, j_z ) = k_mod_sqr
      k_mod                  = DSQRT( k_mod_sqr )

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  P  R  O  J  E  C  T  I  O  N             M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Projection matrix \mathbb{P}_{ij}=\delta_{ij}-\frac{k_ik_j}{k^2}}

      PROJECTION_TENSOR: IF ( Lapla ( j_x, j_y, j_z ) .GT. tol ) THEN
        ! Checking IF k^2 is not too low, to cause NaN (this will happen only for (0,0,0))
        Proj_xx( j_x, j_y, j_z ) = one - ( kx * kx ) / k_mod_sqr
        Proj_yy( j_x, j_y, j_z ) = one - ( ky * ky ) / k_mod_sqr
        Proj_zz( j_x, j_y, j_z ) = one - ( kz * kz ) / k_mod_sqr
        Proj_xy( j_x, j_y, j_z ) = - ( kx * ky ) / k_mod_sqr
        Proj_yz( j_x, j_y, j_z ) = - ( ky * kz ) / k_mod_sqr
        Proj_zx( j_x, j_y, j_z ) = - ( kz * kx ) / k_mod_sqr
      ELSE
        Proj_xx( j_x, j_y, j_z ) = + twothird
        Proj_yy( j_x, j_y, j_z ) = + twothird
        Proj_zz( j_x, j_y, j_z ) = + twothird
        Proj_xy( j_x, j_y, j_z ) = - onethird
        Proj_yz( j_x, j_y, j_z ) = - onethird
        Proj_zx( j_x, j_y, j_z ) = - onethird
      END IF PROJECTION_TENSOR

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  S  H  E  L  L     N  O      M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      dif_ceil                        = DBLE( CEILING( k_mod ) ) - k_mod
      ! figuring the decimal part of the |k|, so that it will be alloted to shell k or k+1

      SHELL_ALLOCATION: IF (dif_ceil .GE. hf) THEN
        Shell( j_x, j_y, j_z )        = FLOOR( k_mod )
      ELSE
        Shell( j_x, j_y, j_z )        = CEILING( k_mod )
      END IF SHELL_ALLOCATION

      Den_k( Shell( j_x, j_y, j_z ) ) = Den_k( Shell( j_x, j_y, j_z ) ) + 1
      ! counts no of grid points that belong to a particular shell, it should go as ~4 pi s^2

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  T  R  U  N  C  A  T  I  O  N  ,   S  H  E  L  L           M  A  T  R  I  X
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! Truncation mask matrix (multiply this with any spectral matrix to DO the truncation)
      TRUNCATION_MASK: IF ( k_mod_sqr .LT. k_tru_sqr ) THEN

        Trunc( j_x, j_y, j_z ) = one
        ! Spherical truncation filter matrix

        num_mod                = num_mod + 1
        ! Total no of active modes inside the truncation sphere.

      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ! D I S S I P A T I O N I N T . F A C T O R M A T R I X
      ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        I_fac( j_x, j_y, j_z ) = DEXP( - vis * k_mod_sqr * dt )

      ELSE

        Trunc( j_x, j_y, j_z ) = zero
        ! Outised the truncation sphere.

      END IF TRUNCATION_MASK

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  F O R C I N G     S H E L L    C O U N T
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      FORCING_CHECK_101: IF ( frc_status .EQ. 1 ) THEN
      FORCING_SHELL_COUNT_101:  IF ( Shell( j_x, j_y, j_z ) .EQ. k_frc )  THEN
        num_mod_frc = num_mod_frc + 1
      END IF FORCING_SHELL_COUNT_101
	    END IF FORCING_CHECK_101

    END DO LOOP_FZ_101
    END DO LOOP_FY_101
    END DO LOOP_FX_101

    Lapla( 0, 0, 0 ) = one
    ! Just to make sure , when something is divided by Lapla , to avoid NaN. Numerator would be zero anyways, ! in most of such cases. So nothing wrong here.

    FORCING_CHECK_201: IF ( frc_status .EQ. 1 ) THEN

      ALLOCATE( F_kx_ind( num_mod_frc ) )
      ALLOCATE( F_ky_ind( num_mod_frc ) )
      ALLOCATE( F_kz_ind( num_mod_frc ) )

      ct = 0
      LOOP_FX_201: DO j_x = 0, Nh
      LOOP_FY_201: DO j_y = -Nh, Nh - 1
      LOOP_FZ_201: DO j_z = -Nh, Nh - 1

        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  F O R C I N G     S H E L L    D A T A
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        FORCING_SHELL_COUNT_201:  IF ( Shell( j_x, j_y, j_z ) .EQ. k_frc )  THEN
          ct        = ct + 1
          F_kx_ind( ct ) = j_x
          F_ky_ind( ct ) = j_y
          F_kz_ind( ct ) = j_z
        END IF FORCING_SHELL_COUNT_201

      END DO LOOP_FZ_201
      END DO LOOP_FY_201
      END DO LOOP_FX_201

    END IF FORCING_CHECK_201

  END

  SUBROUTINE allocate_operators
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate arrays which are constants basically, k2,Trunc, shell no etc.,
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  R  R  A  Y         A  L  L  O  C  A  T  I  O  N .
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(Axis(0:N-1))
    ALLOCATE(Lapla(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(Trunc(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(Shell(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(I_fac(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(K_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),K_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),K_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(Proj_xx(0:Nh,-Nh:Nh-1,-Nh:Nh-1),Proj_yy(0:Nh,-Nh:Nh-1,-Nh:Nh-1),Proj_zz(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(Proj_xy(0:Nh,-Nh:Nh-1,-Nh:Nh-1),Proj_yz(0:Nh,-Nh:Nh-1,-Nh:Nh-1),Proj_zx(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(Den_k(0:k_max))
		IF ( frc_status .EQ. 1 ) THEN
	    ALLOCATE(F_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),F_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),F_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
			F_x = c0
			F_y = c0
			F_z = c0
		END IF

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
    ALLOCATE(U_x(0:N-1,0:N-1,0:N-1),U_y(0:N-1,0:N-1,0:N-1),U_z(0:N-1,0:N-1,0:N-1))
    ALLOCATE(V_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),V_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),V_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(Eng_k(0       :k_max))
    ALLOCATE(Eng_k_avg(0   :k_max))

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
    ALLOCATE(W_x(0:N-1,0:N-1,0:N-1),W_y(0:N-1,0:N-1,0:N-1),W_z(0:N-1,0:N-1,0:N-1))
    ALLOCATE(W_kx(0:Nh,-Nh:Nh-1,-Nh:Nh-1),W_ky(0:Nh,-Nh:Nh-1,-Nh:Nh-1),W_kz(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

  END

  SUBROUTINE allocate_chaos_arrays
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate arrays related to chaos measurement
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		ALLOCATE(Ub_x(0:N-1,0:N-1,0:N-1),Ub_y(0:N-1,0:N-1,0:N-1),Ub_z(0:N-1,0:N-1,0:N-1))
    ALLOCATE(Vb_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),Vb_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),Vb_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

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
		DEALLOCATE(V_x,V_y,V_z)
		DEALLOCATE(U_x,U_y,U_z)
		DEALLOCATE(Eng_k)
    DEALLOCATE(Eng_k_avg)

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
		DEALLOCATE(W_x,W_y,W_z)
		DEALLOCATE(W_kx,W_ky,W_kz)

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
    DEALLOCATE(Axis)
    DEALLOCATE(Lapla,K_x,K_y,K_z)
    DEALLOCATE(Trunc)
    DEALLOCATE(I_fac)
    DEALLOCATE(Proj_xx,Proj_yy,Proj_zz)
    DEALLOCATE(Proj_xy,Proj_yz,Proj_zx)
    DEALLOCATE(Shell,Den_k)
    FORCING_CHECK_202: IF ( frc_status .EQ. 1 ) THEN
      DEALLOCATE(F_kx_ind,F_ky_ind,F_kz_ind)
	    DEALLOCATE(F_x,F_y,F_z)
    END IF FORCING_CHECK_202

	END

	SUBROUTINE deallocate_chaos_arrays
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(Vb_x,Vb_y,Vb_z)
		DEALLOCATE(Ub_x,Ub_y,Ub_z)
	END

END MODULE system_basicdeclaration
