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
! MODULE: system_solver
! LAST MODIFIED: 22 JULY 2021
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER MODULE TO SOLVE 3D NAVIER STOKES EQUATIONS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_solver
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Takes the spectral velocity and moves it one step forward in Navier-Stokes equation using the given algorithm, uses FFTW.
! This has RK4, AB4 algorithm for the Non-linear term
! Exact integration factor for the viscous term
! The spectral equation is
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! dv_i(k)/dt = P_ij . ( u(x) X w(x) )_j(k) - \nu k^2 v_i(k)
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
	!  SUB-MODULES
	!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE system_basicvariables
	USE system_fftw

	IMPLICIT NONE
	! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::uXw_x,uXw_y,uXw_z
	! Real cross product between velocity and vorticity

	! _________________________________________
  ! FOURIER SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::vXw_x, vXw_y, vXw_z
  ! Spectral cross product between velocity and vorticity

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_dot_m1,v_y_dot_m1,v_z_dot_m1
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_dot_m2,v_y_dot_m2,v_z_dot_m2
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_dot_m3,v_y_dot_m3,v_z_dot_m3
	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_dot,v_y_dot,v_z_dot
	! Spectral derivatives for AB4

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_pred,v_y_pred,v_z_pred
  ! Spectral velocity predictor for AB4 matrix

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::dv1_x,dv2_x,dv3_x,dv4_x,dv1_y,dv2_y
  DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::dv3_y,dv4_y,dv1_z,dv2_z,dv3_z,dv4_z
  ! Intermediate matrices for RK4 algorithm

	DOUBLE COMPLEX,DIMENSION(:,:,:),ALLOCATABLE::v_x_temp,v_y_temp,v_z_temp
  ! temporary matrices to store velocities during RK4 algorithm


	CONTAINS

	SUBROUTINE allocate_solver
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays for solver - Commom for both AB4, RK4
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		ALLOCATE(vXw_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),vXw_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),vXw_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
		ALLOCATE(uXw_x(0:N-1,0:N-1,0:N-1),uXw_y(0:N-1,0:N-1,0:N-1),uXw_z(0:N-1,0:N-1,0:N-1))

	END

	SUBROUTINE allocate_solver_RK4
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays for RK4 algoritm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ALLOCATE(dv1_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(dv3_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_x(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(dv1_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(dv3_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_y(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(dv1_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv2_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(dv3_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1),dv4_z(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
    ALLOCATE(v_x_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_temp(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

	END

	SUBROUTINE solver_RK4_algorithm
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
	! Alg: - Runga kutta 4th order
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'v(k)'
		v_x_temp = v_x
		v_y_temp = v_y
		v_z_temp = v_z
		CALL time_increment_RK1() ! This call provides \vec{dv} for the existing \vec{v}
		v_x      = v_x_temp + hf * dv1_x
		v_y      = v_y_temp + hf * dv1_y
		v_z      = v_z_temp + hf * dv1_z
		CALL time_increment_RK2()
		v_x      = v_x_temp + hf * dv2_x
		v_y      = v_y_temp + hf * dv2_y
		v_z      = v_z_temp + hf * dv2_z
		CALL time_increment_RK3()
		v_x      = v_x_temp + dv3_x
		v_y      = v_y_temp + dv3_y
		v_z      = v_z_temp + dv3_z
		CALL time_increment_RK4()
		! Final increment for 'v(k)'

		IF ( forcing_status .EQ. 1 ) THEN
			v_x      = ( v_x_temp + ( dv1_x + two * dv2_x + two * dv3_x + dv4_x ) / six ) * diss_Ifactor + F_k_x
			v_y      = ( v_y_temp + ( dv1_y + two * dv2_y + two * dv3_y + dv4_y ) / six ) * diss_Ifactor + F_k_y
			v_z      = ( v_z_temp + ( dv1_z + two * dv2_z + two * dv3_z + dv4_z ) / six ) * diss_Ifactor + F_k_z
		ELSE
			v_x      = ( v_x_temp + ( dv1_x + two * dv2_x + two * dv3_x + dv4_x ) / six ) * diss_Ifactor
			v_y      = ( v_y_temp + ( dv1_y + two * dv2_y + two * dv3_y + dv4_y ) / six ) * diss_Ifactor
			v_z      = ( v_z_temp + ( dv1_z + two * dv2_z + two * dv3_z + dv4_z ) / six ) * diss_Ifactor
		END IF

		CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
		! Real Velocity

	END

	SUBROUTINE allocate_solver_AB4
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to allocate arrays for AB4 algorithm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  A  L  L  O  C  A  T  I  O  N     F  O  R     S  O  L  V  E  R
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		ALLOCATE(v_x_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_dot(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
		ALLOCATE(v_x_dot_m1(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_x_dot_m2(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_x_dot_m3(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
		ALLOCATE(v_y_dot_m1(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_dot_m2(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_dot_m3(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
		ALLOCATE(v_z_dot_m1(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_dot_m2(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_dot_m3(0:Nh,-Nh:Nh-1,-Nh:Nh-1))
		ALLOCATE(v_x_pred(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_y_pred(0:Nh,-Nh:Nh-1,-Nh:Nh-1),v_z_pred(0:Nh,-Nh:Nh-1,-Nh:Nh-1))

	END

	SUBROUTINE solver_AB4_algorithm
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE AB4 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
	! Alg: - Adam Bashforth 4th Order algorithm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		IF ( t_step .EQ. 0 ) THEN

			! Using initial condition to get -3 step value
			CALL time_derivative_AB()
			v_x_dot_m3 = v_x_dot
			v_y_dot_m3 = v_y_dot
			v_z_dot_m3 = v_z_dot

			CALL allocate_solver_RK4
			! Allocating RK4 for first 3 steps

			CALL solver_RK4_algorithm

		ELSE IF ( t_step .EQ. 1 ) THEN

			! Filling up for -2 step value
			CALL time_derivative_AB()
			v_x_dot_m2 = v_x_dot
			v_y_dot_m2 = v_y_dot
			v_z_dot_m2 = v_z_dot

			CALL solver_RK4_algorithm

		ELSE IF ( t_step .EQ. 2 ) THEN

			! Filling up for -1 step value
			CALL time_derivative_AB()
			v_x_dot_m1 = v_x_dot
			v_y_dot_m1 = v_y_dot
			v_z_dot_m1 = v_z_dot

			CALL solver_RK4_algorithm

			CALL deallocate_solver_RK4
			! No need for RK4 anymore

		ELSE

			CALL time_derivative_AB()

			v_x_pred = ( v_x + dt * ( - 9.0D0 * v_x_dot_m3 + 37.0D0 * v_x_dot_m2 - 59.0D0 * v_x_dot_m1 + &
			 													 55.0D0 * v_x_dot ) / 24.0D0 ) * diss_Ifactor
			v_y_pred = ( v_y + dt * ( - 9.0D0 * v_y_dot_m3 + 37.0D0 * v_y_dot_m2 - 59.0D0 * v_y_dot_m1 + &
			 													 55.0D0 * v_y_dot ) / 24.0D0 ) * diss_Ifactor
			v_z_pred = ( v_z + dt * ( - 9.0D0 * v_z_dot_m3 + 37.0D0 * v_z_dot_m2 - 59.0D0 * v_z_dot_m1 + &
			 													 55.0D0 * v_z_dot ) / 24.0D0 ) * diss_Ifactor

			CALL time_derivative_AB_pred()

			v_x      = ( v_x + dt * ( v_x_dot_m2 - 5.0D0 * v_x_dot_m1 + 19.0D0 * v_x_dot + &
			 								  9.0D0 * v_x_dot_m3 ) / 24.0D0 ) * diss_Ifactor
			v_y      = ( v_y + dt * ( v_y_dot_m2 - 5.0D0 * v_y_dot_m1 + 19.0D0 * v_y_dot + &
			 								  9.0D0 * v_y_dot_m3 ) / 24.0D0 ) * diss_Ifactor
			v_z      = ( v_z + dt * ( v_z_dot_m2 - 5.0D0 * v_z_dot_m1 + 19.0D0 * v_z_dot + &
			 								  9.0D0 * v_z_dot_m3 ) / 24.0D0 ) * diss_Ifactor
			! Predicted 'v_dot' is stored in 'v_dot_m3' - to save space :)

			! Shifting the known velocities for next step
			v_x_dot_m3 = v_x_dot_m2
			v_x_dot_m2 = v_x_dot_m1
			v_x_dot_m1 = v_x_dot

			v_y_dot_m3 = v_y_dot_m2
			v_y_dot_m2 = v_y_dot_m1
			v_y_dot_m1 = v_y_dot

			v_z_dot_m3 = v_z_dot_m2
			v_z_dot_m2 = v_z_dot_m1
			v_z_dot_m1 = v_z_dot

			CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
			! Real Velocity

		END IF

	END

	SUBROUTINE time_derivative_AB()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the Navier-Stokes EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   C  R  O  S  S    P  R  O  D  U  C  T    ( U  X  W  )
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		w_vx = i * ( k_y * v_z - k_z * v_y )
    w_vy = i * ( k_z * v_x - k_x * v_z )
    w_vz = i * ( k_x * v_y - k_y * v_x )
    ! Spectral Vorticity

    CALL fft_c2r( w_vx, w_vy, w_vz, N, Nh, w_ux, w_uy, w_uz )
    ! Real Vorticity

		CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )

		uXw_x = u_y * w_uz - u_z * w_uy
		uXw_y = u_z * w_ux - u_x * w_uz
		uXw_z = u_x * w_uy - u_y * w_ux
		! Cross product between velocity and vorticity

		CALL fft_r2c( uXw_x, uXw_y, uXw_z,  N, Nh, vXw_x, vXw_y, vXw_z )
		! FFT to get spectral cross product

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.

		v_x_dot = truncator * ( proj_xx * vXw_x + proj_xy * vXw_y + proj_zx * vXw_z ) + F_k_x
		v_y_dot = truncator * ( proj_xy * vXw_x + proj_yy * vXw_y + proj_yz * vXw_z ) + F_k_y
		v_z_dot = truncator * ( proj_zx * vXw_x + proj_yz * vXw_y + proj_zz * vXw_z ) + F_k_z

	END

	SUBROUTINE time_derivative_AB_pred()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the Navier-Stokes EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   C  R  O  S  S    P  R  O  D  U  C  T    ( U  X  W  )
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

 		w_vx = i * ( k_y * v_z_pred - k_z * v_y_pred )
    w_vy = i * ( k_z * v_x_pred - k_x * v_z_pred )
    w_vz = i * ( k_x * v_y_pred - k_y * v_x_pred )
    ! Spectral Vorticity

    CALL fft_c2r( w_vx, w_vy, w_vz, N, Nh, w_ux, w_uy, w_uz )
    ! Real Vorticity

		CALL fft_c2r( v_x_pred, v_y_pred, v_z_pred, N, Nh, u_x, u_y, u_z )

		uXw_x = u_y * w_uz - u_z * w_uy
		uXw_y = u_z * w_ux - u_x * w_uz
		uXw_z = u_x * w_uy - u_y * w_ux
		! Cross product between velocity and vorticity

		CALL fft_r2c( uXw_x, uXw_y, uXw_z,  N, Nh, vXw_x, vXw_y, vXw_z )
		! FFT to get spectral cross product

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.

		v_x_dot_m3 = truncator * ( proj_xx * vXw_x + proj_xy * vXw_y + proj_zx * vXw_z ) + F_k_x
		v_y_dot_m3 = truncator * ( proj_xy * vXw_x + proj_yy * vXw_y + proj_yz * vXw_z ) + F_k_y
		v_z_dot_m3 = truncator * ( proj_zx * vXw_x + proj_yz * vXw_y + proj_zz * vXw_z ) + F_k_z

	END

	SUBROUTINE time_increment_RK1()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the Navier-Stokes EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE

		CALL vorticity_X_velocity
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the crossproduce (v x w)_k in spectral space and project it with projection matrix and finally truncate it.
		dv1_x = dt * ( truncator * ( proj_xx * vXw_x + proj_xy * vXw_y + proj_zx * vXw_z ) )
		dv1_y = dt * ( truncator * ( proj_xy * vXw_x + proj_yy * vXw_y + proj_yz * vXw_z ) )
		dv1_z = dt * ( truncator * ( proj_zx * vXw_x + proj_yz * vXw_y + proj_zz * vXw_z ) )

	END

	SUBROUTINE time_increment_RK2()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the Navier-Stokes EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		CALL vorticity_X_velocity
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the crossproduce (v x w)_k in spectral space and project it with projection matrix and finally truncate it.
		dv2_x = dt * ( truncator * ( proj_xx * vXw_x + proj_xy * vXw_y + proj_zx * vXw_z ) )
		dv2_y = dt * ( truncator * ( proj_xy * vXw_x + proj_yy * vXw_y + proj_yz * vXw_z ) )
		dv2_z = dt * ( truncator * ( proj_zx * vXw_x + proj_yz * vXw_y + proj_zz * vXw_z ) )

	END

	SUBROUTINE time_increment_RK3()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the Navier-Stokes EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		CALL vorticity_X_velocity
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the crossproduce (v x w)_k in spectral space and project it with projection matrix and finally truncate it.
		dv3_x = dt * ( truncator * ( proj_xx * vXw_x + proj_xy * vXw_y + proj_zx * vXw_z ) )
		dv3_y = dt * ( truncator * ( proj_xy * vXw_x + proj_yy * vXw_y + proj_yz * vXw_z ) )
		dv3_z = dt * ( truncator * ( proj_zx * vXw_x + proj_yz * vXw_y + proj_zz * vXw_z ) )

	END

	SUBROUTINE time_increment_RK4()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the Navier-Stokes EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE

		CALL vorticity_X_velocity
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the crossproduce (v x w)_k in spectral space and project it with projection matrix and finally truncate it.
		dv4_x = dt * ( truncator * ( proj_xx * vXw_x + proj_xy * vXw_y + proj_zx * vXw_z ) )
		dv4_y = dt * ( truncator * ( proj_xy * vXw_x + proj_yy * vXw_y + proj_yz * vXw_z ) )
		dv4_z = dt * ( truncator * ( proj_zx * vXw_x + proj_yz * vXw_y + proj_zz * vXw_z ) )

	END

	SUBROUTINE vorticity_X_velocity
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to give spectral convection term using v_k.
	! 1. First i(k x v)    --> w  iFFT is done
	! 2. Next (u x w) (x)  --> (u x w)_k  FFT is done
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   C  R  O  S  S    P  R  O  D  U  C  T    ( U  X  W  )
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

		w_vx = i * ( k_y * v_z - k_z * v_y )
    w_vy = i * ( k_z * v_x - k_x * v_z )
    w_vz = i * ( k_x * v_y - k_y * v_x )
    ! Spectral Vorticity

		CALL fft_c2r( w_vx, w_vy, w_vz, N, Nh, w_ux, w_uy, w_uz )
    ! Real Vorticity

		CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
		! Real Velocity

		uXw_x = u_y * w_uz - u_z * w_uy
		uXw_y = u_z * w_ux - u_x * w_uz
		uXw_z = u_x * w_uy - u_y * w_ux
		! Cross product between velocity and vorticity

		CALL fft_r2c( uXw_x, uXw_y, uXw_z,  N, Nh, vXw_x, vXw_y, vXw_z )
		! FFT to get spectral cross product

  END

	SUBROUTINE deallocate_solver_RK4
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		DEALLOCATE(dv1_x,dv2_x)
		DEALLOCATE(dv3_x,dv4_x)
		DEALLOCATE(dv1_y,dv2_y)
		DEALLOCATE(dv3_y,dv4_y)
		DEALLOCATE(dv1_z,dv2_z)
		DEALLOCATE(dv3_z,dv4_z)

	END

	SUBROUTINE deallocate_solver_AB4
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(v_x_pred,v_y_pred,v_z_pred)
		DEALLOCATE(v_x_dot_m1,v_x_dot_m2,v_x_dot_m3)
		DEALLOCATE(v_y_dot_m1,v_y_dot_m2,v_y_dot_m3)
		DEALLOCATE(v_z_dot_m1,v_z_dot_m2,v_z_dot_m3)

	END

	SUBROUTINE deallocate_solver
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		DEALLOCATE(vXw_x,vXw_y,vXw_z)
		DEALLOCATE(uXw_x,uXw_y,uXw_z)

	END

 END MODULE system_solver
