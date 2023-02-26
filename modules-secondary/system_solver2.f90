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
! LAST MODIFIED: 20 FEBRAURY 2023
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER MODULE TO SOLVE 3D NAVIER STOKES EQUATIONS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_solver2
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
	USE system_solver

	IMPLICIT NONE

	CONTAINS

	SUBROUTINE solver2_RK4_algorithm
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
	! Alg: - Runga kutta 4th order
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'v(k)'
		V_x_temp = Vb_x
		V_y_temp = Vb_y
		V_z_temp = Vb_z
		CALL time_increment2_RK1() ! This call provides \vec{dv} for the existing \vec{v}
		Vb_x      = V_x_temp + hf * Dv1_x
		Vb_y      = V_y_temp + hf * Dv1_y
		Vb_z      = V_z_temp + hf * Dv1_z
		CALL time_increment2_RK2()
		Vb_x      = V_x_temp + hf * Dv2_x
		Vb_y      = V_y_temp + hf * Dv2_y
		Vb_z      = V_z_temp + hf * Dv2_z
		CALL time_increment2_RK3()
		Vb_x      = V_x_temp + Dv3_x
		Vb_y      = V_y_temp + Dv3_y
		Vb_z      = V_z_temp + Dv3_z
		CALL time_increment2_RK4()
		! Final increment for 'v(k)'

		Vb_x      = ( V_x_temp + ( Dv1_x + two * Dv2_x + two * Dv3_x + Dv4_x ) / six ) * I_fac
		Vb_y      = ( V_y_temp + ( Dv1_y + two * Dv2_y + two * Dv3_y + Dv4_y ) / six ) * I_fac
		Vb_z      = ( V_z_temp + ( Dv1_z + two * Dv2_z + two * Dv3_z + Dv4_z ) / six ) * I_fac

		CALL fft_c2r( Vb_x, Vb_y, Vb_z, N, Nh, Ub_x, Ub_y, Ub_z )
		! Real Velocity

	END

	SUBROUTINE time_increment2_RK1()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the Navier-Stokes EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE

		CALL vorticity_X_velocity2
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the crossproduce (v x w)_k in spectral space and project it with projection matrix and finally truncate it.
		IF ( frc_status .EQ. 1 ) THEN
			Dv1_x = dt * ( Trunc * ( Proj_xx * VcW_x + Proj_xy * VcW_y + Proj_zx * VcW_z ) + F_x )
			Dv1_y = dt * ( Trunc * ( Proj_xy * VcW_x + Proj_yy * VcW_y + Proj_yz * VcW_z ) + F_y )
			Dv1_z = dt * ( Trunc * ( Proj_zx * VcW_x + Proj_yz * VcW_y + Proj_zz * VcW_z ) + F_z )
		ELSE
			Dv1_x = dt * ( Trunc * ( Proj_xx * VcW_x + Proj_xy * VcW_y + Proj_zx * VcW_z ) )
			Dv1_y = dt * ( Trunc * ( Proj_xy * VcW_x + Proj_yy * VcW_y + Proj_yz * VcW_z ) )
			Dv1_z = dt * ( Trunc * ( Proj_zx * VcW_x + Proj_yz * VcW_y + Proj_zz * VcW_z ) )
		END IF

	END

	SUBROUTINE time_increment2_RK2()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the Navier-Stokes EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		CALL vorticity_X_velocity2
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the crossproduce (v x w)_k in spectral space and project it with projection matrix and finally truncate it.
		IF ( frc_status .EQ. 1 ) THEN
			Dv2_x = dt * ( Trunc * ( Proj_xx * VcW_x + Proj_xy * VcW_y + Proj_zx * VcW_z ) + F_x )
			Dv2_y = dt * ( Trunc * ( Proj_xy * VcW_x + Proj_yy * VcW_y + Proj_yz * VcW_z ) + F_y )
			Dv2_z = dt * ( Trunc * ( Proj_zx * VcW_x + Proj_yz * VcW_y + Proj_zz * VcW_z ) + F_z )
		ELSE
			Dv2_x = dt * ( Trunc * ( Proj_xx * VcW_x + Proj_xy * VcW_y + Proj_zx * VcW_z ) )
			Dv2_y = dt * ( Trunc * ( Proj_xy * VcW_x + Proj_yy * VcW_y + Proj_yz * VcW_z ) )
			Dv2_z = dt * ( Trunc * ( Proj_zx * VcW_x + Proj_yz * VcW_y + Proj_zz * VcW_z ) )
		END IF

	END

	SUBROUTINE time_increment2_RK3()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the Navier-Stokes EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE
		CALL vorticity_X_velocity2
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the crossproduce (v x w)_k in spectral space and project it with projection matrix and finally truncate it.
		IF ( frc_status .EQ. 1 ) THEN
			Dv3_x = dt * ( Trunc * ( Proj_xx * VcW_x + Proj_xy * VcW_y + Proj_zx * VcW_z ) + F_x )
			Dv3_y = dt * ( Trunc * ( Proj_xy * VcW_x + Proj_yy * VcW_y + Proj_yz * VcW_z ) + F_y )
			Dv3_z = dt * ( Trunc * ( Proj_zx * VcW_x + Proj_yz * VcW_y + Proj_zz * VcW_z ) + F_z )
		ELSE
			Dv3_x = dt * ( Trunc * ( Proj_xx * VcW_x + Proj_xy * VcW_y + Proj_zx * VcW_z ) )
			Dv3_y = dt * ( Trunc * ( Proj_xy * VcW_x + Proj_yy * VcW_y + Proj_yz * VcW_z ) )
			Dv3_z = dt * ( Trunc * ( Proj_zx * VcW_x + Proj_yz * VcW_y + Proj_zz * VcW_z ) )
		END IF

	END

	SUBROUTINE time_increment2_RK4()
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to get the time derivative matrix for matrix 'v(k)'
	! This is the Navier-Stokes EQUATION implemented for numerical computation
	! spectral space.
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		IMPLICIT NONE

		CALL vorticity_X_velocity2
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the crossproduce (v x w)_k in spectral space and project it with projection matrix and finally truncate it.
		IF ( frc_status .EQ. 1 ) THEN
			Dv4_x = dt * ( Trunc * ( Proj_xx * VcW_x + Proj_xy * VcW_y + Proj_zx * VcW_z ) + F_x )
			Dv4_y = dt * ( Trunc * ( Proj_xy * VcW_x + Proj_yy * VcW_y + Proj_yz * VcW_z ) + F_y )
			Dv4_z = dt * ( Trunc * ( Proj_zx * VcW_x + Proj_yz * VcW_y + Proj_zz * VcW_z ) + F_z )
		ELSE
			Dv4_x = dt * ( Trunc * ( Proj_xx * VcW_x + Proj_xy * VcW_y + Proj_zx * VcW_z ) )
			Dv4_y = dt * ( Trunc * ( Proj_xy * VcW_x + Proj_yy * VcW_y + Proj_yz * VcW_z ) )
			Dv4_z = dt * ( Trunc * ( Proj_zx * VcW_x + Proj_yz * VcW_y + Proj_zz * VcW_z ) )
		END IF

	END

	SUBROUTINE vorticity_X_velocity2
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

		W_kx = i * ( K_y * Vb_z - K_z * Vb_y )
    W_ky = i * ( K_z * Vb_x - K_x * Vb_z )
    W_kz = i * ( K_x * Vb_y - K_y * Vb_x )
    ! Spectral Vorticity

		CALL fft_c2r( W_kx, W_ky, W_kz, N, Nh, W_x, W_y, W_z )
    ! Real Vorticity

		CALL fft_c2r( Vb_x, Vb_y, Vb_z, N, Nh, Ub_x, Ub_y, Ub_z )
		! Real Velocity

		UcW_x = Ub_y * W_z - Ub_z * W_y
		UcW_y = Ub_z * W_x - Ub_x * W_z
		UcW_z = Ub_x * W_y - Ub_y * W_x
		! Cross product between velocity and vorticity

		CALL fft_r2c( UcW_x, UcW_y, UcW_z,  N, Nh, VcW_x, VcW_y, VcW_z )
		! FFT to get spectral cross product

  END

	SUBROUTINE solver2_AB4_algorithm
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to USE AB4 algorithm to move one step forward in time for the matrix 'v(k,t)-> v(k,t+1)'
	! Alg: - Adam Bashforth 4th Order algorithm
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		IF ( t_step .EQ. 0 ) THEN

			! Using initial condition to get -3 step value
			CALL time_derivative2_AB()
			V_x_dot_m3 = V_x_dot
			V_y_dot_m3 = V_y_dot
			V_z_dot_m3 = V_z_dot

			CALL allocate_solver_RK4
			! Allocating RK4 for first 3 steps

			CALL solver2_RK4_algorithm

		ELSE IF ( t_step .EQ. 1 ) THEN

			! Filling up for -2 step value
			CALL time_derivative2_AB()
			V_x_dot_m2 = V_x_dot
			V_y_dot_m2 = V_y_dot
			V_z_dot_m2 = V_z_dot

			CALL solver2_RK4_algorithm

		ELSE IF ( t_step .EQ. 2 ) THEN

			! Filling up for -1 step value
			CALL time_derivative2_AB()
			V_x_dot_m1 = V_x_dot
			V_y_dot_m1 = V_y_dot
			V_z_dot_m1 = V_z_dot

			CALL solver2_RK4_algorithm

			CALL deallocate_solver_RK4
			! No need for RK4 anymore

		ELSE

			CALL time_derivative2_AB()

			V_x_pred = ( Vb_x + dt * ( - 9.0D0 * V_x_dot_m3 + 37.0D0 * V_x_dot_m2 - 59.0D0 * V_x_dot_m1 + &
			 													 55.0D0 * V_x_dot ) / 24.0D0 ) * I_fac
			V_y_pred = ( Vb_y + dt * ( - 9.0D0 * V_y_dot_m3 + 37.0D0 * V_y_dot_m2 - 59.0D0 * V_y_dot_m1 + &
			 													 55.0D0 * V_y_dot ) / 24.0D0 ) * I_fac
			V_z_pred = ( Vb_z + dt * ( - 9.0D0 * V_z_dot_m3 + 37.0D0 * V_z_dot_m2 - 59.0D0 * V_z_dot_m1 + &
			 													 55.0D0 * V_z_dot ) / 24.0D0 ) * I_fac

			CALL time_derivative2_AB_pred()

			Vb_x      = ( Vb_x + dt * ( V_x_dot_m2 - 5.0D0 * V_x_dot_m1 + 19.0D0 * V_x_dot + &
			 								  9.0D0 * V_x_dot_m3 ) / 24.0D0 ) * I_fac
			Vb_y      = ( Vb_y + dt * ( V_y_dot_m2 - 5.0D0 * V_y_dot_m1 + 19.0D0 * V_y_dot + &
			 								  9.0D0 * V_y_dot_m3 ) / 24.0D0 ) * I_fac
			Vb_z      = ( Vb_z + dt * ( V_z_dot_m2 - 5.0D0 * V_z_dot_m1 + 19.0D0 * V_z_dot + &
			 								  9.0D0 * V_z_dot_m3 ) / 24.0D0 ) * I_fac
			! Predicted 'v_dot' is stored in 'v_dot_m3' - to save space :)

			! Shifting the known velocities for next step
			V_x_dot_m3 = V_x_dot_m2
			V_x_dot_m2 = V_x_dot_m1
			V_x_dot_m1 = V_x_dot

			V_y_dot_m3 = V_y_dot_m2
			V_y_dot_m2 = V_y_dot_m1
			V_y_dot_m1 = V_y_dot

			V_z_dot_m3 = V_z_dot_m2
			V_z_dot_m2 = V_z_dot_m1
			V_z_dot_m1 = V_z_dot

			CALL fft_c2r( Vb_x, Vb_y, Vb_z, N, Nh, Ub_x, Ub_y, Ub_z )
			! Real Velocity

		END IF

	END

	SUBROUTINE time_derivative2_AB()
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

		W_kx = i * ( K_y * Vb_z - K_z * Vb_y )
    W_ky = i * ( K_z * Vb_x - K_x * Vb_z )
    W_kz = i * ( K_x * Vb_y - K_y * Vb_x )
    ! Spectral Vorticity

    CALL fft_c2r( W_kx, W_ky, W_kz, N, Nh, W_x, W_y, W_z )
    ! Real Vorticity

		CALL fft_c2r( Vb_x, Vb_y, Vb_z, N, Nh, Ub_x, Ub_y, Ub_z )

		UcW_x = Ub_y * W_z - Ub_z * W_y
		UcW_y = Ub_z * W_x - Ub_x * W_z
		UcW_z = Ub_x * W_y - Ub_y * W_x
		! Cross product between velocity and vorticity

		CALL fft_r2c( UcW_x, UcW_y, UcW_z,  N, Nh, VcW_x, VcW_y, VcW_z )
		! FFT to get spectral cross product

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.

		V_x_dot = Trunc * ( Proj_xx * VcW_x + Proj_xy * VcW_y + Proj_zx * VcW_z )
		V_y_dot = Trunc * ( Proj_xy * VcW_x + Proj_yy * VcW_y + Proj_yz * VcW_z )
		V_z_dot = Trunc * ( Proj_zx * VcW_x + Proj_yz * VcW_y + Proj_zz * VcW_z )

	END

	SUBROUTINE time_derivative2_AB_pred()
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

 		W_kx = i * ( K_y * V_z_pred - K_z * V_y_pred )
    W_ky = i * ( K_z * V_x_pred - K_x * V_z_pred )
    W_kz = i * ( K_x * V_y_pred - K_y * V_x_pred )
    ! Spectral Vorticity

    CALL fft_c2r( W_kx, W_ky, W_kz, N, Nh, W_x, W_y, W_z )
    ! Real Vorticity

		CALL fft_c2r( V_x_pred, V_y_pred, V_z_pred, N, Nh, Ub_x, Ub_y, Ub_z )

		UcW_x = Ub_y * W_z - Ub_z * W_y
		UcW_y = Ub_z * W_x - Ub_x * W_z
		UcW_z = Ub_x * W_y - Ub_y * W_x
		! Cross product between velocity and vorticity

		CALL fft_r2c( UcW_x, UcW_y, UcW_z,  N, Nh, VcW_x, VcW_y, VcW_z )
		! FFT to get spectral cross product

		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		!   3   D  -   E   U   L   E   R           E   Q   N.
		! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
		! Get the advection term 'u\cdot \Nabla u' in spectral space and project it with projection matrix and finally truncate it.

		V_x_dot_m3 = Trunc * ( Proj_xx * VcW_x + Proj_xy * VcW_y + Proj_zx * VcW_z )
		V_y_dot_m3 = Trunc * ( Proj_xy * VcW_x + Proj_yy * VcW_y + Proj_yz * VcW_z )
		V_z_dot_m3 = Trunc * ( Proj_zx * VcW_x + Proj_yz * VcW_y + Proj_zz * VcW_z )

	END

END MODULE system_solver2
