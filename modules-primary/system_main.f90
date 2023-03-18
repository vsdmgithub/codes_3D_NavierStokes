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
! MODULE: system_main
! LAST MODIFIED: 20 FEBRAURY 2023
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MAIN RUN MODULE TO RUN THE TIME EVOLUTION FOR 3D NSE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_main
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This is the main module. Consists of subroutines
! 1. Pre-analysis
! 2. Time-evolution
! 3. Inter-analysis
! 4. Post-analysis
! All the other modules are sub-modules to this.
! Then nesting of all major sub-modules is as follows

! MAIN-RUN MODULE
!   |
!   ∟ ---> BASIC FUNCTION
!   |
!   ∟ ---> BASIC OUTPUT
!   |
!   ∟ ---> ADV FUNCTION
!   |
!   ∟ ---> ADV OUTPUT
!   |
!   ∟ ---> SOLVER
!   |
!   ∟ ---> INITIAL CONDITION
!   |
!   ∟ ---> BASIC DECLARATION
!   |
!   ∟ ---> ADV DECLARATION
!   |
!   ∟ ---> BASIC VARIABLES
!   |
!   ∟ ---> ADV VARIABLES
!
! There are other modules, which are subsidary modules.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicfunctions
  USE system_advfunctions
  USE system_pvdoutput
  USE system_solver2
  USE system_decorrelator
  IMPLICIT NONE
  ! _________________________
  ! LOCAL VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION:: temp_data
  CONTAINS

  SUBROUTINE pre_analysis
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! Loop of time steps, where at each step the spectral velocities
  ! are updated through any of the algoritm. Meanwhile, analysis and
  ! outputs are printed respectively.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !       T    I    M     E              S    T    E    P              C   H    E   C   K
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    TIME_STEP_MIN_CHECK:IF ( ( dt .LE. dt_max ) ) THEN

      CALL allocate_velocity
      ! REF-> <<< system_basicdeclaration >>>

      CALL allocate_vorticity
      ! REF-> <<< system_basicdeclaration >>>

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  I  N  I  T  I  A  L        C  O  N  D  I  T  I  O  N
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      CALL normalized_initial_condition
      ! REF-> <<< system_basicfunctions >>>

    ELSE

      CALL print_error_timestep()
      ! REF-> <<< system_basicoutput >>>

    END IF TIME_STEP_MIN_CHECK

    PRECHECK_STATUS_101: IF ( pre_status .EQ. 1 ) THEN

      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
      !      S  O  L  V  E  R     A  L  G  O  R  I  T  H  M
      ! ----------------------------------------------------------------
      !      'ab'- ADAMBASHFORTH PRED & CORRECTOR ALG
      !      'rk'- RUNGA KUTTA 4TH ORDER ALG
      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
              sol_algm  = 'rk'
      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      CALL allocate_solver_system
      ! Allocates arrays for solver

      RUN_CODE_CHECK:IF(( run_code .EQ. 'y' ) .AND. ( tst_code .EQ. 'n')) THEN

        CALL prepare_output
        ! Create names, folders to save files, open files in them to write data.
        ! REF-> <<< system_basicoutput >>>

        ! CALL allocate_PVD_subset_arrays
        ! Allocates arrays for PVD output for subset of data
        ! REF-> <<< system_pvdoutput >>>

      END IF RUN_CODE_CHECK

    END IF PRECHECK_STATUS_101

  END

  SUBROUTINE prepare_perturbation
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! Call this to create a perturbation
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

		CALL allocate_lyapunov_arrays
		! REF-> <<< system_basicdeclaration >>>

    Vb_x = V_x
    Vb_y = V_y
    Vb_z = V_z

		CALL create_perturbation_by_forcing

		CALL prepare_output_lyapunov
	  ! REF-> <<< system_basicoutput >>>

		lyp_status = 1

    CALL allocate_decorrelator
    ! REF-> <<< system_decorrelator >>>

    dec_pre = zero + tol_double

    CALL compute_decorrelator
    ! REF-> <<< system_decorrelator >>>

    dec_ini = dec

    WRITE(*,"(A22,ES8.2)")'DECOR_INITIAL:= ',dec_ini

	END

  SUBROUTINE create_perturbation_by_forcing
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! To create the perturbation
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT  NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CALL compute_random_forcing_fixed_dissipation
    ! REF-> <<< system_basicfunctions >>>

    CALL solver_RK4_algorithm
    ! CALL solver_AB4_algorithm
    ! REF-> <<< system_solver >>>

		CALL perturb_forcing
    ! REF-> <<< system_basicfunctions >>>

    CALL solver2_RK4_algorithm
    ! CALL solver2_AB4_algorithm
    ! REF-> <<< system_solver >>>

  END

  SUBROUTINE time_evolution
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! Loop of time steps, where at each step the spectral velocities
  ! are updated through any of the algoritm. Meanwhile, analysis and
  ! outputs are printed respectively.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !                 E   V   O   L   U   T   I   O   N
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    !             S        T         A         R       T
    ! 8888888888888888888888888888888888888888888888888888888888888888

    ! ---------------------------------------------------------------
    !             T   I   M   E       L  O  O  P
    ! _________________________________________________________________
    LOOP_TIME_EVOLUTION: DO t_step = 0, t_step_tot

      CALL inter_analysis ! Does all analysis in between time steps. Including saving data

      DEBUG_ERROR_CHECK:IF ( deb_err .EQ. 1 ) THEN
        EXIT ! Meaning some error in computation.
      END IF DEBUG_ERROR_CHECK

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L  G
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL solver_RK4_algorithm
      ! CALL solver_AB4_algorithm
      ! REF-> <<< system_solver >>>

	  END DO LOOP_TIME_EVOLUTION

	  CALL compute_system_details
		! REF-> <<< system_basicdeclaration >>>

    CALL write_simulation_details('end')
    ! REF-> <<< system_basicoutput >>>
    ! Writes the parameters used in the simulation at end

    ! _________________________________________________________________
    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    !                    E     N     D
    ! 8888888888888888888888888888888888888888888888888888888888888888

	END

  SUBROUTINE time_evolution_lyapunov
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! Loop of time steps, where at each step the spectral velocities
  ! are updated through any of the algoritm. Meanwhile, analysis and
  ! outputs are printed respectively.
	! Done for both systems
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !                 E   V   O   L   U   T   I   O   N
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    !             S        T         A         R       T
    ! 8888888888888888888888888888888888888888888888888888888888888888

    !  ---------------------------------------------------------------
    !             T   I   M   E       L  O  O  P
    ! _________________________________________________________________
		t_step = 0
    LOOP_TIME_EVOLUTION_2: DO WHILE ( dec .LT. 0.8D0 )
      CALL inter_analysis
      ! Does all analysis in between time steps. Including saving data

      DEBUG_ERROR_CHECK:IF ( deb_err .EQ. 1 ) THEN
        EXIT
        ! Meaning some error in computation.
      END IF DEBUG_ERROR_CHECK

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L  G
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL solver_RK4_algorithm
      CALL solver2_RK4_algorithm
      ! REF-> <<< system_solver >>>

      ! CALL solver_AB4_algorithm
      ! CALL solver2_AB4_algorithm
      ! REF-> <<< system_solver >>>

		t_step = t_step + 1
	  END DO LOOP_TIME_EVOLUTION_2

    ! _________________________________________________________________

    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    !                    E     N     D
    ! 8888888888888888888888888888888888888888888888888888888888888888

    dec_fin = dec

    WRITE(*,"(A22,ES8.2)")'DECOR_FINAL:= ',dec_fin

	END

  SUBROUTINE inter_analysis
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! -------------
  ! This does all the inter_analysis, making calls to write output during the evolution, debug and statistics part.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    CALL step_to_time_convert(t_step,t_now,dt)
    ! Converts the 't_step' to actual time 'time_now'
    ! REF-> <<< system_auxilaries >>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  N  A  L  Y  S  I  S       C   A   L   C  .
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

		CALL compute_temporal_data
    ! REF-> <<< system_basicfunctions >>>

		IF ( lyp_status .EQ. 1 ) THEN
	    CALL compute_decorrelator
	    ! REF-> <<< system_decorrelator >>>
		END IF

    FORCING_CHECK_401: IF ( frc_status .EQ. 1 ) THEN

			IF ( lyp_status .EQ. 0 ) THEN
				IF ( t_now .LT. two * t_int ) THEN
		      CALL compute_random_forcing_uncorrelated_velocity
					dis_fix = frc_fac_avg
			    ! IF (MOD(t_step,20) .EQ. 0) THEN
					! 	PRINT*,t_now,dis,frc_fac,eng,frc_fac_avg
					! END IF

				ELSE

		      CALL compute_random_forcing_fixed_dissipation
			    ! IF ( MOD( t_step, 20 ) .EQ. 0 ) THEN
					! 	PRINT*,t_now,dis,eng
					! END IF

				END IF
			ELSE

	      CALL compute_random_forcing_fixed_dissipation
		    ! IF ( MOD( t_step, 20 ) .EQ. 0 ) THEN
				! 	PRINT*,t_now,dis,eng
				! END IF

			END IF

    END IF FORCING_CHECK_401

    ! CALL write_test_data
    ! REF-> <<< system_basicoutput >>>

    SAVE_DATA_CHECK:IF (MOD(t_step,t_step_sav) .EQ. 0) THEN
	    CALL compute_spectral_data
	    ! REF-> <<< system_basicfunctions >>>

	    CALL write_spectral_data
	    ! REF-> <<< system_basicoutput >>>

			IF ( lyp_status .EQ. 1 ) THEN
		    CALL compute_pdf_lyapunov
		    ! REF-> <<< system_decorrelator >>>
		    CALL write_section('dec_fld',Dec_fld( 0, : , : ) )
		    CALL write_section('lyp_fld',Lyp_fld( 0, : , : ) )
		    CALL write_section('dis_fld',Dis_fld( 0, : , : ) )
		    ! REF <<< system_basicoutput >>>
			ELSE
		    CALL write_section('sec_vorZ',W_z( 0 , : , : ) )
		    ! REF-> <<< system_basicoutput >>>
			END IF

    END IF SAVE_DATA_CHECK

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  D  E  B  U  G             F  O  R          N  a   N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEBUG:IF (MOD(t_step,t_step_deb) .EQ. 0) THEN

      CALL perform_debug
      ! REF-> <<< system_basicfunctions >>>

      CALL print_running_status
      ! REF-> <<< system_basicoutput >>>

    END IF DEBUG

  END

  SUBROUTINE post_analysis
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This does all the post analysis, making calls to write output after the evolution, debug and statistics part.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    ! CALL compute_velocity_gradient
    ! REF-> <<< system_advfunctions >>>

    ! CALL fft_c2r( V_x, V_y, V_z, N, Nh, U_x, U_y, U_z )
    ! Making sure, 'v' and 'u' are upto same evolution step

    ! CALL write_spectral_velocity
    ! REF-> <<< system_basicoutput >>>

    ! CALL write_velocity
    ! REF-> <<< system_basicoutput >>>

    ! CALL write_velocity_unformatted
    ! REF-> <<< system_basicoutput >>>

    ! CALL deallocate_PVD_subset_arrays
    ! REF-> <<< system_pvdoutput >>>

    CALL deallocate_solver_system

		IF ( lyp_status .EQ. 1 ) THEN

	    CALL deallocate_decorrelator
	    ! REF-> <<< system_decorrelator >>>

			CALL deallocate_lyapunov_arrays
	    ! REF-> <<< system_basicdeclaration >>>

		END IF

    CALL deallocate_vorticity
    ! REF-> <<< system_basicdeclaration >>>

    CALL deallocate_velocity
    ! REF-> <<< system_basicdeclaration >>>

    CALL deallocate_operators
    ! REF-> <<< system_basicdeclaration >>>

    sim_status = 1
    ! Stating that the simulation has ended.

  END

  SUBROUTINE allocate_solver_system
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This allocates the corresponding solver
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    CALL allocate_solver
      ! REF-> <<< system_solver >>>

    IF ( sol_algm .EQ. 'ab') THEN

      CALL allocate_solver_AB4
      ! REF-> <<< system_solver >>>

    ELSE

      CALL allocate_solver_RK4
      ! REF-> <<< system_solver >>>

    END IF
    ! Allocates arrays for solving

  END

  SUBROUTINE deallocate_solver_system
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This deallocates the corresponding solver
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    CALL deallocate_solver
    ! REF-> <<< system_solver >>>

    IF ( sol_algm .EQ. 'ab') THEN

      CALL deallocate_solver_AB4
      ! REF-> <<< system_solver >>>

    ELSE

      CALL deallocate_solver_RK4
      ! REF-> <<< system_solver >>>

    END IF
    ! Deallocates arrays for solving

  END

 END MODULE system_main
