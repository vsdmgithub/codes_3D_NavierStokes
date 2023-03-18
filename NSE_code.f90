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
! ##################
! PROGRAM THAT SOLVES 3D INCOMPRESSIBLE NAVIER STOKES EQUATION
! LAST MODIFIED: 20 FEBRAURY 2023
! ##################

! ##################
! LIST OF MODULES:
! ___________________
! PRIMARY
! ===================
! system_main
! system_solver
! system_basicfunctions
! system_advfunctions
! system_initialcondition
! system_basicvariables
! system_basicdeclaration
! system_advdeclaration
! system_advvariables
! system_basicoutput
! system_advoutput
! -------------------
! ##################

PROGRAM NSE
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the work is done in the modules.
! Calling the correct modules would finish the simulation
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
!  MODULES
!  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
USE system_main
USE system_test
USE system_timer

	IMPLICIT NONE
	! _________________________
	! LOCAL VARIABLES
	! !!!!!!!!!!!!!!!!!!!!!!!!!

	CALL start_run_timer
	! Clocks the date and time - PROGRAM STARTS
	! REF-> <<< system_timer >>>

  !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !  I  N  I  T  I  A  L  I  Z  A  T  I  O  N
  !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  CALL read_input
	! REF-> <<< system_basicdeclaration >>>

  ! vis = ECCODE ! this will be batched
	vis = 8.0D0

  CALL init_global_variables
	! REF-> <<< system_basicdeclaration >>>

	CALL init_global_arrays
	! REF-> <<< system_basicdeclaration >>>

	tst_code		= 'n'
	! Simple way to ON or OFF the testing of the simulation - Measures the time for evolution

	run_code 		= 'y'
	! Simple way to ON or OFF the running of the simulation.

  CODE_EVALUTION:IF (( run_code .EQ. 'y' ) .AND. ( tst_code .EQ. 'n')) THEN

    CALL pre_analysis
    ! Allocating the evolution arrays, if everything is set, 'pre_status' will be 1.
		! REF-> <<< system_main >>>

    PRECHECK_STATUS_102: IF ( pre_status .EQ. 1 ) THEN

			WRITE(*,'(A60)')	TRIM(ADJUSTL('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'))
			WRITE(*,'(A60)')	TRIM(ADJUSTL('   S  I  M  U  L  A  T  I  O  N        S  T  A  R  T  S '))
			WRITE(*,'(A60)')	TRIM(ADJUSTL('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'))
			WRITE(*,*)

      CALL time_evolution
      ! Solve the 3D NSE equation, in discrete time using pseudospectral method.
			! REF-> <<< system_main >>>

			CALL prepare_perturbation
			! Makes a copy and pertubs by one time step evolution with and without forcing for them.
			! REF-> <<< system_main >>>

      CALL time_evolution_lyapunov
      ! Solve the 3D NSE equation, in discrete time using pseudospectral method.
			! REF-> <<< system_main >>>

			IF ( deb_err .NE. 1 ) THEN

			  CALL post_analysis
	      ! Does the post-analysis, final outputs and deallocation
				! REF-> <<< system_main >>>

			END IF

    END IF PRECHECK_STATUS_102

    IF ( sim_status .EQ. 1 ) THEN

			WRITE(*,*)
			WRITE(*,'(A60)')	TRIM(ADJUSTL('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'))
			WRITE(*,'(A60)')	TRIM(ADJUSTL('    S  I  M  U  L  A  T  I  O  N        E  N  D  S  '))
			WRITE(*,'(A60)')	TRIM(ADJUSTL('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'))

    END IF

  END IF CODE_EVALUTION

	TEST_EVOLUTION:IF ( tst_code .EQ. 'y' ) THEN

		CALL pre_analysis
    ! Allocating the evolution arrays
		! REF-> <<< system_main >>>

		CALL test_fft_time
		! Check forward and inverse fft time
		! REF-> <<< system_test >>>

		CALL test_evolution_time
		! Checks time for one evolution step - predicts total time
		! REF-> <<< system_test >>>

	  CALL post_analysis
		! REF-> <<< system_main >>>

	END IF TEST_EVOLUTION

	CALL end_run_timer
	! Clocks the date and time - PROGRAM STARTS
	! REF-> <<< system_timer >>>


END PROGRAM NSE
