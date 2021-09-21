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
! LAST MODIFIED: 21 JUNE 2021
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
  USE system_solver
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
      ! Allocates velocity arrays for the system to start initialisation
      ! REF-> <<< system_basicdeclaration >>>

      CALL allocate_vorticity
      ! REF-> <<< system_basicdeclaration >>>

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  I  N  I  T  I  A  L        C  O  N  D  I  T  I  O  N
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      CALL normalized_initial_condition
      ! Calls initial condition, checks for NaN and incompressibility
      ! then computes energy, then does FFT
      ! REF-> <<< system_basicfunctions >>>

    ELSE

      CALL print_error_timestep()
      ! REF-> <<< system_basicoutput >>>

    END IF TIME_STEP_MIN_CHECK

    PRECHECK_STATUS_101: IF ( precheck_status .EQ. 1 ) THEN

      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
      !      S  O  L  V  E  R     A  L  G  O  R  I  T  H  M
      ! ----------------------------------------------------------------
      !      'ab'- ADAMBASHFORTH PRED & CORRECTOR ALG
      !      'rk'- RUNGA KUTTA 4TH ORDER ALG
      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
              solver_alg  = 'rk'
      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      CALL allocate_solver_system
      ! Allocates arrays for solver

      RUN_CODE_CHECK:IF ( run_code .EQ. 'y' ) THEN

        CALL prepare_output
        ! Create names, folders to save files, open files in them to write data.
        ! REF-> <<< system_basicoutput >>>

        ! CALL allocate_PVD_subset_arrays
        ! Allocates arrays for PVD output for subset of data
        ! REF-> <<< system_pvdoutput >>>

      END IF RUN_CODE_CHECK

    END IF PRECHECK_STATUS_101

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

    !  ---------------------------------------------------------------
    !             T   I   M   E       L  O  O  P
    ! _________________________________________________________________
    LOOP_TIME_EVOLUTION: DO t_step = 0, t_step_total

      CALL inter_analysis
      ! Does all analysis in between time steps. Including saving data

      DEBUG_ERROR_CHECK:IF ( debug_error .EQ. 1 ) THEN
        EXIT
        ! Meaning some error in computation.
      END IF DEBUG_ERROR_CHECK

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L  G
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      RK4_CHECK:IF ( ( solver_alg .EQ. 'rk') )  THEN

        CALL solver_RK4_algorithm
        ! REF-> <<< system_solver >>>
        GOTO 10101

      END IF RK4_CHECK

      AB4_CHECK:IF ( ( solver_alg .EQ. 'ab') )  THEN

        CALL solver_AB4_algorithm
        ! REF-> <<< system_solver >>>
        GOTO 10101

      END IF AB4_CHECK
      ! Updates v_x,v_y,v_z for next time step

    10101 CONTINUE
    ! Jumps straight out of loop to here.

  END DO LOOP_TIME_EVOLUTION
    ! _________________________________________________________________

    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !      E   U   L   E   R      E   V   O   L   U   T   I   O   N
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    !                    E     N     D
    ! 8888888888888888888888888888888888888888888888888888888888888888

	END

  SUBROUTINE inter_analysis
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! -------------
  ! This does all the inter_analysis, making calls to write output during the evolution, debug and statistics part.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    CALL step_to_time_convert(t_step,time_now,dt)
    ! Converts the 't_step' to actual time 'time_now'
    ! REF-> <<< system_auxilaries >>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  N  A  L  Y  S  I  S       C   A   L   C  .
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    CALL compute_spectral_data
    ! REF-> <<< system_basicfunctions >>>
    FORCING_CHECK_401: IF ( forcing_status .EQ. 1 ) THEN

      CALL compute_forcing_in_modes
      ! REF-> <<< system_basicfunctions >>>

    END IF FORCING_CHECK_401

    CALL write_temporal_data
    ! REF-> <<< system_basicoutput >>>

    CALL write_test_data
    ! REF-> <<< system_basicoutput >>>

    SAVE_DATA_CHECK:IF (MOD(t_step,t_step_save) .EQ. 0) THEN

      CALL write_spectral_data
      ! REF-> <<< system_basicoutput >>>

    END IF SAVE_DATA_CHECK

    SAVE_DATA_CHECK_3D:IF (MOD(t_step,t_step_3D_save) .EQ. 0) THEN

      ! CALL write_PVD_velocity
      ! REF-> <<< system_pvdoutput >>>

      ! CALL compute_vorticity
      ! REF-> <<< system_basicfunctions >>>

      ! CALL write_PVD_vorticity
      ! REF-> <<< system_pvdoutput >>>

      ! CALL write_PVD_vorticity_subset
      ! REF-> <<< system_pvdoutput >>>

    END IF SAVE_DATA_CHECK_3D

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  D  E  B  U  G             F  O  R          N  a   N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEBUG:IF (MOD(t_step,t_step_debug) .EQ. 0) THEN

      ! CALL perform_debug
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

    CALL compute_system_details
    ! REF-> <<< system_basicdeclaration >>>

    CALL write_simulation_end_details
    ! Writes the parameters used in the simulation at end

    ! CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
    ! Making sure, 'v' and 'u' are upto same evolution step

    ! CALL write_spectral_velocity
    ! REF-> <<< system_basicoutput >>>

    CALL write_velocity
    ! REF-> <<< system_basicoutput >>>

    ! CALL write_velocity_unformatted
    ! REF-> <<< system_basicoutput >>>

    ! CALL allocate_dissipation_field
    ! REF-> <<< system_advvariables >>>

    ! CALL compute_dissipation_field
    ! REF-> <<< system_advfunctions >>>

    ! CALL deallocate_dissipation_field
    ! REF-> <<< system_advvariables >>>

    ! CALL deallocate_PVD_subset_arrays
    ! REF-> <<< system_pvdoutput >>>

    CALL deallocate_solver_system

    CALL deallocate_velocity
    ! REF-> <<< system_basicdeclaration >>>

    CALL deallocate_vorticity
    ! REF-> <<< system_basicdeclaration >>>

    CALL deallocate_operators
    ! REF-> <<< system_basicdeclaration >>>

    simulation_status = 1
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

    IF ( solver_alg .EQ. 'ab') THEN

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

    IF ( solver_alg .EQ. 'ab') THEN

      CALL deallocate_solver_AB4
      ! REF-> <<< system_solver >>>

    ELSE

      CALL deallocate_solver_RK4
      ! REF-> <<< system_solver >>>

    END IF
    ! Deallocates arrays for solving

  END

 END MODULE system_main
