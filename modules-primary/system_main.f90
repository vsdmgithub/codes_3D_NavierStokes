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
    IF ( ( dt .LE. dt_max ) ) THEN

      CALL allocate_velocity
      ! Allocates velocity arrays for the system to start initialisation
      ! REF-> <<< system_basicvariables >>>

      CALL allocate_vorticity
      ! REF-> <<< system_basicvariables >>>

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

    END IF

    IF ( check_status .EQ. 1 ) THEN

      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
      !      S  O  L  V  E  R     A  L  G  O  R  I  T  H  M
      ! ----------------------------------------------------------------
      !      'ab'- ADAMBASHFORTH PRED & CORRECTOR ALG
      !      'rk'- RUNGA KUTTA 4TH ORDER ALG
      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH HHHHHHHHHHHHHHHHHHHHHHHHHHHHH
              solver_alg  = 'rk'
      ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

      check_status = 1

      CALL allocate_solver_system
      ! Allocates arrays for solver

      IF ( run_code .EQ. 'y' ) THEN

        CALL prepare_output
        ! Create names, folders to save files, open files in them to write data.
        ! REF-> <<< system_basicoutput >>>

        ! CALL allocate_PVD_subset_arrays
        ! Allocates arrays for PVD output for subset of data
        ! REF-> <<< system_pvdoutput >>>

      END IF

    END IF

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
    DO t_step = 0, t_step_total

      CALL inter_analysis
      ! Does all analysis in between time steps. Including saving data

      IF ( debug_error .EQ. 1 ) THEN
        EXIT
        ! Meaning some error in computation.
      END IF

      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      !  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L  G
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF ( ( solver_alg .EQ. 'rk') )  THEN

        CALL solver_RK4_algorithm
        ! REF-> <<< system_solver >>>
        GOTO 10101

      END IF
      IF ( ( solver_alg .EQ. 'ab') )  THEN

        CALL solver_AB4_algorithm
        ! REF-> <<< system_solver >>>
        GOTO 10101

      END IF
      ! Updates v_x,v_y,v_z for next time step

    10101 CONTINUE
    ! Jumps straight out of loop to here.

    END DO
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
print*,time_now,energy, diss_rate_viscous

    CALL forcing_spectrum
    ! REF-> <<< system_basicfunctions >>>

    CALL write_temporal_data
    ! REF-> <<< system_basicoutput >>>

    ! CALL write_test_data
    ! REF-> <<< system_basicoutput >>>

    IF (MOD(t_step,t_step_save) .EQ. 0) THEN

      CALL write_spectral_data
      ! REF-> <<< system_basicoutput >>>

    END IF

    IF (MOD(t_step,t_step_PVD_save) .EQ. 0) THEN

      ! CALL write_PVD_velocity
      ! REF-> <<< system_pvdoutput >>>

      ! CALL compute_vorticity
      ! REF-> <<< system_basicfunctions >>>

      ! CALL write_PVD_vorticity
      ! REF-> <<< system_pvdoutput >>>

      ! CALL write_PVD_vorticity_subset
      ! REF-> <<< system_pvdoutput >>>

    END IF

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  D  E  B  U  G             F  O  R          N  a   N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF (MOD(t_step,t_step_debug) .EQ. 0) THEN

      CALL perform_debug
      ! REF-> <<< system_basicfunctions >>>

      ! CALL print_running_status
      ! REF-> <<< system_basicoutput >>>

    END IF

   END

  SUBROUTINE post_analysis
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This does all the post analysis, making calls to write output after the evolution, debug and statistics part.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    ! CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
    ! Making sure, 'v' and 'u' are upto same evolution step

    ! CALL write_spectral_velocity
    ! REF-> <<< system_basicoutput >>>

    ! CALL write_velocity
    ! REF-> <<< system_basicoutput >>>

    ! CALL deallocate_PVD_subset_arrays
    ! REF-> <<< system_pvdoutput >>>

    CALL deallocate_solver_system

    CALL deallocate_velocity
    ! REF-> <<< system_basicvariables >>>

    CALL deallocate_vorticity
    ! REF-> <<< system_basicvariables >>>

    CALL deallocate_operators
    ! REF-> <<< system_basicvariables >>>

    state_sim = 1
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
