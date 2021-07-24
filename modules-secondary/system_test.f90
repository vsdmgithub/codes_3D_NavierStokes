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
! MODULE: system_test
! LAST MODIFIED: 3 JUNE 2021
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! TESTING MODULE FOR 3D NSE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_test
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This is the test module. It does the following
! Initializes variables and arrays
! Gets a initial condition
! Calculates the time for one FFT + iFFT operation
! Calculates the time for one time-evolution operation
! Estimates the total time that would be required
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicfunctions
  USE system_solver

  IMPLICIT NONE
  ! _________________________
  !  VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION      :: time_start,time_end
  DOUBLE PRECISION      :: time_micros,time_seconds,time_minutes
  INTEGER(KIND=4)       :: time_estimate
  CHARACTER( LEN = 40 ) :: date_sim
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  CONTAINS

  SUBROUTINE test_fft_time
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! Measures time for FFT and iFFT data
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !                 T I M E R       S T A R T
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CALL CPU_TIME( time_start )

    CALL fft_c2r( v_x, v_y, v_z, N, Nh, u_x, u_y, u_z )
    ! FFT spectral to real velocity

    CALL fft_r2c( u_x, u_y, u_z, N, Nh, v_x, v_y, v_z )
    ! i-FFT real to spectral velocity

    CALL CPU_TIME( time_end )
    time_micros   =    ( time_end - time_start ) * 1000.0D0

    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !                 T I M E R       E N D
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    WRITE(*,*)
		WRITE(*,'(A60)')			TRIM( ADJUSTL( '-_-_-_-_-_-_-_-_-_-_-_--_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-' ) )
		WRITE(*,*)
		WRITE(*,'(A60)')			TRIM( ADJUSTL( ' |||   F F T    T E S T    C O M P L E T E D  ||| :' ) )
    WRITE(*,*)
		WRITE(*,'(A60)')			TRIM( ADJUSTL( ' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH' ) )
		WRITE(*,'(A40,F8.2)')TRIM( ADJUSTL( ' TIME PER FFT (MICROSECONDS) :  ' ) ), time_micros
		WRITE(*,'(A60)')			TRIM( ADJUSTL( ' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH' ) )

	END

  SUBROUTINE test_evolution_time
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! Measures time for increment of one time step
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION :: tol_evolution_time

    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !                 T I M E R       S T A R T
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CALL CPU_TIME( time_start )

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L  G
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( solver_alg .EQ. 'rk') THEN
      CALL solver_RK4_algorithm
    ELSE
      CALL solver_AB4_algorithm
    END IF

    CALL compute_spectral_data

    CALL CPU_TIME( time_end )
    time_seconds   =    time_end - time_start

    ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
    !                 T I M E R       E N D
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    WRITE(*,*)
		WRITE(*,'(A60)')			TRIM( ADJUSTL( '-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-' ) )
		WRITE(*,*)
		WRITE(*,'(A60)')			TRIM( ADJUSTL( ' |||   E V O L U T I O N    T E S T    C O M P L E T E D  ||| :' ) )
    WRITE(*,*)
		WRITE(*,'(A60)')			TRIM( ADJUSTL( ' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH' ) )
		WRITE(*,'(A40,F8.2)')TRIM( ADJUSTL( ' TIME PER EVOLUTION (SECONDS) :  ' ) ), time_seconds
		WRITE(*,'(A60)')			TRIM( ADJUSTL( ' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH' ) )
    WRITE(*,*)

    tol_evolution_time = 0.05

    time_estimate      = CEILING( ( time_end - time_start ) * t_step_total * ( one + tol_evolution_time ) / 3600.0D0 )

    WRITE(*,'(A60)')			TRIM( ADJUSTL( '-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-' ) )
    WRITE(*,*)
    WRITE(*,'(A40,F12.6)')TRIM( ADJUSTL( ' MAXIMUM TIME STEP :  ' ) ), dt_max
    WRITE(*,'(A40,F12.6)')TRIM( ADJUSTL( ' SUGGESTED TIME STEP :  ' ) ), dt
    WRITE(*,'(A40,I8)')TRIM( ADJUSTL( ' TOTAL TIME STEPS :  ' ) ), t_step_total
    WRITE(*,*)
    WRITE(*,'(A60)')			TRIM( ADJUSTL( ' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH' ) )
    WRITE(*,'(A40,I4)')TRIM( ADJUSTL( ' ESTIMATED TIME (HRS) :' ) ), time_estimate
    WRITE(*,'(A60)')			TRIM( ADJUSTL( ' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH' ) )

    CALL deallocate_velocity
    CALL deallocate_operators

	END

 END MODULE system_test
