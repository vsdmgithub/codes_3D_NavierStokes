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
! MODULE: system_timer
! LAST MODIFIED: 10 NOVEMBER 2020
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MAIN TIME KEEPER MODULE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_timer
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This module keeps track of time for the simulation.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	IMPLICIT  NONE
  ! _________________________
  !  VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION      ::time_start,time_end
  DOUBLE PRECISION      ::time_elapsed
  INTEGER(KIND=4)       ::days,hours,mins,secs
  CHARACTER( LEN = 40 ) ::date_sim
  ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

  CONTAINS

	SUBROUTINE start_run_timer
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO: START THE TIMER
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		CALL CPU_TIME( time_start )
		CALL FDATE( date_sim )

		WRITE(*,'(A60)')	TRIM( ADJUSTL( ' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH' ) )
		WRITE(*,'(2A30)') TRIM( ADJUSTL( ' PROGRAM STARTED ON :' ) ),date_sim
		WRITE(*,'(A60)')  TRIM( ADJUSTL( '-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-' ) )
		WRITE(*,*)
		WRITE(*,*)
		WRITE(*,*)
		! Stores time when the program is initiated

	END

	SUBROUTINE end_run_timer
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO: FINISH THE TIMER
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE

		CALL FDATE( date_sim )
		CALL CPU_TIME( time_end )

		time_elapsed = time_end - time_start
		days         = FLOOR( time_elapsed / 86400.0D0 )
		hours        = FLOOR( ( time_elapsed - DBLE( days * 86400 ) ) / 3600.0D0 )
		mins         = FLOOR( ( time_elapsed - DBLE( days * 86400 ) - DBLE( hours * 3600 ) ) / 60.0D0 )
		secs         = FLOOR( time_elapsed - DBLE( days * 86400 ) - DBLE( hours * 3600 ) - DBLE( mins * 60.0D0 ) ) + 1

		WRITE(*,*)
		WRITE(*,*)
		WRITE(*,*)
		WRITE(*,'(A60)')				  			TRIM( ADJUSTL( '-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-' ) )
		WRITE(*,*)
		WRITE(*,'(2A30)')								TRIM( ADJUSTL( ' PROGRAM ENDED ON :' ) ), date_sim
		WRITE(*,'(A60)')						  	TRIM( ADJUSTL( ' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH' ) )
		WRITE(*,'(A30,I2.2,A2,I2.2,A2,I2.2,A2,I2.2,A2)') &
		 																TRIM( ADJUSTL( 'TOTAL RUN TIME :  ' ) ),days, TRIM( ADJUSTL( 'D:' ) ), hours, &
																		TRIM( ADJUSTL( 'H:' ) ), mins,TRIM( ADJUSTL( 'M:' ) ), secs,TRIM( ADJUSTL( 'S:' ) )
		WRITE(*,'(A60)')							  TRIM( ADJUSTL( ' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH' ) )

	END

END MODULE system_timer
