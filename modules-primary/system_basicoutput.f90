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
! MODULE: system_basicoutput
! LAST MODIFIED: 20 FEBRAURY 2023
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! BASIC OUTPUT MODULE - RELATED TO BASIC FUNCTION MODULE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_basicoutput
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all basic outputs produced in the simulation.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicdeclaration

  IMPLICIT NONE
  ! _________________________
  ! OUTPUT VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  CHARACTER(LEN =180)::fil_name
  CHARACTER(LEN =40) ::dir
  CHARACTER(LEN =40) ::fil_time
  CHARACTER(LEN =40) ::pat_sim
  CHARACTER(LEN =60) ::nam_sim
  CHARACTER(LEN =100)::fil_adrs
  CHARACTER(LEN =40) ::dir_sec
  CHARACTER(LEN =40) ::dir_out
  CHARACTER(LEN =40) ::dir_pdf

  CONTAINS

  SUBROUTINE prepare_output
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Name the folders, create them, open files to write system_basicoutput .
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE

    CALL create_output_dir
    ! Creates the directories

    CALL write_simulation_details('start')
    ! Writes the parameters used in the simulation

	END

  SUBROUTINE create_output_dir
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Create the directories.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE

		pat_sim    =   '../data/'
    ! path of the output directory for this simulation relative to this file.

    CALL get_simulation_name(nam_sim)
    ! Creating dated and timed name for the simulation for this particular type
    ! REF:- <<< system_auxilaries >>>

		nam_sim = 'N' // TRIM( ADJUSTL( res_char ) ) // '_' // TRIM( ADJUSTL( nam_sim ) ) // '/'
    ! nam_sim    =   'N256_run_V8'

		dir     = 'tur/'
    dir_sec = 'sec/'
    dir_out = 'out/'
    dir_pdf = 'pdf/'

    fil_adrs =  TRIM( ADJUSTL( pat_sim ) )//TRIM( ADJUSTL( nam_sim ) )//TRIM( ADJUSTL( dir ) )
    ! Address should be added to all file names, if needed sub-dir can be declared later and appended to
    ! this address

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL( pat_sim ) ) )
    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL( pat_sim ) ) // TRIM( ADJUSTL( nam_sim ) ) )
    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( fil_adrs ) ) )
    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( fil_adrs ) ) // TRIM( ADJUSTL( dir_out ) ) )
    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( fil_adrs ) ) // TRIM( ADJUSTL( dir_sec ) ) )

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // 'den_states.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y      S  P  E  C  T  R  U  M
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 1301, FILE = fil_name )
    DO k_ind = 1 , k_tru

      WRITE(1301,f_i8,ADVANCE  ='no')       k_ind
      WRITE(1301,f_i8,ADVANCE  ='no')        Den_k( k_ind )
      WRITE(1301,f_d32p17,ADVANCE ='yes')   two_pi * DBLE( k_ind * k_ind )

    END DO
    CLOSE(1301)

	END

  SUBROUTINE prepare_output_lyapunov
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Name the folders, create them, open files to write system_basicoutput .
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE

		dir   =  'lyp/'

    fil_adrs =  TRIM( ADJUSTL( pat_sim ) )//TRIM( ADJUSTL( nam_sim ) )//TRIM( ADJUSTL( dir ) )
    ! Address should be added to all file names, if needed sub-dir can be declared later and appended to
    ! this address

    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( fil_adrs ) ) )
    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( fil_adrs ) ) // TRIM( ADJUSTL( dir_out ) ) )
    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( fil_adrs ) ) // TRIM( ADJUSTL( dir_sec ) ) )
    CALL SYSTEM('mkdir ' // TRIM( ADJUSTL ( fil_adrs ) ) // TRIM( ADJUSTL( dir_pdf ) ) )

	END

  SUBROUTINE write_test_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the test_data of any new verification/validation conducted in the code
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN=80)::test_fil_name

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   E N E R G Y    V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      test_fil_name = TRIM( ADJUSTL( fil_adrs ) ) // 'test_data.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 9009, file = test_fil_name )
      ! File where energy vs time will be written. With additional data
    END IF

    WRITE(9009,f_d8p4,ADVANCE   ='no')   t_now
    WRITE(9009,f_d32p17,ADVANCE ='yes')  U_x( 2, 2, 2)

    IF ( t_step .EQ. t_step_tot ) THEN
      CLOSE(9009)
    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_simulation_details(out_char)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  !   Write the details of the simulation
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT  NONE
    CHARACTER(LEN=*),INTENT(IN)::out_char

    fil_name = TRIM(ADJUSTL(fil_adrs))//'sys_details_'//TRIM(ADJUSTL(out_char))

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  I  M  U  L  A  T  I  O  N     D  E  T  A  I  L  S
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN(UNIT =233,FILE=TRIM( ADJUSTL( fil_name ) ) // '.dat')

    WRITE(233,"(A50)")TRIM(ADJUSTL('--------------------------------------------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('------3D NSE EQUATION (INCOMPRESSIBLE)------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('-------------------EULERIAN CHAOS----------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('__________________________________________________'))
    WRITE(233,*)
    WRITE(233,"(A50)")TRIM(ADJUSTL('---------------PARAMETERS OF SIMULATION------------'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('--------------------------------------------------'))
    WRITE(233,"(A20,A2,I8)")    'Resolution          ','= ',N
    WRITE(233,"(A20,A2,ES8.2)") 'Time step           ','= ',dt
    WRITE(233,"(A20,A2,F8.4)")  'Total time          ','= ',t_tot
    WRITE(233,"(A20,A2,I8)")    'Total time steps    ','= ',t_step_tot
    WRITE(233,"(A20,A2,I8)")    'No of saves         ','= ',num_sav
    WRITE(233,"(A20,A2,ES8.2)") 'Viscosity           ','= ',vis
    WRITE(233,"(A20,A2,I8)")    'Integr-Reynolds     ','= ',rey_int
    WRITE(233,"(A20,A2,I8)")    'Taylor-Reynolds     ','= ',rey_tay
    WRITE(233,"(A20,A2,I8)")    'Trunc. Mode         ','= ',k_tru
    WRITE(233,"(A20,A2,I8)")    'Kolmo. Mode         ','= ',k_kol
    WRITE(233,"(A20,A2,I8)")    'Taylor. Mode        ','= ',k_tay
    WRITE(233,"(A20,A2,I8)")    'Forcing. Mode       ','= ',k_frc
    WRITE(233,"(A20,A2,I8)")    'Integral. Mode      ','= ',k_int
    WRITE(233,"(A20,A3,I8)")    'Modes Active.       ','= ',num_mod
    WRITE(233,"(A20,A3,I8)")    'Modes Forced.       ','= ',num_mod_frc
    WRITE(233,"(A20,A2,F8.2)")  'Resolving power     ','= ',res_pow
    WRITE(233,"(A20,A2,I8)")    'CFL system          ','= ',cfl
    WRITE(233,"(A20,A2,F8.6)")  'Grid time           ','= ',t_grd
    WRITE(233,"(A20,A2,F8.6)")  'Kolmo. time         ','= ',t_kol
    WRITE(233,"(A20,A2,F8.6)")  'Integral time       ','= ',t_int
    WRITE(233,"(A20,A2,F8.4)")  'Rms Velocity        ','= ',u_int
    WRITE(233,"(A20,A2,F8.4)")  'Kol Velocity        ','= ',u_kol
    WRITE(233,"(A20,A2,F8.4)")  'Energy              ','= ',eng
    WRITE(233,"(A20,A2,F8.4)")  'Enstrophy           ','= ',ens
    ! WRITE(233,"(A20,A2,F8.4)")  'Helicity            ','= ',hel
    WRITE(233,"(A20,A2,F8.4)")  'Dissipation Rate    ','= ',dis
    WRITE(233,"(A20,A2,ES8.2)") 'Compressibility     ','= ',inc
    WRITE(233,"(A20,A2,A8)")    'Initial Condition   ','= ',TRIM( ADJUSTL( icn_type ) )
    WRITE(233,*)
    WRITE(233,"(A50)")TRIM(ADJUSTL('__________________________________________________'))
    WRITE(233,"(A50)")TRIM(ADJUSTL('--------------------------------------------------'))

    CLOSE(233)

    WRITE(*,"(A50)")TRIM(ADJUSTL('------3D NSE EQUATION (INCOMPRESSIBLE)------------'))
    WRITE(*,"(A50)")TRIM(ADJUSTL('-------------------EULERIAN CHAOS----------------'))
    WRITE(*,"(A50)")TRIM(ADJUSTL('__________________________________________________'))
    WRITE(*,*)
    WRITE(*,"(A50)")TRIM(ADJUSTL('---------------PARAMETERS OF SIMULATION------------'))
    WRITE(*,"(A50)")TRIM(ADJUSTL('--------------------------------------------------'))
    WRITE(*,"(A20,A2,I8)")    'Resolution          ','= ',N
    WRITE(*,"(A20,A2,ES8.2)") 'Time step           ','= ',dt
    WRITE(*,"(A20,A2,F8.4)")  'Total time          ','= ',t_tot
    WRITE(*,"(A20,A2,I8)")    'Total time steps    ','= ',t_step_tot
    WRITE(*,"(A20,A2,I8)")    'No of saves         ','= ',num_sav
    WRITE(*,"(A20,A2,ES8.2)") 'Viscosity           ','= ',vis
    WRITE(*,"(A20,A2,I8)")    'Integr-Reynolds     ','= ',rey_int
    WRITE(*,"(A20,A2,I8)")    'Taylor-Reynolds     ','= ',rey_tay
    WRITE(*,"(A20,A2,I8)")    'Trunc. Mode         ','= ',k_tru
    WRITE(*,"(A20,A2,I8)")    'Kolmo. Mode         ','= ',k_kol
    WRITE(*,"(A20,A2,I8)")    'Taylor. Mode        ','= ',k_tay
    WRITE(*,"(A20,A2,I8)")    'Forcing. Mode       ','= ',k_frc
    WRITE(*,"(A20,A2,I8)")    'Integral. Mode      ','= ',k_int
    WRITE(*,"(A20,A3,I8)")    'Modes Active.       ','= ',num_mod
    WRITE(*,"(A20,A3,I8)")    'Modes Forced.       ','= ',num_mod_frc
    WRITE(*,"(A20,A2,F8.2)")  'Resolving power     ','= ',res_pow
    WRITE(*,"(A20,A2,I8)")    'CFL system          ','= ',cfl
    WRITE(*,"(A20,A2,F8.6)")  'Grid time           ','= ',t_grd
    WRITE(*,"(A20,A2,F8.6)")  'Kolmo. time         ','= ',t_kol
    WRITE(*,"(A20,A2,F8.6)")  'Integral time       ','= ',t_int
    WRITE(*,"(A20,A2,F8.4)")  'Rms Velocity        ','= ',u_int
    WRITE(*,"(A20,A2,F8.4)")  'Kol Velocity        ','= ',u_kol
    WRITE(*,"(A20,A2,F8.4)")  'Energy              ','= ',eng
    WRITE(*,"(A20,A2,F8.4)")  'Enstrophy           ','= ',ens
    ! WRITE(*,"(A20,A2,F8.4)")  'Helicity            ','= ',hel
    WRITE(*,"(A20,A2,F8.4)")  'Dissipation Rate    ','= ',dis
    WRITE(*,"(A20,A2,ES8.2)") 'Compressibility     ','= ',inc
    WRITE(*,"(A20,A2,A8)")    'Initial Condition   ','= ',TRIM( ADJUSTL( icn_type ) )
    WRITE(*,*)
    WRITE(*,"(A50)")TRIM(ADJUSTL('__________________________________________________'))
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_temporal_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in time, for every timestep
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   E N E R G Y    V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      fil_name = TRIM( ADJUSTL( fil_adrs ) ) // 'temporal_data.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 4004, file = fil_name )
      ! File where energy vs time will be written. With additional data
    END IF


    WRITE(4004,f_d8p4,ADVANCE   ='no')  t_now
    WRITE(4004,f_d32p17,ADVANCE ='no')  eng
		IF ( lyp_status .EQ. 1 ) THEN
	    eng_B = hf * SUM( Ub_x ** two + Ub_y ** two + Ub_z ** two ) / N3
	    WRITE(4004,f_d32p17,ADVANCE ='no')  eng_B
		END IF
		! WRITE(4004,f_d32p17,ADVANCE ='no')  ens
		! WRITE(4004,f_d32p17,ADVANCE ='no')  hel
		WRITE(4004,f_d32p17,ADVANCE ='no')  dis
		IF ( lyp_status .EQ. 0 ) THEN
			WRITE(4004,f_d32p17,ADVANCE ='no')  frc_fac
			WRITE(4004,f_d32p17,ADVANCE ='no')  frc_fac_avg
		END IF
		WRITE(4004,f_d32p17,ADVANCE ='yes') dis_eng

    IF ( t_step .EQ. t_step_tot ) THEN
      CLOSE(4004)
    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_spectral_data
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Write the spectral energy into a file with time index
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    WRITE (fil_time,f_d8p4) t_now
    ! Writes 't_now' as a CHARACTER

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // TRIM( ADJUSTL( dir_out ) ) &
                // 'Eng_k_t_'//TRIM( ADJUSTL( fil_time ) ) // '.dat'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E  N  E  R  G  Y      S  P  E  C  T  R  U  M
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    OPEN( UNIT = 1001, FILE = fil_name )
    DO k_ind = 1 , k_tru

      WRITE(1001,f_i8,ADVANCE  ='no')       k_ind
      WRITE(1001,f_d32p17,ADVANCE ='no')    Eng_k(k_ind)
      WRITE(1001,f_d32p17,ADVANCE ='yes')   Eng_k_avg(k_ind)

    END DO
    CLOSE(1001)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_section(nam_data,data)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a section data.
  ! -------------

    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::nam_data
    DOUBLE PRECISION,DIMENSION(0:N-1,0:N-1),INTENT(IN)::data

    WRITE (fil_time,f_d8p4) t_now
    ! Writes 't_now' as a CHARACTER

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // TRIM( ADJUSTL ( dir_sec )) // &
	          	 TRIM( ADJUSTL( nam_data ) ) // '_t_' // TRIM( ADJUSTL ( fil_time ) ) // '.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! File where the section is written in readable format

    OPEN( UNIT = 374, FILE = fil_name )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO i_y = 0 , N - 1
    DO i_z = 0 , N - 1

      WRITE(374,f_d32p17,ADVANCE ='YES') data( i_y, i_z)

    END DO
    END DO

    CLOSE(374)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_section_unformatted(nam_data,data)
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a section data unformatted.
  ! -------------

    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)::nam_data
	  DOUBLE PRECISION,DIMENSION(0:N-1,0:N-1),INTENT(IN)::data
	  CHARACTER(LEN=80) :: msg
		INTEGER::istat

    WRITE (fil_time,f_d8p4) t_now
    ! Writes 't_now' as a CHARACTER

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // TRIM( ADJUSTL ( dir_sec )) // &
	          	 TRIM( ADJUSTL( nam_data ) ) // '_unf_t_' // TRIM( ADJUSTL ( fil_time ) ) // '.in'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! File where section is written in binary format

    OPEN( UNIT = 73, FILE = fil_name, FORM = 'unformatted', STATUS = 'new' , IOSTAT=istat, IOMSG=msg )
		! Check for OPEN error
		in_ok: IF ( istat /= 0 ) THEN
			WRITE (*,*) 'Input file OPEN failed: istat = ', istat
			WRITE (*,*) 'Error message = ', msg
		ELSE
	    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
	    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    WRITE(73) ( ( data(i_y,i_z), i_y=0,N-1 ), i_x=0,N-1 )
	    CLOSE(73)
		END IF in_ok
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_velocity_unformatted
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes a real velocity. To read and plot velocity functions
  ! Writtein unformatted.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
	  CHARACTER(LEN=80) :: msg
		INTEGER::istat

    WRITE (fil_time,f_d8p4) t_now
    ! Writes 't_now' as a CHARACTER

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // 'vel_' &
              //TRIM( ADJUSTL( res_char ) ) // '_t_' // TRIM( ADJUSTL ( fil_time ) ) // '.in'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! File where velocity is written in binary file

    OPEN( UNIT = 74, FILE = fil_name, FORM = 'unformatted', STATUS = 'new' , IOSTAT=istat, IOMSG=msg )
		! Check for OPEN error
		in_ok: IF ( istat /= 0 ) THEN
			WRITE (*,*) 'Input file OPEN failed: istat = ', istat
			WRITE (*,*) 'Error message = ', msg
		ELSE
	    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
	    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    WRITE(74) ((( U_x(i_x,i_y,i_z),U_y(i_x,i_y,i_z),U_z(i_x,i_y,i_z), &
	                  i_z = 0 , N - 1 ), i_y = 0 , N - 1 ), i_x = 0 , N - 1 )
	    CLOSE(74)
		END IF in_ok
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE print_running_status
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! print the running status of the program, when called during debug
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

		IF ( lyp_status .EQ. 0 ) THEN
	    IF ( t_step .EQ. 0 ) THEN
	      WRITE(*,'(A68)') &
				TRIM( ADJUSTL( '-----------------------------------------------------------------------------------' ) )
	      WRITE(*,'(A52)') &
				TRIM( ADJUSTL( &
				'| TIME | ENERGY | DISSIPATION | INCOMPRESSIBILITY |'  ) )
	      WRITE(*,'(A68)') &
				TRIM( ADJUSTL( '-----------------------------------------------------------------------------------' ) )
	    END IF
	    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    WRITE(*,'(F6.2,A3,F8.4,A3,F10.6,A5,E10.2)')&
			 t_now,'   ',eng,'    ',dis,'   ',inc
	    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    IF ( t_step .EQ. t_step_tot ) THEN
	      PRINT*,'-----------------------------------------------------------'
	    END IF
		ELSE
	    IF ( t_step .EQ. 0 ) THEN
	      WRITE(*,'(A83)') &
				TRIM( ADJUSTL( '-----------------------------------------------------------------------------------' ) )
	      WRITE(*,'(A82)') TRIM( ADJUSTL( &
				'| TIME | ENERGY-A | ENERGY-B | DECOR | LYP-EXP | DISSIPATION | INCOMPRESSIBILITY |'  ) )
	      WRITE(*,'(A83)') &
				TRIM( ADJUSTL( '-----------------------------------------------------------------------------------' ) )
	    END IF
	    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	    WRITE(*,'(F6.3,A3,F6.4,A4,F6.4,A3,E8.2,A3,F8.4,A5,F8.4,A5,E10.2)')&
			 t_now,'   ',eng,'    ',eng_B,'   ',dec,'   ',lyp_dec,'     ',dis,'     ',inc
	    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		END IF

  END

  SUBROUTINE print_error_timestep()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! print error if time step chosen is large
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,'(A60)')       TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
    WRITE(*,'(A40)')       TRIM( ADJUSTL( 'ERROR: TIME STEP TOO LARGE') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '-----------------------------------------------------------------') )
    WRITE(*,'(A40,F10.6)') TRIM( ADJUSTL( ' RESET THE TIME STEP (AT MAX) AS :') ),dt_max
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '------------------SIMULATION----ABORTED------------------') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '_________________________________________________________________') )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE print_error_nan()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! print error when NaN is encountered
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,'(A60)')       TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
    WRITE(*,'(A40,F8.4)')  TRIM( ADJUSTL( 'ERROR: NAN ENCOUNTERED BEFORE T = ') ), t_now
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '-----------------------------------------------------------------') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '------------------SIMULATION----ABORTED------------------') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '_________________________________________________________________') )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE print_error_inc()
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! print error when NaN is encountered
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE(*,'(A60)')       TRIM( ADJUSTL( 'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT') )
    WRITE(*,'(A40,F8.4)')  TRIM( ADJUSTL( 'ERROR: INCOMPRESSIBILITY LOST BEFORE T = ') ), t_now
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '-----------------------------------------------------------------') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '------------------SIMULATION----ABORTED------------------') )
    WRITE(*,'(A60)')       TRIM( ADJUSTL( '_________________________________________________________________') )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

END MODULE system_basicoutput
