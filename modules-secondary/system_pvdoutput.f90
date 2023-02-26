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
! MODULE: system_pvdoutput
! LAST MODIFIED: 21 JUNE 2021
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! PVD OUTPUT MODULE, WHERE DATA IS WRITTEN IN '.vtr' FORMAT
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_pvdoutput
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major system_pvdoutput for the Code to run.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_VTR
  USE system_VTK
  USE system_advvariables
  USE system_basicoutput

  IMPLICIT NONE
  ! _________________________
  ! OUTPUT VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE(VTR_fil_handle)::fil_pvd
  TYPE(VTK_fil_handle)::fil_pvd
  ! CHARACTER(LEN=40)::pvd_file
  INTEGER(KIND=4)::pvd_N_x,pvd_N_y,pvd_N_z
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::Vec_x,Vec_y,Vec_z,Scal
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE ::Pvd_x,Pvd_y,Pvd_z

  ! creating dataypes to store paraview files for 3d viewing

  CONTAINS

  SUBROUTINE allocate_PVD_subset_arrays
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate the arrays of partial size from the original size for PVD output
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    pvd_N_x = N / 4
    pvd_N_y = N / 4
    pvd_N_z = N / 2

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE( Pvd_x( 0 : pvd_N_x - 1 ) )
    ALLOCATE( Pvd_y( 0 : pvd_N_y - 1 ) )
    ALLOCATE( Pvd_z( 0 : pvd_N_z - 1 ) )
    ALLOCATE( Vec_x( 0 : pvd_N_x - 1, 0 : pvd_N_y - 1, 0 : pvd_N_z - 1 ) )
    ALLOCATE( Vec_y( 0 : pvd_N_x - 1, 0 : pvd_N_y - 1, 0 : pvd_N_z - 1 ) )
    ALLOCATE( Vec_z( 0 : pvd_N_x - 1, 0 : pvd_N_y - 1, 0 : pvd_N_z - 1 ) )
    ALLOCATE(  Scal( 0 : pvd_N_x - 1, 0 : pvd_N_y - 1, 0 : pvd_N_z - 1 ) )

    DO i_x = 0, pvd_N_x - 1

      Pvd_x( i_x ) = i_x * l_grd

    END DO

    DO i_y = 0, pvd_N_y - 1

      Pvd_y( i_y ) = i_y * l_grd

    END DO

    DO i_z = 0, pvd_N_z - 1

      Pvd_z( i_z ) = i_z * l_grd

    END DO


  END

  SUBROUTINE write_PVD_velocity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in PVD format
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // TRIM( ADJUSTL( dir_sec ) ) &
                // 'Vel'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   VORTICITY - PVD FORMAT
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL  VTR_open_file(PREFIX=fil_name,FD=fil_pvd)

    CALL  VTR_write_mesh(FD=fil_pvd,X=Axis,Y=Axis,Z=Axis)

    CALL  VTR_write_var(FD=fil_pvd,NAME="Velocity",VX=U_x,VY=U_y,VZ=U_z )

    CALL  VTR_close_file(FD=fil_pvd)

    ! CALL  VTR_collect_file( FD = fil_pvd )
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_PVD_vorticity
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in PVD format
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    ! WRITE (fil_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // TRIM( ADJUSTL( dir_sec ) ) &
                // 'Vor'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   VORTICITY - PVD FORMAT
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL  VTR_open_file(PREFIX=fil_name,FD=fil_pvd)

    CALL  VTR_write_mesh(FD=fil_pvd,X=Axis,Y=Axis,Z=Axis)

    CALL  VTR_write_var(FD=fil_pvd,NAME="Vorticity",VX=W_x,VY=W_y,VZ=W_z )

    CALL  VTR_close_file(FD=fil_pvd)

    ! CALL  VTR_collect_file( FD = fil_pvd )
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_PVD_vorticity_subset
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in PVD format (only a part of the data to save memory )
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    Vec_x = W_x(0:pvd_N_x-1,0:pvd_N_y-1,0:pvd_N_z-1)
    Vec_y = W_y(0:pvd_N_x-1,0:pvd_N_y-1,0:pvd_N_z-1)
    Vec_z = W_z(0:pvd_N_x-1,0:pvd_N_y-1,0:pvd_N_z-1)
    ! COPYING THE SUBSET DATA

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // TRIM( ADJUSTL( dir_sec ) ) &
                // 'Vor_sub'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   VORTICITY - PVD FORMAT (SUBSET)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL  VTR_open_file(PREFIX=fil_name,FD=fil_pvd)

    CALL  VTR_write_mesh(FD=fil_pvd,X=Pvd_x,Y=Pvd_y,Z=Pvd_z)

    CALL  VTR_write_var(FD=fil_pvd,NAME="Vorticity",VX=Vec_x,VY=Vec_y,VZ=Vec_z )

    CALL  VTR_close_file(FD=fil_pvd)

    ! CALL  VTR_collect_file( FD = fil_pvd )
    ! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

	SUBROUTINE deallocate_PVD_subset_arrays
	! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	! ------------
	! CALL this to deallocate arrays
	! -------------
	! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

		IMPLICIT NONE
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		!  D  E  A  L  L  O  C  A  T  I  O  N
		!  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		DEALLOCATE(Pvd_x,Pvd_y,Pvd_z)
		DEALLOCATE(Vec_x,Vec_y,Vec_z)
		DEALLOCATE(Scal)

	END

END MODULE system_pvdoutput
