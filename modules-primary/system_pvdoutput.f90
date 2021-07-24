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
  USE system_basicvariables
  USE system_advvariables
  USE system_basicoutput

  IMPLICIT NONE
  ! _________________________
  ! OUTPUT VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  ! CHARACTER(LEN=40)::pvd_file
  INTEGER(KIND=4)::pvd_N_x,pvd_N_y,pvd_N_z
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::vec_x,vec_y,vec_z,scalr
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE ::pvd_ax_x,pvd_ax_y,pvd_ax_z

  ! TYPE(VTK_file_handle)::fd
  TYPE(VTR_file_handle)::fd
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
    ALLOCATE( pvd_ax_x( 0 : pvd_N_x - 1 ) )
    ALLOCATE( pvd_ax_y( 0 : pvd_N_y - 1 ) )
    ALLOCATE( pvd_ax_z( 0 : pvd_N_z - 1 ) )
    ALLOCATE( vec_x( 0 : pvd_N_x - 1, 0 : pvd_N_y - 1, 0 : pvd_N_z - 1 ) )
    ALLOCATE( vec_y( 0 : pvd_N_x - 1, 0 : pvd_N_y - 1, 0 : pvd_N_z - 1 ) )
    ALLOCATE( vec_z( 0 : pvd_N_x - 1, 0 : pvd_N_y - 1, 0 : pvd_N_z - 1 ) )
    ALLOCATE( scalr( 0 : pvd_N_x - 1, 0 : pvd_N_y - 1, 0 : pvd_N_z - 1 ) )

    DO i_x = 0, pvd_N_x - 1

      pvd_ax_x( i_x ) = i_x * dx

    END DO

    DO i_y = 0, pvd_N_y - 1

      pvd_ax_y( i_y ) = i_y * dx

    END DO

    DO i_z = 0, pvd_N_z - 1

      pvd_ax_z( i_z ) = i_z * dx

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

    ! WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_3D ) ) &
                // 'V_t'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   VORTICITY - PVD FORMAT
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL  VTR_open_file(PREFIX=file_name,FD=fd)

    CALL  VTR_write_mesh(FD=fd,X=axis,Y=axis,Z=axis)

    CALL  VTR_write_var(FD=fd,NAME="Velocity",VX=u_x,VY=u_y,VZ=u_z )

    CALL  VTR_close_file(FD=fd)

    ! CALL  VTR_collect_file( FD = fd )
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

    ! WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_3D ) ) &
                // 'VX_t'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   VORTICITY - PVD FORMAT
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL  VTR_open_file(PREFIX=file_name,FD=fd)

    CALL  VTR_write_mesh(FD=fd,X=axis,Y=axis,Z=axis)

    CALL  VTR_write_var(FD=fd,NAME="Vorticity",VX=w_ux,VY=w_uy,VZ=w_uz )

    CALL  VTR_close_file(FD=fd)

    ! CALL  VTR_collect_file( FD = fd )
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

    vec_x = w_ux(0:pvd_N_x-1,0:pvd_N_y-1,0:pvd_N_z-1)
    vec_y = w_uy(0:pvd_N_x-1,0:pvd_N_y-1,0:pvd_N_z-1)
    vec_z = w_uz(0:pvd_N_x-1,0:pvd_N_y-1,0:pvd_N_z-1)
    ! COPYING THE SUBSET DATA

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_3D ) ) &
                // 'VX_SUB_t'

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   VORTICITY - PVD FORMAT (SUBSET)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL  VTR_open_file(PREFIX=file_name,FD=fd)

    CALL  VTR_write_mesh(FD=fd,X=pvd_ax_x,Y=pvd_ax_y,Z=pvd_ax_z)

    CALL  VTR_write_var(FD=fd,NAME="Vorticity",VX=vec_x,VY=vec_y,VZ=vec_z )

    CALL  VTR_close_file(FD=fd)

    ! CALL  VTR_collect_file( FD = fd )
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
		DEALLOCATE(pvd_ax_x,pvd_ax_y,pvd_ax_z)
		DEALLOCATE(vec_x,vec_y,vec_z)
		DEALLOCATE(scalr)

	END

END MODULE system_pvdoutput
