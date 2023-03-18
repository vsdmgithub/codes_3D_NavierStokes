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
! MODULE: system_advvariables
! LAST MODIFIED: 20 FEBRAURY 2023
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED VARIABLES TO DO ANALYSIS IN 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advdeclaration
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major advanced functions involving analysis.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_advvariables
  USE system_basicfunctions

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE allocate_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate velocity gradient
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(S_xx(0:N-1,0:N-1,0:N-1),S_yy(0:N-1,0:N-1,0:N-1),S_zz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(S_zx(0:N-1,0:N-1,0:N-1),S_xy(0:N-1,0:N-1,0:N-1),S_yz(0:N-1,0:N-1,0:N-1))
    ALLOCATE(Dis_fld(0:N-1,0:N-1,0:N-1))

    num_bin_dis    = CEILING( ( DBLE(N) / 128.0D0 ) * 100.0D0 )
    print*,num_bin_dis
    ALLOCATE( Dis_val( num_bin_dis ) )
    ALLOCATE( Dis_pdf( num_bin_dis ) )

  END

  SUBROUTINE deallocate_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to deallocate strain tensor
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  D E - A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE(S_xx,S_yy,S_zz)
    DEALLOCATE(S_xy,S_yz,S_zx)
    DEALLOCATE(Dis_fld)
    DEALLOCATE(Dis_val)
    DEALLOCATE(Dis_pdf)

  END

END MODULE system_advdeclaration
