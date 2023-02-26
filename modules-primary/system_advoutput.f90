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
! MODULE: system_advoutput
! LAST MODIFIED: 20 FEBRAURY 2023
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED OUTPUT MODULE - RELATED TO ADVANCED FUNCTION MODULE
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advoutput
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all output calls from advanced functions module.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_advdeclaration

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE write_pdf_dissipation
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes the pdf of dissipation field
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    INTEGER(KIND=4)::bin

    WRITE (fil_time,f_d8p4) t_now
    ! Writes 't_now' as a CHARACTER

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // TRIM( ADJUSTL( dir_pdf ) ) // &
                 'dis_pdf_t_' // TRIM( ADJUSTL ( fil_time ) ) // '.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    OPEN( unit = 759, file = fil_name )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO bin = 1, num_bin_dis

      WRITE(759,f_d12p6,ADVANCE ='NO')   Dis_val( bin )
      WRITE(759,f_d32p17,ADVANCE ='YES') Dis_pdf( bin )

    END DO

    CLOSE(759)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

END MODULE system_advoutput
