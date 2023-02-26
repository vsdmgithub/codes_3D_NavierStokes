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
! MODULE: system_advfunctions
! LAST MODIFIED: 20 FEBRAURY 2023
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED FUNCTIONS MODULE FOR 3D EULER ANALYSIS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advfunctions
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major advanced functions involving analysis.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_basicoutput
  USE system_advoutput
  ! USE matrix_eigenvalues_eigenvectors

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE compute_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the velocity gradient
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! CALL allocate_strain_tensor
    ! REF-> <<< system_advdeclaration >>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  T  R  A  I  N        T  E  N  S  O  R        C  A  L  C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL fft_c2r( i*K_x*V_x, i*K_y*V_y,  i*K_z*V_z , N, Nh, s_xx, s_yy, s_zz)
    CALL fft_c2r( hf*i*( K_y*V_x + K_x*V_y ), hf*i*( K_y*V_z + K_z*V_y ), &
                  hf*i*( K_x*V_z + K_z*V_x ), N, Nh, s_xy, s_yz, s_zx )

    CALL compute_dissipation

    ! CALL write_section('sec_Szz',s_zz(0,:,:))
    ! REF-> <<< system_basicoutput >>>

    ! CALL compute_eigenvalue_distribution

    ! CALL deallocate_strain_tensor
    ! REF-> <<< system_advdeclaration >>>

  END

  SUBROUTINE compute_dissipation
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! -----------------------------------------
    ! DIRECT NOTATION OF DISSIPATION FIELD
    ! -----------------------------------------
    ! CALL fft_c2r( -k_2 * V_x, -k_2 * V_y, -k_2 * V_z, N, Nh, W_x, W_y, W_z )
    ! Dis_fld = vis * ( W_x * U_x + W_y * U_y + W_z * U_z )

    ! -----------------------------------------
    ! POSITIVE DEF NOTATION OF DISSIPATION FIELD
    ! -----------------------------------------
    Dis_fld = s_xx ** two + s_yy ** two + s_zz ** two + two * ( s_xy ** two + s_yz ** two + s_zx ** two )
    Dis_fld = two * vis * Dis_fld

    dis_avg  = SUM( Dis_fld ) / N3
    dis_std  = SUM( Dis_fld ** two ) / N3  - dis_avg ** two

  END

  SUBROUTINE compute_pdf_dissipation
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Find the pdf of dissipation distribution
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4) ::bin,num_pts_pdf
    DOUBLE PRECISION::dis_dum,std_fac

    ! --------------------------------------------------------------------
    ! Following values are for N=128 . Rescale accordingly for different N
    ! Also these are for DIRECT NOTATION of dissipation field
    ! --------------------------------------------------------------------
    ! dis_max = ( N / 128 ) * (+3.0D0)
    ! dis_min = ( N / 128 ) * (-20.0D0)

    std_fac      = 10.0D0
    dis_dum      = dis_avg
    dis_std      = dis_std / ( dis_avg ** two )
    dis_avg      = one
    dis_max      = dis_avg + std_fac * dis_std
    dis_min      = zero
    dis_bin_size = ( dis_max - dis_min ) / num_bin_dis

    DO bin = 1, num_bin_dis
      Dis_val( bin ) = dis_min +  ( DBLE( bin ) - hf ) * dis_bin_size
    END DO

    Dis_pdf     = zero
    num_pts_pdf = 0

    LOOP_RX_702: DO i_x = 0 , N - 1
    LOOP_RY_702: DO i_y = 0 , N - 1
    LOOP_RZ_702: DO i_z = 0 , N - 1

      Dis_fld( i_x, i_y, i_z ) = Dis_fld( i_x, i_y, i_z ) / dis_dum
      bin                      = CEILING( ( Dis_fld( i_x, i_y, i_z ) - dis_min + tol ) / dis_bin_size )
      ! Finding the bin slot

      BIN_CHECK_332: IF ( ( bin .LT. num_bin_dis ) .AND. ( bin .GT. 1 ) ) THEN
        Dis_pdf( bin )         = Dis_pdf( bin ) + one
        num_bin_dis            = num_bin_dis + 1
      ELSEIF (bin .LE. 1 ) THEN
        Dis_pdf( 1 )           = Dis_pdf( 1 ) + one
      ELSEIF (bin .GE. num_bin_dis ) THEN
        Dis_pdf( num_bin_dis ) = Dis_pdf( num_bin_dis ) + one
      END IF BIN_CHECK_332

    END DO LOOP_RZ_702
    END DO LOOP_RY_702
    END DO LOOP_RX_702

    Dis_pdf = Dis_pdf / N3
    ! Normalizing the PDF

    IF  ( DBLE(num_bin_dis) / N3 .LT. 0.95D0 ) THEN
      PRINT*," More than 5 percent of values are falling out of the dissipation histogram bins"
    END IF

    CALL write_pdf_dissipation
    ! REF-> <<< system_advoutput >>>

  END

  ! SUBROUTINE compute_eigenvalue_distribution
  ! ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ! ------------
  ! ! CALL this to compute the distribution of eigenvalues of the strain tensor
  ! ! -------------
  ! ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !
  !   IMPLICIT NONE
  !   ! _________________________
  !   ! LOCAL VARIABLES
  !   ! !!!!!!!!!!!!!!!!!!!!!!!!!
  !   INTEGER(KIND =4)               ::IT_NUM,ROT_NUM
  !   INTEGER(KIND =4)               ::i_pdf
  !   DOUBLE PRECISION,DIMENSION(3,3)::str_mx,eig_vec
  !   DOUBLE PRECISION,DIMENSION(3)  ::eig_val
  !
  !   jump_sz = 4
  !   data_sz = INT( N / jump_sz ) ** 3
  !   i_pdf   = 0
  !
  !   ALLOCATE( eV_mod( data_sz ) )
  !   ALLOCATE( eV_dif( data_sz ) )
  !
  !   !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   !  P  R  I  N   T          O  U  T  P  U  T   -   PDF angles
  !   !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !   LOOP_RX_703: DO i_x = 0 , N - 1, jump_sz
  !   LOOP_RY_703: DO i_y = 0 , N - 1, jump_sz
  !   LOOP_RZ_703: DO i_z = 0 , N - 1, jump_sz
  !   ! Not going over all grid points, jumping over 'jump_sz'. Since we are interested only in the distribution.
  !
  !     i_pdf = i_pdf + 1
  !
  !     str_mx(1,1)=s_xx(i_x,i_y,i_z)
  !     str_mx(1,2)=s_xy(i_x,i_y,i_z)
  !     str_mx(1,3)=s_zx(i_x,i_y,i_z)
  !     str_mx(2,2)=s_yy(i_x,i_y,i_z)
  !     str_mx(2,3)=s_yz(i_x,i_y,i_z)
  !     str_mx(3,3)=s_zz(i_x,i_y,i_z)
  !     str_mx(3,1)=s_zx(i_x,i_y,i_z)
  !     str_mx(2,1)=s_xy(i_x,i_y,i_z)
  !     str_mx(3,2)=s_yz(i_x,i_y,i_z)
  !
  !     ! trace=trace+DABS(str_mx(1,1)+str_mx(2,2)+str_mx(3,3))
  !
  !     CALL jacobi_eigenvalue(3,str_mx,4,eig_vec,eig_val,IT_NUM,ROT_NUM)
  !     ! This finds the eigenvalues and eigenvectors (unit normalized) for the strain tensor
  !
  !     eV_mod( i_pdf ) = DABS( eig_val(3) - eig_val(1) ) / two
  !     eV_dif( i_pdf ) = eig_val(2)
  !
  !    END DO LOOP_RZ_703
  !    END DO LOOP_RY_703
  !    END DO LOOP_RX_703
  !
  !    CALL compute_eigenvalue_joint_pdf
  !
  !    DEALLOCATE(eV_mod,eV_dif)
  !
  ! END

END MODULE system_advfunctions
