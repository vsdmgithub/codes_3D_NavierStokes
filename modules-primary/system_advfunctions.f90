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
! LAST MODIFIED: 21 JUNE 2021
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
  USE matrix_eigenvalues_eigenvectors

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE compute_dissipation_field
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    ! CALL fft_c2r( -k_2 * v_x, -k_2 * v_y, -k_2 * v_z, N, Nh, w_ux, w_uy, w_uz )

    ! ds_rate = w_ux**two + w_uy**two + w_uz**two

    ! CALL write_dissipation_field
    ! REF-> <<< system_advoutput >>>

  END

  SUBROUTINE compute_velocity_gradient
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the velocity gradient
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    CALL allocate_velocity_gradient

    CALL fft_c2r( i * k_x * v_x, i * k_y * v_x, i * k_z * v_x, N, Nh, duxx, duxy, duxz )
    CALL fft_c2r( i * k_x * v_y, i * k_y * v_y, i * k_z * v_y, N, Nh, duyx, duyy, duyz )
    CALL fft_c2r( i * k_x * v_z, i * k_y * v_z, i * k_z * v_z, N, Nh, duzx, duzy, duzz )

    CALL compute_invariant_QR

    CALL write_section('q_inv',q_invar(0,:,:))
    ! REF-> <<< system_basicoutput >>>

    CALL write_section('r_inv',r_invar(0,:,:))
    ! REF-> <<< system_basicoutput >>>

    CALL deallocate_velocity_gradient

  END

  SUBROUTINE compute_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the velocity gradient
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  T  R  A  I  N        T  E  N  S  O  R        C  A  L  C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL fft_c2r( i * k_x * v_x, hf * i * ( k_y * v_x + k_x * v_y ), i * k_z * v_z , N, Nh, s_xx, s_xy, s_zz)
    CALL fft_c2r( i * k_y * v_y, hf * i * ( k_y * v_z + k_z * v_y ), hf * i * ( k_x * v_z + k_z * v_x ), N, Nh, s_yy, s_yz, s_zx)

    ! CALL write_section('sec_Szz',s_zz(0,:,:))
    ! REF-> <<< system_basicoutput >>>

    ! CALL compute_eigenvalue_distribution

  END

  SUBROUTINE compute_invariant_QR
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the Q,R invariants of the velocity field.
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    q_invar =    duxx ** two + duyy ** two + duzz ** two + &
         two * ( duxy * duyx + duyz * duzy + duzx * duxz )
    q_invar = - hf * q_invar

    r_invar = duxx ** thr + two * duxx * ( duxy * duyx + duxz * duzx ) + &
                                  duxy * ( duyx * duyy + duyz * duzx ) + &
                                  duxz * ( duyx * duzy + duzx * duzz ) + &
              duyy ** thr + two * duyy * ( duxy * duyx + duyz * duzy ) + &
                                  duyz * ( duzx * duxy + duzy * duzz ) + &
                                  duyx * ( duzy * duxz + duxy * duxx ) + &
              duzz ** thr + two * duzz * ( duxz * duzx + duyz * duzy ) + &
                                  duzx * ( duxy * duyz + duxx * duxz ) + &
                                  duzy * ( duzx * duyz + duyy * duyz )
    r_invar = - ( 1.0D0/ thr ) * r_invar

    CALL compute_QR_pdf

  END

  SUBROUTINE compute_QR_pdf
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the histogram of Q-R plane
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    FIND_Q_MAX: IF( ABS(MINVAL(q_invar)) .GT. ABS(MAXVAL(q_invar))) THEN
      q_max = ABS(MINVAL(q_invar))
    ELSE
      q_max = ABS(MAXVAL(q_invar))
    END IF FIND_Q_MAX

    FIND_R_MAX: IF( ABS(MINVAL(r_invar)) .GT. ABS(MAXVAL(r_invar))) THEN
      r_max = ABS(MINVAL(r_invar))
    ELSE
      r_max = ABS(MAXVAL(r_invar))
    END IF FIND_R_MAX

    q_bins = 40
    r_bins = 40
    ! No of bins on either side

    ALLOCATE(pdf_QR(q_bins,r_bins))
    ALLOCATE(q_val(q_bins))
    ALLOCATE(r_val(r_bins))

    DO q_b = 1, q_bins
      q_val( q_b ) = -q_max + DBLE( q_b ) * ( two * q_max ) / q_bins
    END DO

    DO r_b = 1, r_bins
      r_val( r_b ) = -r_max + DBLE( r_b ) * ( two * r_max ) / r_bins
    END DO

    pdf_QR  = 0.0D0
    jump_sz = 4
    pdf_sz  = INT( N / jump_sz ) ** 3

    DO i_x = 0 , N - 1, jump_sz
    DO i_y = 0 , N - 1, jump_sz
    DO i_z = 0 , N - 1, jump_sz

      q_b = FLOOR( q_bins * (q_invar( i_x, i_y, i_z ) + q_max ) / ( two * q_max ) )
      r_b = FLOOR( r_bins * (r_invar( i_x, i_y, i_z ) + r_max ) / ( two * r_max ) )

      pdf_QR(q_b,r_b)  = pdf_QR(q_b,r_b) + 1.0D0
      ! Adding to the histogram frequency

    END DO
    END DO
    END DO

    ! pdf_QR = pdf_QR / pdf_sz
    ! Getting the discrete bin pdf.

    CALL write_QR_bins
    ! REF <<< system_advoutput >>>

    CALL write_QR_pdf
    ! REF <<< system_advoutput >>>

    DEALLOCATE(pdf_QR)
    DEALLOCATE(q_val)
    DEALLOCATE(r_val)

  END

  SUBROUTINE compute_eigenvalue_distribution
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the distribution of eigenvalues of the strain tensor
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    INTEGER(KIND =4)               ::IT_NUM,ROT_NUM
    INTEGER(KIND =4)               ::i_pdf
    DOUBLE PRECISION,DIMENSION(3,3)::str_mx,eig_vec
    DOUBLE PRECISION,DIMENSION(3)  ::eig_val

    jump_sz = 4
    pdf_sz  = INT( N / jump_sz ) ** 3
    i_pdf   = 0

    ALLOCATE( ev_avg( pdf_sz ) )
    ALLOCATE( ev_dif( pdf_sz ) )

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   PDF angles
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO i_x = 0 , N - 1 , jump_sz
    DO i_y = 0 , N - 1 , jump_sz
    DO i_z = 0 , N - 1 , jump_sz
    ! Not going over all grid points, jumping over 'jump_sz'. Since we are interested only in the distribution.

      i_pdf = i_pdf + 1

      str_mx(1,1)=s_xx(i_x,i_y,i_z)
      str_mx(1,2)=s_xy(i_x,i_y,i_z)
      str_mx(1,3)=s_zx(i_x,i_y,i_z)
      str_mx(2,2)=s_yy(i_x,i_y,i_z)
      str_mx(2,3)=s_yz(i_x,i_y,i_z)
      str_mx(3,3)=s_zz(i_x,i_y,i_z)
      str_mx(3,1)=s_zx(i_x,i_y,i_z)
      str_mx(2,1)=s_xy(i_x,i_y,i_z)
      str_mx(3,2)=s_yz(i_x,i_y,i_z)

      ! trace=trace+DABS(str_mx(1,1)+str_mx(2,2)+str_mx(3,3))

      CALL jacobi_eigenvalue(3,str_mx,4,eig_vec,eig_val,IT_NUM,ROT_NUM)
      ! This finds the eigenvalues and eigenvectors (unit normalized) for the strain tensor

      ev_avg( i_pdf ) = DABS( eig_val(3) - eig_val(1) )
      ev_dif( i_pdf ) = eig_val(2)

     END DO
     END DO
     END DO

     CALL compute_eigenvalue_pdf

     DEALLOCATE(ev_avg,ev_dif)

  END

  SUBROUTINE compute_eigenvalue_pdf
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the histogram of eigenvalues in 2 dim
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    INTEGER(KIND =4) :: i_pdf
    DOUBLE PRECISION :: ev_avg_max,ev_dif_max
    DOUBLE PRECISION :: avg_bin_size,dif_bin_size

    ev_avg_max = ABS(MAXVAL( ev_avg ))

    FIND_EV_DIF_MAX: IF( ABS(MINVAL(ev_dif)) .GT. ABS(MAXVAL(ev_dif))) THEN
      ev_dif_max = ABS(MINVAL(ev_dif))
    ELSE
      ev_dif_max = ABS(MAXVAL(ev_dif))
    END IF FIND_EV_DIF_MAX

    avg_bins = 60
    dif_bins = 20
    ! No of bins on either side

    avg_bin_size = one * ev_avg_max / avg_bins
    dif_bin_size = two * ev_dif_max / dif_bins

    ALLOCATE(pdf_ev(avg_bins,dif_bins))
    ALLOCATE(ev_avg_val(avg_bins))
    ALLOCATE(ev_dif_val(dif_bins))

    DO avg_b = 1, avg_bins
      ev_avg_val( avg_b ) = DBLE( avg_b ) * avg_bin_size
    END DO

    DO dif_b = 1, dif_bins
      ev_dif_val( dif_b ) = -ev_dif_max + DBLE( dif_b ) * dif_bin_size
    END DO

    pdf_ev = 0.0D0

    LOOP_BINS: DO i_pdf = 1 , pdf_sz
      avg_b               = FLOOR( ev_avg( i_pdf ) / avg_bin_size )
      dif_b               = FLOOR( ( ev_dif( i_pdf ) + ev_dif_max ) / dif_bin_size )
      pdf_ev(avg_b,dif_b) = pdf_ev(avg_b,dif_b) + 1.0D0
      ! Adding to the histogram frequency

    END DO LOOP_BINS

    ! pdf_ev = pdf_ev / pdf_sz
    ! Getting the discrete bin pdf.

    CALL write_ev_bins
    ! REF <<< system_advoutput >>>

    CALL write_ev_pdf
    ! REF <<< system_advoutput >>>

    DEALLOCATE(pdf_ev)
    DEALLOCATE(ev_avg_val)
    DEALLOCATE(ev_dif_val)

  END

END MODULE system_advfunctions
