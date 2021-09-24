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
  USE system_basicfunctions
  USE system_advoutput

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

    CALL deallocate_velocity_gradient

  END

  SUBROUTINE compute_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the velocity gradient
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    CALL allocate_strain_tensor

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  T  R  A  I  N        T  E  N  S  O  R        C  A  L  C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL fft_c2r( i * k_x * v_x, hf * i * ( k_y * v_x + k_x * v_y ), i * k_z * v_z , N, Nh, s_xx, s_xy, s_zz)
    CALL fft_c2r( i * k_y * v_y, hf * i * ( k_y * v_z + k_z * v_y ), hf * i * ( k_x * v_z + k_z * v_x ), N, Nh, s_yy, s_yz, s_zx)

    CALL write_section('sec_Szz',s_zz(0,:,:))
    ! REF-> <<< system_basoutput >>>

    CALL deallocate_strain_tensor

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

    q_bins = 100
    r_bins = 100
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

    pdf_QR = 0.0D0

    DO i_x = 0 , N - 1
    DO i_y = 0 , N - 1
    DO i_z = 0 , N - 1

      q_b = FLOOR( q_bins * (q_invar( i_x, i_y, i_z ) + q_max ) / ( two * q_max ) )
      r_b = FLOOR( r_bins * (r_invar( i_x, i_y, i_z ) + r_max ) / ( two * r_max ) )

      pdf_QR(q_b,r_b)  = pdf_QR(q_b,r_b) + 1.0D0
      ! Adding to the histogram frequency

    END DO
    END DO
    END DO

    pdf_QR = pdf_QR / N3
    ! Getting the discrete bin pdf.

    CALL write_QR_pdf
    ! REF <<< system_advoutput >>>

    DEALLOCATE(pdf_QR)
    DEALLOCATE(q_val)
    DEALLOCATE(r_val)

  END

END MODULE system_advfunctions
