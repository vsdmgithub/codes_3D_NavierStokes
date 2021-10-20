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

  SUBROUTINE compute_dissipation
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    CALL allocate_dissipation
    ! REF-> <<< system_advdeclaration >>>

    ! -----------------------------------------
    ! DIRECT NOTATION OF DISSIPATION FIELD
    ! -----------------------------------------
    CALL fft_c2r( -k_2 * v_x, -k_2 * v_y, -k_2 * v_z, N, Nh, w_ux, w_uy, w_uz )
    ds_rate = viscosity * ( w_ux * u_x + w_uy * u_y + w_uz * u_z )

    ! -----------------------------------------
    ! POSITIVE DEF NOTATION OF DISSIPATION FIELD
    ! -----------------------------------------
    ! ds_rate = s_xx ** two + s_yy ** two + s_zz ** two + two * ( s_xy ** two + s_yz ** two + s_zx ** two )
    ! ds_rate = two * viscosity * ds_rate

    ds_avg  = SUM( ds_rate ) / N3
    ds_std  = DSQRT( SUM( ds_rate ** two ) / N3  - ds_avg ** two )

    ! CALL write_dissipation_field
    ! REF-> <<< system_advoutput >>>

    CALL compute_pdf_dissipation

    CALL deallocate_dissipation
    ! REF-> <<< system_advdeclaration >>>

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
    INTEGER(KIND=4) ::ds_b,bin_count
    DOUBLE PRECISION::ds_max,ds_min,ds_binsize

    ds_max = MAXVAL( ds_rate )
    ds_min = MINVAL( ds_rate )

    print*,"Maximum dissipation is ",ds_max
    print*,"Minimum dissipation is ",ds_min

    ! --------------------------------------------------------------------
    ! Following values are for N=128 . Rescale accordingly for different N
    ! --------------------------------------------------------------------
    !ds_max = +30.0D0

    ds_bins    = ( N / 128 ) * 100
    ds_binsize = ( ds_max - ds_min ) / ds_bins

    ALLOCATE( ds_val( ds_bins ) )
    ALLOCATE( pdf_ds( ds_bins ) )

    DO ds_b = 1, ds_bins

      ds_val( ds_b ) = ds_min + DBLE( ds_b ) * ds_binsize

    END DO

    pdf_ds    = zero
    bin_count = 0

    LOOP_RX_702: DO i_x = 0 , N - 1
    LOOP_RY_702: DO i_y = 0 , N - 1
    LOOP_RZ_702: DO i_z = 0 , N - 1

      ds_b            = CEILING( ( ds_rate( i_x, i_y, i_z ) - ds_min + tol ) / ds_binsize )
      ! Finding the bin slot

      BIN_CHECK: IF( ds_b .LE. ds_bins ) THEN

        pdf_ds( ds_b ) = pdf_ds( ds_b ) + one
        bin_count      = bin_count + 1

      END IF BIN_CHECK

    END DO LOOP_RZ_702
    END DO LOOP_RY_702
    END DO LOOP_RX_702

    pdf_ds = pdf_ds / DBLE( bin_count )
    ! Normalizing the PDF

    CALL write_pdf_dissipation
    ! REF-> <<< system_advoutput >>>

    DEALLOCATE( ds_val )
    DEALLOCATE( pdf_ds )

  END

  SUBROUTINE compute_velocity_gradient
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the velocity gradient
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    CALL allocate_velocity_gradient
    ! REF-> <<< system_advdeclaration >>>

    CALL fft_c2r( i * k_x * v_x, i * k_y * v_x, i * k_z * v_x, N, Nh, duxx, duxy, duxz )
    CALL fft_c2r( i * k_x * v_y, i * k_y * v_y, i * k_z * v_y, N, Nh, duyx, duyy, duyz )
    CALL fft_c2r( i * k_x * v_z, i * k_y * v_z, i * k_z * v_z, N, Nh, duzx, duzy, duzz )

    CALL compute_invariant_qr

    ! CALL write_section('q_inv',q_invar(0,:,:))
    ! REF-> <<< system_basicoutput >>>

    ! CALL write_section('r_inv',r_invar(0,:,:))
    ! REF-> <<< system_basicoutput >>>

    CALL deallocate_velocity_gradient
    ! REF-> <<< system_advdeclaration >>>

  END

  SUBROUTINE compute_strain_tensor
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the velocity gradient
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    CALL allocate_strain_tensor
    ! REF-> <<< system_advdeclaration >>>

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   S  T  R  A  I  N        T  E  N  S  O  R        C  A  L  C.
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL fft_c2r( i * k_x * v_x, hf * i * ( k_y * v_x + k_x * v_y ), i * k_z * v_z , N, Nh, s_xx, s_xy, s_zz)
    CALL fft_c2r( i * k_y * v_y, hf * i * ( k_y * v_z + k_z * v_y ), hf * i * ( k_x * v_z + k_z * v_x ), N, Nh, s_yy, s_yz, s_zx)

    ! CALL write_section('sec_Szz',s_zz(0,:,:))
    ! REF-> <<< system_basicoutput >>>

    CALL compute_eigenvalue_distribution

    CALL deallocate_strain_tensor
    ! REF-> <<< system_advdeclaration >>>

  END

  SUBROUTINE compute_invariant_qr
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
    r_invar = - ( onethird ) * r_invar

    CALL compute_qr_joint_pdf

  END

  SUBROUTINE compute_qr_joint_pdf
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the histogram of Q-R plane
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)  :: q_b, r_b
    INTEGER(KIND=4)  :: bin_count
    DOUBLE PRECISION :: q_max, r_max

    ! FIND_Q_MAX: IF( ABS(MINVAL(q_invar)) .GT. ABS(MAXVAL(q_invar))) THEN
    !   q_max = ABS(MINVAL(q_invar))
    ! ELSE
    !   q_max = ABS(MAXVAL(q_invar))
    ! END IF FIND_Q_MAX
    !
    ! FIND_R_MAX: IF( ABS(MINVAL(r_invar)) .GT. ABS(MAXVAL(r_invar))) THEN
    !   r_max = ABS(MINVAL(r_invar))
    ! ELSE
    !   r_max = ABS(MAXVAL(r_invar))
    ! END IF FIND_R_MAX

    r_max = ( N / 128 ) * 300.0D0
    q_max = ( N / 128 ) * 60.0D0

    q_bins = ( N / 128 ) * 25
    r_bins = ( N / 128 ) * 25
    ! No of bins for each invariant

    ALLOCATE(pdf_qr(q_bins,r_bins))
    ALLOCATE(q_val(q_bins))
    ALLOCATE(r_val(r_bins))

    DO q_b = 1, q_bins
      q_val( q_b ) = -q_max + ( DBLE( q_b ) - hf ) * ( two * q_max ) / q_bins
    END DO

    DO r_b = 1, r_bins
      r_val( r_b ) = -r_max + ( DBLE( r_b ) - hf ) * ( two * r_max ) / r_bins
    END DO

    pdf_qr    = zero
    bin_count = 0

    LOOP_RX_701: DO i_x = 0 , N - 1
    LOOP_RY_701: DO i_y = 0 , N - 1
    LOOP_RZ_701: DO i_z = 0 , N - 1

      q_b = CEILING( q_bins * (q_invar( i_x, i_y, i_z ) + q_max  ) / ( two * q_max ) )
      r_b = CEILING( r_bins * (r_invar( i_x, i_y, i_z ) + r_max  ) / ( two * r_max ) )

      BIN_CHECK_QR: IF ( ( q_b .GE. 1) .AND. ( q_b .LE. q_bins ) .AND. ( r_b .GE. 1) .AND. ( r_b .LE. r_bins ) ) THEN
        pdf_qr(q_b,r_b)  = pdf_qr(q_b,r_b) + 1.0D0
        bin_count        = bin_count + 1
        ! Adding to the histogram frequency
      END IF BIN_CHECK_QR

    END DO LOOP_RZ_701
    END DO LOOP_RY_701
    END DO LOOP_RX_701

    pdf_qr = pdf_qr / DBLE( bin_count )
    ! Getting the discrete bin pdf.

    CALL write_qr_joint_pdf_bins
    ! REF <<< system_advoutput >>>

    CALL write_qr_joint_pdf
    ! REF <<< system_advoutput >>>

    DEALLOCATE(pdf_qr)
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
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND =4)               ::IT_NUM,ROT_NUM
    INTEGER(KIND =4)               ::i_pdf
    DOUBLE PRECISION,DIMENSION(3,3)::str_mx,eig_vec
    DOUBLE PRECISION,DIMENSION(3)  ::eig_val

    jump_sz = 4
    data_sz = INT( N / jump_sz ) ** 3
    i_pdf   = 0

    ALLOCATE( ev_mod( data_sz ) )
    ALLOCATE( ev_dif( data_sz ) )

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   PDF angles
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    LOOP_RX_703: DO i_x = 0 , N - 1, jump_sz
    LOOP_RY_703: DO i_y = 0 , N - 1, jump_sz
    LOOP_RZ_703: DO i_z = 0 , N - 1, jump_sz
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

      ev_mod( i_pdf ) = DABS( eig_val(3) - eig_val(1) ) / two
      ev_dif( i_pdf ) = eig_val(2)

     END DO LOOP_RZ_703
     END DO LOOP_RY_703
     END DO LOOP_RX_703

     CALL compute_eigenvalue_joint_pdf

     DEALLOCATE(ev_mod,ev_dif)

  END

  SUBROUTINE compute_eigenvalue_joint_pdf
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the histogram of eigenvalues in 2 dim
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND =4) :: i_pdf, bin_count
    DOUBLE PRECISION :: ev_mod_max,ev_dif_max
    DOUBLE PRECISION :: ev_mod_bin_size,ev_dif_bin_size
    INTEGER(KIND=4)  :: mod_b,dif_b

    ! ev_mod_max = DABS(MAXVAL( ev_mod ))

    ! FIND_EV_DIF_MAX: IF( ABS(MINVAL(ev_dif)) .GT. ABS(MAXVAL(ev_dif))) THEN
      ! ev_dif_max = DABS(MINVAL(ev_dif))
    ! ELSE
      ! ev_dif_max = DABS(MAXVAL(ev_dif))
    ! END IF FIND_EV_DIF_MAX

    ev_mod_max = ( N / 128 ) * 20.0D0
    ev_dif_max = ( N / 128 ) * 4.0D0

    ev_mod_bins = ( N / 128 ) * 25
    ev_dif_bins = ( N / 128 ) * 25
    ! No of bins on either side

    ev_mod_bin_size = one * ev_mod_max / ev_mod_bins
    ev_dif_bin_size = two * ev_dif_max / ev_dif_bins

    ALLOCATE(pdf_ev(ev_mod_bins,ev_dif_bins))
    ALLOCATE(ev_mod_val(ev_mod_bins))
    ALLOCATE(ev_dif_val(ev_dif_bins))

    DO mod_b = 1, ev_mod_bins
      ev_mod_val( mod_b ) = ( DBLE( mod_b ) - hf ) * ev_mod_bin_size
    END DO

    DO dif_b = 1, ev_dif_bins
      ev_dif_val( dif_b ) = -ev_dif_max + ( DBLE( dif_b ) - hf ) * ev_dif_bin_size
    END DO

    pdf_ev    = zero
    bin_count = 0

    LOOP_BINS: DO i_pdf = 1 , data_sz

      mod_b               = FLOOR(                ev_mod( i_pdf )   / ev_mod_bin_size )
      dif_b               = FLOOR( ( ev_dif_max + ev_dif( i_pdf ) ) / ev_dif_bin_size )

      BIN_CHECK_EV: IF ( ( mod_b .GE. 1) .AND. ( mod_b .LE. ev_mod_bins ) .AND. ( dif_b .GE. 1) &
      .AND. ( dif_b .LE. ev_dif_bins ) ) THEN

        pdf_ev(mod_b,dif_b) = pdf_ev(mod_b,dif_b) + 1.0D0
        bin_count        = bin_count + 1
      ! Adding to the histogram frequency

      END IF BIN_CHECK_EV

    END DO LOOP_BINS

    pdf_ev = pdf_ev / DBLE( bin_count )
    ! Getting the discrete bin pdf.

    CALL write_ev_joint_pdf_bins
    ! REF <<< system_advoutput >>>

    CALL write_ev_joint_pdf
    ! REF <<< system_advoutput >>>

    DEALLOCATE(pdf_ev)
    DEALLOCATE(ev_mod_val)
    DEALLOCATE(ev_dif_val)

  END

END MODULE system_advfunctions
