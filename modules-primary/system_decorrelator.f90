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
! MODULE: system_decorrelator
! LAST MODIFIED: 24 SEPT 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED FUNCTIONS MODULE FOR 3D EULER ANALYSIS
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_decorrelator
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major analysis involving decorrelator
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_advfunctions

  IMPLICIT NONE
  DOUBLE PRECISION::decor,decor_initial,decor_final
  DOUBLE PRECISION::lyp_max,lyp_min,lyp_binsize
  DOUBLE PRECISION::lyp_str_avg,lyp_str_std,lyp_str_bar
  DOUBLE PRECISION::lyp_eta_avg,lyp_eta_std,lyp_eta_bar
  DOUBLE PRECISION::ccrel_str_eta,ccrel_ds_lyp_eta
  DOUBLE PRECISION::ccrel_ds_lyp_str,ccrel_lyp_str_vx_alp,ccrel_ds_vx_alp
  DOUBLE PRECISION::decor_old,lyp_decor
  INTEGER(KIND=4) ::lyp_bins
  ! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::diff_field
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::lyp_str,lyp_eta
  DOUBLE PRECISION,DIMENSION(:)    ,ALLOCATABLE ::lyp_val,pdf_lyp

  CONTAINS

  SUBROUTINE allocate_decorrelator
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate decorrelator
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    lyp_bins = ( N / 128 ) * 100
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(diff_field(0:N-1,0:N-1,0:N-1))
    ALLOCATE(lyp_str(0:N-1,0:N-1,0:N-1))
    ALLOCATE(lyp_eta(0:N-1,0:N-1,0:N-1))
    ALLOCATE(lyp_val(lyp_bins),pdf_lyp(lyp_bins))

  END

  SUBROUTINE compute_decorrelator
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the decorrelators
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  D E C O R R E L A T O R
    !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    diff_field = hf * ( ( u2_x - u_x ) ** two + ( u2_y - u_y ) ** two + ( u2_z - u_z ) ** two )
    ! Matrix of decorelation

    decor      = SUM( diff_field ) / N3
    ! Global decorrelation

    lyp_decor  = ( DLOG( decor ) - DLOG( decor_old ) ) / dt
    ! Finding the exponent at which it grows at 't'

    decor_old  = decor
    ! copying the old decorrelator

    ! CALL write_section('decor',diff_field( 0, : , : ) )
    ! REF <<< system_basicoutput >>>

  END

  SUBROUTINE compute_lyapunov_S
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the decorrelators
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    lyp_str   = -(  ( ( u_x - u2_x )**two ) * s_xx + ( ( u_y - u2_y )**two ) * s_yy + ( ( u_z - u2_z )**two ) * s_zz + &
              two * ( ( u_x - u2_x ) * s_xy * ( u_y - u2_y ) + &
                      ( u_y - u2_y ) * s_yz * ( u_z - u2_z ) + &
                      ( u_z - u2_z ) * s_zx * ( u_x - u2_x ) ) )

    lyp_str_bar = SUM( lyp_str ) / ( N3 * decor )

  END

  SUBROUTINE compute_lyapunov_eta
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the decorrelators
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::ip_x,ip_y,ip_z
    INTEGER(KIND=4)::im_x,im_y,im_z
    DOUBLE PRECISION::df_x, df_y, df_z
    DOUBLE PRECISION::df_x_px, df_y_px, df_z_px
    DOUBLE PRECISION::df_x_py, df_y_py, df_z_py
    DOUBLE PRECISION::df_x_pz, df_y_pz, df_z_pz
    DOUBLE PRECISION::df_x_mx, df_y_mx, df_z_mx
    DOUBLE PRECISION::df_x_my, df_y_my, df_z_my
    DOUBLE PRECISION::df_x_mz, df_y_mz, df_z_mz
    DOUBLE PRECISION::decor_eta_val

    DO i_x = 0 , N - 1
    DO i_y = 0 , N - 1
    DO i_z = 0 , N - 1

      ip_x = i_x + 1
      ip_y = i_y + 1
      ip_z = i_z + 1
      im_x = i_x - 1
      im_y = i_y - 1
      im_z = i_z - 1

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! IMPOSING PERIODICITY DURING FINITE DIFFERENCE FOR LAPLACIAN
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      PBC_PX: IF ( ip_x > N - 1 ) THEN
        ip_x  = ip_x - N
      END IF PBC_PX
      PBC_PY: IF ( ip_y > N - 1 ) THEN
        ip_y  = ip_y - N
      END IF PBC_PY
      PBC_PZ: IF ( ip_z > N - 1 ) THEN
        ip_z  = ip_z - N
      END IF PBC_PZ
      PBC_MX: IF ( im_x < 0 ) THEN
        im_x  = im_x + N
      END IF PBC_MX
      PBC_MY: IF ( im_y < 0 ) THEN
        im_y  = im_y + N
      END IF PBC_MY
      PBC_MZ: IF ( im_z < 0 ) THEN
        im_z  = im_z + N
      END IF PBC_MZ

      df_x    = u2_x( i_x, i_y, i_z ) - u_x( i_x, i_y, i_z )
      df_y    = u2_y( i_x, i_y, i_z ) - u_y( i_x, i_y, i_z )
      df_z    = u2_z( i_x, i_y, i_z ) - u_z( i_x, i_y, i_z )

      df_x_px = u2_x( ip_x, i_y, i_z ) - u_x( ip_x, i_y, i_z )
      df_y_px = u2_y( ip_x, i_y, i_z ) - u_y( ip_x, i_y, i_z )
      df_z_px = u2_z( ip_x, i_y, i_z ) - u_z( ip_x, i_y, i_z )
      df_x_py = u2_x( i_x, ip_y, i_z ) - u_x( i_x, ip_y, i_z )
      df_y_py = u2_y( i_x, ip_y, i_z ) - u_y( i_x, ip_y, i_z )
      df_z_py = u2_z( i_x, ip_y, i_z ) - u_z( i_x, ip_y, i_z )
      df_x_pz = u2_x( i_x, i_y, ip_z ) - u_x( i_x, i_y, ip_z )
      df_y_pz = u2_y( i_x, i_y, ip_z ) - u_y( i_x, i_y, ip_z )
      df_z_pz = u2_z( i_x, i_y, ip_z ) - u_z( i_x, i_y, ip_z )

      df_x_mx = u2_x( im_x, i_y, i_z ) - u_x( im_x, i_y, i_z )
      df_y_mx = u2_y( im_x, i_y, i_z ) - u_y( im_x, i_y, i_z )
      df_z_mx = u2_z( im_x, i_y, i_z ) - u_z( im_x, i_y, i_z )
      df_x_my = u2_x( i_x, im_y, i_z ) - u_x( i_x, im_y, i_z )
      df_y_my = u2_y( i_x, im_y, i_z ) - u_y( i_x, im_y, i_z )
      df_z_my = u2_z( i_x, im_y, i_z ) - u_z( i_x, im_y, i_z )
      df_x_mz = u2_x( i_x, i_y, im_z ) - u_x( i_x, i_y, im_z )
      df_y_mz = u2_y( i_x, i_y, im_z ) - u_y( i_x, i_y, im_z )
      df_z_mz = u2_z( i_x, i_y, im_z ) - u_z( i_x, i_y, im_z )

      decor_eta_val = df_x * ( - 6.0D0 * df_x + df_x_px + df_x_py + df_x_pz + df_x_mx + df_x_my + df_x_mz ) + &
                      df_y * ( - 6.0D0 * df_y + df_y_px + df_y_py + df_y_pz + df_y_mx + df_y_my + df_y_mz ) + &
                      df_z * ( - 6.0D0 * df_z + df_z_px + df_z_py + df_z_pz + df_z_mx + df_z_my + df_z_mz )

      lyp_eta( i_x, i_y, i_z ) = viscosity * decor_eta_val / ( dx * dx )

    END DO
    END DO
    END DO

    lyp_eta_bar = SUM( lyp_eta ) / ( N3 * decor )

  END

  SUBROUTINE deallocate_decorrelator
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to deallocate decorrelator
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  D E - A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DEALLOCATE(diff_field)
    DEALLOCATE(lyp_str)
    DEALLOCATE(lyp_eta)
    ! DEALLOCATE(lyp_val, pdf_lyp)

  END

  SUBROUTINE write_decorrelator_growth
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in time, for every timestep
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   D E C O R    V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      file_name = TRIM( ADJUSTL( file_address ) ) // 'decor_growth.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 4002, file = file_name )
      ! File where energy vs time will be written. With additional data
    END IF

    WRITE(4002,f_d8p4,ADVANCE   ='no') time_now
    WRITE(4002,f_d32p17,ADVANCE ='no') decor
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_decor
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_eta_bar
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_str_bar
    WRITE(4002,f_d32p17,ADVANCE ='yes')lyp_str_bar + lyp_eta_bar

    IF ( t_step .EQ. t_step_total ) THEN
      CLOSE(4002)
    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_lyapunov
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in time, for every timestep
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    file_name = TRIM( ADJUSTL( file_address ) ) // 'lyapunov_t_' &
             // TRIM( ADJUSTL ( file_time ) ) // '.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  Lyapunov Scale
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 4008, file = file_name )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    DO i_x = 0 , N - 1
    DO i_y = 0 , N - 1
    DO i_z = 0 , N - 1

      DIFF_FIELD_CHECK: IF ( diff_field( i_x, i_y, i_z ) .GT. tol ) THEN
        WRITE(4008,f_d32p17,ADVANCE ='no')  lyp_str( i_x, i_y, i_z ) / diff_field( i_x, i_y, i_z )
        WRITE(4008,f_d32p17,ADVANCE ='yes') lyp_eta( i_x, i_y, i_z ) / diff_field( i_x, i_y, i_z )
      END IF DIFF_FIELD_CHECK

    END DO
    END DO
    END DO

    CLOSE(4008)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE compute_pdf_lyapunov_S
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Find the pdf of lyp-S scales distribution
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::l_b,bin_count

    lyp_str     = ( lyp_str ) / ( diff_field )

    lyp_str_avg = SUM( lyp_str ) / N3
    lyp_str_std = DSQRT( SUM( lyp_str ** two ) / N3 - lyp_str_avg ** two )

    ! lyp_max = MAXVAL( lyp_str )
    ! lyp_min = MINVAL( lyp_str )

    ! --------------------------------------------------------------------
    ! Following values are for N=128 . Rescale accordingly for different N
    ! --------------------------------------------------------------------
    !lyp_max = +40.0D0
    !lyp_min = -40.0D0

    lyp_max     = +40.0D0 * ( N / 128 )
    lyp_min     = -40.0D0 * ( N / 128 )

    lyp_binsize = ( lyp_max - lyp_min ) / lyp_bins

    DO l_b = 1, lyp_bins

      lyp_val( l_b ) = lyp_min + DBLE( l_b ) * lyp_binsize

    END DO

    pdf_lyp   = zero
    bin_count = 0

    LOOP_RX_901: DO i_x = 0 , N - 1
    LOOP_RY_901: DO i_y = 0 , N - 1
    LOOP_RZ_901: DO i_z = 0 , N - 1

      l_b            = CEILING( ( lyp_str( i_x, i_y, i_z ) - lyp_min ) / lyp_binsize )
      ! Finding the bin slot

      BIN_CHECK_901: IF( l_b .LE. lyp_bins ) THEN

        pdf_lyp( l_b ) = pdf_lyp( l_b ) + one
        bin_count      = bin_count + 1

      END IF BIN_CHECK_901

    END DO LOOP_RZ_901
    END DO LOOP_RY_901
    END DO LOOP_RX_901

    pdf_lyp = pdf_lyp / DBLE( bin_count )
    ! Normalizing the PDF

    CALL write_pdf_lyapunov_S

  END

  SUBROUTINE compute_pdf_lyapunov_eta
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Find the pdf of lyp-eta scales distribution
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4) ::l_b,bin_count,bin_count_2

    !lyp_eta = lyp_eta / ( diff_field + tol )

    ! lyp_max = MAXVAL( lyp_eta )
    ! lyp_min = MINVAL( lyp_eta )

    ! --------------------------------------------------------------------
    ! Following values are for N=128 . Rescale accordingly for different N
    ! --------------------------------------------------------------------
    !lyp_max   = +75.0D0
    !lyp_min   = -75.0D0

    lyp_max     = +75.0D0 * ( N / 128 )
    lyp_min     = -75.0D0 * ( N / 128 )

    lyp_binsize = ( lyp_max - lyp_min ) / lyp_bins

    DO l_b = 1, lyp_bins

      lyp_val( l_b ) = lyp_min + DBLE( l_b ) * lyp_binsize

    END DO

    pdf_lyp     = zero
    lyp_eta_std = zero
    lyp_eta_avg = zero
    bin_count   = 0
    bin_count_2 = 0

    LOOP_RX_902: DO i_x = 0 , N - 1
    LOOP_RY_902: DO i_y = 0 , N - 1
    LOOP_RZ_902: DO i_z = 0 , N - 1

      DIFF_FIELD_CHECK_2: IF( diff_field( i_x, i_y, i_z ) .GT. tol ) THEN

        lyp_eta( i_x, i_y, i_z ) = lyp_eta( i_x, i_y, i_z ) / diff_field( i_x, i_y, i_z )

        l_b                      = CEILING( ( lyp_eta( i_x, i_y, i_z ) - lyp_min ) / lyp_binsize )
        ! Finding the bin slot

        lyp_eta_std              = lyp_eta_std + lyp_eta( i_x, i_y, i_z ) ** two
        lyp_eta_avg              = lyp_eta_bar + lyp_eta( i_x, i_y, i_z )
        bin_count_2              = bin_count_2 + 1

        BIN_CHECK: IF( l_b .LE. lyp_bins ) THEN

          pdf_lyp( l_b ) = pdf_lyp( l_b ) + one
          bin_count      = bin_count + 1

        END IF BIN_CHECK

      END IF DIFF_FIELD_CHECK_2

    END DO LOOP_RZ_902
    END DO LOOP_RY_902
    END DO LOOP_RX_902

    pdf_lyp     = pdf_lyp / DBLE( bin_count )

    lyp_eta_std = lyp_eta_std / bin_count_2
    lyp_eta_avg = lyp_eta_avg / bin_count_2

    lyp_eta_std = DSQRT( lyp_eta_std - lyp_eta_avg ** two )

    CALL write_pdf_lyapunov_eta
    ! Normalizing the PDF

  END

  SUBROUTINE write_pdf_lyapunov_S
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes the pdf of lyapunov
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::l_b

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_pdf ) ) // 'lyp_str_pdf_' &
              //TRIM( ADJUSTL( N_char ) ) // '_t_' // TRIM( ADJUSTL ( file_time ) ) // '.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    OPEN( unit = 799, file = file_name )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO l_b = 1, lyp_bins

      WRITE(799,f_d12p6,ADVANCE ='NO') lyp_val( l_b )
      WRITE(799,f_d32p17,ADVANCE ='YES') pdf_lyp( l_b )

    END DO

    CLOSE(799)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_pdf_lyapunov_eta
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes the pdf of lyapunov
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::l_b

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_pdf ) ) // 'lyp_eta_pdf' &
              //TRIM( ADJUSTL( N_char ) ) // '_t_' // TRIM( ADJUSTL ( file_time ) ) // '.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    OPEN( unit = 796, file = file_name )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO l_b = 1, lyp_bins

      WRITE(796,f_d12p6,ADVANCE ='NO') lyp_val( l_b )
      WRITE(796,f_d32p17,ADVANCE ='YES') pdf_lyp( l_b )

    END DO

    CLOSE(796)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE compute_cross_correlation
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! To find the cross correlation coefficient between dissipation field and strain lyapunov
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    ccrel_ds_lyp_str = SUM( ds_rate * lyp_str ) / N3
    ccrel_ds_lyp_str = ( ccrel_ds_lyp_str - ds_avg * lyp_str_avg ) / ( ds_std * lyp_str_std )

    ccrel_ds_lyp_eta = SUM( ds_rate * lyp_eta ) / N3
    ccrel_ds_lyp_eta = ( ccrel_ds_lyp_eta - ds_avg * lyp_eta_avg ) / ( ds_std * lyp_eta_std )

    ccrel_str_eta    = SUM( lyp_str * lyp_eta ) / N3
    ccrel_str_eta    = ( ccrel_str_eta - lyp_str_avg * lyp_eta_avg ) / ( lyp_eta_std * lyp_str_std )

    ccrel_lyp_str_vx_alp = SUM( vx_alp * lyp_str ) / N3
    ccrel_lyp_str_vx_alp = ( ccrel_lyp_str_vx_alp - lyp_str_avg * vx_alp_avg ) / ( lyp_str_std * vx_alp_std )

    ccrel_ds_vx_alp  = SUM( ds_rate * vx_alp ) / N3
    ccrel_ds_vx_alp  = ( ccrel_ds_vx_alp - ds_avg * vx_alp_avg ) / ( ds_std * vx_alp_std )

    ! CALL write_extreme_events

    CALL write_cross_correlation

  END

  SUBROUTINE write_extreme_events
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! To find the extreme dissipation points and the lyapunov timescales there
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::event_count
    DOUBLE PRECISION::ds_ratio
    DOUBLE PRECISION::extremity_threshold

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  E X T R E M E    E V E N T    R E C O R D
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    extremity_threshold = 6.0D0
    event_count         = 0
    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // TRIM( ADJUSTL( sub_dir_pdf ) ) // 'extreme_events' &
              //TRIM( ADJUSTL( N_char ) ) // '_t_' // TRIM( ADJUSTL ( file_time ) ) // '.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    OPEN( unit = 866, file = file_name )

    LOOP_RX_908: DO i_x = 0 , N - 1
    LOOP_RY_908: DO i_y = 0 , N - 1
    LOOP_RZ_908: DO i_z = 0 , N - 1

      ds_ratio            = ds_rate( i_x, i_y, i_z ) / ds_std

      EXTREMITY_CHECK: IF( DABS( ds_ratio ) .GT. extremity_threshold ) THEN

        event_count = event_count + 1
        WRITE(866,f_d32p17,ADVANCE ='NO') ds_ratio
        WRITE(866,f_d32p17,ADVANCE ='NO') lyp_str( i_x, i_y, i_z ) / lyp_str_bar
        WRITE(866,f_d32p17,ADVANCE ='NO') lyp_eta( i_x, i_y, i_z ) / lyp_eta_bar
        WRITE(866,f_d32p17,ADVANCE ='YES') vx_alp( i_x, i_y, i_z ) / vx_alp_avg

      END IF EXTREMITY_CHECK

    END DO LOOP_RZ_908
    END DO LOOP_RY_908
    END DO LOOP_RX_908

    CLOSE(866)

  END

  SUBROUTINE write_cross_correlation
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the cross correlation data in time
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !   C R O S S   C O R R E L A T I O N    V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      file_name = TRIM( ADJUSTL( file_address ) ) // 'cross_correlation.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 4008, file = file_name )
      ! File where energy vs time will be written. With additional data
    END IF

    IF ( t_step .GT. 0 ) THEN
      WRITE(4008,f_d8p4,ADVANCE   ='no') time_now
      WRITE(4008,f_d32p17,ADVANCE ='no') ccrel_ds_lyp_str
      WRITE(4008,f_d32p17,ADVANCE ='no') ccrel_ds_lyp_eta
      WRITE(4008,f_d32p17,ADVANCE ='no') ccrel_lyp_str_vx_alp
      WRITE(4008,f_d32p17,ADVANCE ='no') ccrel_ds_vx_alp
      WRITE(4008,f_d32p17,ADVANCE ='yes') ccrel_str_eta

    END IF

    IF ( t_step .EQ. t_step_total ) THEN
      CLOSE(4008)
    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++
    !  S T A T I S T I C S
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IF ( t_step .EQ. 0 ) THEN
      file_name = TRIM( ADJUSTL( file_address ) ) // 'statistics.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 4808, file = file_name )
      ! File where energy vs time will be written. With additional data
    END IF

    WRITE(4808,f_d8p4,ADVANCE   ='no') time_now
    WRITE(4808,f_d32p17,ADVANCE ='no') ds_avg
    WRITE(4808,f_d32p17,ADVANCE ='no') lyp_str_avg
    WRITE(4808,f_d32p17,ADVANCE ='no') lyp_eta_avg
    WRITE(4808,f_d32p17,ADVANCE ='no') ds_std
    WRITE(4808,f_d32p17,ADVANCE ='no') lyp_str_std
    WRITE(4808,f_d32p17,ADVANCE ='yes') lyp_eta_std

    IF ( t_step .EQ. t_step_total ) THEN
      CLOSE(4808)
    END IF
    !  +++++++++
  END

END MODULE system_decorrelator
