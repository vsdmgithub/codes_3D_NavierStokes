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
  DOUBLE PRECISION::lyp_S_avg
  DOUBLE PRECISION::lyp_eta_avg
  DOUBLE PRECISION::decor_old,lyp_decor
  INTEGER(KIND=4),PARAMETER ::lyp_bins=300
  ! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::diff_field
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::lyp_S,lyp_eta
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE ::lyp_val,pdf_lyp

  CONTAINS

  SUBROUTINE allocate_decorrelator
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate decorrelator
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(diff_field(0:N-1,0:N-1,0:N-1))
    ALLOCATE(lyp_S(0:N-1,0:N-1,0:N-1))
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

    CALL write_section('decor',diff_field( 0, : , : ) )
    ! REF <<< system_basicoutput >>>

  END

  SUBROUTINE compute_lyapunov_S
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the decorrelators
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    lyp_S     = -(  ( ( u_x - u2_x )**two ) * s_xx + ( ( u_y - u2_y )**two ) * s_yy + ( ( u_z - u2_z )**two ) * s_zz + &
              two * ( ( u_x - u2_x ) * s_xy * ( u_y - u2_y ) + &
                      ( u_y - u2_y ) * s_yz * ( u_z - u2_z ) + &
                      ( u_z - u2_z ) * s_zx * ( u_x - u2_x ) ) )

    lyp_S_avg = SUM( lyp_S ) / ( N3 * decor )

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

    lyp_eta_avg = SUM( lyp_eta ) / ( N3 * decor )


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
    DEALLOCATE(lyp_S)
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
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_eta_avg
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_S_avg
    WRITE(4002,f_d32p17,ADVANCE ='yes')lyp_S_avg + lyp_eta_avg

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
        WRITE(4008,f_d32p17,ADVANCE ='no')  lyp_S(   i_x, i_y, i_z ) / diff_field( i_x, i_y, i_z )
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
    INTEGER(KIND=4)::l_b,bin_count

    lyp_S     = ( lyp_S ) / ( diff_field )

    ! lyp_max = MAXVAL( lyp_S )
    ! lyp_min = MINVAL( lyp_S )

    ! --------------------------------------------------------------------
    ! Following values are for N=128 . Rescale accordingly for different N
    ! --------------------------------------------------------------------
    !lyp_max = 40.0D0
    !lyp_min = -30.0D0

    lyp_max     = +80.0D0
    lyp_min     = -80.0D0

    lyp_binsize = ( lyp_max - lyp_min ) / lyp_bins

    DO l_b = 1, lyp_bins

      lyp_val( l_b ) = lyp_min + DBLE( l_b ) * lyp_binsize

    END DO

    pdf_lyp   = zero
    bin_count = 0

    LOOP_RX_901: DO i_x = 0 , N - 1
    LOOP_RY_901: DO i_y = 0 , N - 1
    LOOP_RZ_901: DO i_z = 0 , N - 1

      l_b            = CEILING( ( lyp_S( i_x, i_y, i_z ) - lyp_min ) / lyp_binsize )
      ! Finding the bin slot

      BIN_CHECK_901: IF( l_b .LE. lyp_bins ) THEN

        pdf_lyp( l_b ) = pdf_lyp( l_b ) + one
        bin_count      = bin_count + 1

      END IF BIN_CHECK_901

    END DO LOOP_RX_901
    END DO LOOP_RY_901
    END DO LOOP_RZ_901

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
    INTEGER(KIND=4)::l_b,bin_count

    lyp_eta = lyp_eta / ( diff_field + tol )

    ! lyp_max = MAXVAL( lyp_eta )
    ! lyp_min = MINVAL( lyp_eta )

    ! --------------------------------------------------------------------
    ! Following values are for N=128 . Rescale accordingly for different N
    ! --------------------------------------------------------------------
    !lyp_max   = +60.0D0
    !lyp_min   = -80.0D0

    lyp_max     = +150.0D0
    lyp_min     = -150.0D0

    lyp_binsize = ( lyp_max - lyp_min ) / lyp_bins

    DO l_b = 1, lyp_bins

      lyp_val( l_b ) = lyp_min + DBLE( l_b ) * lyp_binsize

    END DO

    pdf_lyp   = zero
    bin_count = 0

    LOOP_RX_902: DO i_x = 0 , N - 1
    LOOP_RY_902: DO i_y = 0 , N - 1
    LOOP_RZ_902: DO i_z = 0 , N - 1

      DIFF_FIELD_CHECK_2: IF( diff_field( i_x, i_y, i_z ) .GT. tol ) THEN

        l_b            = CEILING( ( lyp_eta( i_x, i_y, i_z ) - lyp_min ) / lyp_binsize )
        ! Finding the bin slot

        BIN_CHECK: IF( l_b .LE. lyp_bins ) THEN

          pdf_lyp( l_b ) = pdf_lyp( l_b ) + one
          bin_count      = bin_count + 1

        END IF BIN_CHECK

      END IF DIFF_FIELD_CHECK_2

    END DO LOOP_RX_902
    END DO LOOP_RY_902
    END DO LOOP_RZ_902

    pdf_lyp = pdf_lyp / DBLE( bin_count )

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
    INTEGER(KIND=4)::l_b

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // 'lyp_S_pdf_' &
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
    INTEGER(KIND=4)::l_b

    WRITE (file_time,f_d8p4) time_now
    ! Writes 'time_now' as a CHARACTER

    file_name = TRIM( ADJUSTL( file_address ) ) // 'lyp_eta_pdf' &
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

END MODULE system_decorrelator
