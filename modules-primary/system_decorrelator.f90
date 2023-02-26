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
! LAST MODIFIED: 20 FEBRAURY 2023
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
  DOUBLE PRECISION::dec_old,dec_dot_str,dec_dot_vis
  DOUBLE PRECISION::lyp_max,lyp_min,lyp_bin_size,dec_dot
  DOUBLE PRECISION::lyp_avg,lyp_std
  INTEGER(KIND=4) ::m_ind,m_siz
  INTEGER(KIND=4) ::num_bin_lyp
  ! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::Dec_fld,Lyp_fld
  DOUBLE PRECISION,DIMENSION(:)    ,ALLOCATABLE ::Lyp_val,Lyp_pdf
  DOUBLE PRECISION,DIMENSION(:)    ,ALLOCATABLE ::Mom_val,Lyp_fun,Lyp_mom

  CONTAINS

  SUBROUTINE allocate_decorrelator
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to allocate decorrelator
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE

    num_bin_lyp    = CEILING( ( DBLE(N) / 128.0D0 ) * 100.0D0 )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  A  L  L  O  C  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ALLOCATE(Dec_fld(0:N-1,0:N-1,0:N-1))
    ALLOCATE(Lyp_fld(0:N-1,0:N-1,0:N-1))
    ALLOCATE(Lyp_val(num_bin_lyp),Lyp_pdf(num_bin_lyp))

    CALL allocate_strain_tensor
    ! REF-> <<< system_advdeclaration >>>

    m_siz = 40
    ALLOCATE(Lyp_fun(m_siz),Mom_val(m_siz),Lyp_mom(m_siz))

    DO m_ind = 1,m_siz
      Mom_val( m_ind )            = -4.0D0 + m_ind * 0.2D0
    END DO

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
    CALL compute_lyapunov_field

    lyp_decor   = ( DLOG( dec ) - DLOG( dec_old ) ) / dt
    ! Finding the exponent at which it grows at 't'

    dec_old     = dec
    ! copying the old decorrelator

    lyp_dec_str = dec_dot_str / dec
    lyp_dec_vis = dec_dot_vis / dec
    lyp_dec_cal = dec_dot / dec

    CALL write_decorrelator_growth

  END

  SUBROUTINE compute_generalized_lyapunov
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! Find the generalized exponents L(q)
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  IMPLICIT NONE

    !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  G E N E R A L I Z E D        D E C O R R E L A T O R
    !  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    LOOP_MOMENTS: DO m_ind = 1, m_siz

      Lyp_fun( m_ind ) = DLOG( SUM( DEXP( Lyp_fld * Mom_val( m_ind ) ) ) / N3 )
      Lyp_mom( m_ind ) = SUM( DABS( Lyp_fld ) ** Mom_val( m_ind ) ) / N3

    END DO LOOP_MOMENTS

    CALL write_generalized_lyapunov
    ! Writes all the generalized decorrelators

  END

  SUBROUTINE compute_lyapunov_field
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL this to compute the lyapunov from strain and viscosity
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    IMPLICIT NONE
    INTEGER(KIND=4)::ip_x,ip_y,ip_z
    INTEGER(KIND=4)::im_x,im_y,im_z
    DOUBLE PRECISION::strain_growth,viscous_growth
    DOUBLE PRECISION::df_x, df_y, df_z
    DOUBLE PRECISION::df_x_px, df_y_px, df_z_px
    DOUBLE PRECISION::df_x_py, df_y_py, df_z_py
    DOUBLE PRECISION::df_x_pz, df_y_pz, df_z_pz
    DOUBLE PRECISION::df_x_mx, df_y_mx, df_z_mx
    DOUBLE PRECISION::df_x_my, df_y_my, df_z_my
    DOUBLE PRECISION::df_x_mz, df_y_mz, df_z_mz

    CALL compute_strain_tensor
    ! REF-> <<< system_advfunctions >>>

    dec         = zero
    dec_dot_str = zero
    dec_dot_vis = zero
    lyp_avg     = zero
    lyp_std     = zero

    DO i_x = 0 , N - 1
    DO i_y = 0 , N - 1
    DO i_z = 0 , N - 1

    Dec_fld( i_x, i_y, i_z )= hf*( (Ub_x(i_x,i_y,i_z)-U_x(i_x,i_y,i_z))**two + &
                                   (Ub_y(i_x,i_y,i_z)-U_y(i_x,i_y,i_z))**two + &
                                   (Ub_z(i_x,i_y,i_z)-U_z(i_x,i_y,i_z))**two )
    ! Matrix of decorelation

    DIFF_FIELD_CHECK: IF( Dec_fld( i_x, i_y, i_z ) .GT. zero ) THEN

      dec           = dec + Dec_fld( i_x, i_y, i_z )
      ! Global decorrelation

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  STRAIN BASED TIMESCALES FIELD
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      strain_growth = -( ((U_x(i_x,i_y,i_z)-Ub_x(i_x,i_y,i_z))**two)*s_xx(i_x,i_y,i_z)+&
                         ((U_y(i_x,i_y,i_z)-Ub_y(i_x,i_y,i_z))**two)*s_yy(i_x,i_y,i_z)+&
                         ((U_z(i_x,i_y,i_z)-Ub_z(i_x,i_y,i_z))**two)*s_zz(i_x,i_y,i_z)+&
      two*( (U_x(i_x,i_y,i_z)-Ub_x(i_x,i_y,i_z))*s_xy(i_x,i_y,i_z)*(U_y(i_x,i_y,i_z)-Ub_y(i_x,i_y,i_z))+&
            (U_y(i_x,i_y,i_z)-Ub_y(i_x,i_y,i_z))*s_yz(i_x,i_y,i_z)*(U_z(i_x,i_y,i_z)-Ub_z(i_x,i_y,i_z))+&
            (U_z(i_x,i_y,i_z)-Ub_z(i_x,i_y,i_z))*s_zx(i_x,i_y,i_z)*(U_x(i_x,i_y,i_z)-Ub_x(i_x,i_y,i_z)) ) )
      ! The strain based decorrelation growth

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

      df_x    = Ub_x( i_x, i_y, i_z ) - U_x( i_x, i_y, i_z )
      df_y    = Ub_y( i_x, i_y, i_z ) - U_y( i_x, i_y, i_z )
      df_z    = Ub_z( i_x, i_y, i_z ) - U_z( i_x, i_y, i_z )

      df_x_px = Ub_x( ip_x, i_y, i_z ) - U_x( ip_x, i_y, i_z )
      df_y_px = Ub_y( ip_x, i_y, i_z ) - U_y( ip_x, i_y, i_z )
      df_z_px = Ub_z( ip_x, i_y, i_z ) - U_z( ip_x, i_y, i_z )
      df_x_py = Ub_x( i_x, ip_y, i_z ) - U_x( i_x, ip_y, i_z )
      df_y_py = Ub_y( i_x, ip_y, i_z ) - U_y( i_x, ip_y, i_z )
      df_z_py = Ub_z( i_x, ip_y, i_z ) - U_z( i_x, ip_y, i_z )
      df_x_pz = Ub_x( i_x, i_y, ip_z ) - U_x( i_x, i_y, ip_z )
      df_y_pz = Ub_y( i_x, i_y, ip_z ) - U_y( i_x, i_y, ip_z )
      df_z_pz = Ub_z( i_x, i_y, ip_z ) - U_z( i_x, i_y, ip_z )

      df_x_mx = Ub_x( im_x, i_y, i_z ) - U_x( im_x, i_y, i_z )
      df_y_mx = Ub_y( im_x, i_y, i_z ) - U_y( im_x, i_y, i_z )
      df_z_mx = Ub_z( im_x, i_y, i_z ) - U_z( im_x, i_y, i_z )
      df_x_my = Ub_x( i_x, im_y, i_z ) - U_x( i_x, im_y, i_z )
      df_y_my = Ub_y( i_x, im_y, i_z ) - U_y( i_x, im_y, i_z )
      df_z_my = Ub_z( i_x, im_y, i_z ) - U_z( i_x, im_y, i_z )
      df_x_mz = Ub_x( i_x, i_y, im_z ) - U_x( i_x, i_y, im_z )
      df_y_mz = Ub_y( i_x, i_y, im_z ) - U_y( i_x, i_y, im_z )
      df_z_mz = Ub_z( i_x, i_y, im_z ) - U_z( i_x, i_y, im_z )

      viscous_growth = df_x*(-6.0D0*df_x+df_x_px+df_x_py+df_x_pz+df_x_mx+df_x_my+df_x_mz) + &
                       df_y*(-6.0D0*df_y+df_y_px+df_y_py+df_y_pz+df_y_mx+df_y_my+df_y_mz) + &
                       df_z*(-6.0D0*df_z+df_z_px+df_z_py+df_z_pz+df_z_mx+df_z_my+df_z_mz)
      viscous_growth = viscosity * viscous_growth / ( dx * dx )
      ! The viscous diffusion based decorrelation growth

      Lyp_fld( i_x, i_y, i_z ) = ( strain_growth + viscous_growth ) / Dec_fld( i_x, i_y, i_z )

      dec_dot_str              = dec_dot_str + strain_growth
      dec_dot_vis              = dec_dot_vis + viscous_growth
      lyp_avg                  = lyp_avg + Lyp_fld( i_x, i_y, i_z )
      lyp_std                  = lyp_std + ( Lyp_fld( i_x, i_y, i_z ) ** two )

    ELSE

      Lyp_fld( i_x, i_y, i_z ) = zero

    END IF DIFF_FIELD_CHECK

    END DO
    END DO
    END DO

    dec         = dec         / N3
    dec_dot_str = dec_dot_str / N3
    dec_dot_vis = dec_dot_vis / N3
    dec_dot     = dec_dot_str + dec_dot_vis
    lyp_avg     = lyp_avg     / N3
    lyp_std     = ( lyp_std   / N3 ) - lyp_avg ** two

  END

  SUBROUTINE compute_pdf_lyapunov
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
    INTEGER(KIND=4) ::bin,num_pts_pdf
    DOUBLE PRECISION::lyp_dum,std_fac

    ! --------------------------------------------------------------------
    ! Following values are for N=128 . Rescale accordingly for different N
    ! --------------------------------------------------------------------
    ! lyp_max    = +40.0D0 * ( N / 128 )
    ! lyp_min    = -40.0D0 * ( N / 128 )

    std_fac      = 10.0D0
    lyp_dum      = lyp_avg
    lyp_std      = lyp_std / ( lyp_avg ** two )
    lyp_avg      = one
    lyp_min      = lyp_avg - std_fac * lyp_std
    lyp_max      = lyp_avg + std_fac * lyp_std
    lyp_bin_size = ( lyp_max - lyp_min ) / num_bin_lyp

    DO bin = 1, num_bin_lyp
      Lyp_val( bin ) = lyp_min + ( DBLE( bin ) - hf ) * lyp_bin_size
    END DO

    lyp_pdf     = zero
    num_pts_pdf = 0

    LOOP_RX_902: DO i_x = 0 , N - 1
    LOOP_RY_902: DO i_y = 0 , N - 1
    LOOP_RZ_902: DO i_z = 0 , N - 1

      Lyp_fld( i_x, i_y, i_z ) = Lyp_fld( i_x, i_y, i_z ) / lyp_dum
      bin                      = CEILING( ( Lyp_fld( i_x, i_y, i_z ) - lyp_min ) / lyp_bin_size )

      BIN_CHECK_904: IF( (bin .GT. 1 ) .AND. (bin .LT. num_bin_lyp ) ) THEN
        Lyp_pdf( bin )         = Lyp_pdf( bin ) + one
        num_pts_pdf            = num_pts_pdf + 1
      ELSEIF (bin .LE. 1 ) THEN
        Lyp_pdf( 1 )           = Lyp_pdf( 1 ) + one
      ELSEIF (bin .GE. num_bin_lyp ) THEN
        Lyp_pdf( num_bin_lyp ) = Lyp_pdf( num_bin_lyp ) + one
      END IF BIN_CHECK_904

    END DO LOOP_RZ_902
    END DO LOOP_RY_902
    END DO LOOP_RX_902

    Lyp_pdf   = Lyp_pdf / N3

    IF  ( DBLE(num_pts_pdf) / N3 .LT. 0.95D0 ) THEN
      PRINT*," More than 5 percent of values are falling out of the lyapunov histogram bins"
    END IF

    CALL write_lyp_pdf
    ! Normalizing the PDF

    CALL compute_pdf_dissipation
    ! REF-> <<< system_advfunctions >>>

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
      fil_name = TRIM( ADJUSTL( fil_adrs ) ) // 'decor_growth.dat'
      !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      OPEN(unit = 4002, file = fil_name )
    END IF

    WRITE(4002,f_d8p4,ADVANCE   ='no') t_now
    WRITE(4002,f_d32p17,ADVANCE ='no') dec
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_decor
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_dec_cal
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_dec_str
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_dec_vis
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_avg
    WRITE(4002,f_d32p17,ADVANCE ='no') lyp_std
    WRITE(4002,f_d32p17,ADVANCE ='no') dis_avg
    WRITE(4002,f_d32p17,ADVANCE ='yes')dis_std

    IF ( t_step .EQ. t_step_total ) THEN
      CLOSE(4002)
    END IF
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END

  SUBROUTINE write_generalized_lyapunov
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! CALL THIS SUBROUTINE TO:
  ! write the data in time, for every timestep
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    IMPLICIT NONE

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  G E N .   D E C O R    V S    T I M E
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    WRITE (fil_time,f_d8p4) t_now
    ! Writes 't_now' as a CHARACTER

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // TRIM( ADJUSTL( dir_pdf ) ) // &
               'lyp_fun_t_' // TRIM( ADJUSTL ( fil_time ) ) // '.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    OPEN(unit = 4028, file = fil_name )

    DO m_ind = 1,m_siz
      WRITE(4028,f_d8p4,  ADVANCE   ='no')  Mom_val( m_ind )
      WRITE(4028,f_d32p17,ADVANCE   ='no')  Lyp_fun( m_ind )
      WRITE(4028,f_d32p17,ADVANCE   ='yes') Lyp_mom( m_ind )
    END DO
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    CLOSE(4028)

  END

  SUBROUTINE write_lyp_pdf
  ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  ! ------------
  ! This writes the pdf of lyapunov and timescales
  ! -------------
  ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  IMPLICIT NONE
    ! _________________________
    ! LOCAL VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::bin

    WRITE (fil_time,f_d8p4) t_now
    ! Writes 't_now' as a CHARACTER

    fil_name = TRIM( ADJUSTL( fil_adrs ) ) // TRIM( ADJUSTL( dir_pdf ) ) // &
                  'lyp_pdf_t_' // TRIM( ADJUSTL ( fil_time ) ) // '.dat'
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    OPEN( unit = 799, file = fil_name )
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  P  R  I  N   T          O  U  T  P  U  T   -   DATA FILE
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    DO bin = 1, num_bin_lyp

      WRITE(799,f_d12p6,ADVANCE  ='NO')  Lyp_val( bin )
      WRITE(799,f_d32p17,ADVANCE ='YES') Lyp_pdf( bin )

    END DO

    CLOSE(799)
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
    DEALLOCATE(Dec_fld,Lyp_fld)
    DEALLOCATE(Lyp_val, Lyp_pdf)
    DEALLOCATE(Mom_val,Lyp_fun,Lyp_mom)

    CALL deallocate_strain_tensor
    ! REF-> <<< system_advdeclaration >>>

  END

END MODULE system_decorrelator
