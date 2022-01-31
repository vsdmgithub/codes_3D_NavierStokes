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
! LAST MODIFIED: 21 JUNE 2021
! ##################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! ADVANCED VARIABLES TO DO ANALYSIS IN 3D EULER
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_advvariables
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This contains all major advanced functions involving analysis.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  IMPLICIT NONE
  ! _________________________
  ! VARIABLES
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  :: ev_mod_bins,ev_dif_bins
  INTEGER(KIND=4)  :: q_bins, r_bins
  INTEGER(KIND=4)  :: ds_bins
  INTEGER(KIND=4)  :: vx_bins
  INTEGER(KIND=4)  :: data_sz
  INTEGER(KIND=4)  :: jump_sz
  DOUBLE PRECISION :: ds_avg
  DOUBLE PRECISION :: vx_alp_avg
  ! DOUBLE PRECISION :: ds_std
  ! DOUBLE PRECISION :: vx_alp_std
  ! _________________________
  ! ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  ! DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dummy_ar
  ! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:,:),  ALLOCATABLE ::pdf_qr
  DOUBLE PRECISION,DIMENSION(:,:),  ALLOCATABLE ::pdf_ev
  DOUBLE PRECISION,DIMENSION(:),    ALLOCATABLE ::pdf_ds
  DOUBLE PRECISION,DIMENSION(:),    ALLOCATABLE ::pdf_vx
  DOUBLE PRECISION,DIMENSION(:),    ALLOCATABLE ::ds_val
  DOUBLE PRECISION,DIMENSION(:),    ALLOCATABLE ::vx_val
  DOUBLE PRECISION,DIMENSION(:),    ALLOCATABLE ::q_val,r_val
  DOUBLE PRECISION,DIMENSION(:),    ALLOCATABLE ::ev_mod_val,ev_dif_val
  DOUBLE PRECISION,DIMENSION(:),    ALLOCATABLE ::ev_mod,ev_dif
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::ds_rate                ! Dissipation field
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::vx_alp                 ! Vortex Stretching
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::q_invar,r_invar        ! Invariants of velocity gradient tensor
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::s_xx,s_yy,s_zz
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::s_xy,s_yz,s_zx         ! Strain tensor
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::duxx,duyy,duzz         ! Gradient of velocity
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::duxy,duyx
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::duzy,duyz
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::duxz,duzx              ! Gradient of velocity
  ! _________________________________________
  ! FOURIER SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::dvxx,dvyy,dvzz
  ! DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::dvxy,dvyx
  ! DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::dvzx,dvxz
  ! DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::dvyz,dvzy                 ! Spectral gradient of velocity
  ! DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::h_pos_x,h_pos_y,h_pos_z
  ! DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::h_neg_x,h_neg_y,h_neg_z   ! Helical basis

  ! _________________________________________
  ! SPECTRAL ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE system_advvariables
