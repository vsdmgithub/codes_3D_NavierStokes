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
  DOUBLE PRECISION :: dis_avg,dis_std
  DOUBLE PRECISION::  dis_max,dis_min,dis_bin_size
  ! _________________________
  ! ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  ! DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::dummy_ar
  INTEGER(KIND=4) ::num_bin_dis
  ! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),    ALLOCATABLE ::Dis_pdf,Dis_val
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::Dis_fld                ! Dissipation field
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::s_xx,s_yy,s_zz
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::s_xy,s_yz,s_zx         ! Strain tensor
  ! _________________________________________
  ! FOURIER SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::dvxx,dvyy,dvzz

END MODULE system_advvariables
