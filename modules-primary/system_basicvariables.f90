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

! #########################
! MODULE: system_basicvariables
! LAST MODIFIED: 20 FEBRAURY 2023
! #########################

! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! BASIC VARIABLES AND ARRAYS FOR 3D NSE EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

MODULE system_basicvariables
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the global variables and arrays for the simulation space are declared here.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
  !  SUB-MODULES
  !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
  USE system_constants

  IMPLICIT NONE
  ! _________________________
  ! REAL SPACE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  ::N,Nh
  INTEGER(KIND=4)  ::i_x,i_y,i_z
  DOUBLE PRECISION ::l_sys,l_grd
  DOUBLE PRECISION ::N3
  ! _________________________
  ! FOURIER SPACE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  ::k_tru,k_ind,k_max
  INTEGER(KIND=4)  ::j_x,j_y,j_z
  INTEGER(KIND=4)  ::num_mod
  DOUBLE PRECISION ::k_tru_sqr
  ! _________________________
  ! TIME
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER (KIND=4) ::t_step,t_step_tot,t_step_sav
  INTEGER (KIND=4) ::num_sav,ind_sav
  DOUBLE PRECISION ::dt,dt_max
  DOUBLE PRECISION ::t_now,t_tot,t_sav,t_min
  ! _______________________________
  ! FLUID PARAMETERS AND FUNCTIONS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER (KIND=4) ::frc_status
  INTEGER (KIND=4) ::k_int,k_frc,k_kol,k_tay
  INTEGER (KIND=4) ::num_mod_frc
  INTEGER (KIND=4) ::t_step_kol,t_step_grd,t_step_int
  INTEGER (KIND=4) ::cfl_min,cfl
  INTEGER (KIND=4) ::rey_int,rey_tay,rey_kol
  DOUBLE PRECISION ::t_grd,t_kol,t_int
  DOUBLE PRECISION ::l_int,l_kol,l_tay
  DOUBLE PRECISION ::u_kol,u_int
  DOUBLE PRECISION ::vis,vis_min
  DOUBLE PRECISION ::res_pow
  DOUBLE PRECISION ::frc_fac,frc_fac_avg,nrm_fac
  DOUBLE PRECISION ::eng_0,eng,eng_frc,eng_mod,eng_pre
  DOUBLE PRECISION ::dis_ref,dis,dis_fix,dis_eng
  DOUBLE PRECISION ::hel,ens,inc
  ! _______________________________
  ! CHAOS PARAMETERS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION ::dec,dec_ini,dec_fin,dec_pre
  DOUBLE PRECISION ::lyp_dec,lyp_dec_str,lyp_dec_vis,lyp_dec_cal
  DOUBLE PRECISION ::eng_b
  ! _________________________
  ! CHARACTERS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  CHARACTER( LEN = 6)  :: res_char
  CHARACTER( LEN = 20) :: icn_type
  CHARACTER( LEN = 1 ) :: run_code
  CHARACTER( LEN = 1 ) :: tst_code
  CHARACTER( LEN = 2 ) :: sol_algm
  ! _________________________
  ! DEBUGGING
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  ::num_deb,num_nan
  INTEGER(KIND=4)  ::t_step_deb
  INTEGER(KIND=4)  ::sim_status,pre_status
  INTEGER(KIND=4)  ::tur_status,lyp_status
  INTEGER(KIND=4)  ::inc_err,nan_err,vis_err,deb_err
  ! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::Axis
  ! 1D Grid
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::U_x,U_y,U_z
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::Ub_x,Ub_y,Ub_z
  ! Real velocity  (updated every time step)
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::W_x,W_y,W_z
  ! Real vorticity (updated every time step)
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::Proj_xx,Proj_yy,Proj_zz
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::Proj_xy,Proj_yz,Proj_zx
  ! Projection operators
  ! _________________________________________
  ! FOURIER SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4) ,DIMENSION(:,:,:),ALLOCATABLE ::Shell
  ! Every grid has its modulus, |k|
  INTEGER(KIND=4) ,DIMENSION(:),    ALLOCATABLE ::F_kx_ind,F_ky_ind,F_kz_ind
  ! List of forcing wavevectors
  DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::V_x,V_y,V_z
  ! Spectral velocity (updated every time step)
  DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::Vb_x,Vb_y,Vb_z
  ! Spectral velocity (updated every time step)
  DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::F_x,F_y,F_z
  ! Forcing vector - Spectral
  DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::W_kx,W_ky,W_kz
  ! Spectral vorticity(updated every time step)
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::K_x,K_y,K_z
  ! Wavevector components
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::Lapla
  ! Spectral laplacian factor
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::Trunc
  ! Truncation mask
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::I_fac
  ! Integration factor for the dissipation
  ! _________________________________________
  ! SPECTRAL ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE      ::Den_k
  ! No of modes for each |k|
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::Eng_k,Eng_k_avg
  ! E(k) shell averaged
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::frc_fac_his

  !TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

END MODULE system_basicvariables
