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
! LAST MODIFIED: 16 September 2021
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
  DOUBLE PRECISION ::l_sys,l_int,l_kol
  DOUBLE PRECISION ::dx,dy,dz
  DOUBLE PRECISION ::N3,vol,dxdydz
  ! _________________________
  ! FOURIER SPACE
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  ::k_int,k_for,k_max,k_eta,k_kol
  INTEGER(KIND=4)  ::k_G
  INTEGER(KIND=4)  ::j_x,j_y,j_z
  INTEGER(KIND=4)  ::k_no
  INTEGER(KIND=4)  ::tot_active_modes,tot_modes
  DOUBLE PRECISION ::k_G_2
  ! _________________________
  ! TIME
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER (KIND=4) ::t_step,t_step_total,t_step_save,t_step_3D_save
  INTEGER (KIND=4) ::t_step_kol,t_step_debug
  INTEGER (KIND=4) ::no_of_saves,no_of_3D_saves
  DOUBLE PRECISION ::dt,dt_max
  DOUBLE PRECISION ::time_total,time_now,time_grid,time_save
  DOUBLE PRECISION ::time_kol,time_tur
  ! _________________________
  ! FLUID PARAMETERS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER (KIND=4) ::cfl_min
  INTEGER (KIND=4) ::forcing_status
  INTEGER (KIND=4) ::tot_forced_modes
  DOUBLE PRECISION ::viscosity
  DOUBLE PRECISION ::diss_rate_ref
  DOUBLE PRECISION ::force_factor
  DOUBLE PRECISION ::energy_initial
  ! _________________________
  ! FLUID FUNCTIONS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER (KIND=4) ::cfl_system
  INTEGER (KIND=4) ::rey_no,tay_rey_no
  DOUBLE PRECISION ::u_kol,u_int,u_rms
  DOUBLE PRECISION ::circulation
  DOUBLE PRECISION ::resolving_power
  DOUBLE PRECISION ::energy_forcing_modes
  DOUBLE PRECISION ::diss_rate, diss_rate_viscous
  DOUBLE PRECISION ::energy
  DOUBLE PRECISION ::energy_old
  DOUBLE PRECISION ::helicity
  DOUBLE PRECISION ::k_dot_v_norm
  DOUBLE PRECISION ::enstrophy
  DOUBLE PRECISION ::energy_mode
  DOUBLE PRECISION ::enstrophy_mode
  DOUBLE PRECISION ::helicity_mode
  DOUBLE COMPLEX   ::helicity_mode_complex
  ! _________________________
  ! MISCELLANEOUS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER (KIND=4) ::helicity_comp,velocitygrad_comp
  INTEGER (KIND=4) ::count_forcing_grid
  DOUBLE PRECISION ::norm_factor
  ! _________________________
  ! CHARACTERS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  CHARACTER( LEN = 10) :: N_char
  CHARACTER( LEN = 20) :: ic_type
  CHARACTER( LEN = 3 ) :: run_code
  CHARACTER( LEN = 3 ) :: test_code
  CHARACTER( LEN = 3 ) :: solver_alg
  ! _________________________
  ! DEBUGGING
  ! !!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4)  ::simulation_status,precheck_status,nan_count
  INTEGER(KIND=4)  ::incomp_error,nan_error,viscosity_error,debug_error
  INTEGER(KIND=4)  ::no_of_debug
  ! _________________________________________
  ! REAL SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::axis                   ! 1D Grid
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::u_x,u_y,u_z            ! Real velocity  (updated every time step)
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::u2_x,u2_y,u2_z            ! Real velocity  (updated every time step)
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::w_ux,w_uy,w_uz         ! Real vorticity (updated every time step)
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::proj_xx,proj_yy,proj_zz! Projection operators
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::proj_xy,proj_yz,proj_zx! \mathbb{P}_{ij}=\delta_{ij}-\frac{k_ik_j}{k^2}}
  ! _________________________________________
  ! FOURIER SPACE ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4) ,DIMENSION(:,:,:),ALLOCATABLE ::shell_no                 ! Every grid has its modulus, |k|
  INTEGER(KIND=4) ,DIMENSION(:),    ALLOCATABLE ::fkx,fky,fkz              ! List of forcing wavevectors
  DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::v_x,v_y,v_z              ! Spectral velocity (updated every time step)
  DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::v2_x,v2_y,v2_z              ! Spectral velocity (updated every time step)
  DOUBLE COMPLEX  ,DIMENSION(:,:,:),ALLOCATABLE ::w_vx,w_vy,w_vz           ! Spectral vorticity(updated every time step)
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::k_x,k_y,k_z               ! Wavevector components
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::k_2                       ! Spectral laplacian factor
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::truncator                 ! Truncating mask
  DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE ::integrating_factor        ! Integration factor for the dissipation

  ! _________________________________________
  ! SPECTRAL ARRAYS
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  INTEGER(KIND=4),DIMENSION(:),ALLOCATABLE      ::density_modes                               ! No of modes for each |k|
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::spectral_energy,spectral_energy_avg         ! E(k) shell averaged
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::spectral_enstrophy,spectral_enstrophy_avg   ! Z(k) shell averaged
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE     ::spectral_helicity,spectral_helicity_avg     ! H(k) shell averaged

  !TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

END MODULE system_basicvariables
