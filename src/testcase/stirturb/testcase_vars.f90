MODULE MOD_TestCase_Vars
IMPLICIT NONE
PUBLIC
SAVE

!! ========================================================================= !!
!! Analysis

INTEGER :: nAnalyzeTestCase     = 9999999     !< call AnalyzeTestCase every *th time step. May be adjusted in parameter file
LOGICAL :: doTCSource           = .TRUE.     !< compute source terms for testcase
CHARACTER(LEN=255) :: testcase  = "decaying_turbulence" !< name of testcase

REAL,ALLOCATABLE   :: Time(:)               !< times of log data (nWriteStats)
REAL,ALLOCATABLE   :: writeBuf(:,:)         !< log data (nTGVVars+1,nWriteStats)
INTEGER,PARAMETER  :: nOutputVars   = 12    !< Number of variables to be evaluated for TGV, time not included
!!INTEGER,PARAMETER  :: nOutputVars   = 3    !< Number of variables to be evaluated for TGV, time not included
INTEGER            :: ioCounter     = 0     !< current number of buffer items
INTEGER            :: nWriteStats   = -999  !< Write testcase statistics to file at every n-th AnalyzeTestcase step
CHARACTER(LEN=255) :: Filename              !< filename to store testcase log data

CHARACTER(LEN=255) :: st_outputfp = 'stirring.dat'

!! ========================================================================= !!
!! Initialization

real :: init_dens = 1.
real :: init_velx = 0.
real :: init_vely = 0.
real :: init_velz = 0.
real :: init_pres = 1.

!! ========================================================================= !!
!! Forcing/Stirring

!! Size of the forcing box
logical :: st_verbose
integer :: st_seed

integer :: st_wait_time_steps

real*8 :: st_dVel_max_ampli
real*8 :: st_injection_scale
real*8 :: st_rms_velocity
real*8 :: st_dVel_fraction
real*8 :: st_Lbox
real*8 :: st_tmax
real*8 :: st_Tac
real*8 :: st_kmin, st_kmax
real*8 :: st_zeta
real*8 :: st_F0
real*8 :: st_g0

real*8 :: st_mach
real*8 :: st_force_param
real*8 :: st_force_base

integer :: NperOUstep

LOGICAL :: st_active = .TRUE.
LOGICAL :: st_stop_at_mach = .FALSE.

!! Variables for live OU steps
real*8, dimension(999) :: st_OU_kx_re
real*8, dimension(999) :: st_OU_kx_im
real*8, dimension(999) :: st_OU_ky_re
real*8, dimension(999) :: st_OU_ky_im
real*8, dimension(999) :: st_OU_kz_re
real*8, dimension(999) :: st_OU_kz_im
real*8, dimension(999) :: st_OU_len, st_OU_lenx, st_OU_leny, st_OU_lenz

END MODULE MOD_TestCase_Vars
