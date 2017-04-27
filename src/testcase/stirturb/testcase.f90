#include "flexi.h"

# define DENS 1
# define PRES 2
# define TEMP 3
# define ENER 4
# define EINT 5
# define EKIN 6
# define RMSV 7
# define MASS 8
# define VOLU 9
# define FVV 10
# define PRES_MIN 11
# define PRES_MAX 12
# define MACH_MIN 13
# define MACH_MAX 14
# define RMSV_MAX 15

MODULE MOD_Testcase

IMPLICIT NONE
PRIVATE

INTERFACE DefineParametersTestcase
  MODULE PROCEDURE DefineParametersTestcase
End INTERFACE

INTERFACE InitTestcase
  MODULE PROCEDURE InitTestcase
END INTERFACE

INTERFACE FinalizeTestcase
  MODULE PROCEDURE FinalizeTestcase
END INTERFACE

INTERFACE ExactFuncTestcase
  MODULE PROCEDURE ExactFuncTestcase
END INTERFACE

INTERFACE CalcForcing
  MODULE PROCEDURE DO_NOTHING
END INTERFACE

INTERFACE TestcaseSource
  MODULE PROCEDURE TestcaseSource
END INTERFACE

INTERFACE AnalyzeTestCase
  MODULE PROCEDURE AnalyzeTestCase
END INTERFACE

INTERFACE GetBoundaryFluxTestcase
  MODULE PROCEDURE GetBoundaryFluxTestcase
END INTERFACE

INTERFACE GetBoundaryFVgradientTestcase
  MODULE PROCEDURE GetBoundaryFVgradientTestcase
END INTERFACE

INTERFACE Lifting_GetBoundaryFluxTestcase
  MODULE PROCEDURE Lifting_GetBoundaryFluxTestcase
END INTERFACE

PUBLIC:: DefineParametersTestcase
PUBLIC:: InitTestcase
PUBLIC:: FinalizeTestcase
PUBLIC:: ExactFuncTestcase
PUBLIC:: TestcaseSource
PUBLIC:: CalcForcing
PUBLIC:: AnalyzeTestCase
PUBLIC:: GetBoundaryFluxTestcase
PUBLIC:: GetBoundaryFVgradientTestcase
PUBLIC:: Lifting_GetBoundaryFluxTestcase

CONTAINS

SUBROUTINE DefineParametersTestcase

    USE MOD_Globals
    USE MOD_ReadInTools ,ONLY: prms
    IMPLICIT NONE

    CALL prms%SetSection("Testcase")
    CALL prms%CreateIntOption('nWriteStats', "Write testcase statistics to file at every n-th AnalyzeTestcase step.", '100')
    CALL prms%CreateIntOption('nAnalyzeTestCase', "Call testcase specific analysis routines every n-th timestep. "//&
                                                  "(Note: always called at global analyze level)", '10')

    CALL prms%CreateRealOption('init_dens', "Initial density.",    '1.0')
    CALL prms%CreateRealOption('init_velx', "Initial velocity x.", '0.0')
    CALL prms%CreateRealOption('init_vely', "Initial velocity y.", '0.0')
    CALL prms%CreateRealOption('init_velz', "Initial velocity z.", '0.0')
    CALL prms%CreateRealOption('init_pres', "Initial pressure.",   '1.0')

    CALL  prms%CreateIntOption('st_seed',   "Random Seed",      '713')
    CALL  prms%CreateIntOption('st_wait_time_steps',   "Wait n time steps between forcing", '10')

    CALL prms%CreateRealOption('st_Lbox',   "Box Dimensions",   '1.0')
    CALL prms%CreateRealOption('st_Tac',    "Acceleration Time Scale",  '1.0')
    CALL prms%CreateRealOption('st_kmin',   "Minimum Forcing k number", '1.0')
    CALL prms%CreateRealOption('st_kmax',   "Maximum Forcing k number", '3.0')
    CALL prms%CreateRealOption('st_zeta',   "Ratio compressive to solenoidal forcing.",   '0.5')
    CALL prms%CreateRealOption('st_F0',     "Base Force Value",         '1.0')
    CALL prms%CreateRealOption('st_g0',     "Some g parameter",         '1.0')
    CALL prms%CreateRealOption('st_tmax',   "maximum stirring time",    '99999.0')

    CALL prms%CreateRealOption('st_force_base',     "Base Force",         '1.0')
    CALL prms%CreateRealOption('st_force_param',    "Forcing Parameter",  '1.0')
    CALL prms%CreateRealOption('st_mach',           "Global Mach Number", '1.0')

    CALL prms%CreateLogicalOption('st_active',      "active or not", 'T')
    CALL prms%CreateLogicalOption('st_stop_at_mach',      "stop when mach nr. is reached", 'F')

END SUBROUTINE DefineParametersTestcase

SUBROUTINE InitTestcase

    USE MOD_PreProc,        ONLY: N
    USE MOD_Globals
    USE MOD_ReadInTools,    ONLY: GETINT, GETREAL, GETLOGICAL
    USE MOD_Output_Vars,    ONLY: ProjectName
    USE MOD_TestCase_Vars
    USE MOD_Output,         ONLY: InitOutputToFile
    USE MOD_Analyze_Vars,   ONLY: NAnalyze

    IMPLICIT NONE

    CHARACTER(LEN=31)        :: varnames(nOutputVars)

    integer :: ikxmin, ikxmax, ikymin, ikymax, ikzmin, ikzmax
    integer :: ikx, iky, ikz

    real    :: kx, ky, kz, k, twopi, length
    real, save :: imin, imax, jmin, jmax, kmin, kmax

    REAL*8  :: tif

    integer :: i, iii, it, m, dr_Nstep
    real*8  :: xr8_normal_01, dt, dummy

    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' INIT TESTCASE "STIRTURB" ...'

    nWriteStats      = GETINT( 'nWriteStats','100')
    nAnalyzeTestCase = GETINT( 'nAnalyzeTestCase','10')

    init_dens        = GETREAL('init_dens', '1.0')
    init_velx        = GETREAL('init_velx', '1.0')
    init_vely        = GETREAL('init_vely', '0.0')
    init_velz        = GETREAL('init_velz', '0.0')
    init_pres        = GETREAL('init_pres', '0.0')

    st_seed =  GETINT('st_seed',    '713')
    st_Lbox = GETREAL('st_Lbox',    '1.0')
    st_Tac  = GETREAL('st_Tac',     '1.0')
    st_kmin = GETREAL('st_kmin',    '1.0')
    st_kmax = GETREAL('st_kmax',    '3.0')
    st_zeta = GETREAL('st_zeta',    '0.5')
    st_F0   = GETREAL('st_F0',      '1.0')
    st_g0   = GETREAL('st_g0',      '1.0')
    st_tmax = GETREAL('st_tmax',    '999999.0')
    st_mach = GETREAL('st_mach',    '1.0')

    st_active = GETLOGICAL('st_active',    'T')
    st_stop_at_mach = GETLOGICAL('st_stop_at_mach',    'F')

    st_force_base  = GETREAL('st_force_base',  '1.0')
    st_force_param = GETREAL('st_force_param', '1.0')

    st_wait_time_steps = GETINT('st_wait_time_steps',    '10')

    SWRITE(UNIT_stdOut,'(A)')' INIT TESTCASE "STIRTURB" DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')

    !! ===================================================================== !!
    !! INIT FORCING
 
    ! new simulation, set all vectors to zero
    SWRITE(UNIT_stdOut,'(A)')'   new simulation: initialise OU vectors'

    st_OU_kx_re(:) = 0.0
    st_OU_kx_im(:) = 0.0
    st_OU_ky_re(:) = 0.0
    st_OU_ky_im(:) = 0.0
    st_OU_kz_re(:) = 0.0
    st_OU_kz_im(:) = 0.0

    ! set length of the vectors and the determine, how many modes are needed
    NperOUstep = 0
    do ikz = -CEILING(st_kmax), CEILING(st_kmax) 
        do iky = -CEILING(st_kmax), CEILING(st_kmax) 
            do ikx = -CEILING(st_kmax), CEILING(st_kmax)
                length = SQRT(DBLE(ikx)**2 + DBLE(iky)**2 + DBLE(ikz)**2)
                if((length .le. st_kmax) .and. (length .ge. st_kmin)) then

                    ! valid k vector
                    NperOUstep = NperOUstep + 1
                    st_OU_lenx(NperOUstep) = DBLE(ikx)
                    st_OU_leny(NperOUstep) = DBLE(iky)
                    st_OU_lenz(NperOUstep) = DBLE(ikz)
                    st_OU_len (NperOUstep) = length

                end if
            end do
        end do
    end do
        
END SUBROUTINE InitTestcase

SUBROUTINE ExactFuncTestcase(tIn,x,Resu,Resu_t,Resu_tt)

    !!USE MOD_Globals, ONLY:  Abort
    USE MOD_EOS_Vars, ONLY: kappa
    USE MOD_Preproc, ONLY:  PP_PI
    USE MOD_EOS, ONLY:      PrimToCons
    USE MOD_TestCase_Vars

    IMPLICIT NONE

    REAL,INTENT(IN)     :: x(3)        !< position in physical coordinates
    REAL,INTENT(IN)     :: tIn         !< current simulation time

    REAL,INTENT(OUT)    :: Resu(5)     !< exact fuction evaluated at tIn, returning state in conservative variables
    REAL,INTENT(OUT)    :: Resu_t(5)   !< first time deriv of exact fuction
    REAL,INTENT(OUT)    :: Resu_tt(5)  !< second time deriv of exact fuction

    REAL                :: A,Ms,prim(PP_nVarPrim)

    prim(1) = init_dens
    prim(2) = init_velx
    prim(3) = init_vely 
    prim(4) = init_velz
    prim(5) = init_pres 

    CALL PrimToCons(prim, Resu)

    Resu_t  = 0.
    Resu_tt = 0.

END SUBROUTINE ExactFuncTestcase

SUBROUTINE TestcaseSource(Ut)

    USE MOD_Globals
    USE MOD_PreProc, ONLY: N
    USE MOD_testcase_vars
    USE MOD_TimeDisc_Vars, ONLY: t,dt
    USE MOD_Mesh_Vars, ONLY:Elem_xGP,sJ,nElems

    IMPLICIT NONE

    REAL,INTENT(INOUT) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< solution time derivative

    REAL :: mach_avg, force
    INTEGER :: wait_time_steps = 1
    INTEGER :: time_steps_since_last_forcing = 0

    !! analyze output
    integer :: count = 0
    real*8, dimension(99) :: resBuf

    if (.not. st_active) return

    if (t .gt. st_tmax) then
        st_active = .false.
        return
    end if

    if (wait_time_steps .gt. st_wait_time_steps) then
        call CalcMachAvg(mach_avg)
        if (mach_avg .lt. st_mach) then
            !!force = st_force_base / (1 + exp(st_force_param*(mach_avg - 0.9*st_mach)))
            force = st_force_base

            if (mach_avg .gt. 0.95*st_mach) then
                force = st_force_base * ((st_force_param - 1)/(0.05*st_mach)*(mach_avg-st_mach*0.95) + 1)
            end if

            call OU_time_step()
            call ApplyForcing(force,Ut)

            !!SWRITE(UNIT_stdOut,'(A,5(1X,ES20.8))') 'foc: ', t, dt, force, mach_avg

            time_steps_since_last_forcing = 1
        ELSE
            force = 0.0
            if (st_stop_at_mach) st_active = .FALSE.
            !!SWRITE(UNIT_stdOut,'(A,5(1X,ES20.8))') 'NOF: ', t, dt, force, mach_avg
        end if

        count = 0
        resBuf(incr(count)) = t
        resBuf(incr(count)) = dt
        resBuf(incr(count)) = force
        resBuf(incr(count)) = mach_avg
        call write2file(st_fp_forcing, count, resBuf)
     
        wait_time_steps = 1
        
    end if
        
    time_steps_since_last_forcing = time_steps_since_last_forcing + 1
    wait_time_steps = wait_time_steps + 1

END SUBROUTINE TestcaseSource

SUBROUTINE ApplyForcing(force,Ut)
    USE MOD_Globals
    USE MOD_DG_Vars,        ONLY: U
    USE MOD_PreProc,        ONLY: N
    USE MOD_testcase_vars
    USE MOD_Mesh_Vars,      ONLY: NodeCoords
    USE MOD_Mesh_Vars,      ONLY: Elem_xGP,nElems
    USE MOD_Analyze_Vars,   ONLY: NAnalyze
# if FV_ENABLED
    USE MOD_FV_Vars,        ONLY: FV_Elems
# endif

    IMPLICIT NONE

    REAL,INTENT(IN)    :: force 
    REAL,INTENT(INOUT) :: Ut(PP_nVar,0:PP_N,0:PP_N,0:PP_N,1:nElems) !< solution time derivative

    real,parameter :: twopi = 8. * atan (1.)
    real*8  :: coeff_real, coeff_imag
    real*8  :: accx, accy, accz, norm_acc
    integer :: i, j, k, m, sizeX, mii, iElem

    real :: coords(1:3)

    !! apply forcing
    DO iElem=1,nElems
        DO k=0,PP_N;DO j=0,PP_N;DO i=0,PP_N
            accx = 0.
            accy = 0.
            accz = 0.

# if FV_ENABLED
            IF (FV_Elems(iElem).GT.0) THEN ! FV Element
                coords(1) = NodeCoords(1,0,0,0,iElem) + (NodeCoords(1,1,1,1,iElem)-NodeCoords(1,0,0,0,iElem))/DBLE(PP_N+1) * (DBLE(i) + 0.5)
                coords(2) = NodeCoords(2,0,0,0,iElem) + (NodeCoords(2,1,1,1,iElem)-NodeCoords(2,0,0,0,iElem))/DBLE(PP_N+1) * (DBLE(j) + 0.5)
                coords(3) = NodeCoords(3,0,0,0,iElem) + (NodeCoords(3,1,1,1,iElem)-NodeCoords(3,0,0,0,iElem))/DBLE(PP_N+1) * (DBLE(k) + 0.5)
            ELSE
                coords(1:3) = Elem_xGP(1:3,i,j,k,iElem)        
            END IF
# else
            coords(1:3) = Elem_xGP(1:3,i,j,k,iElem)        
# endif

            DO m = 1, NperOUstep

                coeff_real = cos(twopi/st_Lbox * (&
                            &   st_OU_lenx(m)*coords(1) &
                            & + st_OU_leny(m)*coords(2) &
                            & + st_OU_lenz(m)*coords(3) &
                        &))

                coeff_imag = sin(twopi/st_Lbox * (&
                            &   st_OU_lenx(m)*coords(1) &
                            & + st_OU_leny(m)*coords(2) &
                            & + st_OU_lenz(m)*coords(3) &
                        &))
 
                accx = accx + 2.0 * (st_OU_kx_re(m) * coeff_real + st_OU_kx_im(m) * coeff_imag)
                accy = accy + 2.0 * (st_OU_ky_re(m) * coeff_real + st_OU_ky_im(m) * coeff_imag)
                accz = accz + 2.0 * (st_OU_kz_re(m) * coeff_real + st_OU_kz_im(m) * coeff_imag)

            ENDDO

            Ut(2,i,j,k,iElem) = Ut(2,i,j,k,iElem) + U(1,i,j,k,iElem) * force * accx
            Ut(3,i,j,k,iElem) = Ut(3,i,j,k,iElem) + U(1,i,j,k,iElem) * force * accy
            Ut(4,i,j,k,iElem) = Ut(4,i,j,k,iElem) + U(1,i,j,k,iElem) * force * accz
            Ut(5,i,j,k,iElem) = Ut(5,i,j,k,iElem) &
                            & +  U(2,i,j,k,iElem) * force * accx &
                            & +  U(3,i,j,k,iElem) * force * accy &
                            & +  U(4,i,j,k,iElem) * force * accz

        END DO; END DO; END DO
    END DO
 
END SUBROUTINE ApplyForcing


SUBROUTINE CalcMachAvg(mach_avg)
    USE MOD_PreProc
    USE MOD_Globals
    USE MOD_TestCase_Vars
    USE MOD_DG_Vars,        ONLY: U
    USE MOD_Lifting_Vars,   ONLY: GradUx,GradUy,GradUz
    USE MOD_Analyze_Vars,   ONLY: NAnalyze,Vdm_GaussN_NAnalyze,wGPVolAnalyze
    USE MOD_EOS_Vars,       ONLY: KappaM1,mu0,R,kappa
    USE MOD_Mesh_Vars,      ONLY: sJ
    USE MOD_ChangeBasis,    ONLY: ChangeBasis3D
    USE MOD_Mesh_Vars,      ONLY: nElems,nGlobalElems
    USE MOD_TimeDisc_Vars, ONLY: dt

# if FV_ENABLED
    USE MOD_FV_Vars,        ONLY: FV_Elems
# endif

#   if USE_MPI
    USE MOD_MPI_Vars
#   endif

    IMPLICIT NONE

    REAL, INTENT(OUT)  :: mach_avg !! mass-weighted global mach number

    INTEGER :: ii,i,j,k,count

    REAL    :: sJ_NAnalyze     (1,         0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: U_NAnalyze      (1:PP_nVar, 0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: sJ_N(1, 0:PP_N, 0:PP_N, 0:PP_N) ! local array for sJ

    REAL    :: volu, dens, mass, pres, rmsv, ekin, vels(1:3)

    INTEGER, PARAMETER :: SUMS_LEN = 99 !! plenty enough
    REAL    :: sums(1:SUMS_LEN)
    REAL    :: recv(1:SUMS_LEN)
    sums(:) = 0.0
    recv(:) = 0.0

    DO ii = 1,nElems

# if FV_ENABLED
        IF (FV_Elems(ii).GT.0) THEN ! FV Element
            U_NAnalyze(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze) = U(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze,ii) 
        ELSE
# endif
        ! Interpolate the jacobian to the analyze grid
        sJ_N(1,:,:,:) = sJ(:,:,:,ii,0)
        CALL ChangeBasis3D(1, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, sJ_N(1:1,0:PP_N,0:PP_N,0:PP_N), sJ_NAnalyze(1:1,:,:,:))

        ! Interpolate the solution to the analyze grid
        CALL ChangeBasis3D(PP_nVar, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, U(1:PP_nVar,:,:,:,ii), U_NAnalyze(1:PP_nVar,:,:,:))
# if FV_ENABLED
        END IF
# endif


        DO k=0,NAnalyze
            DO j=0,NAnalyze
                DO i=0,NAnalyze

# if FV_ENABLED
                    IF (FV_Elems(ii).GT.0) THEN ! FV Element
                        !! FIXME: asume unit box
                        volu = 1.0/DBLE(nGlobalElems * (NAnalyze+1)**3) 
                    ELSE
# endif
                    volu = wGPVolAnalyze(i,j,k)/sJ_NAnalyze(1,i,j,k)
# if FV_ENABLED
                    END IF
# endif
                    dens        = U_NAnalyze(1,i,j,k)
                    mass        = volu * dens
                    vels(1:3)   = U_NAnalyze(2:4,i,j,k)/dens
                    rmsv        = SUM(Vels(1:3)*Vels(1:3))
                    !!pres        = KappaM1*(U_NAnalyze(5,i,j,k) - 0.5 * dens * rmsv)
                    !!pres        = MERGE(1.e-6, pres, pres < 1.e-6)

                    sums(VOLU)  = sums(VOLU) + volu
                    sums(MASS)  = sums(MASS) + mass
                    sums(PRES)  = sums(PRES) + pres*volu
                    sums(RMSV)  = sums(RMSV) + rmsv*mass

                END DO
            END DO
        END DO
    END DO

#   if USE_MPI
    CALL MPI_ALLREDUCE(sums,recv,SUMS_LEN,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,iError)
    sums(:) = recv(:)
#   endif
    
    !!mach_avg = SQRT(sums(RMSV)/sums(MASS)/kappa/(sums(PRES)/sums(VOLU)))
    mach_avg = SQRT(sums(RMSV)/sums(MASS))

END SUBROUTINE CalcMachAvg

SUBROUTINE AnalyzeTestcase(simtime)

    USE MOD_PreProc
    USE MOD_Globals
    USE MOD_TestCase_Vars
    USE MOD_TimeDisc_Vars,  ONLY: t,dt
    USE MOD_DG_Vars,        ONLY: U
    USE MOD_Lifting_Vars,   ONLY: GradUx,GradUy,GradUz
    USE MOD_Analyze_Vars,   ONLY: NAnalyze,Vdm_GaussN_NAnalyze,wGPVolAnalyze
    USE MOD_EOS_Vars,       ONLY: KappaM1,mu0,R,kappa
    USE MOD_Mesh_Vars,      ONLY: sJ
    USE MOD_ChangeBasis,    ONLY: ChangeBasis3D
    USE MOD_Mesh_Vars,      ONLY: nElems
    USE MOD_Mesh_Vars,      ONLY: nGlobalElems
    USE Mod_ArtificialViscosity, ONLY: artvisc
    USE MOD_TimeDisc_Vars       ,ONLY: ViscousTimeStep

# if FV_ENABLED
    USE MOD_FV_Vars,        ONLY: FV_Elems
# endif

#   if USE_MPI
    USE MOD_MPI_Vars
#   endif

    IMPLICIT NONE

    REAL,INTENT(IN) :: simtime

    INTEGER :: i,j,k,ii

    REAL    :: volu, mass, ener, ekin, eint, temp
    REAL    :: dens, vels(1:3), pres, rmsv, mach

    REAL    :: pres_min, pres_max, mach_min, mach_max

    REAL    :: sJ_NAnalyze     (1,         0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: U_NAnalyze      (1:PP_nVar, 0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: GradUx_NAnalyze (1:PP_nVar, 0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: GradUy_NAnalyze (1:PP_nVar, 0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: GradUz_NAnalyze (1:PP_nVar, 0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: sJ_N(1, 0:PP_N, 0:PP_N, 0:PP_N) ! local array for sJ

    !! REAL    :: GradVel(1:3,1:3)
    !! REAL    :: vort_(1:3)
    !! REAL    :: strain(1:3,1:3)               ! Strain rate tensor S (symmetric)
    !! REAL    :: strain_devia(1:3,1:3)           ! Deviatoric part of the strain rate tensor S
    !! REAL    :: divVel                        ! Divergence of velocity vector
    !! REAL    :: eps3                          ! Integrand: p*(div u)
    !! REAL    :: u_tens, s_tens, sd_tens       ! matrix : matrix product, integrands of Gradvel, S, Sd
    !! REAL    :: DR_u,DR_S,DR_Sd,DR_p          ! Contributions to dissipation rate

    !! INTEGER :: fv_sum, fv_sum_mpi

    integer :: count = 0
    real*8, dimension(999) :: resBuf

    INTEGER, PARAMETER :: sums_LEN = 99 !! plenty enough
    REAL    :: sums(1:sums_LEN) !! accumulators
    REAL    :: recv(1:sums_LEN) !! accumulators

    sums(:) = 0.0

# if FV_ENABLED
    IF (PP_N .ne. NAnalyze) THEN
        CALL abort(__STAMP__,'Test Case Analyze Points differ from Node Points.',999,999.)
        RETURN
    END IF
# endif

    DO ii = 1,nElems

# if FV_ENABLED
        IF (FV_Elems(ii).GT.0) THEN ! FV Element
            sums(FVV) = sums(FVV) + 1
            U_NAnalyze(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze) = U(1:PP_nVar,0:NAnalyze,0:NAnalyze,0:NAnalyze,ii) 
        ELSE
# endif
        ! Interpolate the gradient of the velocity to the analyze grid
        CALL ChangeBasis3D(PP_nVar, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, GradUx(1:PP_nVar,:,:,:,ii), GradUx_NAnalyze(1:PP_nVar,:,:,:))
        CALL ChangeBasis3D(PP_nVar, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, GradUy(1:PP_nVar,:,:,:,ii), GradUy_NAnalyze(1:PP_nVar,:,:,:))
        CALL ChangeBasis3D(PP_nVar, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, GradUz(1:PP_nVar,:,:,:,ii), GradUz_NAnalyze(1:PP_nVar,:,:,:))

        ! Interpolate the jacobian to the analyze grid
        sJ_N(1,:,:,:) = sJ(:,:,:,ii,0)
        CALL ChangeBasis3D(1, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, sJ_N(1:1,0:PP_N,0:PP_N,0:PP_N), sJ_NAnalyze(1:1,:,:,:))

        ! Interpolate the solution to the analyze grid
        CALL ChangeBasis3D(PP_nVar, PP_N, NAnalyze, Vdm_GaussN_NAnalyze, U(1:PP_nVar,:,:,:,ii), U_NAnalyze(1:PP_nVar,:,:,:))
# if FV_ENABLED
        END IF
# endif

        DO k=0,NAnalyze
            DO j=0,NAnalyze
                DO i=0,NAnalyze

# if FV_ENABLED
                    IF (FV_Elems(ii).GT.0) THEN ! FV Element
                        !! FIXME: asume unit box
                        volu = 1/DBLE(nGlobalElems * (NAnalyze+1)**3) 
                    ELSE
# endif
                    volu = wGPVolAnalyze(i,j,k)/sJ_NAnalyze(1,i,j,k)
# if FV_ENABLED
                    END IF
# endif
                    dens = U_NAnalyze(1,i,j,k)
                    vels(1:3) = U_NAnalyze(2:4,i,j,k)/dens
                    rmsv = SUM(Vels(1:3)*Vels(1:3))
                    ener = U_NAnalyze(5,i,j,k)

                    ekin = 0.5 * dens * rmsv
                    eint = ener - ekin
                    eint = MERGE(1.e-6, eint, eint < 1.e-6)

                    pres = KappaM1*eint
                    temp = pres/R

                    mach = SQRT(rmsv)/SQRT(kappa*pres/dens)

                    sums(VOLU) = sums(VOLU) + volu
                    sums(MASS) = sums(MASS) + dens*volu
                    sums(PRES) = sums(PRES) + pres*volu
                    sums(RMSV) = sums(RMSV) + rmsv*volu*dens

                    sums(ENER) = sums(ENER) + ener*volu
                    sums(EINT) = sums(EINT) + eint*volu
                    sums(EKIN) = sums(EKIN) + ekin*volu
                    sums(TEMP) = sums(TEMP) + temp*volu

                    sums(PRES_MIN) = MIN(sums(PRES_MIN),pres)
                    sums(PRES_MAX) = MAX(sums(PRES_MAX),pres)
                    sums(MACH_MIN) = MIN(sums(MACH_MIN),mach)
                    sums(MACH_MAX) = MAX(sums(MACH_MAX),mach)
                    sums(RMSV_MAX) = MAX(sums(RMSV_MAX),rmsv)

                    !! GradVel(:,1) = 1./(U_NAnalyze(1,i,j,k))*(GradUx_NAnalyze(2:4,i,j,k)-Vels(1:3)*GradUx_NAnalyze(1,i,j,k))
                    !! GradVel(:,2) = 1./(U_NAnalyze(1,i,j,k))*(GradUy_NAnalyze(2:4,i,j,k)-Vels(1:3)*GradUy_NAnalyze(1,i,j,k))
                    !! GradVel(:,3) = 1./(U_NAnalyze(1,i,j,k))*(GradUz_NAnalyze(2:4,i,j,k)-Vels(1:3)*GradUz_NAnalyze(1,i,j,k))

                    !! ! divergence of velocity
                    !! divVel = GradVel(1,1) + GradVel(2,2) + GradVel(3,3)

                    !! ! compute tensor of velocity gradients
                    !! strain = 0.5*(Gradvel + TRANSPOSE(GradVel))

                    !! ! deviatoric part of strain tensor
                    !! strain_devia = strain
                    !! DO p = 1,3
                    !!     strain_devia(p,p) = strain_devia(p,p) - 1./3.*divVel
                    !! END DO

                    !! ! compute vorticity and max(vorticity)
                    !! vort_(1)  = GradVel(3,2) - GradVel(2,3)
                    !! vort_(2)  = GradVel(1,3) - GradVel(3,1)
                    !! vort_(3)  = GradVel(2,1) - GradVel(1,2)

                    !! vort        = SQRT(SUM(vort_(:)*vort_(:)))
                    !! vort_min    = MIN(vort_min,vort)
                    !! vort_max    = MAX(vort_max,vort)

                    !! ! compute enstrophy integrand
                    !! enst        = 0.5*dens * vort*vort 
                    !! enst_sum    = enst_sum + enst
                    !! ! compute integrand for epsilon3, pressure contribution to dissipation (compressiblity effect)
                    !! eps3        = pres*divVel

                    !! ! Matrix : Matrix product for velocity gradient tensor, S:S and Sd:Sd
                    !! u_tens = 0.; s_tens = 0.; sd_tens = 0.
                    !! DO p = 1,3
                    !!     DO q = 1,3
                    !!         u_tens  = u_tens  + GradVel(p,q)*GradVel(p,q)
                    !!         s_tens  = s_tens  + strain(p,q)*strain(p,q)
                    !!         sd_tens = sd_tens + strain_devia(p,q)*strain_devia(p,q)
                    !!     END DO
                    !! END DO

                    !! ! dissipation rate epsilon 1 from deviatoric part of strain rate tensor Sd (compressible)
                    !! DR_SD = DR_SD + sd_tens*volu

                    !! ! dissipation rate epsilon 3 from pressure times div u (compressible)
                    !! DR_p = DR_p + eps3*volu

                END DO
            END DO
        END DO
    END DO

#   if USE_MPI
    CALL MPI_REDUCE(sums, recv, sums_LEN, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    IF (.not. MPIRoot) RETURN
    sums(:) = recv(:)
#   endif

    !!DR_u             = DR_u*mu0/Volume
    !!DR_S             = DR_S*2.*mu0/(rho0*Volume)
    !!DR_SD            = DR_SD*2.*mu0/(rho0*Volume)
    !!DR_p             = -DR_p/(rho0*Volume)

    volu = sums(VOLU)
    mass = sums(MASS)
    pres = sums(PRES)/volu
    rmsv = sums(RMSV)/mass
                      
    ener = sums(ENER)/volu
    eint = sums(EINT)/volu
    ekin = sums(EKIN)/volu
    temp = sums(TEMP)/volu

    pres_min = sums(PRES_MIN)
    pres_max = sums(PRES_MAX)
    mach_min = sums(MACH_MIN)
    mach_max = sums(MACH_MAX)

    count = 0
    resBuf(incr(count)) = simtime
    resBuf(incr(count)) = dt
    resBuf(incr(count)) = REAL(sums(FVV)) / REAL(nGlobalElems)
    !!resBuf(incr(count)) = MERGE(1,0,viscoustimestep)
    !!resBuf(incr(count)) = artvisc%nu_max
    !!resBuf(incr(count)) = artvisc%nu_dg
    resBuf(incr(count)) = ener
    resBuf(incr(count)) = eint
    resBuf(incr(count)) = ekin
    resBuf(incr(count)) = SQRT(sums(RMSV) / sums(MASS)) 
    resBuf(incr(count)) = SQRT(sums(RMSV) / kappa / sums(PRES)) 
    resBuf(incr(count)) = mach_max
    !!resBuf(incr(count)) = sums(RMSV_MAX)**2

    call write2file(st_fp_analyze, count, resBuf)
    call writeTimestamp()

END SUBROUTINE AnalyzeTestCase

subroutine write2file(filepath, count, outBuf)

    USE MOD_Globals

    implicit none

    CHARACTER(LEN=255)              :: filepath
    integer,intent(in)              :: count
    real,DIMENSION(:),intent(in)    :: outBuf

    integer :: ioUnit
    integer :: openStat

    character(range(count)+2) :: countStr

    !!ioUnit = GETFREEUNIT()

    OPEN(NEWUNIT  = ioUnit             , &
         FILE     = TRIM(filepath)     , &
         FORM     = 'FORMATTED'        , &
         STATUS   = 'UNKNOWN'          , &
         POSITION = 'APPEND'           , &
         RECL     = 50000              , &
         IOSTAT   = openStat             )

    IF(openStat.NE.0) THEN
        CALL abort(__STAMP__, 'ERROR: cannot open '// TRIM(filepath))
    END IF

    !! write nr. of columns to string
    write(countStr,'(i0)') count

    WRITE(ioUnit, '('// trim(countStr) //'(1X,ES20.8))') outBuf(1:count)
    CLOSE(ioUnit)

end subroutine write2file

subroutine writeTimeStamp

    USE MOD_Globals
    USE MOD_TimeDisc_Vars,  ONLY: t,dt

    implicit none

    integer :: values(8)
    integer :: ioUnit
    integer :: openStat

    OPEN(NEWUNIT  = ioUnit             , &
         FILE     = 'profiling.dat'    , &
         FORM     = 'FORMATTED'        , &
         STATUS   = 'UNKNOWN'          , &
         POSITION = 'APPEND'           , &
         RECL     = 50000              , &
         IOSTAT   = openStat             )

    IF(openStat.NE.0) THEN
        CALL abort(__STAMP__, 'ERROR: cannot open "profiling.dat"')
    END IF

    call date_and_time(values=values)

    WRITE(ioUnit, '(2(ES16.9),8(I8))') t,dt,values
    CLOSE(ioUnit)

end subroutine writeTimeStamp


function incr(n)

    implicit none

    integer, intent(inout) :: n
    integer :: incr

    n = n + 1
    incr = n

end function incr

SUBROUTINE FinalizeTestcase()

    USE MOD_Globals      ,ONLY:MPIRoot
    USE MOD_TestCase_Vars,ONLY:writeBuf,Time

    IMPLICIT NONE

END SUBROUTINE

SUBROUTINE DO_NOTHING(optionalREAL,optionalREAL2)
    IMPLICIT NONE
    REAL, OPTIONAL,INTENT(IN) :: optionalREAL,optionalREAL2
END SUBROUTINE DO_NOTHING

SUBROUTINE GetBoundaryFluxTestcase(SideID,t,Nloc,Flux,UPrim_master,                   &
#if PARABOLIC
                           gradUx_master,gradUy_master,gradUz_master,&
#endif
                           NormVec,TangVec1,TangVec2,Face_xGP)
! MODULES
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN)   :: SideID  
REAL,INTENT(IN)      :: t       !< current time (provided by time integration scheme)
INTEGER,INTENT(IN)   :: Nloc    !< polynomial degree
REAL,INTENT(IN)      :: UPrim_master( PP_nVarPrim,0:Nloc,0:Nloc) !< inner surface solution
#if PARABOLIC
                                                           !> inner surface solution gradients in x/y/z-direction
REAL,INTENT(IN)      :: gradUx_master(PP_nVarPrim,0:Nloc,0:Nloc)
REAL,INTENT(IN)      :: gradUy_master(PP_nVarPrim,0:Nloc,0:Nloc)
REAL,INTENT(IN)      :: gradUz_master(PP_nVarPrim,0:Nloc,0:Nloc)
#endif /*PARABOLIC*/
                                                           !> normal and tangential vectors on surfaces
REAL,INTENT(IN)      :: NormVec (3,0:Nloc,0:Nloc)
REAL,INTENT(IN)      :: TangVec1(3,0:Nloc,0:Nloc)
REAL,INTENT(IN)      :: TangVec2(3,0:Nloc,0:Nloc)
REAL,INTENT(IN)      :: Face_xGP(3,0:Nloc,0:Nloc)    !< positions of surface flux points
REAL,INTENT(OUT)     :: Flux(PP_nVar,0:Nloc,0:Nloc)  !< resulting boundary fluxes
!==================================================================================================================================
END SUBROUTINE GetBoundaryFluxTestcase


SUBROUTINE GetBoundaryFVgradientTestcase(SideID,t,gradU,UPrim_master)
USE MOD_PreProc
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID
REAL,INTENT(IN)    :: t                                       !< current time (provided by time integration scheme)
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_N) !< primitive solution from the inside
REAL,INTENT(OUT)   :: gradU       (PP_nVarPrim,0:PP_N,0:PP_N) !< FV boundary gradient
!==================================================================================================================================
END SUBROUTINE GetBoundaryFVgradientTestcase


SUBROUTINE Lifting_GetBoundaryFluxTestcase(SideID,t,UPrim_master,Flux)
USE MOD_PreProc
!----------------------------------------------------------------------------------------------------------------------------------
! INPUT / OUTPUT VARIABLES
INTEGER,INTENT(IN) :: SideID
REAL,INTENT(IN)    :: t                                       !< current time (provided by time integration scheme)
REAL,INTENT(IN)    :: UPrim_master(PP_nVarPrim,0:PP_N,0:PP_N) !< primitive solution from the inside
REAL,INTENT(OUT)   :: Flux(        PP_nVarPrim,0:PP_N,0:PP_N) !< lifting boundary flux
!==================================================================================================================================
END SUBROUTINE Lifting_GetBoundaryFluxTestcase

subroutine OU_time_step()

!!****if* source/physics/sourceTerms/Stir/StirGirichidisMain/st_OU_step
!!
!! DESCRIPTION
!!     This routine performs an Ornstein-Uhlenbeck step for the live OU mode
!!
!! PARAMETERS
!!     see the documentation in source/physics/sourceTerms/Stir/StirGirichidisMain/manual
!!
!! AUTHOR
!!     Philipp Girichidis, July 2015
!!
!! CHANGE LOG
!!     2015-07-31: added this live OU computation
!!     2015-08-18: clean up variable names
!!
!!

    USE MOD_TimeDisc_Vars, ONLY: dt
    use MOD_testcase_vars

    implicit none

    integer     :: m
    real        :: q_re, q_im
    real*8      :: sigma, dW

    real        :: dt2 = 0.0

    if (dt .lt. 1.0) then
        dt2 = dt
    end if

    do m = 1, NperOUstep
        !   sigma  =  ( k -            kmin     ) * ( k            - kmax     )
        sigma = (( st_OU_len(m) - st_kmin ) * ( st_OU_len(m) - st_kmax ))**2

# define GGG (3./sqrt(1.-2.*st_zeta+3.*st_zeta*st_zeta))
        dW = r8_normal_01 (st_seed)
        st_OU_kx_re(m) = st_OU_kx_re(m) - (st_g0* GGG *(st_OU_kx_re(m)*dt2/st_Tac + st_F0*sqrt(2.0*sigma/st_Tac*dt2)*dW));
        dW = r8_normal_01 (st_seed)
        st_OU_kx_im(m) = st_OU_kx_im(m) - (st_g0* GGG *(st_OU_kx_im(m)*dt2/st_Tac + st_F0*sqrt(2.0*sigma/st_Tac*dt2)*dW));
        dW = r8_normal_01 (st_seed)
        st_OU_ky_re(m) = st_OU_ky_re(m) - (st_g0* GGG *(st_OU_ky_re(m)*dt2/st_Tac + st_F0*sqrt(2.0*sigma/st_Tac*dt2)*dW));
        dW = r8_normal_01 (st_seed)
        st_OU_ky_im(m) = st_OU_ky_im(m) - (st_g0* GGG *(st_OU_ky_im(m)*dt2/st_Tac + st_F0*sqrt(2.0*sigma/st_Tac*dt2)*dW));
        dW = r8_normal_01 (st_seed)
        st_OU_kz_re(m) = st_OU_kz_re(m) - (st_g0* GGG *(st_OU_kz_re(m)*dt2/st_Tac + st_F0*sqrt(2.0*sigma/st_Tac*dt2)*dW));
        dW = r8_normal_01 (st_seed)
        st_OU_kz_im(m) = st_OU_kz_im(m) - (st_g0* GGG *(st_OU_kz_im(m)*dt2/st_Tac + st_F0*sqrt(2.0*sigma/st_Tac*dt2)*dW));
# undef GGG
    enddo
    
    ! compute symmetric imaginary part for acceleration
    do m = 1, NperOUstep/2
        st_OU_kx_im(m) = -st_OU_kx_im(NperOUstep+1-m)
        st_OU_ky_im(m) = -st_OU_ky_im(NperOUstep+1-m)
        st_OU_kz_im(m) = -st_OU_kz_im(NperOUstep+1-m)
    end do
    
    ! aaply projection for sol / comp
    do m = 1, NperOUstep/2
        q_re = st_OU_lenx(m)*st_OU_kx_re(m) + st_OU_leny(m)*st_OU_ky_re(m) + st_OU_lenz(m)*st_OU_kz_re(m)
        q_im = st_OU_lenx(m)*st_OU_kx_im(m) + st_OU_leny(m)*st_OU_ky_im(m) + st_OU_lenz(m)*st_OU_kz_im(m)
        q_re = q_re / (st_OU_len(m)**2)
        q_im = q_im / (st_OU_len(m)**2)
        
        ! (a+ib)(c+id) = ac-bd + i(cb + ad)
        st_OU_kx_re(m) = st_zeta*st_OU_kx_re(m) + (1.0-2.0*st_zeta)*(st_OU_lenx(m)*q_re);
        st_OU_kx_im(m) = st_zeta*st_OU_kx_im(m) + (1.0-2.0*st_zeta)*(st_OU_lenx(m)*q_im);
        st_OU_ky_re(m) = st_zeta*st_OU_ky_re(m) + (1.0-2.0*st_zeta)*(st_OU_leny(m)*q_re);
        st_OU_ky_im(m) = st_zeta*st_OU_ky_im(m) + (1.0-2.0*st_zeta)*(st_OU_leny(m)*q_im);
        st_OU_kz_re(m) = st_zeta*st_OU_kz_re(m) + (1.0-2.0*st_zeta)*(st_OU_lenz(m)*q_re);
        st_OU_kz_im(m) = st_zeta*st_OU_kz_im(m) + (1.0-2.0*st_zeta)*(st_OU_lenz(m)*q_im);
    end do
  
end subroutine OU_time_step
     
function r8_uniform_01 ( seed )

!*****************************************************************************80
!
!! R8_UNIFORM_01 returns a unit pseudorandom R8.
!
!  Discussion:
!
!    This routine implements the recursion
!
!      seed = 16807 * seed mod ( 2^31 - 1 )
!      r8_uniform_01 = seed / ( 2^31 - 1 )
!
!    The integer arithmetic never requires more than 32 bits,
!    including a sign bit.
!
!    If the initial seed is 12345, then the first three computations are
!
!      Input     Output      R8_UNIFORM_01
!      SEED      SEED
!
!         12345   207482415  0.096616
!     207482415  1790989824  0.833995
!    1790989824  2035175616  0.947702
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    31 May 2007
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Second Edition,
!    Springer, 1987,
!    ISBN: 0387964673,
!    LC: QA76.9.C65.B73.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, December 1986, pages 362-376.
!
!    Pierre L'Ecuyer,
!    Random Number Generation,
!    in Handbook of Simulation,
!    edited by Jerry Banks,
!    Wiley, 1998,
!    ISBN: 0471134031,
!    LC: T57.62.H37.
!
!    Peter Lewis, Allen Goodman, James Miller,
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, 1969, pages 136-143.
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which
!    should NOT be 0.
!    On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R8_UNIFORM_01, a new pseudorandom variate,
!    strictly between 0 and 1.
!
  implicit none

  integer ( kind = 4 ) k
  real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed

  k = seed / 127773

  seed = 16807 * ( seed - k * 127773 ) - k * 2836

  if ( seed < 0 ) then
    seed = seed + 2147483647
  end if
!
!  Although SEED can be represented exactly as a 32 bit integer,
!  it generally cannot be represented exactly as a 32 bit real number!
!
  r8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10

  return
end function r8_uniform_01

function r8_normal_01 ( seed )

!*****************************************************************************80
!
!! R8_NORMAL_01 returns a unit pseudonormal R8.
!
!  Discussion:
!
!    The standard normal probability distribution function (PDF) has
!    mean 0 and standard deviation 1.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 August 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random
!    number generator.
!
!    Output, real ( kind = 8 ) R8_NORMAL_01, a normally distributed
!    random value.
!
  implicit none

  real ( kind = 8 ) r1
  real ( kind = 8 ) r2
  real ( kind = 8 ) r8_normal_01
  real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
  !real ( kind = 8 ) r8_uniform_01
  integer ( kind = 4 ) seed
  real ( kind = 8 ) x

  r1 = r8_uniform_01 ( seed )
  r2 = r8_uniform_01 ( seed )
  x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )

  r8_normal_01 = x

  return

end function r8_normal_01

!! function xr8_uniform_01 ( seed )
!!   implicit none
!!   integer ( kind = 4 ) k
!!   real ( kind = 8 ) xr8_uniform_01
!!   integer ( kind = 4 ) seed
!!   k = seed / 127773
!!   seed = 16807 * ( seed - k * 127773 ) - k * 2836
!!   if ( seed < 0 ) then
!!     seed = seed + 2147483647
!!   end if
!!   xr8_uniform_01 = real ( seed, kind = 8 ) * 4.656612875D-10
!!   return
!! end function xr8_uniform_01
!! 
!! function xr8_normal_01 ( seed )
!!   implicit none
!!   real ( kind = 8 ) r1
!!   real ( kind = 8 ) r2
!!   real ( kind = 8 ) xr8_normal_01
!!   real ( kind = 8 ), parameter :: r8_pi = 3.141592653589793D+00
!!   real ( kind = 8 ) xr8_uniform_01
!!   integer ( kind = 4 ) seed
!!   real ( kind = 8 ) x
!!   r1 = xr8_uniform_01 ( seed )
!!   r2 = xr8_uniform_01 ( seed )
!!   x = sqrt ( - 2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * r8_pi * r2 )
!!   xr8_normal_01 = x
!!   return
!! end function xr8_normal_01

END MODULE MOD_Testcase
