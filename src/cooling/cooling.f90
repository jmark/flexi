#include "flexi.h"

MODULE MOD_Cooling

IMPLICIT NONE

TYPE, PRIVATE :: PARAMS

    LOGICAL :: enabled          = .FALSE.
    REAL    :: polytropicConst  = 0.6

    CHARACTER(LEN=255) :: outputfp = 'cooling.dat'
END TYPE PARAMS

TYPE(PARAMS) cooling

CONTAINS

SUBROUTINE DefineParametersCooling

    USE MOD_Globals
    USE MOD_ReadInTools ,ONLY: prms

    CALL prms%SetSection("Cooling")
    CALL prms%CreateLogicalOption('cooling%enabled', "dito", value='F')
    CALL prms%CreateRealOption('cooling%polytropicConst', "dito", value='0.6')

END SUBROUTINE

SUBROUTINE InitCooling

    USE MOD_Globals
    USE MOD_PreProc
    USE MOD_ReadInTools
    !!USE MOD_Cooling_vars

    IMPLICIT NONE

    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' INIT Cooling ...'

    cooling%enabled = GETLOGICAL('cooling%enabled','F')
    cooling%polytropicConst = GETREAL('cooling%polytropicConst','0.6')

    SWRITE(UNIT_stdOut,'(A)')' INIT Cooling DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE

SUBROUTINE ApplyCooling

    USE MOD_PreProc
    USE MOD_Globals
    USE MOD_TestCase_Vars
    USE MOD_TimeDisc_Vars,  ONLY: t,dt
    USE MOD_DG_Vars,        ONLY: U
    USE MOD_Lifting_Vars,   ONLY: GradUx,GradUy,GradUz
    USE MOD_Analyze_Vars,   ONLY: NAnalyze,Vdm_GaussN_NAnalyze,wGPVolAnalyze
    USE MOD_EOS_Vars,       ONLY: KappaM1,mu0,R,kappa,sKappaM1
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

    INTEGER :: i,j,k,ii

    REAL    :: volu, mass, ener, ekin, eint, temp
    REAL    :: dens, vels(1:3), pres, rmsv, mach

    REAL    :: sJ_NAnalyze     (1,         0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: U_NAnalyze      (1:PP_nVar, 0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: GradUx_NAnalyze (1:PP_nVar, 0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: GradUy_NAnalyze (1:PP_nVar, 0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: GradUz_NAnalyze (1:PP_nVar, 0:NAnalyze, 0:NAnalyze, 0:NAnalyze)
    REAL    :: sJ_N(1, 0:PP_N, 0:PP_N, 0:PP_N) ! local array for sJ

    integer :: count = 0
    real*8, dimension(999) :: resBuf

    INTEGER, PARAMETER :: IDX_DENS      = 1
    INTEGER, PARAMETER :: IDX_PRES      = 2
    INTEGER, PARAMETER :: IDX_TEMP      = 3
    INTEGER, PARAMETER :: IDX_ENER      = 4
    INTEGER, PARAMETER :: IDX_EINT      = 5
    INTEGER, PARAMETER :: IDX_EKIN      = 6
    INTEGER, PARAMETER :: IDX_RMSV      = 7
    INTEGER, PARAMETER :: IDX_MASS      = 8
    INTEGER, PARAMETER :: IDX_VOLU      = 9
    INTEGER, PARAMETER :: IDX_FVV       = 10
    INTEGER, PARAMETER :: IDX_PRES_MIN  = 11
    INTEGER, PARAMETER :: IDX_PRES_MAX  = 12
    INTEGER, PARAMETER :: IDX_MACH_MIN  = 13
    INTEGER, PARAMETER :: IDX_MACH_MAX  = 14
    INTEGER, PARAMETER :: IDX_RMSV_MAX  = 15
    INTEGER, PARAMETER :: SUMS_LEN  = 20

    REAL :: sums(1:SUMS_LEN) !! accumulators
    REAL :: recv(1:SUMS_LEN) !! accumulators

# if FV_ENABLED
    IF (PP_N .ne. NAnalyze) THEN
        CALL abort(__STAMP__,'Cooling Analyze Points differ from Node Points.',999,999.)
        RETURN
    END IF
# endif

    !! ========================================================================= !!
    !! Initialize Result Buffer

    count = 0
    resBuf(incr(count)) = t
    resBuf(incr(count)) = dt

    !! ========================================================================= !!
    !! Energy before Cooling

    sums(:) = 0.0
    recv(:) = 0.0

    DO ii = 1,nElems

# if FV_ENABLED
        IF (FV_Elems(ii).GT.0) THEN ! FV Element
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
                        volu = 1/DBLE(nElems * (NAnalyze+1)**3) 
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

                    sums(IDX_VOLU) = sums(IDX_VOLU) + volu
                    sums(IDX_MASS) = sums(IDX_MASS) + dens*volu

                    sums(IDX_ENER) = sums(IDX_ENER) + ener*volu
                    sums(IDX_EINT) = sums(IDX_EINT) + eint*volu
                    sums(IDX_EKIN) = sums(IDX_EKIN) + ekin*volu

                END DO
            END DO
        END DO
    END DO

#   if USE_MPI
    CALL MPI_REDUCE(sums, recv, sums_LEN, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    sums(:) = recv(:)
#   endif

    volu = sums(IDX_VOLU)
    mass = sums(IDX_MASS)
    ener = sums(IDX_ENER)/volu
    eint = sums(IDX_EINT)/volu
    ekin = sums(IDX_EKIN)/volu

    resBuf(incr(count)) = ener
    resBuf(incr(count)) = eint
    resBuf(incr(count)) = ekin

    !! ========================================================================= !!
    !! Apply isothermal polytropic process: Cooling

    DO ii = 1,nElems
        DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
            U(5,i,j,k,ii) = sKappaM1 * cooling%polytropicConst * U(1,i,j,k,ii) + 0.5 * SUM(U(2:4,i,j,k,ii)*U(2:4,i,j,k,ii)/U(1,i,j,k,ii))
        END DO; END DO; END DO
    END DO

    !! ========================================================================= !!
    !! Energy after Cooling

    sums(:) = 0
    recv(:) = 0

    DO ii = 1,nElems

# if FV_ENABLED
        IF (FV_Elems(ii).GT.0) THEN ! FV Element
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
                        volu = 1/DBLE(nElems * (NAnalyze+1)**3) 
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

                    sums(IDX_VOLU) = sums(IDX_VOLU) + volu
                    sums(IDX_MASS) = sums(IDX_MASS) + dens*volu
                    sums(IDX_ENER) = sums(IDX_ENER) + ener*volu
                    sums(IDX_EINT) = sums(IDX_EINT) + eint*volu
                    sums(IDX_EKIN) = sums(IDX_EKIN) + ekin*volu

                END DO
            END DO
        END DO
    END DO

#   if USE_MPI
    CALL MPI_REDUCE(sums, recv, sums_LEN, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
    IF (.not. MPIRoot) RETURN
    sums(:) = recv(:)
#   endif

    volu = sums(IDX_VOLU)
    mass = sums(IDX_MASS)
    ener = sums(IDX_ENER)/volu
    eint = sums(IDX_EINT)/volu
    ekin = sums(IDX_EKIN)/volu

    resBuf(incr(count)) = ener
    resBuf(incr(count)) = eint
    resBuf(incr(count)) = ekin

    call write2file(count, resBuf)

END SUBROUTINE

subroutine write2file(count, outBuf)

    USE MOD_Globals

    implicit none

    integer,intent(in)              :: count
    real,DIMENSION(:),intent(in)    :: outBuf

    integer :: ioUnit
    integer :: openStat

    character(range(count)+2) :: countStr

    !!ioUnit = GETFREEUNIT()

    OPEN(NEWUNIT  = ioUnit             , &
         FILE     = TRIM(cooling%outputfp)  , &
         FORM     = 'FORMATTED'        , &
         STATUS   = 'UNKNOWN'          , &
         POSITION = 'APPEND'           , &
         RECL     = 50000              , &
         IOSTAT   = openStat             )

    IF(openStat.NE.0) THEN
        CALL abort(__STAMP__, 'ERROR: cannot open '// TRIM(cooling%outputfp))
    END IF

    !! write nr. of columns to string
    write(countStr,'(i0)') count

    WRITE(ioUnit, '('// trim(countStr) //'(1X,ES20.8))') outBuf(1:count)
    CLOSE(ioUnit)

end subroutine write2file

function incr(n)

    implicit none

    integer, intent(inout) :: n
    integer :: incr

    n = n + 1
    incr = n

end function incr


SUBROUTINE FinalizeCooling

    IMPLICIT NONE

END SUBROUTINE

END MODULE MOD_Cooling
