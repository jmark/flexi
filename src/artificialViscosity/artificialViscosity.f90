#include "flexi.h"

MODULE MOD_ArtificialViscosity

IMPLICIT NONE

TYPE, PRIVATE :: VARS

    INTEGER :: elemid, nbelemid, indvar
    LOGICAL :: enabled = .FALSE.
    REAL    :: nu_max, nu_dg, nu_fv
    REAL    :: eps0, eta_min, eta_max
    REAL    :: mu_fv

    REAL, ALLOCATABLE, PUBLIC  :: nu(:)         !! artificial viscosity for each element

END TYPE VARS 

TYPE(VARS), PUBLIC :: artvisc

REAL, ALLOCATABLE, PRIVATE :: sVdm_Leg(:,:)     !! 1D inverse Vandermondematrix to Legendre polynomials

CONTAINS

SUBROUTINE DefineParametersArtificialViscosity

    USE MOD_Globals
    USE MOD_ReadInTools ,ONLY: prms

    IMPLICIT NONE

    CALL prms%SetSection("Artificial Viscosity")

    CALL prms%CreateLogicalOption('artvisc%enabled', "Apply shock capturing", 'F')
    CALL prms%CreateIntOption('artvisc%indvar',   "Indicator Variable",    '5')
    CALL prms%CreateRealOption('artvisc%eps0',    "Shock Capturing Base Foo",    '0.1')
    CALL prms%CreateRealOption('artvisc%eta_min', "Shock Capturing Foo",    '-8.0')
    CALL prms%CreateRealOption('artvisc%eta_max', "Shock Capturing Foo",    '-3.0')
    CALL prms%CreateRealOption('artvisc%mu_fv',   "Artificial Viscosity in FV elements",    '0.1')

END SUBROUTINE DefineParametersArtificialViscosity


SUBROUTINE InitArtificialViscosity()

    USE MOD_Globals
    USE MOD_PreProc
    USE MOD_ReadInTools
    USE MOD_Mesh_Vars, ONLY: nElems
    USE MOD_Basis, ONLY: buildLegendreVdm
    USE MOD_Interpolation_Vars, ONLY: xGP, InterpolationInitIsDone

    IMPLICIT NONE

    REAL, DIMENSION(0:PP_N,0:PP_N) :: Vdm_Leg !! dummy array

    SWRITE(UNIT_StdOut,'(132("-"))')
    SWRITE(UNIT_stdOut,'(A)') ' INIT ArtificialViscosity...'

    artvisc%enabled   = GETLOGICAL('artvisc%enabled','.FALSE.') 
    artvisc%indvar    = GETINT('artvisc%indvar','5') 
    artvisc%eps0      = GETREAL('artvisc%eps0', '0.1')
    artvisc%eta_min   = GETREAL('artvisc%eta_min', '-8.0')
    artvisc%eta_max   = GETREAL('artvisc%eta_max', '-3.0')
    artvisc%mu_fv     = GETREAL('artvisc%mu_fv', '0.1')

    IF (PP_N.LT.2.AND.artvisc%enabled) THEN
      CALL abort(__STAMP__,'Polynomial Degree too small for Shock Capturing!',999,999.)
      RETURN
    END IF

    ALLOCATE(sVdm_Leg(0:PP_N,0:PP_N))
    CALL buildLegendreVdm(PP_N,xGP,Vdm_Leg,sVdm_Leg)

    ALLOCATE(artvisc%nu(nElems))

    artvisc%nu(:)  = 0.0
    artvisc%nu_max = 0.0

    SWRITE(UNIT_stdOut,'(A)')' INIT ArtificialViscosity DONE!'
    SWRITE(UNIT_StdOut,'(132("-"))')

END SUBROUTINE InitArtificialViscosity


SUBROUTINE CalcArtificialViscosity(U)
    !! Use framework of Persson and Peraire to measure shocks with DOF energy
    !! indicator and calculate artificial viscosity, if necessary

    USE MOD_Globals
    USE MOD_PreProc
    USE MOD_EOS,                ONLY: ConsToPrim
    USE MOD_EOS_Vars,           ONLY: KappaM1,kappa
    USE MOD_ChangeBasis,        ONLY: ChangeBasis3D
    USE MOD_Mesh_Vars,          ONLY: sJ,Metrics_fTilde,Metrics_gTilde,Metrics_hTilde,nElems
# if FV_ENABLED
    USE MOD_FV_Vars,        ONLY: FV_Elems
# endif

    IMPLICIT NONE

    REAL, DIMENSION(PP_nVar,0:PP_N,0:PP_N,0:PP_N,nElems), INTENT(IN)  :: U

    REAL, DIMENSION(0:PP_N,0:PP_N,0:PP_N) :: Uind,Umod
    REAL                                  :: LU,LUM1,LUM2,LU_N,LU_NM1,eta_dof,eta_min,eta_max,eps0!,c2,ca2,va2,astar,cf,lambda
    REAL                                  :: v(3),Prim(PP_nVarPrim),cf,Max_Lambda(6),lambda_max,h,lambda_max2
    INTEGER                               :: ii,i,j,k
    REAL*8 :: rmsv

    artvisc%nu_fv   = 0.0
    artvisc%nu_dg   = 0.0
    artvisc%nu_max  = 0.0
    artvisc%nu(:)   = 0.0

    !! artificial Viscosity
    eps0    = artvisc%eps0
    eta_min = artvisc%eta_min
    eta_max = artvisc%eta_max

    DO ii = 1, nElems

# if FV_ENABLED
        IF (FV_Elems(ii).GT.0) then !! FV Element

            rmsv = 0.0
            DO k=1,PP_N; DO j=1,PP_N; DO i=1,PP_N
                rmsv = MAX(rmsv, SUM((U(2:4,i,j,k,ii)/U(1,i,j,k,ii))**2))
            END DO; END DO; END DO

            artvisc%nu(ii) = artvisc%mu_fv * rmsv
            artvisc%nu_fv  = MAX(artvisc%nu_fv, artvisc%nu(ii))
        
        ELSE
# endif
        !! Choose indicator (density: ind=1, pressure: ind=2)
        SELECT CASE(artvisc%indvar)
            CASE(1:5)
                Uind = U(artvisc%indvar,:,:,:,ii)
            CASE(6)
                Uind = KappaM1*(U(5,:,:,:,ii) - 0.5*(SUM(U(2:4,:,:,:,ii)*U(2:4,:,:,:,ii))/U(1,:,:,:,ii))) !! -s2mu_0*SUM(U(6:8,:,:,:,ii)*U(6:8,:,:,:,ii)))
        END SELECT
      
        !! Transform Uind into modal Legendre interpolant Umod
        CALL ChangeBasis3D(PP_N,PP_N,sVdm_Leg,Uind,Umod)

        !! Compute (truncated) error norms
        LU      = SUM(Umod(:,:,:)**2)
        LUM1    = SUM(Umod(0:PP_N-1,0:PP_N-1,0:PP_N-1)**2)
        LUM2    = SUM(Umod(0:PP_N-2,0:PP_N-2,0:PP_N-2)**2)
        LU_N    = LU-LUM1
        LU_NM1  = LUM1-LUM2

        !! DOF energy indicator
        eta_dof = LOG10(MAX( LU_N/(LU + EPSILON(0.)), LU_NM1/(LUM1 + EPSILON(0.)) ) + EPSILON(0.))

        IF (eta_dof.GE.eta_max) THEN
            artvisc%nu(ii) = eps0
        ELSE IF (eta_dof.LE.eta_min) THEN
            artvisc%nu(ii) = 0.
        ELSE
            artvisc%nu(ii) = 0.5*eps0*(1.0+SIN(PP_Pi*(eta_dof-0.5*(eta_max+eta_min))/(eta_max-eta_min)))
        END IF
       
        !! Get (transformed) max eigenvalue:
        !! Max_Lambda = 0.0
        !! DO i=1,PP_N; DO j=1,PP_N; DO k=1,PP_N

        !!     CALL ConsToPrim(prim,U(:,i,j,k,ii))
        !!     v(:) = prim(2:4) 
        !!     cf   = SQRT(kappa*prim(5)/prim(1)) !! NOTE: only for euler equation

        !!     Max_Lambda(1) = MAX(Max_Lambda(1),sJ(i,j,k,ii,0)*(ABS(SUM(Metrics_fTilde(:,i,j,k,ii,0)*v)) + &
        !!                      cf*SQRT(SUM(Metrics_fTilde(:,i,j,k,ii,0)*Metrics_fTilde(:,i,j,k,ii,0)))))
        !!     Max_Lambda(2) = MAX(Max_Lambda(2),sJ(i,j,k,ii,0)*(ABS(SUM(Metrics_gTilde(:,i,j,k,ii,0)*v)) + &
        !!                      cf*SQRT(SUM(Metrics_gTilde(:,i,j,k,ii,0)*Metrics_gTilde(:,i,j,k,ii,0)))))
        !!     Max_Lambda(3) = MAX(Max_Lambda(3),sJ(i,j,k,ii,0)*(ABS(SUM(Metrics_hTilde(:,i,j,k,ii,0)*v)) + &
        !!                      cf*SQRT(SUM(Metrics_hTilde(:,i,j,k,ii,0)*Metrics_hTilde(:,i,j,k,ii,0)))))

        !!     Max_Lambda(4) = MAX(Max_Lambda(4),ABS(v(1)) + cf)
        !!     Max_Lambda(5) = MAX(Max_Lambda(5),ABS(v(2)) + cf)
        !!     Max_Lambda(6) = MAX(Max_Lambda(6),ABS(v(3)) + cf)

        !! END DO; END DO; END DO

        !! v(1) = MAXVAL(Metrics_fTilde(:,:,:,:,ii,0))
        !! v(2) = MAXVAL(Metrics_gTilde(:,:,:,:,ii,0))
        !! v(3) = MAXVAL(Metrics_hTilde(:,:,:,:,ii,0))

        !! !! eps0 = 1.0/SQRT(MAXVAL(v))
        !! !! h    = 2.0*eps0/MINVAL(sJ(:,:,:,ii,0))
        !! !! lambda_max = MAXVAL(Max_Lambda(1:3))*h

        !! h    = (8.0/MINVAL(sJ(:,:,:,ii,0)))**(1.0/3.0)
        !! lambda_max2 = MAXVAL(Max_Lambda(4:6))*h
      
        !! ! print*,lambda_max
        !! ! print*,lambda_max2

        !! ! c2=kappa*Uind(i,j,k)/U(1,i,j,k,l)
        !! ! ca2=U(6,i,j,k,l)*U(6,i,j,k,l)/U(1,i,j,k,l) !Alfen wave speed
        !! ! va2=SUM(U(7:8,i,j,k,l)*U(7:8,i,j,k,l))/U(1,i,j,k,l)+ca2
        !! ! astar=SQRT((c2+va2)*(c2+va2)-4.*c2*ca2)
        !! ! cf=SQRT(0.5*(c2+va2+astar))
        !! ! lambda = max(lambda,abs(U(2,i,j,k,l)/U(1,i,j,k,l))+cf)
        !! ! nu(l) = 0.05*lambda/(REAL(PP_N))*nu(l)

        !! !! ! Scaling of artificial viscosity
        !! nu(ii) = nu(ii)*lambda_max2/(REAL(PP_N))
        !! nu_dg  = nu(ii)*lambda_max2/(REAL(PP_N))

        !! !print*,l
        !! !print*,h
        !! !print*,lambda_max2
        !! !print*,nu(l)

        !! nu_dg = MAX(nu_dg, nu(ii))

# if FV_ENABLED
        END IF
# endif

        !! Save max artificial viscosity for DFL timestepping
        artvisc%nu_max = MAX(artvisc%nu_max,artvisc%nu(ii))

    END DO !! loop over elements
    
    !!SWRITE(UNIT_stdOut,'(A,5(1X,ES20.8))') 'ArtificialViscosity: ', t, dt, nu_max 

END SUBROUTINE CalcArtificialViscosity


SUBROUTINE FinalizeArtificialViscosity()

    IMPLICIT NONE
END SUBROUTINE FinalizeArtificialViscosity

END MODULE MOD_ArtificialViscosity
