#include "flexi.h"

MODULE MOD_Cooling

IMPLICIT NONE

TYPE, PRIVATE :: PARAMS

    LOGICAL :: enabled          = .FALSE.
    REAL    :: polytropicConst  = 0.6

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
   
    USE MOD_DG_Vars,        ONLY: U
    USE MOD_Preproc,        ONLY: N
    USE MOD_Mesh_Vars,      ONLY: nElems
    USE MOD_EOS_Vars,       ONLY: sKappaM1

    IMPLICIT NONE

    INTEGER :: ii,i,j,k

    !! apply isothermal polytropic process
    DO ii = 1,nElems
        DO k=0,PP_N; DO j=0,PP_N; DO i=0,PP_N
            U(5,i,j,k,ii) = sKappaM1 * cooling%polytropicConst * U(1,i,j,k,ii) + 0.5 * SUM(U(2:4,i,j,k,ii)*U(2:4,i,j,k,ii)/U(1,i,j,k,ii))
        END DO; END DO; END DO
    END DO

END SUBROUTINE

SUBROUTINE FinalizeCooling

    IMPLICIT NONE

END SUBROUTINE

END MODULE MOD_Cooling
