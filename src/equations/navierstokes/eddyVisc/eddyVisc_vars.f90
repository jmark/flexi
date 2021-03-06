!=================================================================================================================================
! Copyright (c) 2010-2016  Prof. Claus-Dieter Munz 
! This file is part of FLEXI, a high-order accurate framework for numerically solving PDEs with discontinuous Galerkin methods.
! For more information see https://www.flexi-project.org and https://nrg.iag.uni-stuttgart.de/
!
! FLEXI is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License 
! as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! FLEXI is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
! of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License v3.0 for more details.
!
! You should have received a copy of the GNU General Public License along with FLEXI. If not, see <http://www.gnu.org/licenses/>.
!=================================================================================================================================

!===================================================================================================================================
!> Contains the parameters needed for eddy viscosity models
!===================================================================================================================================
MODULE MOD_EddyVisc_Vars
! MODULES
IMPLICIT NONE
PUBLIC
SAVE

ABSTRACT INTERFACE
  SUBROUTINE EddyViscInt(grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho,iElem,i,j,k,muSGS)
  INTEGER,INTENT(IN)  :: iElem  !< index of current element
  !> indices of the c
  INTEGER,INTENT(IN)  :: i,j,k
  !> gradients of the directions
  REAL,INTENT(IN)     :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32
  REAL,INTENT(IN)     :: rho    !< Density
  REAL,INTENT(INOUT)  :: muSGS  !< local SGS viscosity
  END SUBROUTINE
END INTERFACE

ABSTRACT INTERFACE
  SUBROUTINE EddyVisc_surfInt(grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32,rho,DeltaSS,SGS_Ind,muSGS,Face_xGP)
  !> gradients of the velocities w.r.t. all directions
  REAL,INTENT(IN)   :: grad11,grad22,grad33,grad12,grad13,grad21,grad23,grad31,grad32
  REAL,INTENT(IN)   :: rho      !< Density
  REAL,INTENT(IN)   :: DeltaSS  !< Filter width
  REAL,INTENT(IN)   :: SGS_Ind  !< Indicator for SGS model
  REAL,INTENT(IN)   :: Face_xGP !< Coordinate for van-Driest damping
  REAL,INTENT(OUT)  :: muSGS    !< local SGS viscosity
  END SUBROUTINE
END INTERFACE

!-----------------------------------------------------------------------------------------------------------------------------------
! GLOBAL VARIABLES 
!-----------------------------------------------------------------------------------------------------------------------------------
INTEGER                             :: eddyViscType           !< type of eddy viscosity
PROCEDURE(EddyViscInt)     ,POINTER :: eddyViscosity          !< pointer to routine for computing volume eddy viscosity
PROCEDURE(EddyVisc_surfInt),POINTER :: eddyViscosity_surf     !< pointer to routine for computing surface eddy viscosity

!Smagosinsky Standard
REAL,ALLOCATABLE  :: DeltaS(:)         !< filter width, used by Smagorinsky modell
REAL,ALLOCATABLE  :: DeltaS_master(:)
REAL,ALLOCATABLE  :: DeltaS_slave(:)
REAL,ALLOCATABLE  :: muSGS(:,:,:,:,:)  !< Viscosity for the sub-grid
REAL,ALLOCATABLE  :: muSGSmax(:)       !< Viscosity for the sub-grid
REAL              :: CS                !< Smagorinsky constant, LES
REAL              :: PrSGS             !< Prandtl number for the sub-grid scales
!TKECUbed
REAL,ALLOCATABLE  :: SGS_Ind(:,:,:,:,:) 
REAL,ALLOCATABLE  :: SGS_Ind_master(:,:,:,:) 
REAL,ALLOCATABLE  :: SGS_Ind_slave(:,:,:,:) 

LOGICAL           :: VanDriest=.FALSE.
LOGICAL           :: SmagorinskyInitIsDone=.FALSE.

END MODULE MOD_EddyVisc_Vars
