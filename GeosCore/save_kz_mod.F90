! ---------------------------------------------------------------
! Compute and save out Kz values 
! Karen Yu
! 21 June 2016
! ----------------------------------------------------------------
  MODULE SAVE_KZ_MOD
! !USES:
!
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC :: WRITE_WC
!
! !REVISION HISTORY:
! 21 June 2016 - K. Yu - Initial version
!EOP
! --------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
  INTEGER, PARAMETER    :: dp = kind(1.0d0)
  INTEGER, PARAMETER    :: sp = kind(1.0)
! 
CONTAINS
!EOC
! ---------------------------------------------------------------
!BOP
! !ROUTINE: WRITE_WC
!
! !DESCRIPTION:
!
! !INTERFACE:
!
  SUBROUTINE WRITE_WC( am_I_Root, Input_Opt, State_Chm, RC)
!
! !USES:
!
  USE CMN_SIZE_MOD
  USE GIGC_ErrCode_Mod
  USE GIGC_Input_Opt_Mod, ONLY : OptInput
  USE GIGC_State_Chm_Mod, ONLY : ChmState
  USE Ncdf_Mod,           ONLY : NC_Create, NC_Close
  USE Ncdf_Mod,           ONLY : Nc_Var_Def, Nc_Var_Write
  USE DIAG_MOD,           ONLY : OMEGAINST
  USE HCO_State_MOd,      ONLY : HCO_State
  USE HCOI_GC_MAIN_MOD,   ONLY : GetHcoState
  USE ERROR_MOD,          ONLY : ERROR_STOP
  USE TIME_MOD,           ONLY : GET_NYMD, GET_NHMS
! !INPUT PARAMETERS:
!
  LOGICAL,          INTENT(IN)    :: am_I_Root
  TYPE(OPtInput),   INTENT(IN)    :: Input_Opt
  TYPE(ChmState),   INTENT(INOUT) :: State_Chm
! !OUTPUT PARAMETERS
  INTEGER,          INTENT(INOUT) :: RC
! !REVISION HISTORY:
! 21 June 2016 - K. Yu - Initial version
!
! !LOCAL VARIABLES
! 
  REAL(sp), POINTER     :: Arr1D(:)
  INTEGER,  POINTER     :: Int1D(:)
  REAL(sp), POINTER     :: Arr3D(:,:,:)
  INTEGER               :: fId, lonId, latId, levId, timeId
  INTEGER               :: VarCt
  INTEGER               :: nTime
  TYPE(HCO_STATE), POINTER :: HcoState => NULL()
  LOGICAL, SAVE         :: FIRST = .TRUE.
  INTEGER               :: I
  REAL(dp)              :: xmid(IIPAR)
  REAL(dp)              :: ymid(JJPAR)
  CHARACTER(LEN=255)    :: fname
  CHARACTER(LEN=255)    :: ffname
  INTEGER               :: DDATE
!
!====================================================================
! Code starts here
!====================================================================

  nTime = 1

  ! Get HcoState
  IF (FIRST) THEN
    CALL GetHcoState( HcoState )
    IF ( .NOT. ASSOCIATED( HcoState )) &
        CALL ERROR_STOP ( 'Cannot get HcoState', 'Save Kz' )
  
    xmid = HcoState%Grid%Xmid%Val(:,1)
    ymid = HcoState%Grid%ymid%Val(:,1)
    HcoState => NULL()
    FIRST = .FALSE.
  ENDIF

  ! construct name of output file
  DDATE = GET_NYMD()
  WRITE(fname, '(i8)') DDATE
  DDATE = GET_NHMS()
  WRITE(ffname, '(i6.6)') DDATE
  fname = TRIM(fname) // TRIM(ffname)
  fname = TRIM(fname) // '.nc'
  fname = '/n/regal/jacob_lab/kyu/gc_wc/' // TRIM(fname)
  ! Create output file
  CALL NC_CREATE(TRIM(fname), IIPAR, JJPAR, LLPAR, nTime, &
                 fId,    lonId, latId, levId, timeId, VarCt )

  ! Add logitude
  CALL NC_VAR_DEF( fId, lonId, -1, -1, -1, &
                   'lon', 'Longitude', 'degrees_east', 4, VarCt )
  ALLOCATE( Arr1d( IIPAR ) )
  Arr1D = xmid 
  CALL NC_VAR_WRITE( fId, 'lon', Arr1D=Arr1D )
  DEALLOCATE( Arr1D )

  ! Add latitude
  CALL NC_VAR_DEF( fId, -1, latId, -1, -1, &
                   'lat', 'Latitude', 'degrees_north', 4, VarCt )
  ALLOCATE( Arr1D( JJPAR ) )
  Arr1D = ymid
  CALL NC_VAR_WRITE( fId, 'lat', Arr1D=Arr1D )
  DEALLOCATE( Arr1D )

  ! Add level
  CALL NC_VAR_DEF( fId, -1, levId, -1, -1, &
                   'lev', 'GEOS-Chem level', 'unitless', 1, VarCt )
  ALLOCATE(Int1D(LLPAR))
  DO I = 1, LLPAR 
    INT1D(I) = I
  ENDDO
  CALL NC_VAR_WRITE (fId, 'lev', Arr1D=Int1D )
  DEALLOCATE(Int1D)

  ! Write out Rn
  ALLOCATE( Arr3D(IIPAR, JJPAR, LLPAR) )
  Arr3D = State_Chm%Tracers(:,:,:,1)
  CALL NC_VAR_DEF( fId, lonId, latId, levId, -1, &
                   'Rn', 'Rn222 mixing ratio', 'v/v', 4, VarCt )
  CALL NC_VAR_WRITE( fId, 'Rn', Arr3D=Arr3D)
  DEALLOCATE( Arr3D )

  ! Write out Pb
  ALLOCATE( Arr3D( IIPAR, JJPAR, LLPAR ) )
  Arr3D = State_Chm%Tracers(:,:,:,2)
  CALL NC_VAR_DEF( fId, lonId, latId, levId, -1, &
                   'Pb', 'Pb210 mixing ratio', 'v/v', 4, VarCt )
  CALL NC_VAR_WRITE( fId, 'Pb', Arr3D=Arr3D )
  DEALLOCATE( Arr3D )

  ! Write out Be
  ALLOCATE( Arr3D(IIPAR, JJPAR, LLPAR ) )
  Arr3D = State_Chm%Tracers(:,:,:,3)
  CALL NC_VAR_DEF( fId, lonId, latId, levId, -1, &
                   'Be', 'Be7 mixing ratio', 'v/v', 4, VarCt )
  CALL NC_VAR_WRITE( fId, 'Be', Arr3D=Arr3D )
  DEALLOCATE( Arr3D )

  ! Write out omega
  ALLOCATE( Arr3D( IIPAR, JJPAR, LLPAR ) )
  Arr3D = OMEGAINST
  CALL NC_VAR_DEF( fId, lonId, latId, levId, -1, &
                   'w', 'tpcore omega', 'Pa/s', 4, VarCt )
  CALL NC_VAR_WRITE( fId, 'w', Arr3D=Arr3D ) 
  DEALLOCATE( Arr3D )

  ! CLose file
  CALL NC_CLOSE(fId )

  END SUBROUTINE WRITE_WC
!EOC
  END MODULE SAVE_KZ_MOD
