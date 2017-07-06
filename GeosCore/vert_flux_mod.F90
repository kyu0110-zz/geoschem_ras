!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !MODULE: vert_flux_mod
!
! !DESCRIPTION: \subsection*{Overview}
!  Module VERT\_FLUX\_MOD contains routines to increase vertical
!  transport in the global simulations to make up for the reduced
!  convective mass flux in the 0.25x0.3125 simulation. 
!
! !INTERFACE: 
!
MODULE VERT_FLUX_MOD
!
! !USES:
! 
  IMPLICIT NONE
  PRIVATE
!
! !PUBLIC MEMBER FUNCTIONS:
!
  PUBLIC  ::  GEOSFP_VERT_FLUX 
!
! !PRIVATE MEMBER FUNCTIONS:
!
  PRiVATE :: GET_AIR_MASS 
!
! !PRIVATE DATA MEMBERS:
!  
!
! !AUTHOR:
! !REVISION HISTORY:
! 07 March 2014 - K. Yu  - Initial version
!EOP
!------------------------------------------------------------------------------
!BOC
CONTAINS
!EOC
!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GEOSFP_VERT_FLUX
!
! !DESCRIPTION: Adds vertical mass flux to correct for reduced convective
! mass flux in GEOS-FP 
!\\
!\\
! !INTERFACE:
!
  SUBROUTINE GEOSFP_VERT_FLUX ( D_DYN, N_TRACERS, TPAUSE,  & 
                                Q, OMEGA, NEGOMEGA, P_SURF, Ap, Bp )
!
! !USES:
!
    USE GIGC_State_Met_Mod, ONLY : MetState
    USE CMN_SIZE_MOD
    USE GRID_MOD,           ONLY : GET_AREA_M2

    ! Include file w/ physical constants
    USE CMN_GCTM_MOD
    
!
! !INPUT PARAMETERS: 
!
    ! Number of tracers
    INTEGER, INTENT(IN)   :: N_TRACERS
    
    ! Dynamic time stemp
    REAL*8,  INTENT(IN)   :: D_DYN

    ! Vertical mass fluxes beign read in [Pa/s]
    REAL*8,  INTENT(IN)   :: TPAUSE(:,:)
    REAL*8,  INTENT(IN)   :: OMEGA(:,:,:)
    REAL*8,  INTENT(IN)   :: NEGOMEGA(:,:,:)
    REAL*8,  INTENT(IN)   :: P_SURF(:,:)
    REAL*8,  INTENT(IN)   :: Ap(:), Bp(:)

! !INPUT/OUTPUT PARAMETERS: 
!
    ! Tracer mixing ratios [v/v]
    REAL*8,  INTENT(INOUT), TARGET :: Q(:,:,:,:)

! !OUTPUT PARAMETERS:
!
!
! !AUTHOR:
! 
! !REVISION HISTORY: 
!  07 March 2014 - K. Yu     - Initial version                              
!EOP
!------------------------------------------------------------------------------
!BOC
!
! !DEFINED PARAMETERS:
!
!
! !LOCAL VARIABLES:
!
    INTEGER             :: I_TRACER
    INTEGER             :: I
    INTEGER             :: J
    INTEGER             :: L
    REAL*8              :: ZMASS(IIPAR,JJPAR,LLPAR)
    REAL*8              :: Q_NEW(IIPAR,JJPAR,LLPAR)
    REAL*8              :: TRACERFLUX(IIPAR,JJPAR,LLPAR)
    ! Area of each grid box
    REAL*8              :: AREA_M2(IIPAR,JJPAR,LLPAR)
    REAL*8              :: AIRMASS(IIPAR,JJPAR,LLPAR)


    !     ----------------
    !     Begin execution.
    !     ----------------

    ! Print info
    WRITE( 6, '(a)' ) 'Increasing vertical fluxes for global simulations' 

    ! Compute surface area of grid box 
    DO I = 1, IIPAR
      DO J = 1, JJPAR
        DO L = 1, LLPAR
          AREA_M2(I,J,L) = GET_AREA_M2(I,J,L)
          AIRMASS(I,J,L) = GET_AIR_MASS(I,J,L,P_SURF(I,J), Ap, Bp)
        ENDDO
      ENDDO
    ENDDO

    print*, 'AIRMASS', SUM(AIRMASS)

    ! UP/DOWN mass flux [kg air / m^2]
    ! If the net omega is in the same direction as negomega, the additional
    ! flux is the difference between the two. If the net flux is in the 
    ! opposite direction, the additional flux is just negomega
    WHERE ( OMEGA < 0d0 ) &       
         ZMASS = ((NEGOMEGA / -9.81) * D_DYN) - ((OMEGA / -9.81) * D_DYN)
    WHERE ( OMEGA .ge. 0d0) &
       ZMASS = (NEGOMEGA / -9.81) * D_DYN

    print*, 'OMEGA', SUM(OMEGA)
    print*, 'NEGOMEGA', SUM(NEGOMEGA)  
    ! UP/DOWNW mass flux [kg air]
    DO L = 1, LLPAR
      DO J = 1, JJPAR
        DO I = 1, IIPAR
          ZMASS(I,J,L) = ZMASS(I,J,L) * AREA_M2(I,J,L)
        ENDDO
      ENDDO
    ENDDO

    ! Can't more more than half the air of the box each box
    DO L = 1, LLPAR
    DO J = 1, JJPAR
    DO I = 1, IIPAR
    ZMASS(I,J,L) = min(0.5 * AIRMASS(I,J,L), ZMASS(I,J,L))
    ENDDO
    ENDDO
    ENDDO

    print*, 'LLTROP', LLTROP
    print*, 'ZMASS', SUM(ZMASS)
    print*, 'AIRMASS', SUM(AIRMASS)
    ! Loop through each tracer
    DO I_TRACER = 1, N_TRACERS

    print*, 'Q in', SUM(Q(:,:,:,I_TRACER))
      ! Compute the additional mass flux of each tracer, which is the mass
      ! flux of air multiplied by the mixing ratio of the tracer divided by
      ! ratio of molec weight of tracer to molec weight of air
      ! [kg tracer] = [kg air] * [kg tracer/kg air] 
      TRACERFLUX(:,:,:) = ZMASS(:,:,:) * Q(:,:,:,I_TRACER)

      ! Make tracerflux zero above tropopause
      TRACERFLUX(:,:,1:LLPAR-LLTROP) = 0.0

    print*, 'TRACERFLUX', SUM(TRACERFLUX)

      ! Bottom grid box first: new mass = mass currently in cell minus 
      ! mass out + mass in from box above
      Q_NEW(:,:,LLPAR) = Q(:,:,LLPAR,I_TRACER) * AIRMASS(:,:,LLPAR) - TRACERFLUX(:,:,LLPAR) + TRACERFLUX(:,:,LLPAR-1)
   
      ! Do transpot layer by layer
      DO L = LLPAR-1, LLPAR-LLTROP+2, -1

        ! transport tracers out of grid box (both up and down)
        ! New mass = mass currently in cell - 2 * mass out (one for up and one
        ! for down) + mass in from cell above + mass in from cell below
        Q_NEW(:,:,L) = Q(:,:,L,I_TRACER) * AIRMASS(:,:,L) - (2.0 * TRACERFLUX(:,:,L)) + & 
                       TRACERFLUX(:,:,L-1) + TRACERFLUX(:,:,L+1)

      ENDDO   ! L
   
      ! Top layer
      ! New mass = mass currently in cell - mass out (down) + mass in 
      ! from cell below
      Q_NEW(:,:,LLPAR-LLTROP+1) = Q(:,:,LLPAR-LLTROP+1,I_TRACER) * AIRMASS(:,:,LLPAR-LLTROP+1) - TRACERFLUX(:,:,LLPAR-LLTROP+1) + TRACERFLUX(:,:,LLPAR-LLTROP+2)

      ! Reset original array with the new values, converting from [kg] to
      ! [v/v]
      Q(:,:,LLPAR-LLTROP+1:LLPAR,I_TRACER) = Q_NEW(:,:,LLPAR-LLTROP+1:LLPAR) / AIRMASS(:,:,LLPAR-LLTROP+1:LLPAR)

        print*, 'Q_NEW', SUM(Q(:,:,:,I_TRACER))
    ENDDO     ! I_TRACER
    
    WRITE( 6, '(a)' ) 'Done vertical fluxes for global simulations' 

  END SUBROUTINE GEOSFP_VERT_FLUX
!EOC
!-------------------------------------------------------------------------------  
!BOP
    FUNCTION GET_AIR_MASS( I, J, L, P_SURF, Ap, Bp )    RESULT( AIR_MASS )
! !USES:
    USE CMN_SIZE_MOD
    USE CMN_GCTM_MOD
    !USE PRESSURE_MOD,   ONLY    : AP, BP, GET_AP, GET_BP
    USE GRID_MOD,       ONLY    : GET_AREA_M2
! !INPUT PARAMETERS:
    INTEGER, INTENT(IN)         :: I, J, L
    REAL(fp), INTENT(IN)        :: P_SURF
    REAL(fp), INTENT(IN)        :: Ap(:), Bp(:)
!EOP
!BOC
!   !LOCAL VARIABLES:
    REAL(fp) :: P_BOT, P_TOP, AIR_MASS
    INTEGER  :: L2

    L2 = L + 1

    P_BOT = AP(L2) + ( BP(L2) * P_SURF )
    P_TOP = AP(L2-1) + ( BP(L2-1) * P_SURF )
    AIR_MASS = (P_BOT - P_TOP) * G0_100 * GET_AREA_M2(I,J,L)

    END FUNCTION GET_AIR_MASS
!EOC
END MODULE VERT_FLUX_MOD
!EOC
