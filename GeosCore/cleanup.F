!------------------------------------------------------------------------------
!                  GEOS-Chem Global Chemical Transport Model                  !
!------------------------------------------------------------------------------
!BOP
!
! !ROUTINE: cleanup
!
! !DESCRIPTION: Subroutine CLEANUP deallocates the memory assigned to 
!  dynamically allocatable arrays just before exiting a GEOS-Chem simulation.
!\\
!\\
! !INTERFACE:
!
      SUBROUTINE CLEANUP( am_I_Root, Input_Opt, ERROR, RC )
!
! !USES:
!
      USE AEROSOL_MOD,             ONLY : CLEANUP_AEROSOL
      USE C2H6_MOD,                ONLY : CLEANUP_C2H6
      USE CARBON_MOD,              ONLY : CLEANUP_CARBON
      USE CHEMGRID_MOD,            ONLY : CLEANUP_CHEMGRID
      USE CO2_MOD,                 ONLY : CLEANUP_CO2
      USE COMODE_MOD,              ONLY : CLEANUP_COMODE
      USE GCKPP_COMODE_MOD,        ONLY : CLEANUP_GCKPP_COMODE
      USE DIAG_MOD,                ONLY : CLEANUP_DIAG
      USE DIAG03_MOD,              ONLY : CLEANUP_DIAG03
      USE DIAG04_MOD,              ONLY : CLEANUP_DIAG04
      USE DIAG41_MOD,              ONLY : CLEANUP_DIAG41
      USE DIAG50_MOD,              ONLY : CLEANUP_DIAG50
      USE DIAG51_MOD,              ONLY : CLEANUP_DIAG51
      USE DIAG53_MOD,              ONLY : CLEANUP_DIAG53 !(clf, 3/11/11)
      USE DIAG_OH_MOD,             ONLY : CLEANUP_DIAG_OH
      USE DIAG_PL_MOD,             ONLY : CLEANUP_DIAG_PL
      USE DRYDEP_MOD,              ONLY : CLEANUP_DRYDEP
      USE DUST_MOD,                ONLY : CLEANUP_DUST
      USE ERROR_MOD,               ONLY : DEBUG_MSG
      USE GET_NDEP_MOD,            ONLY : CLEANUP_GET_NDEP
      USE GLOBAL_CH4_MOD,          ONLY : CLEANUP_GLOBAL_CH4
      USE GRID_MOD,                ONLY : CLEANUP_GRID
      USE HDF_MOD,                 ONLY : CLEANUP_HDF
      USE ISOROPIAII_MOD,          ONLY : CLEANUP_ISOROPIAII
      USE LINOZ_MOD,               ONLY : CLEANUP_LINOZ
      USE MERCURY_MOD,             ONLY : CLEANUP_MERCURY
      USE MODIS_LAI_MOD,           ONLY : CLEANUP_MODIS_LAI
      USE OCEAN_MERCURY_MOD,       ONLY : CLEANUP_OCEAN_MERCURY
      USE DEPO_MERCURY_MOD,        ONLY : CLEANUP_DEPO_MERCURY
      USE LAND_MERCURY_MOD,        ONLY : CLEANUP_LAND_MERCURY
      USE PBL_MIX_MOD,             ONLY : CLEANUP_PBL_MIX
      USE PJC_PFIX_MOD,            ONLY : CLEANUP_PJC_PFIX
      USE PLANEFLIGHT_MOD,         ONLY : CLEANUP_PLANEFLIGHT
      USE PRESSURE_MOD,            ONLY : CLEANUP_PRESSURE
      USE SEASALT_MOD,             ONLY : CLEANUP_SEASALT
      USE SULFATE_MOD,             ONLY : CLEANUP_SULFATE
      USE STRAT_CHEM_MOD,          ONLY : CLEANUP_STRAT_CHEM
      USE TAGGED_CO_MOD,           ONLY : CLEANUP_TAGGED_CO
#if   defined( TOMAS )
      USE TOMAS_MOD,               ONLY : CLEANUP_TOMAS  !sfarina, 1/16/13
#endif
      USE TOMS_MOD,                ONLY : CLEANUP_TOMS
      USE TPCORE_FVDAS_MOD,        ONLY : EXIT_TPCORE
      USE TPCORE_GEOS5_WINDOW_MOD, ONLY : EXIT_GEOS5_TPCORE_WINDOW
      USE TPCORE_GEOSFP_WINDOW_MOD,ONLY : EXIT_GEOSFP_TPCORE_WINDOW !clb
      USE UCX_MOD,                 ONLY : CLEANUP_UCX
      USE WETSCAV_MOD,             ONLY : CLEANUP_WETSCAV
#if   defined( APM )
      USE APM_WETS_MOD,            ONLY : CLEANUP_WETSCAVBIN
#endif

#if ! defined( ESMF_ )
      USE TRANSPORT_MOD,           ONLY : CLEANUP_TRANSPORT
#endif
      USE CMN_O3_Mod,              ONLY : Cleanup_CMN_O3
      USE CMN_SIZE_Mod,            ONLY : Cleanup_CMN_SIZE
      USE COMODE_LOOP_Mod,         ONLY : Cleanup_COMODE_LOOP
      USE VDIFF_PRE_Mod,           ONLY : Cleanup_VDIFF_PRE
      USE CMN_FJX_Mod,             ONLY : Cleanup_CMN_FJX
      USE GIGC_Input_Opt_Mod,      ONLY : OptInput
      USE Regrid_A2A_Mod,          ONLY : Cleanup_Map_A2a
      
      ! HEMCO
!#if !defined(ESMF_) 
      USE EMISSIONS_MOD,           ONLY : EMISSIONS_FINAL
!#endif

#if defined( DEVEL )
      ! Cleanup routine for netCDF diagnostics structure
      USE DIAGNOSTICS_MOD,         ONLY : DIAGNOSTICS_FINAL
#endif

#if   defined( RRTMG )
      USE RRTMG_RAD_TRANSFER_MOD,  ONLY : CLEANUP_SURFACE_RAD
      USE RRTMG_RAD_TRANSFER_MOD,  ONLY : CLEANUP_STRAT_CLIM
      USE RRTMG_RAD_TRANSFER_MOD,  ONLY : CLEANUP_MCICA_CLOUDS
#endif


      IMPLICIT NONE
!
! !INPUT PARAMETERS: 
!
      LOGICAL,        INTENT(IN)  :: am_I_Root    ! Are we on the root CPU?
      TYPE(OptInput), INTENT(IN)  :: Input_Opt    ! Input Options object 
      LOGICAL,        INTENT(IN)  :: ERROR        ! Cleanup after error? 
!
! !OUTPUT PARAMETERS:
!
      INTEGER,        INTENT(OUT) :: RC           ! Success or failure
!
! !REVISION HISTORY: 
!  29 Nov 1999 - R. Yantosca - Initial version
!  (1 ) CLEANUP is written in Fixed-Format F90.
!  (2 ) Now calls CLEANUP_WETSCAV, which deallocates arrays from 
!        "wetscav_mod.f". (bmy, 3/9/00)
!  (3 ) Add call to CLEANUP_SULFATE, which deallocates arrays from
!        "sulfate_mod.f".  Also now deallocate ND32 arrays. (bmy, 6/6/00)
!  (4 ) Add call to CLEANUP_DAO, which deallocates arrays from "dao_mod.f".  
!        (bmy, 6/26/00)
!  (5 ) Add call to CLEANUP_TAGGED_CO and CLEANUP_COMODE, which deallocates 
!        arrays from and "comode_mod.f". (bmy, 7/19/00)
!  (6 ) Add call to CLEANUP_GLOBAL_OH and CLEANUP_COMODE, which deallocates 
!        arrays from "global_oh_mod.f". (bmy, 7/28/00)
!  (7 ) Add calls to CLEANUP_BIOMASS and CLEANUP_BIOFUEL, which deallocates 
!        arrays from "biomass_mod.f" and "biofuel_mod.f".  Also deallocate
!        the AD32_bf array for the biofuel NOx diagnostic. (bmy, 9/12/00)
!  (8 ) Add call to CLEANUP_DIAG51, to deallocate module arrays from
!        "diag51_mod.f" (bmy, 11/29/00)
!  (9 ) Removed obsolete code from 11/29/00 (bmy, 12/21/00)
!  (10) Add call to CLEANUP_CH4, to deallocate module arrays from
!        "global_ch4_mod.f" (bmy, 1/16/01)
!  (11) Now deallocate the AD34 array.  Also updated comments and
!        made some cosmetic changes. (bmy, 3/15/01)
!  (12) Now deallocate the AD12 array (bdf, bmy, 6/15/01)
!  (13) Add call to CLEANUP_ACETONE, to deallocate module arrays from 
!        "acetone_mod.f"  Also deallocate AD11 array.  Also deallocate 
!        variables from dao_mod.f last, to try to avoid bus error on 
!        SGI (bmy, 8/3/01) 
!  (14) Added call to CLEANUP_UVALBEDO from "uvalbedo_mod.f".  Also removed
!        obsolete code from 9/01.  Also only include references to CLEANUP_* 
!        subroutines in other modules for clarity. (bmy, 1/15/02)
!  (15) Added call to CLEANUP_C2H6 from "c2h6_mod.f" (bmy, 1/25/02)
!  (16) Added call to CLEANUP_AIRCRAFT_NOX from "aircraft_nox_mod.f" 
!        (bmy, 2/14/02)
!  (17) Now deallocate CTNO2, CTHO2, LTNO2, LTHO2 arrays (rvm, bmy, 2/27/02)
!  (18) Now reference CLEANUP_PLANEFLIGHT from "planeflight_mod.f".
!        Now also deallocate AD01 and AD02 arrays. (mje, bmy, 8/7/02)
!  (19) Now reference cleanup routines from "global_nox_mod.f", 
!        "global_hno3_mod.f", "global_no3_mod.f", "drydep_mod.f", and
!        "rpmares_mod.f". (bmy, 12/16/02)
!  (20) Now reference cleanup routine from "transport_mod.f" (bmy, 2/10/03)
!  (21) Now reference cleanup routine from "pjc_pfix_mod.f" and 
!        "tpcore_fvdas_mod.f90". (bmy, 5/9/03)
!  (22) Now reference cleanup routine from "toms_mod.f" (bmy, 7/14/03)
!  (23) Now reference cleanup routine from "carbon_mod.f", "dust_mod.f", and
!        "dust_dead_mod.f". (bmy, 7/14/03)
!  (23) Now references cleanup routine from "lightning__nox_mod.f" 
!        (bmy, 4/14/04)
!  (24) Now references cleanup routine from "seasalt_mod.f" (bmy, 4/26/04)
!  (25) Now references cleanup routines from new modules (bmy, 7/20/04)
!  (26) Now calls cleanup routine from "epa_nei_mod.f" (bmy, 11/5/04)
!  (27) Now call CLEANUP_MERCURY from "mercury_mod.f" (eck, bmy, 12/7/04)
!  (28) Now call CLEANUP_OCEAN_MERCURY from "ocean_mercury_mod.f".  Also
!        reordered the calling sequence. (sas, bmy, 1/21/05)
!  (29) Now call CLEANUP_PBL_MIX from "pbl_mix_mod.f".  Now call CLEANUP_DIAG41
!        from "diag41_mod.f". (bmy, 2/17/05)
!  (30) Now calls CLEANUP_HCN_CH3CN from "hcn_ch3cn_mod.f (bmy, 6/23/05)
!  (31) Now calls CLEANUP_DIAG04, CLEANUP_CO2, and CLEANUP_TROPOPAUSE 
!         (bmy, 8/15/05) 
!  (32) Now calls CLEANUP_LAI from "lai_mod.f", CLEANUP_MEGAN from
!        "megan_mod.f" and CLEANUP_REGRID_1x1 from "regrid_1x1_mod.f"
!        (tmf, bdf, bmy, 10/24/05)
!  (33) Now calls CLEANUP_EMEP from "emep_mod.f" (bdf, bmy, 11/1/05)
!  (34) Now calls CLEANUP_GC_BIOMASS and CLEANUP_GFED2_BIOMASS (bmy, 4/5/06)
!  (35) Now calls CLEANUP_DIAG56 from "diag56_mod.f" and
!        CLEANUP_LIGHTNING_NOX_NL from "lightning_nox_nl_mod.f"
!        (ltm, bmy, 5/5/06)
!  (36) Now references CLEANUP_BRAVO from "bravo_mod.f" and CLEANUP_EDGAR
!        from "edgar_mod.f" (bmy, 7/6/06)
!  (37) Now calls CLEANUP_H2_HD from "h2_hd_mod.f" and CLEANUP_GLOBAL_O1D
!        from "global_o1d_mod.f".  Remove call to CLEANUP_LIGHTNING_NOx_NL 
!        from "lightning_nox_nl_mod.f (hup, phs, bmy, 10/2/07)
!  (38) Now calls GEOS5_EXIT_TPCORE_WINDOW to finalize the TPCORE for
!        GEOS-5 nested window simulations (yxw, dan, bmy, 11/6/08)
!  (39) Now references CLEANUP_CAC_ANTHRO (amv, phs, 3/10/08)
!  (40) Now references CLEANUP_ARCTAS_SHIP (phs, 3/10/08)
!  (41) Now references CLEANUP_VISTAS_ANTHRO (phs, 3/10/08)
!  (41) Now references CLEANUP_LINOZ (phs, 10/16/09)
!  (42) Now references CLEANUP_HDF (amv, bmy, 12/21/09)
!  (43) Now references CLEANUP_ISOROPIAII (ccc, bmy, 1/29/09)
!  (44) Now references CLEANUP_DEPO_MERCURY and CLEANUP_LAND_MERCURY
!       (ccc, 5/6/10)
!  (45) Added call to CLEANUP_GLOBAL_OC, which deallocates arrays from
!       "global_oc_mod.f"  (clf, 2/28/2011)
!  (46) Added call to CLEANUP_GLOBAL_BC, which deallocates arrays from
!       "global_bc_mod.f"  (clf, 2/28/2011) 
!  13 Aug 2010 - R. Yantosca - Added ProTeX headers
!  28 Feb 2011 - C. Friedman - Added call to CLEANUP_GLOBAL_OC, which 
!                              deallocates arrays from "global_oc_mod.f"
!  22 Aug 2011 - R. Yantosca - Add modifications for APM microphysics (G. Luo)
!  07 Sep 2011 - P. Kasibhatla - Add modifications for GFED3
!  19 Mar 2012 - M. Payer    - Remove call to CLEANUP_ACETONE. It is no longer
!                              needed after removal of JO1D and RESP routines.
!  05 Apr 2012 - R. Yantosca - Now call CLEANUP_MODIS_LAI
!  11 Apr 2012 - R. Yantosca - Remove reference to obsolete lai_mod.F
!  01 May 2012 - M. Payer    - Added calls for CLEANUP_GLOBAL_GRID and
!                              CLEANUP_GRID
!  19 Nov 2012 - R. Yantosca - Add cleanup calls for modules in Headers;
!                              these were omitted during development
!  28 Nov 2012 - R. Yantosca - Remove reference to CLEANUP_DAO; we have now
!                              removed all allocatable arrays from dao_mod.F
!  03 Dec 2012 - R. Yantosca - Now call Cleanup_CMN_SIZE
!  13 Dec 2012 - R. Yantosca - Remove reference to obsolete CMN_DEP_mod.F
!  06 Mar 2013 - H. Amos     - Merge C. Friedman's PAH code
!  14 Mar 2013 - M. Payer    - Restore reference to CLEANUP_DAO. Some arrays
!                              have been restored to dao_mod.F for use in
!                              the vertical regridding of OH for offline
!                              simulations.
!  24 Mar 2013 - S.D. Eastham- Switched TROPOPAUSE_MOD to CHEMGRID_MOD
!  04 Apr 2013 - S.D. Eastham- Added call for CLEANUP_UCX
!  20 Aug 2013 - R. Yantosca - Removed "define.h", this is now obsolete
!  26 Sep 2013 - R. Yantosca - Renamed GEOS_57 Cpp switch to GEOS_FP
!  26 Sep 2013 - R. Yantosca - Now reference EXIT_GEOSFP_TPCORE_WINDOW
!  29 Oct 2013 - R. Yantosca - Now call CLEANUP_DAO for all simulations (if
!                              arrays are not allocated it'll exit gracefully)
!  03 Jun 2014 - R. Yantosca - Don't finalize HEMCO if emissions are turned
!                              off.  This prevents a seg fault error.
!  03 Jun 2014 - R. Yantosca - Now pass Input_Opt via the argument list
!  03 Jun 2014 - R. Yantosca - Remove references to emissions modules that
!                              HEMCO renders obsolete
!  24 Jun 2014 - R. Yantosca - Remove CLEANUP_BROMOCARB; HEMCO replaces this
!  25 Jun 2014 - R. Yantosca - Remove reference to tracer_mod.F
!  21 Jul 2014 - R. Yantosca - Remove reference to regrid_1x1_mod.F
!  23 Jul 2014 - R. Yantosca - Remove reference to obsolete CMN_mod.F
!  23 Jul 2014 - R. Yantosca - Remove reference to obsolete global_grid_mod.F
!  23 Jul 2014 - R. Yantosca - Remove reference to obsolete CMN_NOX_mod.F
!  25 Jul 2014 - R. Yantosca - Remove reference to commsoil_mod.F90
!  25 Jul 2014 - R. Yantosca - Now call CLEANUP_GET_NDEP
!  15 Aug 2014 - R. Yantosca - Remove reference to biofuel_mod.F
!  15 Aug 2014 - R. Yantosca - Remove reference to biomass_mod.F
!  10 Sep 2014 - M. Sulprizio- Remove references to global_oc_mod.F and
!                              global_bc_mod.F
!  23 Sep 2014 - M. Sulprizio- Remove references to global_hno3_mod.F and
!                              global_no3_mod.F
!  12 Jan 2015 - R. Yantosca - Remove CLEANUP_UVALBEDO routine
!  11 Mar 2015 - R. Yantosca - Remove call to CLEANUP_GLOBAL_OH
!  11 Mar 2015 - R. Yantosca - Remove call to CLEANUP_GLOBAL_NOX
!  16 Mar 2015 - R. Yantosca - Remove call to CLEANUP_HCN_CH3CN
!  25 Mar 2015 - C. Keller   - Now cleanup UCX before emissisions.
!  16 Apr 2015 - R. Yantosca - Remove call to CLEANUP_DAO; it's obsolete
!EOP
!------------------------------------------------------------------------------
!BOC
      !=================================================================
      ! CLEANUP begins here!
      !=================================================================

      ! Echo info
      WRITE( 6, 100 ) 
 100  FORMAT( '     - CLEANUP: deallocating arrays now...' )

      ! UCX needs to be cleaned up before emissions, because the UCX 
      ! restart variables needs to be passed to the HEMCO diagnostics
      ! first.
      CALL CLEANUP_UCX( am_I_Root )

#if defined( DEVEL )
      ! Cleanup the new netCDF diagnostics structure
      CALL DIAGNOSTICS_FINAL( am_I_Root, Input_Opt, RC )
#endif

      ! Cleanup HEMCO
      IF ( Input_Opt%LEMIS ) THEN
         CALL EMISSIONS_FINAL ( am_I_Root, ERROR )
      ENDIF

      ! Call cleanup routines from individual F90 modules
      CALL CLEANUP_AEROSOL()
      CALL CLEANUP_C2H6()
      CALL CLEANUP_CARBON()
      CALL CLEANUP_CHEMGRID()
      CALL CLEANUP_CO2()
      CALL CLEANUP_COMODE()
      CALL CLEANUP_GCKPP_COMODE()
      CALL CLEANUP_DIAG()
      CALL CLEANUP_DIAG03()
      CALL CLEANUP_DIAG04()
      CALL CLEANUP_DIAG41()
      CALL CLEANUP_DIAG50()
      CALL CLEANUP_DIAG51()
      CALL CLEANUP_DIAG53()
      CALL CLEANUP_DIAG_OH()
      CALL CLEANUP_DIAG_PL()
      CALL CLEANUP_DRYDEP()
      CALL CLEANUP_DUST()
      CALL CLEANUP_GLOBAL_CH4()
      CALL CLEANUP_GRID()
      CALL CLEANUP_HDF()
      CALL CLEANUP_ISOROPIAII()
      CALL CLEANUP_LINOZ()
      CALL CLEANUP_MAP_A2A()
      CALL CLEANUP_MERCURY()
      CALL CLEANUP_MODIS_LAI()
      CALL CLEANUP_OCEAN_MERCURY()
      CALL CLEANUP_DEPO_MERCURY()
      CALL CLEANUP_LAND_MERCURY()
      CALL CLEANUP_PBL_MIX()
      CALL CLEANUP_PJC_PFIX()
      CALL CLEANUP_PLANEFLIGHT()
      CALL CLEANUP_PRESSURE()
      CALL CLEANUP_SEASALT()
      CALL CLEANUP_SULFATE()
      CALL CLEANUP_STRAT_CHEM()
      CALL CLEANUP_TAGGED_CO()
#if ! defined( ESMF_ )
      CALL CLEANUP_TRANSPORT()
#endif
      CALL CLEANUP_TOMS()
      CALL CLEANUP_WETSCAV()

#if   defined( APM )
      CALL CLEANUP_WETSCAVBIN()  ! Cleanup arrays for APM wet scavenging
#endif

#if   defined( TOMAS )
      CALL CLEANUP_TOMAS()
#endif

#if   defined( RRTMG )
      CALL CLEANUP_SURFACE_RAD  ! Cleanup arrays for surface properties
      CALL CLEANUP_STRAT_CLIM
      CALL CLEANUP_MCICA_CLOUDS
#endif

#if   defined(NEST_CH) || defined(NEST_NA) || defined(NEST_EU)
#if   defined( GEOS_5 ) && defined( GRID05x0666 )
      CALL EXIT_GEOS5_TPCORE_WINDOW() 
#elif defined( GEOS_FP) && defined( GRID025x03125)
      ! clb
      CALL EXIT_GEOSFP_TPCORE_WINDOW()
#endif
#else
      CALL EXIT_TPCORE()
#endif

      ! Call extra cleanup routines, mostly from modules in Headers/
      CALL Cleanup_CMN_O3     ( am_I_Root, RC )
      CALL Cleanup_CMN_SIZE   ( am_I_Root, RC )
      CALL Cleanup_Get_Ndep   ( am_I_Root, RC )
      CALL Cleanup_COMODE_LOOP( am_I_Root, RC )
      CALL Cleanup_VDIFF_PRE  ( am_I_Root, RC )
      CALL Cleanup_CMN_FJX    ( am_I_Root, RC )

      END SUBROUTINE CLEANUP
!EOC
