#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  #
#------------------------------------------------------------------------------
#BOP
#
# !IROUTINE: Makefile_MakeDoc.mk (in doc subdirectory)
#
# !DESCRIPTION: This Makefile fragment contains commands to build the 
#  documentation for the GEOS-Chem Makefiles  It is inlined into
#  the Makefile (in the doc subdirectory) by an "include" command.
#\\
#\\
# !REMARKS:
# To build the documentation, call "make" with the following syntax:
#                                                                             .
#   make TARGET [ OPTIONAL-FLAGS ]
#                                                                             .
# To display a complete list of options, type "make help".
#                                                                             .
# You must have the LaTeX utilities (latex, dvips, dvipdf) installed
# on your system in order to build the documentation.
#
# !REVISION HISTORY: 
#  14 Sep 2010 - R. Yantosca - Initial version, split off from Makefile
#  16 Dec 2010 - R. Yantosca - Renamed output files to "GC_Ref_Vol_1.*"
#  15 Jan 2014 - R. Yantosca - Now only create *.pdf output
#  15 Jan 2014 - R. Yantosca - Now only prints prologues, avoids printing code
#  10 Jul 2015 - R. Yantosca - Use ./protex to avoid problems on some systems
#EOP
#------------------------------------------------------------------------------
#BOC


# List of source code files (order is important)
SRC2 :=                          \
./intro.make                     \
$(ROOT)/Makefile                 \
$(ROOT)/Makefile_header.mk       \
$(HCO)/Makefile                  \
$(HCOI)/Makefile                 \
$(HCOX)/Makefile                 \
$(UTIL)/Makefile                 \
$(ISO)/Makefile                  \
$(CORE)/Makefile                 \
$(GTMM)/Makefile                 \
$(GCRT)/Makefile                 \
$(KPP)/Makefile                  \
$(KPP)/benchmark/Makefile        \
$(KPP)/NOx_Ox_HC_Aer_Br/Makefile \
$(KPP)/SOA/Makefile              \
$(KPP)/SOA_SVPOA/Makefile        \
$(KPP)/UCX/Makefile              \
$(DOC)/Makefile                  \
$(DOC)/Makefile_SrcDoc.mk        \
$(DOC)/Makefile_UtilDoc.mk       \
$(DOC)/Makefile_GtmmDoc.mk       \
$(DOC)/Makefile_MakeDoc.mk       \
$(HELP)/Makefile


# Output file names
TEX2 := GC_Ref_Vol_1.tex
DVI2 := GC_Ref_Vol_1.dvi
PDF2 := GC_Ref_Vol_1.pdf


# Make command
makedoc: 
	rm -f $(TEX2)
	./protex -sfS $(SRC2) > $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	latex $(TEX2)
	dvipdf $(DVI2) $(PDF2)
	rm -f *.aux *.dvi *.log *.toc

#EOC
