# $Id: GNUmakefile,v 1.9 2007/06/22 12:44:42 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := EventGeneratorROOT
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL :=/usr/local/geant/geant4.9.0.p01
endif
CLHEP_BASE_DIR    := /usr/local/CLHEP
include $(G4INSTALL)/config/architecture.gmk
LDFLAGS += -L/usr/X11R6/lib64 -lXaw
LDFLAGS += -L/usr/lib64 


#===============================================================
# ROOT configuration

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --libs)

CXXFLAGS += $(ROOT_CFLAGS)
CPPFLAGS += $(ROOT_CFLAGS)
LDFLAGS  += $(ROOT_LIBS)
LDFLAGS  += -L/usr/lib -lfreetype
#===============================================================


.PHONY: all
all: lib bin

#### G4ANALYSIS_USE := true

#ifdef G4ANALYSIS_USE
 # CPPFLAGS += -DG4ANALYSIS_USE
#endif

include $(G4INSTALL)/config/architecture.gmk

#ifdef G4ANALYSIS_USE
#  # for the aida-config command see the README file
#  CPPFLAGS += `aida-config --include`
#  LOADLIBS += `aida-config --lib`
#endif

include $(G4INSTALL)/config/binmake.gmk

#visclean:
#	rm -f g4*.prim g4*.eps g4*.wrl
#	rm -f .DAWN_*
#
#histclean:
#	rm ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/RunAction.o
