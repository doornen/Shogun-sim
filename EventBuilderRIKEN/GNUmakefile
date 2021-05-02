# $Id: GNUmakefile,v 1.9 2007/06/22 12:44:42 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := EventBuilder
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL=/usr/local/geant4/geant4.9.4.p04
endif

include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/binmake.gmk
#LDFLAGS += -L/usr/X11R6/lib64 -lXaw
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
