name := ShogunReconstructor
#g++ -o ShogunReconstructor ShogunReconstructor.C `root-config --cflags --glibs`

#===============================================================
# ROOT configuration

ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS   := $(shell root-config --glibs)

CXXFLAGS += $(ROOT_CFLAGS)
#CPPFLAGS += $(ROOT_CFLAGS)
#LDFLAGS  += $(ROOT_LIBS)
#LDFLAGS  += -L/usr/lib -lfreetype
#===============================================================

.PHONY: all
all: $(name)





