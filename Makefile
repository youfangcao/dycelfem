# SPDX-License-Identifier: BSD-3-Clause
# SPDX-FileCopyrightText: 2026 The Board of Trustees of the University of Illinois
# =============================================================================
#  DyCelFEM - build for macOS and Linux
#
#      make                 build bin/dycelfem
#      make -j8             parallel build
#      make deps            fetch and build libSBML into third_party/install
#      make clean           remove objects and the binary
#      make distclean       also remove third_party/
#      make check           quick smoke test of the example
#
#  Dependencies
#      libSBML 5.x   'make deps' builds it locally; a system copy is used if
#                    pkg-config or the usual prefixes can find one.
#      GLUT + OpenGL macOS: the system frameworks. No install needed.
#                    Linux: freeglut3-dev libglu1-mesa-dev (Debian/Ubuntu) or
#                           freeglut-devel mesa-libGLU-devel (Fedora/RHEL).
# =============================================================================

ROOT      := $(CURDIR)
BINDIR    := $(ROOT)/bin
OBJDIR    := $(ROOT)/build/obj
TP        := $(ROOT)/third_party
LOCALSBML := $(TP)/install

UNAME_S := $(shell uname -s)

CXX      ?= c++
CXXFLAGS ?= -O2 -g
# The code predates C++11 and relies on the older dialect in places.
STD      := -std=c++03
WARN     := -Wno-deprecated-declarations -Wno-writable-strings -Wno-format-security \
            -Wno-dangling-else -Wno-logical-op-parentheses -Wno-parentheses \
            -Wno-return-type -Wno-unused-value -Wno-unknown-warning-option
FORCEINC :=

# ---- locate libSBML ---------------------------------------------------------
# Priority:
#   1. SBML_PREFIX=/path given on the command line
#   2. the copy built by 'make deps' in third_party/install
#   3. pkg-config
#   4. common system prefixes
# 'make deps' is preferred because a system libSBML may be old, or (seen in the
# wild) built with a relative install_name that cannot be loaded at run time.
SBML_CFLAGS :=
SBML_LIBS   :=
SBML_RPATH  :=
SBML_SOURCE :=

ifneq ($(SBML_PREFIX),)
  SBML_CFLAGS := -I$(SBML_PREFIX)/include
  SBML_LIBS   := -L$(SBML_PREFIX)/lib -lsbml
  SBML_RPATH  := -Wl,-rpath,$(SBML_PREFIX)/lib
  SBML_SOURCE := SBML_PREFIX=$(SBML_PREFIX)
else ifneq ($(wildcard $(LOCALSBML)/include/sbml/SBMLTypes.h),)
  SBML_CFLAGS := -I$(LOCALSBML)/include
  SBML_LIBS   := -L$(LOCALSBML)/lib -lsbml
  SBML_RPATH  := -Wl,-rpath,$(LOCALSBML)/lib
  SBML_SOURCE := third_party/install (make deps)
else
  PKGCONFIG := $(shell command -v pkg-config 2>/dev/null)
  ifneq ($(PKGCONFIG),)
    ifeq ($(shell pkg-config --exists libsbml && echo yes),yes)
      SBML_CFLAGS := $(shell pkg-config --cflags libsbml)
      SBML_LIBS   := $(shell pkg-config --libs libsbml)
      SBML_RPATH  := $(addprefix -Wl$(comma)-rpath$(comma),$(patsubst -L%,%,$(filter -L%,$(SBML_LIBS))))
      SBML_SOURCE := pkg-config (system)
    endif
  endif
  ifeq ($(SBML_LIBS),)
    SYSPREFIX := $(firstword $(foreach p,/usr/local /opt/homebrew /usr, \
                   $(if $(wildcard $(p)/include/sbml/SBMLTypes.h),$(p))))
    ifneq ($(SYSPREFIX),)
      SBML_CFLAGS := -I$(SYSPREFIX)/include
      SBML_LIBS   := -L$(SYSPREFIX)/lib -lsbml
      SBML_RPATH  := -Wl,-rpath,$(SYSPREFIX)/lib
      SBML_SOURCE := $(SYSPREFIX) (system)
    endif
  endif
endif
comma := ,

# ---- platform: GL/GLUT ------------------------------------------------------
ifeq ($(UNAME_S),Darwin)
  GL_LIBS := -framework GLUT -framework OpenGL
else
  GL_LIBS := -lglut -lGLU -lGL
  # Debian/Ubuntu multiarch and the usual local prefixes
  GL_CFLAGS := $(if $(wildcard /usr/include/GL/glut.h),,-I/usr/local/include)
endif

INCLUDES := -I$(ROOT) $(SBML_CFLAGS) $(GL_CFLAGS)
ALLFLAGS := $(STD) $(CXXFLAGS) $(WARN) $(FORCEINC)
LDLIBS   := $(SBML_LIBS) $(GL_LIBS) -lm
LDFLAGS  := $(SBML_RPATH)

# ---- sources ----------------------------------------------------------------
# libBioModel here is the subset the simulator actually links; see
# doc/libBioModel-pruning.md for how it was derived and what was removed.
# RK4A.cpp is an unrelated standalone tool with its own main().
BIO_SRCS    := $(wildcard $(ROOT)/libBioModel/*.cpp)
COMMON_SRCS := $(wildcard $(ROOT)/common/*.cpp)
SIM_SRCS    := $(filter-out $(ROOT)/src/RK4A.cpp,$(wildcard $(ROOT)/src/*.cpp))
SRCS        := $(SIM_SRCS) $(BIO_SRCS) $(COMMON_SRCS)
OBJS        := $(patsubst $(ROOT)/%.cpp,$(OBJDIR)/%.o,$(SRCS))
DEPS        := $(OBJS:.o=.d)
TARGET      := $(BINDIR)/dycelfem

.PHONY: all deps clean distclean check help check-tree

all: check-tree $(TARGET)

# Running make from the wrong directory used to produce a baffling
# "Undefined symbols: _main" from an empty object list. Say what is wrong.
check-tree:
ifeq ($(strip $(SIM_SRCS)),)
	@echo "ERROR: no sources found under $(ROOT)/src."
	@echo ""
	@echo "  You are running make from '$(ROOT)',"
	@echo "  which is not the DyCelFEM package root."
	@echo ""
	@echo "  The package root is the directory containing src/, libBioModel/,"
	@echo "  common/ and examples/ - for example:"
	@echo ""
	@echo "      cd dycelfem-1.0"
	@echo "      make deps && make -j8 && make check"
	@echo ""
	@false
endif

$(TARGET): check-sbml $(OBJS) | $(BINDIR)
	$(CXX) $(ALLFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)
	@echo ""
	@echo "  Built $@   (libSBML from $(SBML_SOURCE))"
	@echo "  Try:  ./bin/dycelfem --write-default-config my.cfg"
	@echo "        make check"
	@echo ""

# Objects depend on the libSBML configuration, so switching between a system
# copy and 'make deps' forces a rebuild instead of silently linking objects
# compiled against different headers.
SBMLSTAMP := $(OBJDIR)/.sbml-config
.PHONY: FORCE
FORCE:
$(SBMLSTAMP): FORCE
	@mkdir -p $(dir $@)
	@printf '%s\n' '$(SBML_CFLAGS) $(SBML_LIBS)' | cmp -s - $@ 2>/dev/null || \
	  { printf '%s\n' '$(SBML_CFLAGS) $(SBML_LIBS)' > $@; }

$(OBJDIR)/%.o: $(ROOT)/%.cpp $(SBMLSTAMP)
	@mkdir -p $(dir $@)
	$(CXX) $(ALLFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

$(BINDIR):
	@mkdir -p $(BINDIR)

check-sbml:
ifeq ($(SBML_LIBS),)
	@echo "ERROR: libSBML not found."
	@echo "  Build a local copy with:   make deps"
	@echo "  or point at an existing one: make SBML_PREFIX=/path/to/libsbml"
	@false
endif
	@echo "libSBML: $(SBML_SOURCE)"
	@echo "         $(SBML_CFLAGS) $(SBML_LIBS)"

deps:
	@bash $(ROOT)/build-libsbml.sh

clean:
	rm -rf $(OBJDIR) $(BINDIR)

distclean: clean
	rm -rf $(TP)

check: $(TARGET)
	@bash $(ROOT)/examples/smoke-test.sh

help:
	@echo "make        - build bin/dycelfem"
	@echo "make deps   - fetch and build libSBML locally"
	@echo "make check  - build, then run a 2-step smoke test"
	@echo "make clean / distclean"

-include $(DEPS)
