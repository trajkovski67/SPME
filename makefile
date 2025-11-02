#!/bin/make -f
# ============================================================================
#  Makefile for PME program with Intel Fortran + Intel MKL (FFTW wrapper)
# ============================================================================

# --------------------------------------------------------------------------
# Program name
# --------------------------------------------------------------------------
PROG := pme

# --------------------------------------------------------------------------
# Source files (relative to vpath)
# --------------------------------------------------------------------------
SRCS := main.f90 \
        integrals.f90 \
        slater.f90 \
        prog.f90 \
        print_matrix.f90 \
        linear_algebra.f90 \
        io_tools.f90

# --------------------------------------------------------------------------
# Build directories
# --------------------------------------------------------------------------
BUILD := build
OBJS  := $(patsubst %.f90,$(BUILD)/%.o,$(SRCS))

# --------------------------------------------------------------------------
# Compiler & flags
# --------------------------------------------------------------------------
FC     := ifort
LD     := $(FC)
RM     := rm -f

# Default flags (optimized build)
FLAGS := -O3 -xHost -fp-model fast=2 -qopenmp -ipo -module $(BUILD)
FFLAGS += -module $(BUILD)

# Search for sources in these directories
vpath %.f90 src app

# --------------------------------------------------------------------------
# Intel MKL configuration (with FFTW3 interface)
# --------------------------------------------------------------------------
# Make sure Intel oneAPI environment is loaded
MKLROOT ?= /opt/intel/mkl

# Include FFTW headers
FFLAGS  += -I$(MKLROOT)/include/fftw

# Link MKL libraries (Intel threading + FFTW Fortran interface)
LIBS := -L$(MKLROOT)/lib/intel64 \
        -Wl,--start-group \
        -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
        -lmkl_cdft_core \
        -Wl,--end-group \
        -liomp5 -lpthread -lm -ldl

# --------------------------------------------------------------------------
# Targets
# --------------------------------------------------------------------------
.PHONY: all setup clean veryclean debug release

all: setup $(PROG)

# Link executable
$(PROG): $(OBJS)
	$(LD) $(OBJS) $(LIBS) -o $@

# Compile objects
$(BUILD)/%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Explicit dependencies (optional)
$(BUILD)/main.f90.o: $(BUILD)/integrals.f90.o \
                     $(BUILD)/slater.f90.o \
                     $(BUILD)/print_matrix.f90.o \
                     $(BUILD)/linear_algebra.f90.o \
                     $(BUILD)/io_tools.f90.o
$(BUILD)/prog.f90.o: $(BUILD)/main.f90.o

# Create build directories
setup: $(BUILD) $(BUILD)/lib
$(BUILD) $(BUILD)/lib:
	mkdir -p $@

# Clean object files
clean:
	$(RM) $(OBJS)

# Remove everything, including modules and binary
veryclean:
	$(RM) $(wildcard $(BUILD)/*.o) $(wildcard $(BUILD)/lib/*.o) \
	      $(wildcard $(BUILD)/*.mod) $(PROG)

# --------------------------------------------------------------------------
# Extra build modes
# --------------------------------------------------------------------------
debug: FFLAGS = -O0 -g -warn all -check all -fpe0 -traceback -module $(BUILD)
debug: veryclean all

release: FFLAGS = -O3 -xHost -fp-model fast=2 -qopenmp -ipo -module $(BUILD)
release: veryclean all


