# Makefile for the EZmock program in C

TOPDIR = ..
include $(TOPDIR)/options.mk

LIBS = -lm

# Settings for FFTW
ifneq ($(strip $(FFTW_DIR)),)
  LIBS += -L$(strip $(FFTW_DIR))/lib
  INCL += -I$(strip $(FFTW_DIR))/include
endif

# Settings for CFITSIO
ifeq ($(strip $(WITH_FITS)), T)
  CFLAGS += -DWITH_CFITSIO
  LIBS += -lcfitsio
  ifneq ($(strip $(CFITSIO_DIR)),)
    LIBS += -L$(strip $(CFITSIO_DIR))/lib
    INCL += -I$(strip $(CFITSIO_DIR))/include
  endif
endif

# Settings for the precision
ifeq ($(strip $(SINGLE_PREC)), T)
  CFLAGS += -DSINGLE_PREC
  LIBS += -lfftw3f
  ifeq ($(strip $(USE_OMP)), T)
    LIBS += -DOMP -fopenmp -lfftw3f_omp
  endif
else
  LIBS += -lfftw3
  ifeq ($(strip $(USE_OMP)), T)
    LIBS += -DOMP -fopenmp -lfftw3_omp
  endif
endif

INCL += -I$(TOPDIR)/libEZmock/EZmock -I$(TOPDIR)/libEZmock/math -I$(TOPDIR)/libEZmock/prand/src/header -Imain -Iio -Ilib
EXEC = $(TOPDIR)/EZmock
SRCS = $(wildcard main/*.c lib/*.c io/*.c $(TOPDIR)/libEZmock/*/*.c $(TOPDIR)/libEZmock/prand/src/*.c)

all:
	$(CC) $(CFLAGS) -o $(EXEC) $(SRCS) $(LIBS) $(INCL)

clean:
ifneq ($(strip $(wildcard $(EXEC))),)
	rm $(EXEC)
endif

