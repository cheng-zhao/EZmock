CC = gcc
CFLAGS = -std=c99 -O3 -Wall -flto

SINGLE_PREC = T
USE_OMP = T
WITH_FITS = F  # T for enabling FITS-format outputs

# Available targets: EZmock, libEZmock.so, libEZmock.a
TARGETS = EZmock libEZmock.so

# Directory for installing the header and library files
INSTALL_DIR = ./build

# Directory for FFTW-3 and CFITSIO (>=4.2) libraries
# The corresponding header files should be in
#     $(FFTW_DIR)/include    and    $(CFITSIO_DIR)/include
# The library files should be in
#     $(FFTW_DIR)/lib        and    $(CFITSIO_DIR)/lib
FFTW_DIR = 
CFITSIO_DIR = 
