# Makefile for libEZmock

TOPDIR = ..
include $(TOPDIR)/options.mk

LIBS = -lm

# Settings for FFTW
ifneq ($(strip $(FFTW_DIR)),)
  LIBS += -L$(strip $(FFTW_DIR))/lib
  INCL += -I$(strip $(FFTW_DIR))/include
endif

# Settings for the precision
ifeq ($(strip $(SINGLE_PREC)), T)
  LIBS += -DSINGLE_PREC -lfftw3f
  ifeq ($(strip $(USE_OMP)), T)
    LIBS += -DOMP -fopenmp -lfftw3f_omp
  endif
else
  LIBS += -lfftw3
  ifeq ($(strip $(USE_OMP)), T)
    LIBS += -DOMP -fopenmp -lfftw3_omp
  endif
endif

ifeq ($(TARGET), libEZmock.so)
  CFLAGS += -fPIC
endif

INCL += -IEZmock -Imath -Iprand/src/header
SRCS = $(wildcard EZmock/*.c math/*.c prand/src/*.c)
OBJS = $(patsubst %.c, %.o, $(SRCS))

# List the compiled files
FILES = $(strip $(wildcard $(OBJS) $(addprefix libEZmock,.so .a)))

# Set the installation path
IS_ROOT := $(if $(patsubst /%,,$(INSTALL_DIR)),,yes)
IS_HOME := $(if $(patsubst ~%,,$(INSTALL_DIR)),,yes)
ifeq ($(strip $(IS_ROOT)$(IS_HOME)$(IS_NETWORK)$(IS_DRIVE)),yes)
  PREFIX := $(INSTALL_DIR)
else
  PREFIX := $(TOPDIR)/$(INSTALL_DIR)
endif

all: $(TARGET)

libEZmock.so: $(OBJS)
	$(CC) $(CFLAGS) -shared -o $@ $^ $(INCL) $(LIBS)

libEZmock.a: $(OBJS)
	ar rcs $@ $^

%.o: %.c
	$(CC) $(CFLAGS) -c $^ -o $@ $(INCL) $(LIBS)

install:
ifneq ($(strip $(wildcard $(addprefix libEZmock,.so .a))),)
	install -d $(PREFIX)/lib/
  ifneq ($(strip $(wildcard libEZmock.so)),)
	install -m 755 libEZmock.so $(PREFIX)/lib/
  endif
  ifneq ($(strip $(wildcard libEZmock.a)),)
	install -m 644 libEZmock.a $(PREFIX)/lib/
  endif
	install -d $(PREFIX)/include/
	install -m 644 EZmock/EZmock.h $(PREFIX)/include/
endif

clean:
ifneq ($(FILES),)
	rm $(FILES)
endif

