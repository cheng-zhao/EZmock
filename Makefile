# Top level Makefile for EZmock

TOPDIR = .
include $(TOPDIR)/options.mk

all: $(TARGETS)

EZmock:
	$(MAKE) -C src

libEZmock.so:
	$(MAKE) -C libEZmock TARGET=$@

libEZmock.a:
	$(MAKE) -C libEZmock TARGET=$@

install:
	$(MAKE) -C libEZmock $@

clean:
	$(MAKE) -C src clean

cleanlib:
	$(MAKE) -C libEZmock clean

cleanall: clean cleanlib
