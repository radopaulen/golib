# This is the main makefile for GOLIB

#####

export rootpath = $(shell cd .. ; pwd)
export srcpath  = $(rootpath)/src
export libpath  = $(rootpath)/lib
export incpath  = $(rootpath)/include
export fppath   = $(srcpath)/polyhedral
export gopath   = $(srcpath)/global

#####

lib: dirs fp_lib go_lib

fp_lib:
	(cd $(fppath); make lib)

go_lib:
	(cd $(gopath); make lib)

#####

install: dirs fp_install go_install

dirs :
	@if test ! -e $(incpath)/; then \
		echo creating header files directory; \
		mkdir $(incpath); \
	fi
	@if test ! -e $(libpath)/; then \
		echo creating header files directory; \
		mkdir $(libpath); \
	fi

fp_install:
	(cd $(fppath); make install)

go_install:
	(cd $(gopath); make install)

#####

clean: fp_clean go_clean

fp_clean:
	-(cd $(fppath); make clean)

go_clean:
	-(cd $(gopath); make clean)

#####

cleaninstall: fp_cleaninstall go_cleaninstall

fp_cleaninstall:
	-(cd $(fppath); make cleaninstall)

go_cleaninstall:
	-(cd $(gopath); make cleaninstall)
