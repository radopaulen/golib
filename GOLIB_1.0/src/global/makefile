# This makefile creates symbolic links to the header files in $(incpath)

version   = 1.0

include $(srcpath)/makeoptions.mk

#####

libobjs = 
incobjs = nlpstruct.hpp mcipopt.hpp dostruct.hpp mccvodes.hpp doseq.hpp \
          sbb.hpp nlpsbb.hpp dosbb.hpp setinv.hpp pestruct.hpp odegpe.hpp \
	  structure.hpp
incdeps = 

#####

lib:

#####

install : dispInstall
	@for INC in $(incobjs); do \
		if test ! -e $(incpath)/$$INC; then \
			echo creating symbolic link to header file $$INC; \
			cd $(incpath); ln -s $(gopath)/$$INC $$INC; \
		fi; \
	done
	@if test ! -e $(libpath)/$(libtarget); then \
		echo creating symolic link to library $(libtarget); \
		cd $(libpath); ln -s $(gopath)/$(libtarget) $(libtarget); \
	fi
	@echo

dispInstall :
	@echo
	@(echo '***Installing global optimization library (ver.' $(version) ')***')
	@echo

#####

clean :

#####

cleaninstall : dispCleanInstall
	cd $(incpath) ; rm -i $(incobjs)
	cd $(libpath) ; rm -i $(libtarget)

dispCleanInstall :
	@echo
	@(echo '***Uninstalling global optimization library (ver.' $(version) ')***')
	@echo
