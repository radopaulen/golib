# This makefile creates symbolic links to the header files in $(incpath)

include $(srcpath)/makeoptions.mk

#####

libobjs = 
incobjs = fprelax.hpp fprgurobi.hpp fprcplex.hpp

#####

lib:

#####

install: dispInstall
	@for INC in $(incobjs); do \
		if test ! -e $(incpath)/$$INC; then \
			echo creating symbolic link to header file $$INC; \
			cd $(incpath); ln -s $(fppath)/$$INC $$INC; \
		fi; \
	done
	@echo

dispInstall:
	@echo
	@(echo '***Installing polyhedral relaxation components***')
	@echo

#####

clean:

#####

cleaninstall: dispCleanInstall
	cd $(incpath) ; rm -i $(incobjs)

dispCleanInstall:
	@echo
	@(echo '***Uninstalling polyhedral relaxation components***')
	@echo
