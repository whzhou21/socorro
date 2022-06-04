# ----- Socorro Makefile --------------------------------------------- #

SHELL	= /bin/bash
PYTHON	= python

# ----- Definitions -------------------------------------------------- #

CSRC	= ctof_io.c fft_3d.c

EXE	= socorro

SRC	= $(shell ls *.F90)
OBJ	= $(SRC:.F90=.o) $(CSRC:.c=.o)

OBJ_DIR	= Obj_build
LIB_DIR	= $(@:lib-%=%)

# ----- Cleanup ------------------------------------------------------ #

clean	:
	@ echo 'Removing the object files in src ...' ; sleep 0.5s
	@ rm -rf $(OBJ_DIR)

cleanall:
	@ echo 'Removing all build related files in src ...' ; sleep 0.5s
	@ rm -rf $(OBJ_DIR) $(EXE) Makefile.depends

# ----- Library installation ----------------------------------------- #

lib-%   :
	@ if [ -e ../lib/$(LIB_DIR)/Install.py ] ; then \
	  echo -e "Installing ../lib/$(@:lib-%=%)\n" ; sleep 0.5s ; \
	  ( cd ../lib/$(LIB_DIR) ; $(PYTHON) Install.py $(args) ) ; \
	else \
	  echo "Install.py file for lib $(@:lib-%=%) does not exist" ; \
	fi ; touch socorro.F90

# ----- Socorro compilation ------------------------------------------ #

f90mods	:
	@ rm -rf Makefile.depends ; MAKE/TOOLS/f90makedep.pl $(SRC)

depends	:
	@ echo 'Creating the module dependency list ...' ; sleep 0.5s
	@ make f90mods

.DEFAULT:
	@ test -f MAKE/MINE/Makefile.$@ -o -f MAKE/OPTIONS/Makefile.$@
	@ if [ ! -d $(OBJ_DIR) ] ; then mkdir $(OBJ_DIR) ; fi
	@ echo 'Gathering build information ...' ; sleep 0.5s
	@ if [ -f MAKE/MINE/Makefile.$@ ]; \
	  then cp MAKE/MINE/Makefile.$@ $(OBJ_DIR)/Makefile; fi
	@ if [ -f MAKE/OPTIONS/Makefile.$@ ]; \
	  then cp MAKE/OPTIONS/Makefile.$@ $(OBJ_DIR)/Makefile; fi
	@ make f90mods ; mv Makefile.depends $(OBJ_DIR)
	@ echo 'Compiling socorro with Makefile.$@ settings ...' ; sleep 0.5s
	@ cd $(OBJ_DIR) ; $(MAKE) "OBJ = $(OBJ)" "EXE = ../$(EXE)"

# ----- Create a tarball of the src dir ------------------------------ #

tar	:
	@ echo 'Creating socorro_src.tar.gz ...'
	@ make cleanall && cd ../ && tar -czf socorro_src.tar.gz src
	@ echo 'Finished socorro_src.tar.gz'

# ----- Help message ------------------------------------------------- #

help	:
	@ echo ''
	@ echo 'Execute any of the following make commands from the src dir:'
	@ echo ''
	@ echo 'make doc                     Create the documentation files for socorro'
	@ echo 'make tar                     Create socorro_src.tar.gz for the src dir'
	@ echo ''
	@ echo 'make clean                   Delete the object files in src'
	@ echo 'make cleanall                Delete all build related files in src'
	@ echo ''
	@ echo 'make lib-package             Print help message to download/build a package in ../lib'
	@ echo 'make lib-package args="..."  Download/build a package in ../lib'
	@ echo ''
	@ echo 'make depends                 Build the Fortran module dependency list'
	@ echo 'make machine                 Build socorro using Makefile.machine settings'
	@ echo ''
	@ echo 'where Makefile.machine is one of these files from src/MAKE/OPTIONS:'
	@ echo ''
	@ ( cd MAKE/OPTIONS/ && ls -1 Makefile.* && cd ../ )
	@ echo ''
	@ echo '... or one of these files from src/MAKE/MINE'
	@ echo ''
	@ files="`ls MAKE/MINE/Makefile.* 2>/dev/null`"; \
	  for file in $$files; do head -1 $$file; done
	@ echo ''
	@ echo '... or a new Makefile.machine can be created and placed in src/MAKE/MINE for use'
	@ echo ''

# ----- End of the Makefile ------------------------------------------ #
