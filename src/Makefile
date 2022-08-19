# ----- Socorro multi-machine makefile ------------------------------- #

SHELL	= /bin/bash
PYTHON	= python

# ----- Definitions -------------------------------------------------- #

FSRC	= $(shell ls *.F90)
CSRC	= $(shell ls *.c)

EXE	= socorro
OBJ	= $(FSRC:.F90=.o) $(CSRC:.c=.o)

OBJ_DIR	= Obj_build
LIB_DIR	= $(@:lib-%=%)

# ----- Help message ------------------------------------------------- #

help	:
	@ echo ''
	@ echo 'Execute any of the following make commands from the src dir:'
	@ echo ''
	@ echo 'make doc                     Create the documentation files for socorro'
	@ echo 'make tar                     Create socorro_src.tar.gz for the src dir'
	@ echo ''
	@ echo 'make clean                   Delete the object files in src'
	@ echo 'make purge                   Delete all build related files in src'
	@ echo ''
	@ echo 'make lib-package             Print help message to download/build a package in ../lib'
	@ echo 'make lib-package args="..."  Download/build a package in ../lib'
	@ echo ''
	@ echo 'make depends                 Build the Fortran module dependency list'
	@ echo 'make machine                 Build socorro using Makefile.machine settings'
	@ echo ''
	@ echo 'where Makefile.machine is one of these files from src/MAKE/OPTIONS:'
	@ echo ''
	@ ( cd MAKE/OPTIONS/ ; ls -1 Makefile.* 2>/dev/null ; cd ../ )
	@ echo ''
	@ echo '... or one of these files from src/MAKE/MACHINES'
	@ echo ''
	@ ( cd MAKE/MACHINES/ ; ls -1 Makefile.* 2>/dev/null ; cd ../ )
	@ echo ''
	@ echo '... or one of these files from src/MAKE/MINE'
	@ echo ''
	@ ( cd MAKE/MINE/ ; ls -1 Makefile.* 2>/dev/null ; cd ../ )
	@ echo ''
	@ echo '... or a new Makefile.machine can be created and placed in src/MAKE/MINE for use'
	@ echo ''

# ----- Create a tar file of the src dir ----------------------------- #

tar	:
	@ echo 'Creating socorro_src.tar.gz ...'
	@ make  purge && cd ../ && tar -czf socorro_src.tar.gz src
	@ echo 'Finished socorro_src.tar.gz'

# ----- Cleanup build files ------------------------------------------ #

clean	:
	@ echo 'Removing the object files in src ...'
	@ rm -rf $(OBJ_DIR)

purge	:
	@ echo 'Removing all build related files in src ...'
	@ rm -rf $(OBJ_DIR) $(EXE) Makefile.depends

# ----- Library compilation ------------------------------------------ #

lib-%   :
	@ if [ -e ../lib/$(LIB_DIR)/Install.py ] ; then \
	  echo -e "Installing ../lib/$(@:lib-%=%)\n" ; \
	  ( cd ../lib/$(LIB_DIR) ; $(PYTHON) Install.py $(args) ) ; \
	else \
	  echo "Install.py file for lib $(@:lib-%=%) does not exist" ; \
	fi ; touch main.F90

# ----- Socorro compilation ------------------------------------------ #

f90mods	:
	@ rm -rf Makefile.depends ; MAKE/TOOLS/f90makedep.pl $(FSRC)

depends	:
	@ echo 'Creating the module dependency list ...'
	@ make f90mods

.DEFAULT:
	@ test -f MAKE/MACHINES/Makefile.$@ -o \
	       -f MAKE/MINE/Makefile.$@ -o \
	       -f MAKE/OPTIONS/Makefile.$@
	@ if [ ! -d $(OBJ_DIR) ] ; then mkdir $(OBJ_DIR) ; fi
	@ echo 'Gathering build information ...'
	@ if [ -f MAKE/MACHINES/Makefile.$@ ] ; \
	  then cp MAKE/MACHINES/Makefile.$@ $(OBJ_DIR)/Makefile ; fi
	@ if [ -f MAKE/MINE/Makefile.$@ ] ; \
	  then cp MAKE/MINE/Makefile.$@ $(OBJ_DIR)/Makefile ; fi
	@ if [ -f MAKE/OPTIONS/Makefile.$@ ] ; \
	  then cp MAKE/OPTIONS/Makefile.$@ $(OBJ_DIR)/Makefile ; fi
	@ make f90mods ; mv Makefile.depends $(OBJ_DIR)
	@ echo 'Compiling socorro with Makefile.$@ settings ...'
	@ cd $(OBJ_DIR) ; $(MAKE) "OBJ = $(OBJ)" "EXE = ../$(EXE)"

# ----- End of the makefile ------------------------------------------ #
