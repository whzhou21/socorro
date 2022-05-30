# ----- Socorro Makefile --------------------------------------------- #

SHELL	= /bin/bash
PYTHON	= python

# ----- Definitions -------------------------------------------------- #

CSRC	= ctof_io.c fft_3d.c remap_2d.c

EXE	= socorro

SRC	= $(shell ls *.F90)
OBJ	= $(SRC:.F90=.o) $(CSRC:.c=.o)

OBJ_DIR	= Obj_build
LIB_DIR	= $(@:lib-%=%)

# ----- Cleanup ------------------------------------------------------ #

purge	:
	@( echo 'Removing all build related files in src and lib...' ; sleep 1s )
	@( rm -rf $(OBJ_DIR) $(EXE) Makefile.depends )

clean	:
	@( echo 'Removing object files in src ...' ; sleep 1s )
	@( rm -rf $(OBJ_DIR) )

cleanall:
	@( echo 'Removing all build related files in src ...' ; sleep 1s )
	@( rm -rf $(OBJ_DIR) $(EXE) Makefile.depends )

# ----- Compilation -------------------------------------------------- #

depends	:
	@( echo 'Creating the module dependency list ...' ; sleep 1s )
	@( rm -rf Makefile.depends ; ../tools/bin/f90makedep.pl $(SRC) )

.DEFAULT:
	@( echo 'Gathering build information ...' ; sleep 1s )
	@( if [ ! -d $(OBJ_DIR) ] ; then mkdir $(OBJ_DIR) ; fi )
	@( if [ ! -f MAKE/Makefile.$@ ] ; then echo 'Error: Makefile.$@ was not found' ; exit ; fi )
	@( cp MAKE/Makefile.$@ $(OBJ_DIR)/Makefile )
	@( make depends ; mv Makefile.depends $(OBJ_DIR) )
	@( echo 'Compiling socorro with Makefile.$@ settings ...' ; sleep 1s )
	@( cd $(OBJ_DIR) ; $(MAKE) "OBJ = $(OBJ)" "EXE = ../$(EXE)" )

# ----- Create a tarball of the src dir ------------------------------ #

tar	:
	@( echo 'Creating socorro_src.tar.gz ...' )
	@( make cleanall && cd ../ && tar -czf socorro_src.tar.gz src )
	@( echo 'Finished socorro_src.tar.gz')

# ----- Help message ------------------------------------------------- #

help	:
	@echo ''
	@echo 'Execute any of the following make commands:'
	@echo ''
	@echo 'make doc                     Create the documentation files for socorro'
	@echo 'make tar                     Create socorro_src.tar.gz for the src dir'
	@echo ''
	@echo 'make purge                   Delete all build related files in src and lib'
	@echo 'make clean                   Delete the object files in src'
	@echo 'make cleanall                Delete all build related files in src'
	@echo ''
	@echo 'make lib-package             Help for download/build/install a package in lib'
	@echo 'make lib-package args="..."  Download/build/install a package in lib'
	@echo ''
	@echo 'make depends                 Build the Fortran module dependency list'
	@echo 'make machine                 Build socorro using Makefile.machine settings'
	@echo ''
	@echo 'where Makefile.machine is one of these files from src/MAKE:'
	@echo ''
	@( cd MAKE/ && ls -1 Makefile.* && cd ../ )
	@echo ''
	@echo '... or a new Makefile.machine can be created and placed in src/MAKE for use'
	@echo ''

# ----- End of the Makefile ------------------------------------------ #
