
include ../make.conf

all:
	perl makedocs ../src/*_mod.f;\
        perl makedocs2 ../src/*_mod.f90; \
	mv ../src/*.doc .
	(cd LATEX; make all;)

clean:
	rm -f *.doc
	(cd LATEX; make clean)
