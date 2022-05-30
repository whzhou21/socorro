include ../../make.conf

LATEX=latex
DVIPS=dvips
FIGS= figs/biglc_hierarchy.eps \
       figs/complex.eps \
       figs/complexfork.eps \
       figs/complexfork2.eps \
       figs/constraint.eps \
       figs/datamap.eps \
       figs/fork.eps \
       figs/grid.eps \
       figs/handle.eps \
       figs/lc_hierarchy.eps \
       figs/share.eps

all: tutorial.dvi tutorial.ps high_level.dvi high_level.ps

tutorial.ps:  tutorial.dvi
	$(DVIPS) -t letter -o tutorial.ps tutorial.dvi

tutorial.dvi: $(FIGS) modularity.tex memory.tex socorro.tex tutorial.tex 
	$(LATEX) tutorial.tex
	$(LATEX) tutorial.tex

high_level.ps:  high_level.dvi
	$(DVIPS) -t letter -o high_level.ps high_level.dvi

high_level.dvi: high_level.tex
	$(LATEX) high_level.tex
	$(LATEX) high_level.tex

clean: 
	rm -f *.ps *.dvi *.aux *.log *~ *.bak
