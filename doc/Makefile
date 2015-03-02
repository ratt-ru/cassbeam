TARGETS = cassbeam.ps

all: $(TARGETS)

cassbeam.dvi: cassbeam.tex
	latex cassbeam.tex
	latex cassbeam.tex

cassbeam.ps: cassbeam.dvi
	dvips cassbeam.dvi -o cassbeam.ps

clean:
	-rm -f *.aux *.dvi *.log $(TARGETS)
