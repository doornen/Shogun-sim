# the - sign means to not worry about the return value
# the rules

# latex
%.dvi : %.tex
	- @latex -interaction=batchmode $<
	- @bibtex $(basename $<)
	- @latex -interaction=batchmode $<
	- @bibtex $(basename $<)
	- @latex -interaction=batchmode $<
	- @latex -interaction=batchmode $<
	- @latex -interaction=scrollmode $<

# ps from dvi
%.ps : %.dvi
	- @dvips -o $@.tmp $< >/dev/null
	- mv $@.tmp $@

# ps from gnuplot
%.ps : %.gp
	gnuplot $< >$@

# potrace from jpg
%.potrace.eps : %.jpg
	djpeg $< | potrace >$@

# epsi from ps
%.epsi: %.ps
	pstoepsi $< $@

# eps from ps
%.eps: %.ps
	ps2eps -q -l -P -B -s b0 -n <$< >$@

# pdf from ps
%.pdf: %.ps
	ps2pdf $< $@

# ps.gz from ps
%.ps.gz : %.ps
	gzip -c $< > $@

# from tgif.obj to tgif.ps
%.ps: %.obj
	tgif -print -stdout -ps -color $< > $@

# from tgif.obj to tgif.ps
%.eps: %.obj
	tgif -print -stdout -eps -color $< > $@

# from tgif.obj to tgif.txt
%.txt: %.obj
	tgif -print -stdout -text $< > $@

%.jpeg.ps: %.jpeg
	jpeg2ps $< > $@

%.jpg.ps: %.jpg
	jpeg2ps $< > $@

%.jpeg.eps: %.jpeg
	jpeg2ps $< > $@

%.jpg.eps: %.jpg
	jpeg2ps $< > $@

%.gif.eps: %.gif
	giftopnm $< | sam2p - EPS: $@

%.eps: %
	sam2p $< EPS: $@

# the standard targets
ALLTEX=$(patsubst %.tex,%.ps,$(wildcard *.tex))
ALLTEXDVI= $(ALLTEX) $(patsubst %.tex,%.dvi,$(wildcard *.tex))  # also generate AND KEEP dvi files
ALLGP=$(patsubst %.gp,%.eps,$(wildcard *.gp))  # gnuplot (gp) files

ALLALL=$(ALLTEX)   # make all tex files into ps files

.PHONY: all clean force
all: $(ALLALL)

clean:
	- rm -f $(ALLALL) *.dvi *.log *.aux *.bbl *.blg *.toc *.end *.bm

cb.bib:
	wget -O - http://www.mpi-hd.mpg.de/cb/bibtex/$@ >$@

dirname.texinp:
	- echo creating file $<
	echo -n "`hostname`:`pwd`\relax " >$@

talk.dvi: powerdot-hs2.sty

toc:
	@gawk '{print} /^\\input/ {system("cat "$$2);}' <talk.tex | \
           sed -n '/part/ {s/^.*part{\(.*\)}.*/\1/p}; s/^.*begin.slide}{\(.*\)}.*/    \1/p'

#$(foreach f,$(wildcard *.tex),$(shell ./mk_tex_depend.sh $(f)))
depend.mk: $(wildcard *.tex)
	./mk_tex_depend.sh > $@

include depend.mk

