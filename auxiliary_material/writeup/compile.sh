f='main'
latex $f.tex
latex $f.tex
latex $f.tex
bibtex $f
latex $f.tex
latex $f.tex
latex $f.tex
dvips $f.dvi
ps2pdf $f.ps

rm $f.aux $f.bbl $f.blg $f.log $f.out $f.synctex* ${f}Notes.bib $f.dvi $f.ps

