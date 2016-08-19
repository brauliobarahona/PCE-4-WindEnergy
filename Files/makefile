all: PCE-for-WE_Presentation.md
		pandoc PCE-for-WE_Presentation.md -V geometry:margin=0.7in -V fontsize=10pt  -t beamer -f markdown -s -S --bibliography=bib.bib --csl=CitationStyle.cls -H head.tex -o PCE-for-WE_Presentation.pdf
		open PCE-for-WE_Presentation.pdf