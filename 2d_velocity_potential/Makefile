LATEX = lualatex

all: article.pdf

article.pdf: article.tex NACA2412.pdf
	latexmk -lualatex $<

NACA2412.pdf: plot_NACA2412.jl
	julia --project=./ $<

.PHONY:clean-temp
clean-temp:
	latexmk -c article.tex
	rm -f article.run.xml
	rm -f article.bbl
	rm -rf _minted-article/

.PHONY: clean
clean: clean-temp
	latexmk -C article.tex