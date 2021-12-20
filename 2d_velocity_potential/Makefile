LATEX = lualatex

all: article.pdf

article.pdf: article.tex
	latexmk -lualatex $<

.PHONY:clean-temp
clean-temp:
	latexmk -c article.tex
	rm -f article.run.xml
	rm -f article.bbl
	rm -rf _minted-article/

.PHONY: clean
clean: clean-temp
	rm article.pdf