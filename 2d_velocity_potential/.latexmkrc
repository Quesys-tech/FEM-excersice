$lualatex = 'lualatex %O -synctex=1 -interaction=nonstopmode --shell-escape %S';
$biber = 'biber %O --bblencoding=utf8 -u -U --output_safechars %B';
$bibtex = 'upbibtex %O %B';
$makeindex = 'upmendex %O -o %D %S';