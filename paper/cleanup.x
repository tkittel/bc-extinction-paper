#!/usr/bin/env bash
NAME=bcextnupdate
#Step into paper directory:
cd "$( dirname "${BASH_SOURCE[0]}" )"
rm -f x.log ${NAME}.aux ${NAME}.log ${NAME}.toc ${NAME}.out ${NAME}.bbl ${NAME}.dvi ${NAME}.blg ${NAME}.spl output.log texput.log
