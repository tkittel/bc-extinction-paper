#!/bin/bash
NAME=bcextnupdate
PREFIX="===> COMPILE: "
#always run once + bibtex + one more time. After that, only run as long as it says "Rerun" in the output.
BAD=0
RUN=1
pdflatex -halt-on-error ${NAME}.tex 2>&1 | tee output.log
if [ ${PIPESTATUS[0]} != 0 ]; then
    echo
    echo "${PREFIX}Errors encountered but continuing anyway (in case it was bibtex errors)."
    echo
fi
bibtex ${NAME}
if [ $? != 0 ]; then
    RUN=0
    BAD=1
    echo
    echo "${PREFIX}Errors encountered during bibtex processing! Aborting!"
    echo
fi
N=1
while [  $RUN == 1 ]; do
    if [ $N == 10 ]; then
        echo
        echo "${PREFIX}Infinite loop detected. Aborting!!"
        BAD=1
        break
    fi
    pdflatex -halt-on-error ${NAME}.tex 2>&1 | tee output.log
    if [ ${PIPESTATUS[0]} != 0 ]; then
        RUN=0
        echo "${PREFIX}Errors encountered! Aborting!"
        BAD=1
        break
    fi
    egrep -q 'LaTeX.*Rerun|may have changed' output.log
    if [ $? == 1 ]; then
        echo
        echo "${PREFIX}Done normally after "$((N+1))" compilations!"
        break
    fi
    N=$((N+1))
done
echo
NERRWARN=`egrep -i 'warn|error|No file |^!' output.log|wc -l`
if [ $NERRWARN == 0 ]; then
    echo "${PREFIX}No errors or warnings detected during the last run!"
else
    echo "${PREFIX}Errors or warnings detected during last run:"
    echo
    egrep -i 'warn|error|No file |^!|^l\.' output.log
fi
return $BAD 2>/dev/null || exit $BAD
