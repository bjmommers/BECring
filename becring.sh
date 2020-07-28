#!/bin/bash
# run code
./becring "input.txt" "output.txt";
# plot summary stats
#gnuplot -c statplot.gpi "statplot.eps" "summary stats";
./statplot.gpi "output.txt" "output" "Overlap with t=0 state"
# create movies based on arguments
if [ "$1" = "-m" ]; then
	./makemovies.sh; #auto-makes imag time animation if datafiles exist
fi
#delete plots
if [ -e "frames/frame00000.dat" ]; then
    rm frames/frame*.dat;
#    sleep 1
fi
if [ -e "iframes/iframe00000.dat" ]; then
    rm iframes/iframe*.dat;
fi
# all done
exit 0
