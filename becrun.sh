#!/bin/bash

# Script to run BECring 1D simulations

#set default parameters (arguments will override if present)

infile="input.txt"
outfile="output.txt"
outdir="data/output-$(date +%Y%m%d-%H%M%S)"
animate=0
animaterealtime=0
animateimagtime=0
denphase=0
realimag=0
video=0
keepframes=0
batch=0

# check arguments for input, output, movies, batch mode
################################################################################
# This block adapted from Drew Stokes' blog
# https://medium.com/@Drew_Stokes/bash-argument-parsing-54f3b81a6a8f

PARAMS=""

while (( "$#" )); do
  case "$1" in
    -i|--input)
      infile=$2
      shift 2
      ;;
    -o|--output)
      outfile=$2
      shift 2
      ;;
    -d|--dir)
      outdir="data/$2"
      shift 2
      ;;
    -D|--DIR)
      outdir="$2"
      shift 2
      ;;
    -m|--movies)
      animate=1
      animaterealtime=1
      animateimagtime=1
      shift
      ;;
    -ma|--fullmovies)
      animate=1
      denphase=1
      shift
      ;;
    -mr|--realtime)
      animate=1
      animaterealtime=1
      denphase=1
      shift
      ;;
    -mi|--imagtime)
      animate=1
      animateimagtime=1
      denphase=1
      shift
      ;;
    -mri|--realimag)
      animate=1
      animaterealtime=1
      animateimagtime=1
      realimag=1
      shift
      ;;
    -v|--video)
      video=1
      shift
      ;;
    -f|--frames)
      keepframes=1
      shift
      ;;
    -b|batch)
      batch=1
      shift
    #batch mode currently unused
      ;;
    --) # end argument parsing
      shift
      break
      ;;
    -*|--*=) # unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # preserve positional arguments
      PARAMS="$PARAMS $1"
      shift
      ;;
  esac
done

# set positional arguments in their proper place
# not really needed in this implementation, but kept for completeness
eval set -- "$PARAMS"

################################################################################

#RUN

#Create output directories
if [ -d "$outdir" ]; then
    rm -rf "$outdir"
fi
mkdir "$outdir"
mkdir "$outdir/frames"
mkdir "$outdir/iframes"

#copy input file to output dir (for reference)
cp $infile "$outdir/input.txt"

#run the code
#./src/becring "input.txt" "output.txt" "$outdir";
./src/becring "$infile" "$outfile" "$outdir" >& "$outdir/stdout.txt" 2> "$outdir/stderr.txt" ;

#> "$outdir/stdout.txt"

################################################################################
#PLOTTING

#copy plotting scripts across
cp ./src/histplot.gpi ./src/histcompare.gpi ./src/phaseplot-single-rotating.gpi ./src/densitytime.gpi ./src/makemovies.sh ./src/frameplot.gpi ./src/reimplot.gpi ./src/fullplot.gpi ./src/modeplot.gpi ./src/phaseplot.gpi ./src/phaseplot-single.gpi ./src/statplot.gpi ./src/trajectory-plot.gpi "$outdir/"

#switch to output directory for plotting
cd $outdir

# produce statplot
gnuplot -c statplot.gpi $outfile "output" "1D simulation"
gnuplot -c trajectory-plot.gpi $outfile

# produce phaseplot
if [ -f "phasext.dat" ]; then
    gnuplot -c phaseplot.gpi
    gnuplot -c phaseplot-single.gpi
fi

if [ -f "modes.dat" ]; then
    gnuplot -c modeplot.gpi
fi

if [ -f "modeampmatrix-abs.dat" ]; then
    gnuplot -c histplot.gpi "modeampmatrix-abs.dat"
fi

# generate animations if -m flag set
if [ $animate -eq 1 ]; then
    case "$denphase$realimag$video$animaterealtime$animateimagtime" in
        "00000")
            ./makemovies.sh
            ;;
        "00100")
            ./makemovies.sh -v
            ;;
        "10000")
            ./makemovies.sh -f
            ;;
        "10000")
            ./makemovies.sh -r
            ;;
        "10100")
            ./makemovies.sh -f -v
            ;;
        "01100")
            ./makemovies.sh -r -v
            ;;
        "00001")
            ./makemovies.sh -it
            ;;
        "00101")
            ./makemovies.sh -v -it
            ;;
        "10001")
            ./makemovies.sh -f -it
            ;;
        "10001")
            ./makemovies.sh -r -it
            ;;
        "10101")
            ./makemovies.sh -f -v -it
            ;;
        "01101")
            ./makemovies.sh -r -v -it
            ;;
        "00010")
            ./makemovies.sh -rt
            ;;
        "00110")
            ./makemovies.sh -v -rt
            ;;
        "10010")
            ./makemovies.sh -f -rt
            ;;
        "10010")
            ./makemovies.sh -r -rt
            ;;
        "10110")
            ./makemovies.sh -f -v -rt
            ;;
        "01110")
            ./makemovies.sh -r -v -rt
            ;;
        "00011")
            ./makemovies.sh -it -rt
            ;;
        "00111")
            ./makemovies.sh -v -it -rt
            ;;
        "10011")
            ./makemovies.sh -f -it -rt
            ;;
        "10011")
            ./makemovies.sh -r -it -rt
            ;;
        "10111")
            ./makemovies.sh -f -v -it -rt
            ;;
        "01111")
            ./makemovies.sh -r -v -it -rt
            ;;
        "*")
            echo "invalid arguments to makemovies.sh"
            exit 1
            ;;
    esac
fi

#delete plots unless -f flag set
if [ $keepframes -eq 0 ]; then
    if [ -e "frames/frame00000.dat" ]; then
        rm -rf frames;
    fi
    if [ -e "iframes/iframe00000.dat" ]; then
        rm -rf iframes;
    fi
fi
# revert back to home folder
if [ -e "makemovies.sh" ]; then
    cd ..
fi

# all done
#notify-send "BECring: calculation complete"
exit 0
