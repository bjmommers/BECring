#!/bin/bash

# Run 3 BECring simulations and calculate Fisher information

# The rotation rate is taken as the central value, and two
#   extra sims are run with slightly higher and lower
#   rotation rates in order to allow central finite-differencing
#   for calculation of Fisher information


# Set default parameters (arguments will override if present)
# Animation-related arguments apply only to the simulation
#   with the central rotation value

infile="temp-input.txt"
outfile="output.txt"
#outdir should be defined here, but is replaced later
outdir="data/fisher-output-error"
rootdir="data/fisher-test"
animate=0
animaterealtime=0
animateimagtime=0
denphase=0
realimag=0
video=0
keepframes=0
batch=0
# add the -qfi flag to compute the quantum Fisher info (comp. expensive)
computeqfi=0

# Default delta: the value to offset the rotation rate by
#   for finite-differencing
delta=0.01

# check arguments for input, output, movies, batch mode
################################################################################
# This block adapted from Drew Stokes' blog
# https://medium.com/@Drew_Stokes/bash-argument-parsing-54f3b81a6a8f

PARAMS=""

while (( "$#" )); do
  case "$1" in
    -i|--input)
      # do nothing, infile must be 'temp-input.txt'
      #infile=$2
      shift 2
      ;;
    -o|--output)
      outfile=$2
      shift 2
      ;;
    -d|--dir)
      rootdir="data/$2"
      shift 2
      ;;
    -D|--delta|--Delta)
      delta="$2"
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
      keepframes=0
      shift
      ;;
    -b|batch)
      batch=1
      shift
    #batch mode currently unused
      ;;
    -qfi)
      computeqfi=1
      shift
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
#LOOP
# By looping over -1,0,1; we can use the iterator $j to modify out rotation rate
for j in {-1..1}
do
    # Find the original rotation rate, omega
    oldomega=`sed -n '/omegareal/{n;p;}' fisher-input.txt`
    # Calculate omega for this iteration
    newomega=`bc <<< "$oldomega + ($delta * $j)"`
    # bc doesn't put in a leading zero, so we manually include it in our sed cmd
    # Set the appropriate output directory
    case $j in
        -1)
            outdir="$rootdir/fisher-output-minus"
            ;;
        0)
            outdir="$rootdir/fisher-output-omega"
            ;;
        1)
            outdir="$rootdir/fisher-output-plus"
            ;;
        *)
            echo "ERROR: unexpected Fisher loop value!"
            exit
    esac



    #RUN

    #Create output directories
    if [ -d "$outdir" ]; then
        rm -rf "$outdir"
    fi
    mkdir -p "$outdir"
    mkdir -p "$outdir/frames"
    mkdir -p "$outdir/iframes"
    # root directory 
    cd "$rootdir"

    # Create a temporary input file from the provided fisher-input.txt
    touch temp-input.txt
    cat ../../fisher-input.txt > temp-input.txt
    # Replace the rotation rate with the calculated value
    sed -i "/omegareal/{n;s/.*/0$newomega/;}" temp-input.txt
    sed -i "/omegaimag/{n;s/.*/0$newomega/;}" temp-input.txt

    # Store delta in a file for reference
    echo "$delta" > "fisher-delta.dat"


    cd ../..
    #copy input file to output dir (for reference)
    cp fisher-input.txt "$outdir/input.txt"

    #run the code
    ./src/becring "$rootdir/$infile" "$outfile" "$outdir" >& "$outdir/stdout.txt" 2> "$outdir/stderr.txt" ;



    ################################################################################
    #PLOTTING

    #copy plotting scripts across
    cp ./src/histplot.gpi ./src/histcompare.gpi ./src/phaseplot-single-rotating.gpi ./src/densitytime.gpi ./src/makemovies.sh ./src/frameplot.gpi ./src/reimplot.gpi ./src/fullplot.gpi ./src/modeplot.gpi ./src/phaseplot.gpi ./src/phaseplot-single.gpi ./src/statplot.gpi ./src/trajectory-plot.gpi "$outdir/"

    #switch to output directory for plotting
    cd $outdir
    

    # produce statplot
    gnuplot -c statplot.gpi $outfile "output" "1D simulation" 2> /dev/null
    gnuplot -c trajectory-plot.gpi $outfile 2> /dev/null

    # produce phaseplot
    if [ -f "phasext.dat" ]; then
        gnuplot -c phaseplot.gpi 2> /dev/null
        gnuplot -c phaseplot-single.gpi 2> /dev/null
    fi

    if [ -f "modes.dat" ]; then
        gnuplot -c modeplot.gpi 2> /dev/null
    fi

    if [ -f "modeampmatrix-abs.dat" ]; then
        gnuplot -c histplot.gpi "modeampmatrix-abs.dat" 2> /dev/null
    fi

    # generate animations if -m flag set ONLY for central value of omega
    # This case statement is obscure as hell, and I apologise
    if [ $animate -eq 1 ] && [ $j -eq 0 ]; then
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

    # revert back to home folder
    if [ -e "makemovies.sh" ]; then
        cd ..
    fi

    # Return to BECring directory to start loop again
    cd ../..

done

# Calculate Fisher info using Python script
python3 src/fisher-info.py "$rootdir" "$computeqfi"

# Plot classical fisher information
gnuplot -c src/cfi-plot.gpi "$rootdir"

if [ $computeqfi -eq 1 ]; then
    gnuplot -c src/qfi-plot.gpi "$rootdir";
fi


#delete plots unless -f flag set
# This has to wait until the very end since frames are needed
# to calculate the quantum Fisher information
if [ $keepframes -eq 0 ]; then
    if [ -e "$rootdir/fisher-output-minus/frames/frame00000.dat" ]; then
        rm -rf "$rootdir/fisher-output-minus/frames";
    fi
    if [ -e "$rootdir/fisher-output-minus/iframes/iframe00000.dat" ]; then
        rm -rf "$rootdir/fisher-output-minus/iframes";
    fi
    if [ -e "$rootdir/fisher-output-omega/frames/frame00000.dat" ]; then
        rm -rf "$rootdir/fisher-output-minus/frames";
    fi
    if [ -e "$rootdir/fisher-output-omega/iframes/iframe00000.dat" ]; then
        rm -rf "$rootdir/fisher-output-minus/iframes";
    fi
    if [ -e "$rootdir/fisher-output-plus/frames/frame00000.dat" ]; then
        rm -rf "$rootdir/fisher-output-minus/frames";
    fi
    if [ -e "$rootdir/fisher-output-plus/iframes/iframe00000.dat" ]; then
        rm -rf "$rootdir/fisher-output-minus/iframes";
    fi
fi
    
# all done
# Uncomment the following line to send a desktop notification on completion
#notify-send "BECring: calculation complete"
exit 0
