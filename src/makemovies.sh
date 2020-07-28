#!/bin/bash
# real time wavefunction plots

# set parameters
fullplot=0
reim=0
video=0
realtime=0
imagtime=0

# check arguments
################################################################################
# This block adapted from Drew Stokes' blog
# https://medium.com/@Drew_Stokes/bash-argument-parsing-54f3b81a6a8f

PARAMS=""

while (( "$#" )); do
  case "$1" in
    -f|--full)
      fullplot=1
      shift 1
      ;;
    -r|--reim)
      reim=1
      shift 1
      ;;
    -v|--video)
      video=1
      shift 1
      ;;
    -rt|--realtime)
      realtime=1
      shift 1
      ;;
    -it|--imagtime)
      imagtime=1
      shift 1
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
echo "Making movies..."

# set which type of plot to do
if [ $fullplot -eq 1 ]; then
    script="../fullplot.gpi" # includes density & phase
elif [ $reim -eq 1 ]; then
    script="../reimplot.gpi" # includes real & imaginary parts
else
    script="../frameplot.gpi" # density and potential only [default]
fi

if [ $realtime -eq 1 ]; then
    cd frames/;
    for i in frame*.dat; do cp $i data; gnuplot -c $script > $i.png 2>/dev/null ; done
    ffmpeg -y -loglevel quiet -r 60 -f image2 -s 1920x1080 -i frame%05d.dat.png -vcodec libx264 -crf 25  -pix_fmt yuv420p animation.mp4
    rm frame*.png
    rm data
    mv animation.mp4 "../animation.mp4"
    cd ..;
fi


#make imaginary time animation if data files exist
if [ $imagtime -eq 1 ]; then
	cd iframes/;
	for i in iframe*.dat; do cp $i data; gnuplot $script > $i.png 2>/dev/null ; done
    ffmpeg -y -loglevel quiet -r 60 -f image2 -s 1920x1080 -i iframe%05d.dat.png -vcodec libx264 -crf 25  -pix_fmt yuv420p ianimation.mp4
	rm iframe*.png
	rm data
	mv ianimation.mp4 "../ianimation.mp4"
    cd ..
fi

#if [ $realtime -eq 1 ]; then
#    if [ $imagtime -eq 1 ]; then
#        convert ianimation.gif animation.gif combined.gif
#        if [ $video -eq 1 ]; then
#            ffmpeg -v quiet -i combined.gif combined.mp4
#        fi
#    fi
#fi
exit 0
