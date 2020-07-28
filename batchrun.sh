#!/bin/bash
# Script to run batches of BECring input files

#set default parameters for animation
animate=0
denphase=0
realimag=0
video=0
keepframes=0

#array for directly passing arguments to ./becrun.sh
args=()

# check arguments for input, output, movies, batch mode
################################################################################
# This block adapted from Drew Stokes' blog
# https://medium.com/@Drew_Stokes/bash-argument-parsing-54f3b81a6a8f

PARAMS=""

while (( "$#" )); do
  case "$1" in
    -m|--movies)
      animate=1
      args+=" -m"
      shift
      ;;
    -ma|--densityphase)
      animate=1
      denphase=1
      args+=" -ma"
      shift
      ;;
    -mr|--realimag)
      animate=1
      realimag=1
      args+=" -mr"
      shift
      ;;
    -v|--video)
      video=1
      args+=" -v"
      shift
      ;;
    -f|--frames)
      keepframes=1
      args+=" -f"
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

# check input files supplied
#if [ ! -e batch/*.in ]; then
    #echo "Error: no input files in batch/ directory!"
    #exit
#fi

# create folder
runtime=$(date +%Y%m%d-%H%M%S)
mkdir "batch/data/$runtime"
cp batch/*.in "batch/data/$runtime/"


for j in batch/data/$runtime/*.in
do
    echo "Now running ./becrun.sh -i $j -D $j-output $args"
    ./becrun.sh -i $j -D $j-output $args
done

exit 0
