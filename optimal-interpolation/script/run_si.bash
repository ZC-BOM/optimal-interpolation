#!/bin/bash

##------------------------------------------------------------------------------
## Source environment. Set this to where env.sh is stored.
##------------------------------------------------------------------------------

source PATH_TO_DIR/env.sh

##------------------------------------------------------------------------------
## Source Configurations. Set this to where config.conf is stored.
##------------------------------------------------------------------------------

source PATH_TO_DIR/config/config.conf

export PYTHONPATH=$PYTHONPATH:$script_dir/utils

##------------------------------------------------------------------------------
## For logging.
##------------------------------------------------------------------------------##

LOGTIME=$(date '+%Y%m%d-%H%M%S')
LOGFILE=
LOGLIST=''
function log
{
    printf "$(date --rfc-3339=seconds)  :  $*\n"
    [[ "$LOGFILE" != "" ]] && printf "$(date --rfc-3339=seconds)  :  $*\n" >> $LOGFILE
}

##------------------------------------------------------------------------------
## Run si_satellite_grid.
##------------------------------------------------------------------------------

function si_satellite_grid
{
echo "Running: si_satellite_grid"
LOGFILE=${log_dir}/si_satellite_grid_${LOGTIME}.log
LOGLIST="${LOGLIST} ${LOGFILE}"
for extent in $SI_EXTENT_LIST
do
    export EXTENT=$extent
    for ACCUM_MON in '1'
    do
        for MAP_TYPE in 'gsmap'
            do
                export ACCUM_MON
                export MAP_TYPE
                python ${script_dir}/si_satellite_grid.py -d ${INPUT_DATE} >> $LOGFILE 2>&1 || echo "Error in si_satellite.py"
            done
    done
    unset EXTENT
done
}

##------------------------------------------------------------------------------
## Commands
##------------------------------------------------------------------------------

# Default command list
export cmdlist=""

while getopts hd: c
do
    case $c in
        h)     usage;;
        d)     INPUT_DATE="${OPTARG:-""}";;
        \?)    usage;;
    esac
done

shift `expr $OPTIND - 1`

if [[ $INPUT_DATE != "" ]]; then
    echo "INPUT_DATE: " $INPUT_DATE
    export INPUT_DATE
else
    echo "NOTE: INPUT_DATE is not specified"
fi

# Custom command list
[[ $# > 0 ]] && cmdlist="$*"

for cmd in $cmdlist
do
    case $cmd in
        si_satellite_grid)                  si_satellite_grid;;
        *)                                  printf "Unknown cmd $cmd\n"
                                            usage;;
    esac
done
exit




