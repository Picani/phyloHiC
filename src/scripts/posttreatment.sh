#!/bin/sh

show_help(){
    cat <<EOF
Usage: ${0##/*} [-RSf] [-e EXPNAME] N ROOT ROOT_SPECIES SPECIES_TREE
Walk throught ROOT, making trees from .phylip files, then root them on
ROOT_SPECIES and compute their support against SPECIES_TREE.
The .phylip files are expected to contains N matrices.

Options:
  -h          Display this help and exit.
  -e EXPNAME  Work only on this experiment.
  -f          Overwrite already existing files.
  -R          Do not perform the rooting and support values steps.
  -S          Do not perform the support values step.
EOF
}

force=false
no_rooting=false
no_support=false
expname="*"

while getopts :hfRSe: opt
do
    case $opt in
        h)
            show_help
            exit 0
            ;;
        f)
            force=true
            ;;
        R)
            no_rooting=true
            no_support=true
            ;;
        S)
            no_support=true
            ;;
        e)
            expname=$OPTARG
            ;;
        \?)
            echo "Unknown option: -$OPTARG"
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument."
            exit 1
            ;;
    esac
done
shift "$((OPTIND-1))"   # Discard the options and sentinel --

# Are all the arguments here?
if [ "$no_support" = true ]; then
    if [ "$no_rooting" = true ]; then
        if [ $# -lt 2 ]; then
            echo "Error: not enough argument."
            show_help
            exit 1
        fi
    else
        if [ $# -lt 3 ]; then
            echo "Error: not enough argument."
            show_help
            exit 1
        fi
    fi
else
    if [ $# -lt 4 ]; then
        echo "Error: not enough argument."
        show_help
        exit 1
    fi
fi

# We check for the needed executables
if [ ! $(command -v fastme) ]
then
    echo "Error: fastme is not available. See http://www.atgc-montpellier.fr/fastme/"
    exit 1
fi

if [ ! $(command -v nw_reroot) ]
then
    echo "Error: nw_reroot is not available. See http://cegg.unige.ch/newick_utils"
    exit 1
fi

if [ ! $(command -v nw_support) ]
then
    echo "Error: nw_support is not available. See http://cegg.unige.ch/newick_utils"
    exit 1
fi


args=("$@")
nb=${args[0]}
root=${args[1]}
species_root=${args[2]}
species_tree=${args[3]}

for d in ${root}/**/${expname}/**/*_${nb}/
do
    for phylip in $d/*.phylip
    do

        # Tree construction
        outname=$(basename $phylip .phylip).nwk

        if [ -f $d/$outname ]
        then
            if [ "$force" = false ]
            then
                echo "Skipping $phylip because its .nwk counterpart already exists."
                continue
            fi
        fi

        nlines=$(grep -c -P '\S' $phylip)
        nspecies=$(head -n1 $phylip)
        nb_datasets=$(echo "scale = 0; $nlines / ($nspecies + 1)" | bc)
        fastme -mI -nB -s -D$nb_datasets -i $phylip -o $d/$outname

        # Re-rooting
        if [ "$no_rooting" = true ]
        then
            continue
        fi

        re_rooutname=$(basename $outname .nwk)_rerooted.nwk
        if [ -f $d/$re_rooutname ]
        then
            if [ "$force" = false ]
            then
                echo "Skipping $outname because its re_rooted counterpart already exists."
                continue
            fi
        fi

        nw_reroot $d/$outname $species_root > $d/$re_rooutname

        # Support
        if [ "$no_support" = true ]
        then
            continue
        fi

        support=$(basename $re_rooutname .nwk)_support.nwk
        if [ -f $support ]
        then
            if [ "$force" = false ]
            then
                echo "Skipping $re_rooutname because its support counterpart already exists."
                continue
            fi
        fi

        nw_support -p $species_tree $d/$re_rooutname > $d/$support

    done
done
