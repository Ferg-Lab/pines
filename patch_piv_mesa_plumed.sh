#!/bin/bash

if [ "$#" -lt 1 ]; then
	echo "Please specify PLUMED path for patching it with MESA PIV codes."
	#exit 1
else

    # Check if PLUMED directory exists and extract absolute path to src
    if [ -d $1 ]; then
    
        plumed_absolute_location=$( realpath $1 )
        plumed_src_directory="${plumed_absolute_location}/src"

        ## PRINT TO SCREEN
        echo -e "\n"
        echo $plumed_absolute_location
        echo -e "PLUMED src: $plumed_src_directory \n"
        echo -e "Copying MESA PIV codes to $plumed_src_directory \n"

        # copy MESA PIV codes to src directory
        # remove existing ANN code and copy updated ANN, and add PIV codes
        rm -rf $plumed_src_directory/annfunc
        cp -rv src/annfunc $plumed_src_directory/.
        cp -rv src/pines $plumed_src_directory/.

        echo -e "Copied.\n"
    else
    
        echo "Given PLUMED directory not found."
    fi
fi
