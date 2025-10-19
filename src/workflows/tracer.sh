#!/bin/bash

# ==============================================================================
#? Multiparametric Calcium Trace Analyzer
# ==============================================================================
 
# --- General settings and variables -------------------------------------------
source ./src/bash_commons.sh

in_path="./data/in/"
out_path="./data/out/"

# --- The pipeline starts here -------------------------------------------------
echo -e "\n${mag}STARTING TRaceR${end}"

# Check if target directory exists
if [ ! -d "$out_path" ]; then
  mkdir -p "$out_path"
fi

echo -e "\n${grn}STEP 1: TRaceR${end}"
# To loop through files with spaces in their names or paths, change the default
# IFS `$' \n\t'` with the less eager `$'\n'` to properly parse `find` output.
OIFS="$IFS"
IFS=$'\n'
for sub_folder in $(find "${in_path}" -mindepth 1 -maxdepth 3 -type d | sort)
do
	# Mirror the input filesystem into the output directory
	out_sub_folder="${sub_folder/\/data\/in\//\/data\/out\/}"
	if [ ! -d "$out_path" ]; then
		mkdir -p "$out_path"
	fi

	# Run the TRaceR
	Rscript --vanilla "./src/tracer.R" \
		"$sub_folder" \
		"$out_sub_folder"
done

echo -e "\n${grn}STEP 2: Stat_Maker${end}"
for sub_folder in $(find "${out_path}" -mindepth 1 -maxdepth 1 -type d | sort)
do
	# Run the TRaceR
	Rscript --vanilla "./src/stat_maker.R" \
		"${sub_folder/\/data\/out\//\/data\/in\/}" \
		"$sub_folder"
done
IFS="$OIFS"

# --- The pipeline ends here ---------------------------------------------------
if [[ $? -eq 0 ]]; then
    echo -e "\n${mag}PIPELINE COMPLETED SUCCESSFULLY${end}"
else
    echo -e "\n${red}PIPELINE FAILED${end}"
fi
