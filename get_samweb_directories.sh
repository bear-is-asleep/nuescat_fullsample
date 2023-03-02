#!/bin/bash

DEFNAME="$1"
INPUT_DIR="$2"

# create output file in input directory
output_file="$INPUT_DIR/output_new.txt"

# loop over each line in the input file
while read -r RUN SUBRUN EVENT; do

  # list files matching DEFNAME and RUN.SUBRUN
  output=$(samweb -e sbnd list-files "defname:$DEFNAME and run_number $RUN.$SUBRUN")

  # extract file name from output
  filename=$(echo "$output" | awk '{print $1}')

  # get file access url
  url=$(samweb -e sbnd get-file-access-url "$filename" --schema root)

  # extract desired substring
  substring=$(echo "$url" | awk -F 'reco2_caf/' '{print $2}')

  # prepend prefix to substring
  substring="/pnfs/sbn/data/sbn_nd/poms_production/official/MCP2022A/v09_37_02_04/prodoverlay_corsika_cosmics_proton_genie_rockbox_sce/reco2_caf/$substring"

  # append substring as a new column to input file
  echo "$RUN $SUBRUN $EVENT $substring" >> "$output_file"

done < "$INPUT_DIR/selected_events.txt"
