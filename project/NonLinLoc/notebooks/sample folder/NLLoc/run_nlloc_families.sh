#!/bin/bash

# Path to your runfile
runfile="run/nlloc_svi_P.in"

# Directory to your pick files
pick_files_dir="/projects/amt/lpapin/SVI/NLLoc/obs/picks_family"
output_loc_dir="/projects/amt/lpapin/SVI/NLLoc/loc/SVI_Bostock"

# List of family numbers (customize if needed)
family_numbers=$(seq -f "%03g" 1 300)

# Backup the original runfile to restore later
cp $runfile "${runfile}.bak"

# Loop through each family number
for fam in $family_numbers; do
    echo "Processing family ${fam}..."

    # Create the folder for the current family if it doesn't exist
    mkdir -p "${output_loc_dir}/${fam}"

    # Update the LOCFILES line with the current family number and other paths
    sed -i "s|LOCFILES .*|LOCFILES ${pick_files_dir}/*${fam}* NLLOC_OBS /projects/amt/lpapin/SVI/NLLoc/time/layer ${output_loc_dir}/${fam}/loc #LP|g" "$runfile"
    
    # Run NonLinLoc with the updated runfile
    NLLoc "$runfile"
    
    echo "Finished processing family ${fam}"
done

# Restore the original runfile after the loop is done
mv "${runfile}.bak" "$runfile"

echo "All families processed, original runfile restored."
