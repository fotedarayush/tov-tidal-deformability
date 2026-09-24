#!/bin/bash

# ============================================================
# Batch runner for TOV + tidal deformability calculations
# ============================================================

# Folder containing the EoS tables
EOS_DIR="eos/EoS_dd2_npY_T0_beta_eq_eta_D-eta_V-Mg"

# Text file containing one EoS filename per line
EOS_LIST="eos/lists/Eta_D_Eta_V_EoS_list.txt"

# Output directory
OUTPUT_DIR="results/tov"

# Numerical sequence parameters
RHO_START="1.e14"
RHO_STEP="0.3e14"
NMODELS="100"

# Number of EoSs to run
MAX_EOS = 10

# Create output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# Clear previous failed-run list
> "$OUTPUT_DIR/failed_runs.txt"

count=0

while IFS= read -r eos_name; do

    # Skip empty lines
    [ -z "$eos_name" ] && continue

    count=$((count + 1))
    if ["$count" -gt "MAX_EOS"]; then
    	break
    fi
    eos_path="$EOS_DIR/$eos_name"

    echo "=============================================="
    echo "Running model $count"
    echo "EoS: $eos_name"
    echo "=============================================="

    # Check that the EoS actually exists
    if [ ! -f "$eos_path" ]; then
        echo "ERROR: EoS file not found:"
        echo "$eos_path"

        echo "$eos_name" >> "$OUTPUT_DIR/failed_runs.txt"
        continue
    fi

    # Rewrite infile for this EoS
    echo "$RHO_START   $RHO_STEP   $NMODELS" > infile
    echo "$eos_path" >> infile

    # Run the Fortran code
    ./logtov_seq_geom_tidal.out

    # Check whether output was produced
    if [ -f "logtov_seq_geom_tidal.dat" ]; then

        # Remove .dat from EoS name
        stem="${eos_name%.dat}"

        # Save result under a unique name
        mv logtov_seq_geom_tidal.dat \
           "$OUTPUT_DIR/${stem}_TOV.dat"

        echo "Saved:"
        echo "$OUTPUT_DIR/${stem}_TOV.dat"

    else

        echo "ERROR: No output produced for $eos_name"
        echo "$eos_name" >> "$OUTPUT_DIR/failed_runs.txt"

    fi

done < "$EOS_LIST"

echo
echo "=============================================="
echo "Batch run complete"
echo "Processed $count EoS files"
echo "Results stored in: $OUTPUT_DIR"
echo "=============================================="
