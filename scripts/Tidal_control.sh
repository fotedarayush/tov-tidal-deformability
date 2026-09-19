#!/bin/bash
# declare STRING variable
#file="/home/sculkaputz/Documents/Tidal_control_bash/EoS_RMF_2023.txt"
#destdir="/home/sculkaputz/Documents/Tidal_control_bash/infile"
file="/home/dalvarez/Documents/Tidal_deformability_code/Eta_D_Eta_V_EoS_list.txt"
destdir="/home/dalvarez/Documents/Tidal_deformability_code/infile"
while IFS= read -r varname; do
    printf '%s\n' "$varname"
if [ -f "$destdir" ]
then
    echo "1.e14   0.3e14   100" > "$destdir"
    echo "$varname" >> "$destdir"
fi
./logtov_seq_geom_tidal.out
mv logtov_seq_geom_tidal.dat "TOV_$varname"
done < "$file"
#print variable on a screen
