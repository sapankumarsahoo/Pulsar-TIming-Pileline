#!/bin/bash

# gedit GHRSS_fil_file ; chmod a+rwx GHRSS_fil_file

source /home/ssahoo/.bashrc   #bhima
location="/home/ssahoo/soft/gptool_ver4.4.4/gptool"   #bhima 

#source /data/bhaswati/sapan/.bashrc   #fs9
#location="/data/bhaswati/sapan/software/gptool_ver4.4.5.PSR/gptool"   #fs9


#Cration of gptool.in file for GHRSS search
# Define the filename
FILE="gptool.in"

# Write the content to the file
cat > "$FILE" <<EOF
#*#*#gptool input file v2.0#*#*#
-------------------------------------------------
#****Mode of observation****#
IA		: Beam mode
0		: Polarization mode (0-> intensity, 1-> stokes data)
1		: Sample size of data (in bytes, usually 2)
-------------------------------------------------
#****Observation Parameters****#
300		: Frequency band (lowest value in MHz)
200		: Bandwidth(in MHz)
-1		: Sideband flag (-1-> decreasing +1-> increasing)
4096		: Number of channels
0.08192		: Sampling Interval (in ms)
-------------------------------------------------
#****Pulsar Parameters****#
J1207-5050	: Pulsar name
4.84275732	: Pulsar period (in milliseconds)
50.640		: DM (in pc/cc)
-------------------------------------------------
#****Dedispersion & Folding parameters****#
1		: Number of bins in folded profile (-1 for native resolution)
0		: Phase offset for folding
12		: Number of coefficients for each polyco span (nCoeff)
60		: Validity of each span (nSpan in mins)
12		: Maximum hour angle
-------------------------------------------------
#****Display Parameters****#
0		: Polarization channel to display (0-3 or -1 for all four)
2		: Display window size (seconds, 0-> pulsar period)
0		: Update mode    	(0-> automatic, 1-> manual)
0		: Time delay between window updates (0-> no delay, 1-> emulate real time)
-------------------------------------------------
#****Spectral line RFI mitigation options****#
200		: Number of channels to flag at band beginning
100		: Number of channels to flag at band end
1		: Frequency flagging options (0-> no flagging, 1-> real-time calculation)
3		: Bandshape to use for frequency flagging (1-> normalized bandshape, 2-> mean-to-rms bandshape, 3-> Both)
3.0		: Threshold for frequency flagging (in units of RMS deviation)
-------------------------------------------------
#****Time domain impulsive RFI mitigation options****#
1		: Time flagging options 	(0-> no flagging, 1-> real-time calculation)
1		: Data normalization before filtering (0-> no, 1-> yes)
1		: Time flagging algorithm	(1-> histogram-based, 2-> MAD-based)
3.0		: Threshold for time flagging (in units of RMS deviation)
-------------------------------------------------
#****Other options****#
20		: Smoothing window size for bandshape normalization (in number of channels)
1		: Normalization procedure (1-> cumulative smooth bandshape, 2-> externally supplied bandshape.dat)
1		: Replace by median values (0-> Ignore flagged samples, 1-> Replace flagged samples by window median, 2-> Replace by smooth bandshape)
-------------------------------------------------
#****I/O options****#
0		: Write channel flag file (0-> no,1-> yes)
0		: Write time flag file (0-> no, 1-> yes)
1		: Write out filtered 2D raw data (0-> no, 1-> yes)
0		: Write out fullDM.raw	(0-> no, 1-> yes)
-------------------------------------------------
#****manual flagging options****#
0		: Number of bad channel blocks
List		: #in next line, example: [200,400],[1200,1400]
EOF

# Notify user
echo "File '$FILE' has been created successfully."



# Step 1: Convert raw files to .gpt files

# Specify the directory to iterate through
# cd /data/ssahoo/data/42_079_gwbh7_500_200_4K_81us.4jul2k22

# Outer loop to iterate over directories starting with "HR"
#folders=(
#    "HR_2134-35_ia_500_200_4096_4_1_8_04jul2022"
#    "HR_2351-35_ia_500_200_4096_4_1_8_04jul2022"
#)


#for file in "${folders[@]}"; do
for file in *raw; do
    base_folder=$(basename "$file" .raw)
    # Get first 10 and last 9 characters
    ra="${base_folder:3:2}:${base_folder:5:2}:00.00"
    dec="${base_folder:7:3}:00:00.00"
    mkdir "$base_folder"
    echo "$base_folder"
    cp gptool.in "$base_folder"/
    mv "$base_folder"* "$base_folder"/ 
    cd "$base_folder" || continue  # Move into the directory, skip if cd fails
#    echo "$base_folder"
    pwd
    ls
    
    # Inner loop to process each .raw file
    for raw_file in *.raw; do
        # Get the base name of the raw file
        base_name=$(basename "$raw_file" .raw)
        # Step 2: Create the header file
        timestamp_file="$raw_file.timestamp"
        echo "$raw_file.timestamp"
        head -n 1 "$raw_file.timestamp"
        
        # Extract date and time information from the first row of the timestamp file
        first_row=$(head -n 1 "$timestamp_file")
        year=$(echo "$first_row" | awk '{print $1}')
        #month=$(echo "$first_row" | awk '{printf "%02d", $2}') # Ensure two-digit month with leading zero
        #day=$(echo "$first_row" | awk '{printf "%02d", $3}')   # Ensure two-digit day with leading zero
	    month=$(echo "$first_row" | awk '{printf "%d", $2}')
        month=$(echo "$month" | sed 's/^-0*/-/; s/^0*//')  #remove first zeros
	    day=$(echo "$first_row" | awk '{printf "%d", $3}')
        day=$(echo "$day" | sed 's/^-0*/-/; s/^0*//')
        hour=$(echo "$first_row" | awk '{print $4}')
        hour=$(echo "$hour" | sed 's/^-0*/-/; s/^0*//')
        minute=$(echo "$first_row" | awk '{print $5}')
        minute=$(echo "$minute" | sed 's/^-0*/-/; s/^0*//')
        sec=$(echo "$first_row" | awk '{print $6}')
        sub=$(echo "$first_row" | awk '{print $7}')
        second=$(echo "$sec + $sub" | bc)
        
        # Subtract 5 hours and 30 minutes from the time
        new_hour=$(($hour - 5))
        new_hour=$(echo "$new_hour" | sed 's/^-0*/-/; s/^0*//')
        new_minute=$(($minute - 30))
        new_minute=$(echo "$new_minute" | sed 's/^-0*/-/; s/^0*//')
        
        # Initialize utc_day to day
        utc_day="$day"
        utc_day=$(echo "$utc_day" | sed 's/^-0*/-/; s/^0*//')
        
        # Handle negative minute values
        if [ $new_minute -lt 0 ] && [ $new_hour -gt 0 ]; then
            utc_hour=$(($new_hour - 1))
            utc_minute=$(($new_minute + 60))
        elif [ $new_hour -lt 0 ] && [ $new_minute -gt 0 ]; then
            utc_hour=$(($new_hour + 24))
            utc_minute=$new_minute
            utc_day=$(($day - 1))
        elif [ $new_hour -lt 0 ] && [ $new_minute -lt 0 ]; then
            utc_hour=$(($new_hour + 23))
            utc_minute=$(($new_minute + 60))
            utc_day=$(($day - 1))
        elif [ $new_hour -ge 0 ] && [ $new_minute -ge 0 ]; then
            utc_hour="$new_hour"
            utc_minute="$new_minute"
        elif [ $new_hour -eq 0 ] && [ $new_minute -lt 0 ]; then
            utc_hour=$(($new_hour + 23))
            utc_minute=$(($new_minute + 60))
            utc_day=$(($day - 1))
        fi
        
        # If utc_day is less than or equal to zero, adjust the day and month
        if [ $utc_day -le 0 ]; then
            # Adjust the month
            month=$(($month - 1))

            # Handle month wrap-around if necessary
            if [ $month -le 0 ]; then
                month=12
                year=$(($year - 1))
            fi

            # Determine the number of days in the previous month
            case $month in
                1|3|5|7|8|10|12) days_in_month=31 ;;
                4|6|9|11) days_in_month=30 ;;
                2)
                    # Check for leap year
                    if [ $((year % 4)) -eq 0 ] && ([ $((year % 100)) -ne 0 ] || [ $((year % 400)) -eq 0 ]); then
                        days_in_month=29
                    else
                        days_in_month=28
                    fi
                    ;;
            esac

            utc_day=$days_in_month
        fi
        
        # Convert date and time to UTC format with leading zeros for day and month
        # Convert UTC date and time to MJD
        utc_mjd=$(cal2mjd $year $month $utc_day | awk '{print $NF}')
        utc_time="$utc_hour:$utc_minute:$second"
        utc_date="$year/$month/$utc_day"
        
        # Write header file
        echo " update the template with the updated details for ----->  "$base_name".raw "
        echo "# DATA FILE HEADER #" > "$base_name.raw.hdr"
        echo "Site            : GMRT" >> "$base_name.raw.hdr"
        echo "Observer        : ssahoo" >> "$base_name.raw.hdr"
        echo "Proposal        : 47_113" >> "$base_name.raw.hdr"
        echo "Array Mode      : IA" >> "$base_name.raw.hdr"
        echo "Observing Mode  : GHRSS" >> "$base_name.raw.hdr"
        echo "Date            : "$utc_date"" >> "$base_name.raw.hdr"
        echo "Num Antennas    : 30" >> "$base_name.raw.hdr"
        echo "Antenna List    : C00 C01 C02 C03 C04 C05 C06 C08 C09 C10 C11 C13 C14" >> "$base_name.raw.hdr"
        echo "Num Channels    : 4096" >> "$base_name.raw.hdr"
        echo "Channel width   : -0.04882" >> "$base_name.raw.hdr"
        echo "Frequency Ch.1  : 500.000000" >> "$base_name.raw.hdr"
        echo "Sampling Time   : 81.92" >> "$base_name.raw.hdr"
        echo "Num bits/sample : 8" >> "$base_name.raw.hdr"
        echo "Data Format     : integer binary, little endian" >> "$base_name.raw.hdr"
        echo "Polarizations   : Total I" >> "$base_name.raw.hdr"
        echo "MJD             :"$utc_mjd"" >> "$base_name.raw.hdr"
        echo "UTC             : "$utc_time"" >> "$base_name.raw.hdr"
        echo "Source          : ${file: 0:10}" >> "$base_name.raw.hdr"
        echo "Coordinates     : $ra, $dec" >> "$base_name.raw.hdr"
        echo "Coordinate Sys  : J2000" >> "$base_name.raw.hdr"
        echo "Drift Rate      : 0.0, 0.0" >> "$base_name.raw.hdr"
        echo "Obs. Length     : 7200" >> "$base_name.raw.hdr"
        echo "Bad Channels    : 2:1-400,3700-4096" >> "$base_name.raw.hdr"
        echo "Bit shift value : 0" >> "$base_name.raw.hdr"
        
        cat "$base_name.raw.hdr"

        # Run gptool command to convert raw file to .gpt file
        #$location -f "$raw_file" -nodedisp -m 32 -t 4 -o .
        # Generating filterbank file
        #filterbank "$base_name.raw.gpt" > "$base_name.fil"
        #rm -r "$base_name.raw.gpt"
#        # rm -r "$base_name.raw"
        #prepfold -timing J1207-5050.par -n 128 -nsub 64 -npart 60 "$base_name.fil"

    done  # End of inner for loop
    cd ..  # Move back to the parent directory
done  # End of outer for loop
