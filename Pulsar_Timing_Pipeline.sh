#!/bin/bash
echo -e "\e[38;5;208m                             
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------
||    _____   _    _   _    __________________ _        _________   _   __    __    _   _    _   ________ _   ||
||   |  __ \ | |  | | | |  /  ________________(_)      |___   ___| | | |   \/   |  | | | \  | | |  ______(_)  ||
||   | |__) || |  | | | | (  (_____     ___  _  ____       | |     | | | |\  /| |  | | |  \ | | | |    ___    ||
||   |  ___/ | |  | | | |  \_____  \   /   \) || ___(      | |     | | | | \/ | |  | | | \ \| | | |   |_  |   ||
||   | |     | |__| | | |  _ ____)  ) | ( )   || |         | |     | | | |    | |  | | | |\   | |  \____| |   ||
||   |_|      \____/  |_| (_)______/   \___/(_||_|         |_|     |_| |_|    |_|  |_| |_| \__|  \________|   ||
||                                                                                                            ||
----------------------------------------------------------------------------------------------------------------
----------------------------------------------------------------------------------------------------------------                            
\e[0m"
#xterm  &
source /home/jroy/.bashrc
source /home/jroy/.bashrc_anaconda
source /home/jroy/.bashrc_anaconda2
source /data/jroy/.bashrc
source /data/jroy/.bashrc_anaconda
source /data/jroy/.bashrc_anaconda2


while true; do
    # Display options
    echo "Choose an option:"
    echo -e " \e[38;5;208m   a) Option A: Default Procedure For One Epoch \n    b) Option B: Filterbank File Generation For Multiple Epochs \n    c) Option C : Step Wise Timing Procedure For One Epoch \n    Press Enter to exit \e[0m"


         
    # Read user input
    read -p "Enter your choice: " choice

    # Check user input
    case $choice in
#----------------------------Option A -----------------------------------------------------------
        a|A)
#For other Pulsar data analysis, take input of pulsar name.

 #            echo "Enter the date in the below format"
	echo "For PA data date format --> 01jan2023"
	echo "CD data date format --> 01Jan2K23"
#	read dt
	echo "Enter location of the raw and timestamp file"
#	read location
	echo "Enter the frequency(MHz)"
#	read freq
	echo "Enter number of Channels (512 / 2048 / 4096) "
#	read bw
	echo "I or C  --> I for PA (Phased array Mode) & C (Coherent De-dispersed Phase Array Mode)"
#	read mode
#             echo "$dt" "$location" "$freq" "$mode" "$bw" 0 > file_list

#             read -p "Enter the pulsar name :-  " pulsar_name
# Extracting values from file_details
result_location=$(grep 'Result location' file_details | cut -d '"' -f 2)
timing_file=$(grep 'Previouis Timing file Location' file_details | cut -d '"' -f 2)
par_file=$(grep 'Parameter file location' file_details | cut -d '"' -f 2)

# Output the details for further analysis
echo "$result_location"
echo "$timing_file"
echo "$par_file"

#             read -p "Enter the location, where you want to save the results:-" result_location
#result_location="/data/Sapan_grad_proj/pipeline_test/results"
#             read -p "Enter the file with location of previous timing file" timing_file
#timing_file="/data/Sapan_grad_proj/pipeline_test/results/backup_J0248+4230_GSB.tim"
#             read -p "Enter the location and parameter file details" par_file
#par_file="/data/Sapan_grad_proj/pipeline_test/results/J0248+4230.par"

            
             
if [ "$mode" == I ]; then #For automated ia or pa in the last word of the command
echo "Please check if pulsar name is correct or not."
echo "=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>= Starting the process =<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>="
echo "=== /data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh $pulsar_name file_list cd ==="
xterm -e /data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh "$pulsar_name" file_list pa

else  #For automated ia or pa in the last word of the command
echo "Please check if pulsar name is correct or not."
echo "=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>= Starting the process =<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>="
echo "=== /data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh $pulsar_name file_list cd ==="
#/data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh "$pulsar_name" file_list cd
fi

#move the result for ToAs analysis.
#            cd "$result_location"/"$dt"/J*/
#            ls -ltr

#            mv -r /data/Sapan_grad_proj/analysis2/"$dt"/J*/*.pfd /data/Sapan_grad_proj/analysis2/"$dt"/J*/*.polycos /data/Sapan_grad_proj/analysis2/"$dt"/J*/*.bestprof "$result_location"


#to get the SNR with pdmp command:-
cd "$result_location"
echo "$result_location"
for pfd in *pfd; do                # xy
SNR=$(pdmp "$pfd" | grep "Best S/N" | awk '{print $4}')
echo "Number of ToA is calculated on the basis of SNR (SNR / 10):   $SNR"
# Calculate the number of ToAs
n_toa_close=$(echo "scale=0; $SNR / 10" | bc)  # Perform floating-point division
# Remove decimal values from n_toa
n_toa_close=${n_toa_close%.*}  # Remove everything after the decimal point
echo "n_toa_close value: $n_toa_close"  # Print the modified value of n_toa_close

#Avoid for SNR value 11to19, with 1 ToA, timing plot won't run
if  [[ $n_toa_close -eq 1 ]]; then
    n_toa_close=$((n_toa_close + 1))
    echo "n_toa_close: $n_toa_close"
else
    echo "n_toa_close: $n_toa_close"
fi
# List of numbers
numbers=(1 2 3 4 5 6 10 12 15 30 60)
n_toa=0
# Find nearest smaller value
for num in "${numbers[@]}"; do
    if (( num <=  n_toa )); then
        n_toa=$num
    elif (( num <= n_toa_close && num > n_toa )) ; then 
    n_toa=$num
    fi
done
echo "The nearest number of ToA which is a multiple of 60:  $n_toa"
cd "$result_location"

for file in *pfd; do   #xx
source /home/jroy/.bashrc
get_TOAs.py -n 2 -t /data/Sapan_grad_proj/TOAs/all_files/band_3/4096/J0248+4230_pa_500_200_4096_4_1_8_27nov2017.raw.gptool_PSR_0248+4230.pfd.bestprof "$file" > "$result_location"/"$file"_band3_TOAs"$n_toa"_SNR_"$SNR".tim
more "$result_location"/"$file"_band3_TOAs"$n_toa"_SNR_"$SNR".tim # To show the results

#Updating the tim file with the cat. Before this copy, your existing timing study tim file to the same location (/data/Sapan_grad_proj/TOAs/tim_files/tim_files/) & rename it as "final_tim_file.tim"  
cat "$timing_file" "$result_location"/"$file"_band3_TOAs"$n_toa"_SNR_"$SNR".tim > "$result_location"/temporary_tim_file"$n_toa".tim #to generate temporfary tim file
more "$result_location"/temporary_tim_file"$n_toa".tim
#showing plots for old + new added timing residuals
echo  "tempo2 -f "$par_file" /data/Sapan_grad_proj/analysis2/timing_result/temporary_tim_file"$n_toa".tim -gr plk"
echo "showing plots for old and newly added timing residuals"
#cd /data/Sapan_grad_proj/TOAs/tim_files/
source /home/jroy/.bashrc
#tempo2 -f /data/Sapan_grad_proj/TOAs/tim_files/J0248+4230.par /data/Sapan_grad_proj/TOAs/tim_files/temporary_tim_file"$i".tim -gr plk
tempo2 -f "$par_file" "$result_location"/"$file"_band3_TOAs"$n_toa"_SNR_"$SNR".tim -gr plk
tempo2 -f "$par_file" "$result_location"/temporary_tim_file"$n_toa".tim -gr plk

done #xx
done   # xyz            
# xy
;;          
#----------------------------Option B -----------------------------------------------------------            
        b|B)
            # Execute command B
            echo "Executing command B"
            # Add your command B here
            # Function to read inputs and write to file
            read -p "Enter the location, where you want to save the results:-" result_location
            cd "$result_location"                  
#             rm -r file_list            
write_to_file() {
             echo "Enter the date in the below format (--> For PA data date format --> 01jan2023, --> CD data date format --> 01Jan2K23):-"
	read date
	echo "Enter location of the raw and timestamp file"
	read location
	echo "Enter the frequency(MHz)"
	read freq
	echo "Enter band-width"
	read bandwidth
	echo "I or C  --> I for PA (Phased array Mode) & C (Coherent De-dispersed Phase Array Mode)"
	read mode

    # Write details to file
    echo "$date $location $freq $bandwidth $resolution $mode" 0 >> file_list
}

# Main loop
while true; do
    write_to_file

    read -p "Press Enter to continue for adding more epochs or Enter 'stop' or 'STOP' to Proceed Further:- " choice

    if [[ "$choice" == "stop" || "$choice" == "STOP" ]]; then
        break
    fi
done

# After loop, execute command for Processing the file_list.
echo "Executing command Folding"
cd "$result_location"

if ["$mode" == I]; then #For automated ia or pa in the last word of the command
echo "Please check if pulsar name is correct or not."
echo "=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>= Starting the process =<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>="
echo "=== /data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh J0248+4230 file_list cd ==="
source /home/jroy/.bashrc
source /data/jroy/.bashrc
/data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh J0248+4230 file_list pa

else  #For automated ia or pa in the last word of the command
echo "Please check if pulsar name is correct or not."
echo "=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>= Starting the process =<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>="
echo "=== /data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh J0248+4230 file_list cd ==="
source /home/jroy/.bashrc
source /data/jroy/.bashrc
/data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh J0248+4230 file_list cd
fi

echo "Folding of all the epochs are done and results are saved in the given result location. Please enter to processes futher"
read xyz
#Giving loop command line for correcting any error during running the script
while true; do #this will continue untill user enter stop.
    # Prompt the user to enter a command line
    echo "Go to the location folder by 'cd' and then change or modify the timestamp folder as require"
    read -p "Enter a command line. Keep repeating it, untill you find your desire files. Then write (or 'stop' to exit) to exit from it: " error_script

    # Check if the user wants to stop
    if [ "$error_script" == "stop" ]; then
        echo "Exiting the loop. Moving forward to the next steps"
        break
    fi
    eval "$error_script"  #This will execute the above command line. It is already there.
done


#time to do prepfold
echo "Jumping to the folder --> /data/Sapan_grad_proj/analysis2/"$date"/J*/ "
cd /data/Sapan_grad_proj/analysis2/"$date"/J*/
ls -ltr

echo "Do you want change the parameters and do the presfold again ?"
read -p "Enter Y or N or stop: " choice
if [[ "$choice" =~ ^[Yy]$ ]]; then #3
source /home/jroy/.bashrc
source /data/jroy/.bashrc
 
while true; do #this will continue untill userd enter stop, for doing prefold at different profile bins.
    # Prompt the user to enter a command line
    echo "Enter your command in the format --> prepfold -timing ___.par ___.fil -n timechunks_nos(8/16/32/64/128/254/512)"
    read -p "Enter a command line. Keep repeating it, untill you find your desire profile bins and pulse profile. Then write (or 'stop' to exit) to exit to the next step: " prepfold_command

    # Check if the user wants to stop
    if [ "$prepfold_command" == "stop" ]; then
        echo "Exiting the loop. Goodbye!"
        break
    fi
    eval "$prepfold_command"  #This will execute the above command line untill you type 'stop'. It is already there.
done

echo "After producing multiples ps files, open same machine from another tab & only kept the data files (ps / pfd / fil / polycos) which will be used for timing , others files you can move into a new folder"

elif [[ "$choice" =~ ^[Nn]$ ]]; then #3
echo "You chose NO. Moving forward to the next steps."
else #3
echo "Invalid choice."
fi  #3

            
            
            
         
            ;;           
#----------------------------Option C -----------------------------------------------------------
        c|C)
echo "Enter the number to choose the corresponding process."
echo "1: To generate Filterbank file, bestprofile, pfd, polycos etc"
echo "2: To generate profile plot"
echo "3: To find the ToAs"
echo "4: To account for in-band sub-banding and correction for DM variation with in a bandwidth"
echo "5: To view the timing residual plot"
read -p "Enter 1 or 2 or 3 or 4:- " choice

if [[ "$choice" -eq 1 ]]; then  #1
cd /data/Sapan_grad_proj/analysis2/
source /home/jroy/.bashrc
source /data/jroy/.bashrc
	echo "Enter the date in the below format"
	echo "For PA data date format --> 01jan2023 , CD data date format --> 01Jan2K23"
	read dt
	echo "Enter location of the raw and timestamp file"
	read location
	echo "Enter the frequency(MHz)"
	read freq
	echo "Enter number of Channel "
	read bw
	echo "I or C  --> I for PA (Phased array Mode) & C (Coherent De-dispersed Phase Array Mode)"
	read mode
	
	cd /data/Sapan_grad_proj/analysis2/
	rm -r "$dt"
	echo "===  rm -r "$dt"  ==="
	echo "=== cd /data/Sapan_grad_proj/analysis2/ ==="
	echo "$dt" "$location" "$freq" "$mode" "$bw" 0 > file_list
	more file_list
if ["$mode" == I]; then #For automated ia or pa in the last word of the command
echo "Please check if pulsar name is correct or not."
echo "=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>= Starting the process =<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>="
echo "=== /data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh J0248+4230 file_list cd ==="
xterm -e /data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh J0248+4230 file_list pa

else  #For automated ia or pa in the last word of the command
echo "Please check if pulsar name is correct or not."
echo "=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>= Starting the process =<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>=<>="
echo "=== /data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh J0248+4230 file_list cd ==="
xterm -e /data/Sapan_grad_proj/scripts/code2/auto_fold_GWB_GPTOOL_FREQ_OFFSET.csh J0248+4230 file_list cd
fi

#Giving loop command line for correcting any error during running the script
echo "If it is showing error that so such file (timestamp file) or anything else, you repeate the step. Click Y to repeate the above steps or N for continue further."
read -p "Enter Y or N: " choice
if [[ "$choice" =~ ^[Yy]$ ]]; then #2
source /home/jroy/.bashrc
 
while true; do #this will continue untill user enter stop.
    # Prompt the user to enter a command line
    echo "Go to the location folder by 'cd' and then change or modify the timestamp folder as require"
    read -p "Enter a command line. Keep repeating it, untill you find your desire files. Then write (or 'stop' to exit) to exit from it: " error_script

    # Check if the user wants to stop
    if [ "$error_script" == "stop" ]; then
        echo "Exiting the loop. Moving forward to the next steps"
        break
    fi
    eval "$error_script"  #This will execute the above command line. It is already there.
done
elif [[ "$choice" =~ ^[Nn]$ ]]; then #2
echo "You chose NO. Moving forward to the next steps."
else #2
echo "Invalid choice. Please enter Y or N."
fi #2


#time to do prepfold
cd /data/Sapan_grad_proj/analysis2/"$dt"/J*/
ls -ltr

echo "Do you want change the parameters and do the presfold again ?"
read -p "Enter Y or N or stop: " choice
if [[ "$choice" =~ ^[Yy]$ ]]; then #3
source /home/jroy/.bashrc
source /data/jroy/.bashrc
 
while true; do #this will continue untill userd enter stop, for doing prefold at different profile bins.
    # Prompt the user to enter a command line
    echo "Enter your command in the format --> prepfold -timing ___.par ___.fil -n timechunks_nos(8/16/32/64/128/254/512)"
    read -p "Enter a command line. Keep repeating it, untill you find your desire profile bins and pulse profile. Then write (or 'stop' to exit) to exit to the next step: " prepfold_command

    # Check if the user wants to stop
    if [ "$prepfold_command" == "stop" ]; then
        echo "Exiting the loop. Goodbye!"
        break
    fi
    eval "$prepfold_command"  #This will execute the above command line untill you type 'stop'. It is already there.
done

echo "After producing multiples ps files, open same machine from another tab & only kept the data files (ps / pfd / fil / polycos) which will be used for timing , others files you can move into a new folder"

elif [[ "$choice" =~ ^[Nn]$ ]]; then #3
echo "You chose NO. Moving forward to the next steps."
else #3
echo "Invalid choice. Please enter Y or N."
fi  #3

#Running the script again
echo "Do you want to retun to the main menu ?"
read -p "Enter yes ot no " main_menu_3
if [[ "$main_menu_3" =~ ^[Yy]$ ]]; then
/data/Sapan_grad_proj/code/working.sh
elif [[ "$main_menu_3" =~ ^[Nn]$ ]]; then
echo "Continue..."
else
echo "Continue..."
fi

#STEP-2================To generate the profile plot==============
elif [[ "$choice" -eq 2 ]]; then #**

# Navigate to the directory
cd /data/Sapan_grad_proj/TOAs/all_files/band_3/4096

# Source the necessary bashrc files
source /home/jroy/.bashrc_anaconda2
source /home/jroy/.bashrc_anaconda

# Loop through each profile
for profile in *pfd; do

# Prompt user to enter the value of number of channels
#    read -p "Enter the value of number of channels: " ch_num 
echo "Continue Processing for --> $profile"
read -p "Enter whether it is PA (4098) data or CD (512) data (PA / CD): " PA_CD

if [ "$PA_CD" == "PA" ]; then
    ch_num="4096"
    echo "$ch_num"
elif [ "$PA_CD" == "CD" ]; then
    ch_num="512"
    echo "$ch_num"
else
    echo "Enter the correct option (PA or CD)"
fi
    # Perform frequency and time averaging with pam
    echo "Begining the processing for --> $profile"
    pam -a PRESTO -T -F "$profile" -e ar
    #phase bin info
    phase_bin=$(psrstat "$profile" | grep "Number of pulse phase bins" | awk '{print $7}')
    # Extract central frequency
    freq=$(psrstat "$profile" | grep "Centre frequency" | awk '{print $5}')
    # Extract bandwidth
    bw=$(psrstat "$profile" | grep "Bandwidth" | awk '{print $4}')
    # Extract DM
    dm=$(psrstat "$profile" | grep "Dispersion measure" | awk '{print $5}') 
    scat_time_scale=$(echo "scale=6; 10^(-6.46 + 0.154 * (l($dm) / l(10)) + 1.07 * (l($dm)*l($dm) / (l(10)*l(10)) ) - 3.86 * l($freq /1000)/l(10))" | bc -l)
    # Calculate frequency resolution in Hz
    freq_reso=$(echo "scale=6; $bw / $ch_num" | bc -l)
    #scateering time scale (ms) from the DM vs freq emperical ralation:
    inra_ch_dispersion=$(echo "scale=10; 4.148808 * 10^6 * $dm * ( (1/ ($freq)^2) - (1/ ($freq + $freq_reso)^2) )" | bc -l)
        
    echo "Phase bin number (nbin) is $phase_bin"
    echo "Central frequency (MHz) is: $freq"
    echo "Bandwidth (MHz) is: $bw"
    echo "Frequency resolution (MHz) is: $freq_reso "
    echo "DM (in pc/cm^3) is: $dm"
    echo "scattering time scale from empericla relation is (ms): $scat_time_scale" #Check the values again with the VS Code values.
    echo "Intrachannel dispersion: $inra_ch_dispersion ms"  
done #1

#To produce the profle plotting with ipython
for ar_file in *ar; do
echo "Processing $ar_file..."
source /home/jroy/.bashrc_anaconda2
source /home/jroy/.bashrc_anaconda
ipython <<EOF
cd /data/Sapan_grad_proj/TOAs/all_files/band_3/4096
# Loop through each profile and perform pam command
#for profile in *pfd; do
#    pam -a PRESTO -T -F "\$profile" -e ar
#done

import scipy
import astropy.io
from numpy import *
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
    
# Extract data
hdul = fits.open('$ar_file') 
column0 = hdul[3].data[0][10][0][0] 
column1 = [hdul[3].data[0][10][0][0][i] for i in range(len(column0))]   # Relative flux densities values
# Find the shift distance to bring the maximum value to the center, maxima of the profile at the center
max_index = column1.index(max(column1))
shift_distance = len(column1) // 2 - max_index

# Shift the values to the right
column1_shifted = column1[-shift_distance:] + column1[:-shift_distance]

value = np.min(column1_shifted)
#To avoid negative relative flux values and increase all by adding the abs of minimum value of the series.
abs_value = np.abs(value)
relative_flux_den = column1_shifted + abs_value



#=======================Mualti-Gaussian Fit==========================
#adding relative flux density information and phase bin into a temporary file name "sapan.txt"
print("Generating the process to multi-componenet Gaussian fit into the profile")
phase_bins = np.arange(len(column1))
# Write phase bins and intensity data to sapan.txt
with open("sapan.txt", "w") as f:
    for phase_bin, intensity in zip(phase_bins, relative_flux_den):
        f.write(f"{phase_bin} {intensity}\n")


#for multi-gaussian fitting :- 
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

file = np.loadtxt('sapan.txt')
x = file[:, 0]
y = file[:, 1]
peaks, _ = find_peaks(y, height=1)  # Adjust threshold as needed

# Define Gaussian function
def gaussian(x, *params):
    A = params[::3]
    mu = params[1::3]
    sigma = params[2::3]
    result = np.zeros_like(x)
    for i in range(len(A)):
        result += A[i] * np.exp(-(x - mu[i])**2 / (2 * sigma[i]**2))
    return result

def fwhm(x, y):
    half_max = max(y) / 2
    pos = np.where(y >= half_max)[0]
    fwhm_value = x[pos[-1]] - x[pos[0]]
    return fwhm_value

def plot_peaks(x, y, peaks, params=None):
    plt.plot(x, y, label='Profile', color='black')
    if params is not None:
        plt.plot(x, gaussian(x, *params), label='Multi-component Gaussian Fit', linestyle='dashed', color='r')
        for i, peak in enumerate(peaks):
            # Calculate FWHM for each peak
            peak_height = params[3*i]
            peak_center = params[3*i + 1]
            peak_sigma = params[3*i + 2]
            fwhm_value = fwhm(x, gaussian(x, *params))
            print(f"FWHM for Peak {i+1}: {fwhm_value}")
        print("Fitting Parameters:")
        print("Amplitudes:", params[::3])
        print("Centers:", params[1::3])
        print("Sigmas:", params[2::3])

    plt.plot(x, y, label='Profile', color='black')
    plt.plot(x[peaks], y[peaks], "x", label="Peaks", color='r')
    plt.legend()
    plt.title("Profile with Peaks")
    plt.xlabel("Profile Bin")
    plt.ylabel("Relative Flux Density")
    plt.show()

# Initial plot with peaks and Gaussian fit
plot_peaks(x, y, peaks)

# Fit Gaussian to remaining peaks
initial_guess = []
for peak in peaks:
    initial_guess.extend([y[peak], x[peak], 1.0])

# Adjust initial guess for better convergence
popt, _ = curve_fit(gaussian, x, y, p0=initial_guess, maxfev=1000000)
initial_guess = [initial_guess[i] if i < len(popt) else 1.0 for i in range(len(initial_guess))]

# Final plot with remaining peaks and Gaussian fit
plot_peaks(x, y, peaks, popt)

plt.title('Profile with Multi-component Gaussian Fit to Peaks')
plt.xlabel('Profile Bin')
plt.ylabel('Relative Flux Density')
plt.legend()
plt.show()





print("Profile plotting is processing for the file -->","$ar_file")
#print("Shift distance:", shift_distance)
#print("Min value:", abs_value)
print("Original Relative Flux Densities Values:", column1)
print("Shifted Relative Flux Densities Values:", column1_shifted)
print("Relative Flux Densities, after adding the abs(minimum value), to bring the profile in +ve Y-asix: (Shifted by bin numbers --> ",shift_distance,")" , relative_flux_den)
print("Max value of relative Flux density:", np.max(column1_shifted))
#plt.plot(relative_flux_den / np.max(relative_flux_den), label='Profile')

print("========================Profile Plotting done======================="
EOF
done #**

#Running the script again
echo "Do you want to retun to the main menu ?"
read -p "Enter yes ot no " main_menu_3
if [[ "$main_menu_3" =~ ^[Yy]$ ]]; then
/data/Sapan_grad_proj/code/working.sh
elif [[ "$main_menu_3" =~ ^[Nn]$ ]]; then
echo "Continue..."
else
echo "Continue..."
fi

#STEP-3================Calculating the TOAs & error in the TOAs==============

elif [[ "$choice" -eq 3 ]]; then #1
    echo "You chose option 3 for generating ToAs"
    echo "Begining the process for calculating TOAs and their uncertainities"

    read -p "Enter the datails of the template profile with location (e.g. /Location/of/Template_profile/template_profile.bestprof ):-  " template_profile
    #  /data/Sapan_grad_proj/TOAs/all_files/band_3/4096/J0248+4230_pa_500_200_4096_4_1_8_27nov2017.raw.gptool_PSR_0248+4230.pfd.bestprof
    read -p "Enter the location of pfd file, bestprof, polycos of the epoch to calculate ToAs:-   " location_pfd
    #cd /data/Sapan_grad_proj/TOAs/all_files/band_3/4096 
    read -p "Enter the location where you want to save the results:-   " result_location
    read -p "Enter the previous tim location:-  " tim_location
    #            /data/Sapan_grad_proj/pipeline_test/results

    cd "$location_pfd"
    source /home/jroy/.bashrc
    source /data/jroy/.bashrc
     
# Define the array of file locations
    for file in "$location_pfd"/*pfd #1
    do
        # Define the array of file locations
        file_locations=( "$file"
            # Add more file locations here if needed
        )

        # Iterate over each file location and extract the last filename
        for location in "${file_locations[@]}"; do #temp
            last_filename=$(basename "$location")
            echo "$last_filename"
            echo "Begining the process of ToA generation for -->  $file"
            evince "$file".ps &		

            # SNR from pdmp , taken as m
            m=$(pdmp "$file" | grep "Best S/N" | awk '{print $4}')
            echo "The SNR of "$file" is --> $m"

            for i in 1 2 3 4 5 6 10 12 15 30
            do #when ever you are doing, please insert your __.bestprof file name in the below name.
                get_TOAs.py -n "$i" -t "$template_profile" "$file"> "$result_location"/"$last_filename"_band3_TOAs"$i"_SNR"$m".tim
                echo "$i"
                more "$result_location"/"$last_filename"_band3_TOAs"$i"_SNR"$m".tim #to show the uncertainty in the calculated TOAs

                # Updating the tim file with the cat. Before this copy, your existing timing study tim file to the same location (/data/Sapan_grad_proj/TOAs/tim_files/tim_files/) & rename it as "final_tim_file.tim"  
                cat "$tim_location"/final_tim.tim "$result_location"/"$last_filename"_band3_TOAs"$i"_SNR"$m".tim > "$result_location"/temporary_tim_file"$i".tim #to generate temporary tim file
            done #2

            # this is to delete extra TOAs that we don't want	
            echo "Enter the value of no of TOAs wants to retain:"	
            read k	#to delete the undesire created TOAs
            { echo "#""$last_filename"_band3_TOAs"$k"_SNR"$m" ; cat "$result_location"/"$last_filename"_band3_TOAs"$k"_SNR"$m".tim ;} >> "$tim_location"/final_tim.tim  #adding the new TOAs into the old tim file i,e final_tim.tim
            rm -r "$result_location"/temporary_tim_file"$k".tim #deleting the temporary tim file
            source /home/jroy/.bashrc
            tempo2 -f "$result_location"/"J0248+4230.par" "$result_location"/"$last_filename"_band3_TOAs"$k"_SNR"$m".tim -gr plk
            tempo2 -f "$result_location"/"J0248+4230.par" "$tim_location"/final_tim.tim -gr plk

            for j in 1 2 3 4 5 6 10 12 15 30 ; do  #5
                if [ "$j" -ne "$k" ]; then
                    echo Deleting file: TOAs -n "$j" 
                    rm -r "$result_location"/"$last_filename"_band3_TOAs"$j"_SNR"$m".tim 
                    echo 'files deleted' "$last_filename"_band3_TOAs"$j"_SNR"$m".tim 'successfully'
                    rm -r "$result_location"/temporary_tim_file"$j".tim
                    echo 'files deleted' temporary_tim_file"$j".tim 'successfully'     
                else
                    echo "You have entered wronge number, please try again."    
                fi
            done #5
        done #temp		
    done #1

    # Running the script again
    echo "Do you want to return to the main menu ?"
    read -p "Enter yes or no " main_menu_3

    if [[ "$main_menu_3" =~ ^[Yy]$ ]]; then
        /data/Sapan_grad_proj/code/working.sh
    elif [[ "$main_menu_3" =~ ^[Nn]$ ]]; then
        echo "Continue..."
    else
        echo "Continue..."
    fi

#============================(-.-)==================================
#STEP-4: IN-Band Sub-banding and calculation for DM.
elif [[ "$choice" -eq 4 ]]; then  
    read -p "Enter the location where you want to save the result:-  " result_location
    #result_location="/data/Sapan_grad_proj/pipeline_test/results"
    read -p "Enter the datails of the template profile with location (e.g. /Location/of/Template_profile/template_profile.bestprof ):-  " template_profile
    #template_profile="/data/Sapan_grad_proj/TOAs/all_files/band_3/4096/J0248+4230_pa_500_200_4096_4_1_8_27nov2017.raw.gptool_PSR_0248+4230.pfd.bestprof"
    read -p "Enter the location of pfd file, bestprof, polycos of the epoch to calculate ToAs:-   " location_pfd
    #location_pfd="/data/Sapan_grad_proj/TOAs/all_files/band_3/4096"

    echo "FORMAT 1" > "$result_location"/your_file.tim
    echo "MODE 1" >> "$result_location"/your_file.tim
    echo "EFAC 1.05000" >> "$result_location"/your_file.tim

    for files in "$location_pfd"/*pfd; do 
        source /home/jroy/.bashrc
        for file in "$location_pfd"/*pfd; do
            # Define the array of file locations
            file_locations=( "$file"
                # Add more file locations here if needed
            )
            # Iterate over each file location and extract the last filename
            for location in "${file_locations[@]}"; do 
                last_filename=$(basename "$location")
                echo "$last_filename"
                
                echo "#" "$last_filename"
                for i in 1; do 
                    get_TOAs.py -s 8 -n "$i" -t "$template_profile" "$files" > "$result_location"/DM_Subband_fit.tim

                    { echo "#" "$last_filename"; cat "$result_location"/DM_Subband_fit.tim; } >> "$result_location"/DM_fit.tim

                    cat "$result_location"/your_file.tim "$result_location"/DM_fit.tim > "$result_location"/DM_fit_subband.tim
                    rm -r "$result_location"/DM_fit.tim

                    # adding 'r' to the third column for tim file format:
                    awk 'NR <= 4 {print $0; next} {print $0, "r"}' "$result_location"/DM_fit_subband.tim > "$result_location"/DM_final_fit_subband.tim
                    rm -r "$result_location"/DM_fit_subband.tim
                    more "$result_location"/DM_final_fit_subband.tim
                    cd "$result_location"

                    source /home/jroy/.bashrc
                    tempo2 -f "$result_location"/J0248+4230.par "$result_location"/DM_final_fit_subband.tim -gr plk

                    echo "Enter the DM value after doing fitting"
                    read j
                    awk 'NR <= 4 {print $0; next} {print $0, "-dm '$j'"}' "$result_location"/DM_final_fit_subband.tim > "$result_location"/DM_fiting_"$last_filename".tim
                    more "$result_location"/DM_fiting_"$last_filename".tim
                    echo "Please press 'Enter' to show the Timing Residual Plot with updated DM values:-"
                    read xyz
                    source /home/jroy/.bashrc
                    tempo2 -f "$result_location"/J0248+4230.par "$result_location"/DM_fiting_"$last_filename".tim -gr plk
                done
            done
        done
    done
else
    echo "Invalid choice. Please enter 1, 2, 3, or 4."
fi


#Running the script again
echo "Do you want to retun to the main menu ?"
read -p "Enter yes or no " main_menu_3
if [[ "$main_menu_3" =~ ^[Yy]$ ]]; then
/data/Sapan_grad_proj/code/working.sh
elif [[ "$main_menu_3" =~ ^[Nn]$ ]]; then
echo "Continue..."
else
echo "Continue..."
fi


#	rm -r /data/Sapan_grad_proj/TOAs/tim_files/"$file"_band3_TOAs"$j".tim
#replace band$
#replcae *500*pfd
#replace bestprof file name
#You can add i and j range 1 2 3 4 5 6 10 12 15 20 30 60 (should be a multiple of 60 as instructed by get_TOAs.py)









            ;;
        "")
            # Exit the script if user presses Enter
            echo "Exiting..."
            exit 0
            ;;
        *)
            # Handle invalid input
            echo "Invalid option. Please choose again."
            ;;
    esac
done
