#! /bin/csh -f

# Set MATH alias
    alias MATH 'set \!:1 = `echo "\!:3-$" | bc `'
    alias MATH1 'set \!:1 = `echo "\!:3-$" | bc -l`'

set code=/data/ankita/scripts/code2
set output = `pwd`

#set output=/data/ankita/analysis/J2101-4802/

set filename=filelist
set numfile=`wc -w $filename | awk '{print $1}'`
set resdir=`pwd`
set date=`echo $resdir | awk -F/ '{print $4}'`

if($#argv<12) then
    echo " "
    echo "Usage: timing.csh frequency_(MHz) [322 or 607] full_data_dir [/data/data/31_036_10dec2k16.2048.322.GHRSS] DM [38.000] Period [4.84] bin[32] backend[GSB/GHRSS_2K/GHRSS_1K/GWB] mode [ia/pa] scratchdir [/data/ankita/analysis/] zmax [0 or 200] numharm [8 or 16] dspsr [0 for both prepso-dspsr, 1 for only dspsr] raw_id[0-5]"
    echo "Exiting..... "
    exit
  endif

# set the variables from the command line arguments :
set date = $argv[1]
set freq = $argv[2]
set rawdir = $argv[3]
set dm = $argv[4]
set period = $argv[5]
set bin = $argv[6]
set dedisp_type = $argv[7]
set mode =  $argv[8]
set scratchdir = $argv[9]
set zmax = $argv[10]
set numharm = $argv[11]
set dspsr = $argv[12]
set nchan = $argv[13]
set raw_id = $argv[14]

set scale = 4

set period_1 = `echo $period | awk -F. '{print $1}'`
set period_2 = `echo $period | awk -F. '{print $2}'`

echo $period_1 $period_2

$code/gmrt_hdr.timing_GWB.csh $date $freq $rawdir $dedisp_type $mode $nchan $raw_id

@ i = 1
while ($i <= $numfile)
	set scan=`cat $filename | head -$i | tail -1`
	cd $scan
	rm -rf ${date}_possible_candidates
	#rm -rf *fil *par *pfd* *.dat *.inf *.fft *ACCEL*

	#Setting the value of the bandwidth and intergration factor from the filename
        if ("$dedisp_type" == "I")then
                cd ${rawdir}
                if ("$raw_id" == 0)then
			 set rawfile = "`find . -name ${scan}*raw`"
                        set src_id_temp = `echo $rawfile | awk -F/ '{print $2}'`
                        set src_id = `echo $src_id_temp | awk -F_ '{print $1}'`
                        if ("$src_id" == "HR") then
                                set bandwidth = `echo $rawfile | awk -F_ '{print $5}'`
                                set int_factor = `echo $rawfile | awk -F_ '{print $7}'`
                        else
                                set bandwidth = `echo $rawfile | awk -F_ '{print $4}'`
                                set int_factor = `echo $rawfile | awk -F_ '{print $6}'`
                        endif
                else
                        set rawfile = "`find . -name ${scan}*raw.${raw_id}`"
                        set src_id_temp = `echo $rawfile | awk -F/ '{print $2}'`
                        set src_id = `echo $src_id_temp | awk -F_ '{print $1}'`
                        if ("$src_id" == "HR") then
                                set bandwidth = `echo $rawfile | awk -F_ '{print $5}'`
                                set int_factor = `echo $rawfile | awk -F_ '{print $7}'`
                        else
                                set bandwidth = `echo $rawfile | awk -F_ '{print $4}'`
                                set int_factor = `echo $rawfile | awk -F_ '{print $6}'`
                        endif



                endif
                cd -
        else
                cd ${rawdir}
                set rawfile = "`find . -name ${scan}*raw${raw_id}*dat`"
                set bandwidth = `echo $rawfile | awk -F_ '{print $3}'`
                set int_factor = `echo $rawfile | awk -F '[_.]' '{print $6}'`
                cd -
        endif
	
	if ("$dedisp_type" == "I")then
		if ("$raw_id" == 0)then
	        	set tres=`cat ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_hdr | grep "Sampling Time" | awk -F: '{print $2}'`
        		MATH1 tres = $tres/1000000
        		set duration=`wc -l ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp |awk '{print $1}'`
        		MATH duration = $duration * 0.251658240
        		set duration=`echo $duration | awk -F. '{print $1}'`

	       		echo "$nchan $tres $duration"
                        
                        #GPTool
                        cp /data/ankita/scripts/header_files/gptool.${freq}MHz_${bandwidth}MHzBW_${nchan}channel.in gptool.in
                        /data/jroy/softwares/gptool_ver4.2.1/gptool -f ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw -o . -m 32 -nodedisp  -t 4
                        ln -s ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gpt ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_dat

 
			#GPTool
			#cp ../../header_files/gptool.${freq}MHz_${bandwidth}MHzBW_${nchan}channel.in gptool.in
			#/data/jroy/softwares/gptool_ver4.2.1/gptool -f ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw -o . -m 32 -nodedisp  -t 4
			#ln -s ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gpt ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_dat
	
			#Filterbank Conversion
	                filterbank ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_dat >${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.fil  

			echo "Numout generation for $scan" >> Time.log
        		prepfold -dm 0.0 -p 0.1 -nosearch -noxwin ${scan}*.fil  -o ${scan}.numout
			set numout=`cat ${scan}.numout*bestprof | grep "Data Folded" | awk -F= '{print $2}'`
		        echo "Calulation of numout $numout is complete for $scan" >> Time.log 
	
			echo "Barycentric for $scan" >> Time.log 
			prepsubband -lodm $dm  -dmstep 0 -numdms 1 -numout $numout ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.fil -o ${scan}_search
			realfft -fwd ${scan}_search*.dat
			accelsearch -zmax $zmax -numharm $numharm ${scan}_search*.fft 
	
			#Candidate selection
        		set cand = `cat ${scan}_search*_ACCEL_${zmax}  | awk '{print $1 " " $6}' | grep " ${period_1}\.${period_2}" | awk '{print $1}'`
			echo "candidate number is $cand"

			set length = `echo "$cand" | wc -w`

			echo "$length candidate(s) found"

			if ( "$cand" == "" )then
			echo "No candidate matched"
	
			else 
				"In case multiple candidates are matching your search, these candidates along with their periods will be saved here. 21nov2k16 epoch is an example of such multiple matches." >> ${date}_possible_candidates

				cat ${scan}_search*_ACCEL_${zmax}  | awk '{print $1 " " $6}' | grep " ${period_1}\.${period_2}" | awk '{print $1 " " $2}' >> ${date}_possible_candidates
 
				set cand=`cat ${date}_possible_candidates | head -$i | tail -1 | awk '{print $1}'`

				echo "Candidate number is $cand"

				prepfold -accelcand $cand -accelfile ${scan}_search*_ACCEL_${zmax}.cand -dm $dm -nodmsearch -nosearch -noxwin ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.fil -o ${scan}_search_bary.${date}.${freq} 
				set topo_p=`cat ${scan}_search_bary*bestprof | grep P_topo | awk -F= '{print $2}' | awk -F+ '{print $1}'`
				set topo_pd=`cat ${scan}_search_bary*bestprof | grep "P'_topo" | awk -F= '{print $2}' | awk -F+ '{print $1}'`
				set topo_pdd=`cat ${scan}_search_bary*bestprof | grep "P''_topo" | awk -F= '{print $2}' | awk -F+ '{print $1}'`
				MATH1 topo_ps = $topo_p/1000.0
				echo "Period, Pd and Pdd are $topo_ps $topo_p $topo_pd $topo_pdd"
				prepfold -dm $dm -nodmsearch -n $bin -p $topo_ps -pd $topo_pd -pdd $topo_pdd -nosearch -topo -noxwin ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.fil -o ${scan}_search_topo.${date}.${freq} 

				echo "Folding is complete for $scan" >> Time.log  
        			echo `pwd`
				cp -p ${scan}_search_bary*pfd* ${output}/${scan}/
				cp -p ${scan}_search_topo*pfd* ${output}/${scan}/

			endif

			##rm -rf ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.fil
			#rm -rf ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gpt
			#rm -rf ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_dat
			#ln -s ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_dat

		else

			set tres=`cat ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_hdr | grep "Sampling Time" | awk -F: '{print $2}'`
                MATH1 tres = $tres/1000000
                set duration=`wc -l ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp.${raw_id} |awk '{print $1}'`
                MATH duration = $duration * 0.251658240
                set duration=`echo $duration | awk -F. '{print $1}'`

                echo "$nchan $tres $duration"
		
		#GPTool
                cp data/ankita/scripts/header_files/gptool.${freq}MHz_${bandwidth}MHzBW_${nchan}channel.in gptool.in
                /data/jroy/softwares/gptool_ver4.2.1/gptool -f ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id} -o . -m 32 -nodedisp  -t 4
                ln -s ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${int_factor}.gpt ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_dat

                #Filterbank Conversion
		filterbank ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_dat >${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.fil

                echo "Numout generation for $scan" >> Time.log
                prepfold -dm 0.0 -p 0.1 -nosearch -noxwin ${scan}*.fil -o ${scan}.numout
                set numout=`cat ${scan}.numout*bestprof | grep "Data Folded" | awk -F= '{print $2}'`
                echo "Calulation of numout $numout is complete for $scan" >> Time.log

                echo "Barycentric for $scan" >> Time.log
                prepsubband -lodm $dm  -dmstep 0 -numdms 1 -numout $numout ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.fil -o ${scan}_search
                realfft -fwd ${scan}_search*.dat
                accelsearch -zmax $zmax -numharm $numharm ${scan}_search*.fft

		#Candidate selection
		set cand = `cat ${scan}_search*_ACCEL_${zmax}  | awk '{print $1 " " $6}' | grep " ${period_1}\.${period_2}" | awk '{print $1}'`
                echo "candidate number is $cand"

                set length = `echo "$cand" | wc -w`

                echo "$length candidate(s) found"

                if ( "$cand" == "" ) then
                echo "No candidate matched"

                else
			"In case multiple candidates are matching your search, these candidates along with their periods will be saved here. 21nov2k16 epoch is an example of such multiple matches." >> ${date}_possible_candidates

			cat ${scan}_search*_ACCEL_${zmax}  | awk '{print $1 " " $6}' | grep " ${period_1}\.${period_2}" | awk '{print $1 " " $2}' >> ${date}_possible_candidates

                        set cand=`cat ${date}_possible_candidates | head -$i | tail -1 | awk '{print $1}'`

                        echo "Candidate number is $cand"

                        prepfold -accelcand $cand -accelfile ${scan}_search*_ACCEL_${zmax}.cand -dm $dm -nodmsearch -nosearch -noxwin ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.fil -o ${scan}_search_bary.${date}.${freq}
                        set topo_p=`cat ${scan}_search_bary*bestprof | grep P_topo | awk -F= '{print $2}' | awk -F+ '{print $1}'`
                        set topo_pd=`cat ${scan}_search_bary*bestprof | grep "P'_topo" | awk -F= '{print $2}' | awk -F+ '{print $1}'`
                        set topo_pdd=`cat ${scan}_search_bary*bestprof | grep "P''_topo" | awk -F= '{print $2}' | awk -F+ '{print $1}'`
                        MATH1 topo_ps = $topo_p/1000.0
                        echo "Period, Pd and Pdd are $topo_ps $topo_p $topo_pd $topo_pdd"
                        prepfold -dm $dm -nodmsearch -n $bin -p $topo_ps -pd $topo_pd -pdd $topo_pdd -nosearch -topo -noxwin ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.fil -o ${scan}_search_topo.${date}.${freq}

                        echo "Folding is complete for $scan" >> Time.log
                        echo `pwd`
                        cp -p ${scan}_search_bary*pfd* ${output}/${scan}/
                        cp -p ${scan}_search_topo*pfd* ${output}/${scan}/
		endif
		#rm -rf ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.fil
                #rm -rf ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gpt
                #rm -rf ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_dat
                #ln -s ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id} ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_dat
	endif
	else

		#GPTool
		#cp ../../header_files/gptool.${freq}MHz_${bandwidth}MHzBW_${nchan}channel.in gptool.in
                #/data/jroy/softwares/gptool_ver4.2.1/gptool -f ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.dat -o . -m 32 -nodedisp  -t 4
                #ln -s ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.dat.gpt ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.gmrt_dat
	
		#Fiterbank Conversion
		filterbank ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.gmrt_dat > ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.fil 

                echo "Numout generation for $scan" >> Time.log
                prepfold -dm 0.0 -p 0.1 -nosearch -noxwin ${scan}*.fil -o ${scan}.numout
                set numout=`cat ${scan}.numout*bestprof | grep "Data Folded" | awk -F= '{print $2}'`
                echo "Calulation of numout $numout is complete for $scan" >> Time.log

                echo "Barycentric for $scan" >> Time.log
                prepsubband -nsub 128 -lodm $dm  -dmstep 0 -numdms 1 -numout $numout ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.fil -o ${scan}_search
                realfft -fwd ${scan}_search*.dat
                accelsearch -zmax $zmax -numharm $numharm ${scan}_search*.fft
		
		#Candidate selection
		set cand = `cat ${scan}_search*_ACCEL_${zmax}  | awk '{print $1 " " $6}' | grep " ${period_1}\.${period_2}" | awk '{print $1}'`
                echo "candidate number is $cand"

                set length = `echo "$cand" | wc -w`

                echo "$length candidate(s) found"

                if ( "$cand" == "" ) then
                echo "No candidate matched"

                else
			"In case multiple candidates are matching your search, these candidates along with their periods will be saved here. 21nov2k16 epoch is an example of such multiple matches." >> ${date}_possible_candidates

			 cat ${scan}_search*_ACCEL_${zmax}  | awk '{print $1 " " $6}' | grep " ${period_1}\.${period_2}" | awk '{print $1 " " $2}' >> ${date}_possible_candidates

                        set cand=`cat ${date}_possible_candidates | head -$i | tail -1 | awk '{print $1}'`

                        echo "Candidate number is $cand"

                        prepfold -n 128 -nsub 128 -accelcand $cand -accelfile ${scan}_search*_ACCEL_${zmax}.cand -dm $dm -nodmsearch -nosearch -noxwin ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.fil -o ${scan}_search_bary.${date}.${freq}
                        set topo_p=`cat ${scan}_search_bary*bestprof | grep P_topo | awk -F= '{print $2}' | awk -F+ '{print $1}'`
                        set topo_pd=`cat ${scan}_search_bary*bestprof | grep "P'_topo" | awk -F= '{print $2}' | awk -F+ '{print $1}'`
                        set topo_pdd=`cat ${scan}_search_bary*bestprof | grep "P''_topo" | awk -F= '{print $2}' | awk -F+ '{print $1}'`
                        MATH1 topo_ps = $topo_p/1000.0
                        echo "Period, Pd and Pdd are $topo_ps $topo_p $topo_pd $topo_pdd"
                        prepfold -dm $dm -nodmsearch -n $bin -nsub 128 -p $topo_ps -pd $topo_pd -pdd $topo_pdd -nosearch -topo -noxwin ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.fil -o ${scan}_search_topo.${date}.${freq}

                        echo "Folding is complete for $scan" >> Time.log
                        echo `pwd`
                        cp -p ${scan}_search_bary*pfd* ${output}/${scan}/
                        cp -p ${scan}_search_topo*pfd* ${output}/${scan}/
		endif

		rm -rf ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.fil
                #rm -rf ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.dat.gpt
                #rm -rf ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.gmrt_dat
                #ln -s ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.dat ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.gmrt_dat
	endif
	cd ..
        @ i++ 
end
