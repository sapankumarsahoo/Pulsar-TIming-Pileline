#! /bin/csh -f

# Set MATH alias
   alias MATH 'set \!:1 = `echo "\!:3-$" | bc `'
   alias MATH1 'set \!:1 = `echo "\!:3-$" | bc -l`'


set code=/data/Sapan_grad_proj/scripts/code2 
set output = `pwd`

#set output=/data/Sapan_grad_proj/analysis2/J0248+4230/

set filename=filelist
set numfile=`wc -w $filename | awk '{print $1}'`


if($#argv<5) then
    echo " "
    echo "Usage: timing.csh date [01nov2019/1nov2k19] frequency_(MHz) [322 or 607] full_data_dir [/data/data/31_036_10dec2k16.2048.322.GHRSS] dedisp_type [I/C] mode [ia/pa] number_of_channels[512/1024/2048/4096] raw_id[0-5]"
    echo "Exiting..... "
    exit
  endif

# set the variables from the command line arguments :
set date = $argv[1]
set freq = $argv[2]
set rawdir = $argv[3]
set dedisp_type = $argv[4]
set mode = $argv[5]
set nchan = $argv[6]
set raw_id = $argv[7]

echo "$freq $rawdir"

$code/gmrt_hdr.timing_GWB_FREQ_OFFSET.csh $date $freq $rawdir $dedisp_type $mode $nchan $raw_id

@ i = 1
while ($i <= $numfile)
	set scan=`cat $filename | head -$i | tail -1`
	cd $scan
	r#m -rf *par *pfd*
	cp /data/Sapan_grad_proj/analysis2/J0248-4230/${scan}*par .

	#Set the value of bandwidth and integration factor from the filname
        if ("$dedisp_type" == "I")then 
		cd ${rawdir}       
		if ("$raw_id" == 0)then
			set rawfile = "`find . -name ${scan}*raw`"
                        set bandwidth = `echo $rawfile | awk -F_ '{print $4}'`
                        set int_factor = `echo $rawfile | awk -F_ '{print $6}'`
                else
                        set rawfile = "`find . -name ${scan}*raw.${raw_id}`"
                        set bandwidth = `echo $rawfile | awk -F_ '{print $4}'`
                        set int_factor = `echo $rawfile | awk -F_ '{print $6}'`
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
	                set tres=`cat ${scan}_${mode}_${freq}_200_4096_${int_factor}_1_8_${date}.raw.gmrt_hdr | grep "Sampling Time" | awk -F: '{print $2}'`
	                MATH1 tres = $tres/1000000
	                set duration=`wc -l ${rawdir}/${scan}_${mode}_${freq}_200_4096_${int_factor}_1_8_${date}.raw.timestamp |awk '{print $1}'`
	                MATH duration = $duration * 0.251658240
	                set duration=`echo $duration | awk -F. '{print $1}'`
			
			#GPTool
	                #cp /data/analysis/header_files/gptool.${freq}MHz_${bandwidth}MHzBW_${nchan}channel.in gptool.in
        	        #/home/jroy/softwares/gptool_ver4.2.1/gptool -f ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw -o . -m 32 -nodedisp  -t 4
                	#ln -s ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gpt ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_dat

		        # filterbank covnverion
		        filterbank ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_dat > ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.fil

	                prepfold -timing *par ${scan}*fil -noxwin -o ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.new_freq
        	        echo "Folding is complete for $scan" >> Time.log
	
        	       # rm -rf *.fil
			#rm -rf ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gpt
	                #rm -rf *.gmrt_dat
        	        #ln -s ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_dat
			
			echo `pwd`
	                cp -p *pfd* $output/${date}/${scan}
			cd ..

		else
			set tres=`cat ${scan}_${mode}_${freq}_200_4096_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_hdr | grep "Sampling Time" | awk -F: '{print $2}'`
                        MATH1 tres = $tres/1000000
                        set duration=`wc -l ${rawdir}/${scan}_${mode}_${freq}_200_4096_${int_factor}_1_8_${date}.raw.timestamp.${raw_id} |awk '{print $1}'`
                        MATH duration = $duration * 0.251658240
                        set duration=`echo $duration | awk -F. '{print $1}'`
	
       	 		#GPTool
			#cp /data/analysis/header_files/gptool.${freq}MHz_${bandwidth}MHzBW_${nchan}channel.in gptool.in
			#/home/jroy/softwares/gptool_ver4.2.1/gptool -f ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id} -o . -m 32 -nodedisp  -t 4
	        	#ln -s ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gpt ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_dat
	
			# filterbank covnverion
	        	filterbank ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_dat > ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.fil
	
			prepfold -timing *par ${scan}*fil -noxwin -o ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.new_freq
	       		echo "Folding is complete for $scan" >> Time.log
        
			#rm -rf *.fil
			#rm -rf ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gpt
			#rm -rf *.gmrt_dat
			#ln -s ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id} ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_dat

			echo `pwd`
                        cp -p *pfd* $output/${date}/${scan}
                        cd ..
		endif
	
	 else

		#GPTool
		#cp /data/analysis/header_files/gptool.${freq}MHz_${bandwidth}MHzBW_${nchan}channel.in gptool.in
                #/home/jroy/softwares/gptool_ver4.2.1/gptool -f ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.dat -o . -m 32 -nodedisp  -t 4
                #ln -s ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.dat.gpt ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.gmrt_dat

		#Filterbank conversion
		filterbank ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.gmrt_dat > ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.fil

                prepfold -n 128 -nsub 128 -timing ${scan}*par ${scan}*fil -noxwin -o ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.new_freq

                echo "Folding is complete for $scan" >> Time.log

                #rm -rf *.fil
                #rm -rf ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.dat.gpt
                #rm -rf *.gmrt_dat
                #ln -s ln -s ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.dat ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.gmrt_dat
                echo `pwd`
                cp -p *pfd*  $output/${date}/${scan}
	endif

        cd .. 
        @ i++
end
