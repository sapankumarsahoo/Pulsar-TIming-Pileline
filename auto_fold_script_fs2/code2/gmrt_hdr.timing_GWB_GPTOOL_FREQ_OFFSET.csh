#! /bin/csh -f

# Set MATH alias
  alias MATH 'set \!:1 = `echo "\!:3-$" | bc `'
  alias MATH1 'set \!:1 = `echo "scale=2; \!:3-$" | bc -l`'
  alias MATH2 'set \!:1 = `echo "scale=9; \!:3-$" | bc -l`'

set srcpath = "/data/ssahoo/scripts/code2/sources"
set code = "/data/ssahoo/scripts/code2/"

set filename=filelist
set numfile=`wc -w $filename | awk '{print $1}'`

if($#argv<5) then
    echo " "
    echo "Usage: gmrt_hdr.timing.csh date [01nov2019/1nov2k19] frequency_(MHz) [322 or 607] full_data_dir [/mnt/node104a/data/27_041_PA_5feb2k15] dedisp_type [I/C] mode [ia/pa] number_of_channels[512/1024/2048/4096] raw_id[0-5]" 
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

@ i = 1
while ($i <= $numfile)
	set out = "param"           # this file is used by template_GWB.c, contains all the parameters needed by gmrt_hdr
	set scan=`cat $filename | head -$i | tail -1`

	#Set the value of bandwidth and integration factor from the filname
        if ("$dedisp_type" == "I")then
		cd ${rawdir}
		if ("$raw_id" == 0)then
                        set rawfile = "`find . -name ${scan}*raw`"
                        set bandwidth = `echo $rawfile | awk -F_ '{print $4}'`
                        echo "Bandwidth is $bandwidth"
                        set int_factor = `echo $rawfile | awk -F_ '{print $6}'`
                else
                        set rawfile = "`find . -name ${scan}*raw.${raw_id}`"
                        set bandwidth = `echo $rawfile | awk -F_ '{print $4}'`
                        echo "Bandwidth is $bandwidth"
                        set int_factor = `echo $rawfile | awk -F_ '{print $6}'`
                endif
                cd -
        else
		cd ${rawdir}
                set rawfile = "`find . -name ${scan}*raw${raw_id}*dat`"
                set bandwidth = `echo $rawfile | awk -F_ '{print $3}'`
                echo "Bandwidth is $bandwidth"
                set int_factor = `echo $rawfile | awk -F '[_.]' '{print $6}'`
                cd -
        endif

	 if ("$dedisp_type" == "I")then
		if ("$raw_id" == 0)then
			echo "$scan of $numfile with ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw is getting used for gmrt_hdr generation....  "
		else
			echo "$scan of $numfile with ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id} is getting used for gmrt_hdr generation....  "
		endif
	else
		echo "$scan of $numfile with ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.dat is getting used for gmrt_hdr generation....  "
	endif

	mkdir ${scan}
	cp ${code}/template.${freq}.${bandwidth}.${nchan}.GWB.gmrt_hdr template.gmrt_hdr
	cp template.gmrt_hdr $scan 
	cd $scan

	 if ("$dedisp_type" == "I")then
		# obtaining time corordinates
		if ("$raw_id" == 0)then
			rm -rf ${scan}*raw.gmrt_dat ${scan}*raw.gmrt_hdr
	        	set tstmp=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp | head -1 | awk '{print $14}'`
			set year=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp | head -1 | awk '{print $8}'`
	        	set month=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp | head -1 | awk '{print $9}'`
	        	set day=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp | head -1 | awk '{print $10}'`
	        	set hour=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp | head -1 | awk '{print $11}'`
	        	set min=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp | head -1 | awk '{print $12}'`
	        	set sec=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp | head -1 | awk '{print $13}'`

		else
			rm -rf ${scan}*raw.${raw_id}.gmrt_dat ${scan}*raw.${raw_id}}.gmrt_hdr
			set tstmp=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp.${raw_id} | head -1 | awk '{print $14}'`
		        set year=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp.${raw_id} | head -1 | awk '{print $8}'`
                        set month=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp.${raw_id} | head -1 | awk '{print $9}'`
                        set day=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp.${raw_id} | head -1 | awk '{print $10}'`
                        set hour=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp.${raw_id} | head -1 | awk '{print $11}'`
                        set min=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp.${raw_id} | head -1 | awk '{print $12}'`
                        set sec=`cat ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.timestamp.${raw_id} | head -1 | awk '{print $13}'`
		endif
                
		set odate=`echo "$day":"$month":"$year"`
                set ist=`echo "$hour":"$min":"$sec"`  

		# Calculating the time resolution
		set TIME_RES
                MATH1 TIME_RES = (($int_factor / ( 2 * $bandwidth)) * (2 * $nchan))
                echo "Sampling Time = " $TIME_RES

		#Value of Frequency in *gmrt_hdr
		set FREQ_hdr
		MATH2 FREQ_hdr = $freq.000000

		${code}/conv_time $tstmp $odate $ist 0
		set MJD = `grep "MJD" cal2mjd.out |awk '{print $3}'`
		set dayflg = `cat day_flg`
        	if ($dayflg == 1)then
		MATH MJD = $MJD -1 
		endif
  		set Date = `grep "Date" time_date |awk '{print $2}'`
  		set UTC = `grep "UTC" time_date |awk '{print $2}'`
		echo "MJD for $scan = " $MJD
  		echo "Date for $scan = " $Date
  		echo "UTC for $scan = " $UTC

	# obtaining source coordinates
        	set r1 = `cat $srcpath/source* | grep $scan | awk '{print $2}' | cut -c1-2`
        	set r2 = `cat $srcpath/source* | grep $scan | awk '{print $2}' | cut -c4-5`
        	set r3 = `cat $srcpath/source* | grep $scan | awk '{print $2}' | cut -c7-11`
        	set RAJ = `echo "$r1":"$r2":"$r3"`
       		echo "RAJ = "$RAJ
        	set d1 = `cat $srcpath/source* | grep $scan | awk '{print $3}' | cut -c1-3` 
        	set d2 = `cat $srcpath/source* | grep $scan | awk '{print $3}' | cut -c5-6` 
        	set d3 = `cat $srcpath/source* | grep $scan | awk '{print $3}' | cut -c8-11` 
        	set DECJ = `echo "$d1":"$d2":"$d3"`
        	echo "DECJ = "$DECJ

	# putting all the parameters in the output file  param....
  		echo $Date >$out
		echo $MJD >>$out
		echo $UTC >>$out
		echo $scan >>$out
		echo $RAJ >>$out
        	echo $DECJ >>$out
		echo $TIME_RES >>$out
		echo $FREQ_hdr >> $out
  		echo "Running ${code}/template_GWB $Date $MJD $UTC $scan $RAJ $DECJ $TIME_RES $FREQ_hdr"
  		${code}/template_GWB
		if ("$raw_id" == 0)then
			cp gmrt_hdr ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_hdr
			#ln -s ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.gmrt_dat
		else
			cp gmrt_hdr ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_hdr
                        #ln -s ${rawdir}/${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id} ${scan}_${mode}_${freq}_${bandwidth}_${nchan}_${int_factor}_1_8_${date}.raw.${raw_id}.gmrt_dat
		endif
	
	else
		#Obtaining time coordinates
		rm -rf ${scan}*raw${raw_id}.gmrt_dat ${scan}*raw${raw_id}.gmrt_hdr
		set tstmp=0.`cat ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.timestamp | head -2 | awk -F '[:.]' '{print $5}'`
                set year=`cat ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.timestamp | head -3 | tail -1 | awk -F: '{print $4}'`
                set month=`cat ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.timestamp | head -3 | tail -1 | awk -F: '{print $3}'`
                set day=`cat ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.timestamp | head -3 | tail -1 | awk -F: '{print $2}'`
                set hour=`cat ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.timestamp | head -2 | awk -F '[:.]' '{print $2}'`
                set min=`cat ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.timestamp | head -2 | awk -F '[:.]' '{print $3}'`
                set sec=`cat ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.timestamp | head -2 | awk -F '[:.]' '{print $4}'`
		
                set odate=`echo "$day":"$month":"$year"`
                set ist=`echo "$hour":"$min":"$sec"`

		# Calculate the time resolution
		set TIME_RES
		MATH1 TIME_RES = (((2 *$int_factor) / ( 2 * $bandwidth)) * (2 * $nchan))
		echo "Sampling Time = " $TIME_RES

		#Value of Frequency in *gmrt_hdr
		set FREQ_hdr
		MATH2 FREQ_hdr = ($freq - ($bandwidth / ($nchan * 2)))
		
                ${code}/conv_time $tstmp $odate $ist 0
                set MJD = `grep "MJD" cal2mjd.out |awk '{print $3}'`
                set dayflg = `cat day_flg`
                if ($dayflg == 1)then
                MATH MJD = $MJD -1
                endif
                set Date = `grep "Date" time_date |awk '{print $2}'`
		set UTC = `grep "UTC" time_date |awk '{print $2}'`
                echo "MJD for $scan = " $MJD
                echo "Date for $scan = " $Date
                echo "UTC for $scan = " $UTC
		
	#Obtaining source coordinates
		set r1 = `cat $srcpath/source* | grep $scan | awk '{print $2}' | cut -c1-2`
                set r2 = `cat $srcpath/source* | grep $scan | awk '{print $2}' | cut -c4-5`
                set r3 = `cat $srcpath/source* | grep $scan | awk '{print $2}' | cut -c7-11`
                set RAJ = `echo "$r1":"$r2":"$r3"`
                echo "RAJ = "$RAJ
                set d1 = `cat $srcpath/source* | grep $scan | awk '{print $3}' | cut -c1-3`
                set d2 = `cat $srcpath/source* | grep $scan | awk '{print $3}' | cut -c5-6`
                set d3 = `cat $srcpath/source* | grep $scan | awk '{print $3}' | cut -c8-11`
                set DECJ = `echo "$d1":"$d2":"$d3"`
                echo "DECJ = "$DECJ

	#puting all the parameters in the outout file param.... 
		echo $Date >$out
                echo $MJD >>$out
                echo $UTC >>$out
                echo $scan >>$out
                echo $RAJ >>$out
                echo $DECJ >>$out
		echo $TIME_RES >>$out
		echo $FREQ_hdr >>$out
                echo "Running ${code}/template_GWB $Date $MJD $UTC $scan $RAJ $DECJ $TIME_RES $FREQ_hdr"
                ${code}/template_GWB
                cp gmrt_hdr ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.gmrt_hdr
		#ln -s ${rawdir}/${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.dat ${scan}_${freq}_${bandwidth}_${nchan}_${int_factor}.${date}.raw${raw_id}.gmrt_dat
	endif

# cleaning up...
  	rm -f param
  	rm -f gmrt_hdr
  	rm -f time_date
  	rm -f cal2mjd.out
  	rm -f template.gmrt_hdr

        echo `pwd`
        cd .. 
	@ i++
end
