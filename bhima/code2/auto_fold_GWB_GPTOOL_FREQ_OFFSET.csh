#! /bin/csh -f

# Set MATH alias
  alias MATH 'set \!:1 = `echo "\!:3-$" | bc `'
  alias MATH1 'set \!:1 = `echo "\!:3-$" | bc -l`'

if($#argv<2) then
    echo " "
    echo "Usage: autofold.csh pulsarname filename mode [ia/pa]"
    echo "Exiting..... "
    exit
  endif

# set the variables from the command line arguments :
set psrname = $argv[1]
set list = $argv[2]
set mode = $argv[3]

set code=/home/bhaswati/soft/scripts/code2
set numlist=`wc -l $list | awk '{print $1}'`
set w_dir=`pwd`
echo "$list $numlist"
cp $list $w_dir
cd $w_dir

@ i = 1
while ($i <= $numlist)
        set date=`cat $list | head -$i | tail -1| awk '{print $1}'`
        set rawdir=`cat $list | head -$i | tail -1| awk '{print $2}'`
        set freq=`cat $list | head -$i | tail -1| awk '{print $3}'`
        set dedisp_type=`cat $list | head -$i | tail -1| awk '{print $4}'`
	set nchan=`cat $list | head -$i | tail -1| awk '{print $5}'`
        set raw_id=`cat $list | head -$i | tail -1| awk '{print $6}'`
echo "$date $rawdir $freq $dedisp_type $nchan $raw_id"

mkdir $date
cd $date

echo "$psrname" >filelist
$code/timing_GWB_GPTOOL_FREQ_OFFSET.csh $date $freq $rawdir $dedisp_type $mode $nchan $raw_id
cd ../
@ i++

end
