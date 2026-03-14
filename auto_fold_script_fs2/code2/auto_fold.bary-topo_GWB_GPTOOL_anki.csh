#! /bin/csh -f

# Set MATH alias
  alias MATH 'set \!:1 = `echo "\!:3-$" | bc `'
  alias MATH1 'set \!:1 = `echo "\!:3-$" | bc -l`'

if($#argv<9) then
    echo " "
    echo "Usage: autofold.csh pulsarname filename DM [38.000] Period [4.84] bin[32] mode [ia/pa] zmax [0 or 200] numharm [8 or 16] dspsr [0 for both prepso-dspsr, 1 for only dspsr]"
    echo "Exiting..... "
    exit
  endif

# set the variables from the command line arguments :
 set psrname = $argv[1]
 set list = $argv[2]
 set dm = $argv[3]
 set period = $argv[4]
 set bin = $argv[5]
 set mode = $argv[6]
 set zmax = $argv[7]
 set numharm = $argv[8]
 set dspsr = $argv[9]
#

set code=/data/ankita/scripts/code2
set numlist=`wc -l $list | awk '{print $1}'`
set scratchdir=`pwd`
echo "$list $numlist"
cp $list $scratchdir
cd $scratchdir


@ i = 1
while ($i <= $numlist)
        set date=`cat $list | head -$i | tail -1| awk '{print $1}'`
        set rawdir=`cat $list | head -$i | tail -1| awk '{print $2}'`
        set freq=`cat $list | head -$i | tail -1| awk '{print $3}'`
        set dedisp_type=`cat $list | head -$i | tail -1| awk '{print $4}'`
	set nchan=`cat $list | head -$i | tail -1| awk '{print $5}'`
        set raw_id=`cat $list | head -$i | tail -1| awk '{print $6}'`

mkdir $date
cd $date

echo "$date $rawdir $freq $dedisp_type $nchan $raw_id"
echo "$psrname" >filelist
$code/timing_topo_bary_GWB_GPTOOL_anki.csh $date $freq $rawdir $dm $period $bin $dedisp_type $mode $scratchdir $zmax $numharm $dspsr $nchan $raw_id

cd ../
@ i++

end
