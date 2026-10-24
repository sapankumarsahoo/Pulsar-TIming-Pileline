import sys, os
import time
import datetime
sys.path.insert(0,'/data/AllPulsars/Denis/astroutils/') # for met, mjd conversions and also coords (sepangle)
import numpy
from coords import ec2eq, ra_parse, dec_parse
import conv_dates_lib as conv

'''
24 Jan 2018, INAF-Milano.
  35 of Simon's 135 low Edot pulsars he sent you have ELONG instead of RAJ.
  So add ec2eq.

12 March 2021, LAB
  Finish porting to centos7 on laptop.
  
 4 May 2023, LAB
gtselect takes around 40 minutes. Why is it so much faster at SLAC? Anyway...
Some said maybe it'll be much faster to do gtselect on individual weekly files,
and then concatenate the little ones into a big one.
Give it a try.
 11 May 2023 -- did line command speed tests, and individual gtselects on individual weekly files takes the same time.
 ON THE OTHER HAND, directing gtselect's output to the same disk as the input files, /mnt/hgfs/echanges/, gets you a factor of two.


'''

if len(sys.argv) < 2:
    print("Usage: FoldSome.py <parfile> [True/refold/False]")
    print("       Reads a .par file, does gtselect tempo2 mv ln plot_phaso unless you say not to.")
    print("        <parfile> can have wildcards, to do several in a row.")
    print(" python /data/TEMPO/DavesStuff/FoldSome.py parfiles/*.par False     ")
    sys.exit(1)

mypars = []
for i in range(len(sys.argv)):
   if i>0 and i<(len(sys.argv)-1):
      mypars.append(sys.argv[i])

ExecCommands = True
ReFold = False
if sys.argv[-1] == 'False':
   ExecCommands = False
if sys.argv[-1] == 'refold':
   ExecCommands = False
   ReFold = True
print(sys.argv)
print("ExecCommands", ExecCommands)

tmin = 54682.9 # launch (mission start, actually) is 54682.7
tmax = 63595. # 
startMJD = tmin
finishMJD = tmax

gtpath = '/data/Fermitools/miniconda3/envs/fermi/bin/fermitools/'
ft1name='/mnt/hgfs/echanges/p305/FT1s.lst'  # 21 Dec 2018

# Syntax examples:
#gtselect infile=J1817-1742_gtsrcprob.fits outfile=J1813-1749.fits ra=273.396525 dec=-17.832633 emin=100. emax=999999. rad=2. zmax=180. tmin=INDEF tmax=INDEF
#tempo2 -gr fermi -ft1  J1813-1749.fits  -ft2 ../ft2-30s.fits -f J1813-1749_hgc12.par  -phase -bary -graph 0 -cacheft2
#plot_phaso J1813-1749.fits ra=273.396525 dec=-17.832633 title=J1813-1749 wbkgd=1 wcol=3FGLJ1814.0-1757c start=54918 finish=55969  wmin=0.01 tmin=54918 tmax=55969
#SearchPulsation psfpar0=5.445 psfpar1=0.848 psfpar2=0.084 fitsfile=a.fits ra=273.811133 dec=-17.634167 psrtitle=J1815-1738 -b tmin=292896000. tmax=422496000.


for parname in mypars:
   elong = -999.
   pf = open(parname,"r")
   for line in pf.readlines():
      if line.startswith('PSR'):
        psrname = line.split()[1]
      if line.startswith('START'):
        startMJD  = float(line.split()[1])
      if line.startswith('FINISH'):
        finishMJD = float(line.split()[1])

      if line.startswith('RAJ'):
        ra  = float(ra_parse (line.split()[1]))
      if line.startswith('DECJ'):
        dec  = float(dec_parse (line.split()[1]))

      if line.startswith('ELONG') or line.startswith('LAMBDA'):
        elong  = float(line.split()[1])
      if line.startswith('ELAT') or line.startswith('BETA'):
        elat  = float(line.split()[1])

   if elong > -888.: ra, dec =  ec2eq(elong, elat)
   print(' ')

# 15 Nov 2019, changed zmax from 100 to 105, seems to really improve J0653+04706
   thisft1 = '%s.fits ' %(psrname.split('.')[0])
   
   gtselcmd='%sgtselect infile=%s outfile=/mnt/hgfs/echanges/%s ra=%0.3f dec=%0.3f emin=100. emax=1.e6 rad=2. zmax=105. tmin=%d tmax=%d' %(gtpath, ft1name, thisft1, ra, dec, conv.mjd2met(tmin), conv.mjd2met(tmax))
   print(gtselcmd)
   bob = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
   calen_date  = datetime.datetime.fromtimestamp(time.mktime(time.strptime(bob,"%Y-%m-%d %H:%M:%S")))
   print ("Started at ", calen_date)
   mjdnow = conv.calen2mjd(calen_date)
   beginMET = conv.mjd2met(mjdnow)
   if ExecCommands: os.system(gtselcmd) 
   bob = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime())
   calen_date  = datetime.datetime.fromtimestamp(time.mktime(time.strptime(bob,"%Y-%m-%d %H:%M:%S")))
   mjdnow = conv.calen2mjd(calen_date)
   endMET = conv.mjd2met(mjdnow)
   print ("gtselect took %.2f minutes"%((endMET-beginMET)/60.))
   
   cpcmd = 'cp --remove-destination /mnt/hgfs/echanges/%s ./' %(thisft1)
   print(cpcmd)
   if ExecCommands or ReFold: 
      os.system(cpcmd) 
# It may bitch and moan that it can't copy because the link exists. 
# I THINK IT DOESN'T LIKE FT2 ON OTHER DISK EITHER. So... cd /mnt/hgfs/echanges for now...
      
   tempo2cmd = 'tempo2 -gr fermi -ft1 %s  -ft2 /mnt/hgfs/echanges/lat_spacecraft_merged.fits -f %s -phase -bary -graph 0 -cacheft2 -nofit' %(thisft1, parname)
   print(tempo2cmd)
   if ExecCommands or ReFold: os.system(tempo2cmd) 

   mvcmd = 'mv %s  /mnt/hgfs/echanges/' %(thisft1)
   print(mvcmd)
   if ExecCommands or ReFold: os.system(mvcmd)

   lncmd = 'ln -s  /mnt/hgfs/echanges/%s %s ' %(thisft1, thisft1)
   print(lncmd)
   if ExecCommands or ReFold: os.system(lncmd)

   plotcmd = 'plot_phaso /mnt/hgfs/echanges/%s ra=%0.3f dec=%0.3f title=%s wbkgd=1 zerosup=1 logesig=0.5 logeref=3.6 start=%d finish=%d grdev=%s.png/png' %(thisft1, ra, dec, psrname.split('.')[0], startMJD, finishMJD, psrname.split('.')[0])
   print(plotcmd)
   if ExecCommands or ReFold: os.system(plotcmd) 

   eogcmd = 'eog %s.png &' %(psrname.split('.')[0])
   print(eogcmd)
   if ExecCommands or ReFold: os.system(eogcmd) 

   addwtscmd = 'add_weights --par %s  %s.fits' %(parname, psrname.split('.')[0])
   print(addwtscmd)

