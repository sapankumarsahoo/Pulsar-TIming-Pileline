### Routine written by Lucas Guillemot
### Laboratoire de Physique et Chimie de l'Environnement et de l'Espace -- Universit\'e d'Orl\'eans / CNRS, F-45071 Orl\'eans Cedex 02, France
### Station de radioastronomie de Nan\c{c}ay, Observatoire de Paris, CNRS/INSU, F-18330 Nan\c{c}ay, France

import datetime
def aide():
	print('')
	print('Aide conv_dates.py')
	print('==================')
	print('')
	print('Se lance par : python conv_dates.py XXX=XXX')
	print('')
	print('- met=XXX.XXX')
	print('- mjd=XXX.XXX')
	print('- date=\"YYYY-MM-JJ hh:mm:ss\"')
	print('- annee=XXX.XXX')
	print('- doy=XXX-YYY')
	
	return

def mjd2met(mjd):
	MJDref_int  = 51910.
	MJDref_frac = 0.0007428703703703703  # Lucas uses this in plot_phaso... does not improve J2022+3842.
#	MJDref_frac = 0.  # das 20 June 2017, it DOES help Crab/HistoFromPar.py

	MJDint	= int(mjd)
	MJDfrac = mjd - MJDint
	
	met	= float(MJDint - MJDref_int) + (MJDfrac - MJDref_frac)
	met    *= 86400.
	
	return met

def met2mjd(met):
	MJDref_int  = 51910.
	MJDref_frac = 0.0007428703703703703
#	MJDref_frac = 0.

	MJDint	= int(met / 86400.)
	MJDfrac = met / 86400. - MJDint + MJDref_frac
	MJDint += MJDref_int
	
	if (MJDfrac >= 1.):
		MJDfrac -= 1.
		MJDint  += 1
	elif (MJDfrac < 0.):
		MJDfrac += 1.
		MJDint  -= 1
	
	mjd = float(MJDint) + MJDfrac
	
	return mjd

def RXTEmjd2met(mjd):
	MJDref_int  = 49353.
	MJDref_frac = 6.965740740000000E-04  # Lucas uses this in plot_phaso... does not improve J2022+3842.

	MJDint	= int(mjd)
	MJDfrac = mjd - MJDint
	
	met	= float(MJDint - MJDref_int) + (MJDfrac - MJDref_frac)
	met    *= 86400.
	
	return met

def RXTEmet2mjd(met):
	MJDref_int  = 49353.
	MJDref_frac = 6.965740740000000E-04 # from header of Gotthelf pca J2022+3842 file

	MJDint	= int(met / 86400.)
	MJDfrac = met / 86400. - MJDint + MJDref_frac
	MJDint += MJDref_int
	
	if (MJDfrac >= 1.):
		MJDfrac -= 1.
		MJDint  += 1
	elif (MJDfrac < 0.):
		MJDfrac += 1.
		MJDint  -= 1
	
	mjd = float(MJDint) + MJDfrac
	
	return mjd

def mjd2calen(mjd):
	ref_date     = datetime.datetime(1995, 10, 10, 0, 0, 0)
	ref_date_mjd = 50000.
	
	diffdate    = mjd - ref_date_mjd
	delta       = datetime.timedelta(days=diffdate)
	calen  	    = ref_date + delta
	
	return calen

def calen2mjd(calen):
	ref_date     = datetime.datetime(1995, 10, 10, 0, 0, 0)
	ref_date_mjd = 50000.
	
	delta	     = calen - ref_date
	mjd	     = ref_date_mjd + delta.days + delta.seconds / 86400. + delta.microseconds / (1000000. * 86400.)
	
	return mjd
	
def calen2year(calen):
	newyear_day = datetime.datetime(calen.year, 1, 1, 0, 0, 0)
	diffdate    = calen - newyear_day
	
	year        = float(calen.year)
	year       += (diffdate.days + (diffdate.seconds + diffdate.microseconds / 1000000.) / 86400.) / 365.
	
	return year

def year2calen(year):
	year_int    = int(year)
	
	diffdate    = (year - year_int) * 365.
	delta       = datetime.timedelta(days=diffdate)
	
	newyear_day = datetime.datetime(year_int, 1, 1, 0, 0, 0)
	calen	    = newyear_day + delta
	
	return calen

def mjd2doy(mjd):
#	mjd_ref     = [48987,49352,49717,50082,50448,50813,51178,51543,51909,52274,52639,53004,53370,53736,54101,54466,54832,55197,55562,55927,56293,56658,57023,57387,57753,58118,58483,58848,59214]
#	year_ref    = [1993, 1994, 1995, 1996, 1997,  1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021]

	mjd_ref     = [48987,49352,49717,50082,50448,50813,51178,51543,51909,52274,52639,53004,53370,53736,54101,54466,54832,55197,55562,55927,56293,56658,57023,57387,57753,58118,58483,58848, 59214, 59579, 59944, 60309, 60675, 61040, 61405, 61770]
	year_ref    = [1993, 1994, 1995, 1996, 1997,  1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028]
	
	mjd_int     = int(mjd)
	
	for i in range(len(mjd_ref) - 1):
		mjd_inf = mjd_ref[i]
		mjd_sup = mjd_ref[i+1]
		year    = year_ref[i]
		
		if (mjd_int >= mjd_inf) and (mjd_int < mjd_sup):
			doy = str(mjd_int - mjd_inf + 1)
			
			if len(doy) == 1:
				doy = '00' + doy
			elif len(doy) == 2:
				doy = '0'  + doy
		
			doy = str(year) + '-' + doy
			
			break
			
	return doy
	
def doy2mjd(doy):
	yy = int(doy.split('-')[0])
	dd = int(doy.split('-')[1])

#	mjd_ref     = [48987,49352,49717,50082,50448,50813,51178,51543,51909,52274,52639,53004,53370,53736,54101,54466,54832,55197,55562,55927,56293,56658,57023,57387,57753,58118,58483,58848,59214]
#	year_ref    = [1993, 1994, 1995, 1996, 1997,  1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021]
	
	mjd_ref     = [48987,49352,49717,50082,50448,50813,51178,51543,51909,52274,52639,53004,53370,53736,54101,54466,54832,55197,55562,55927,56293,56658,57023,57387,57753,58118,58483,58848, 59214, 59579, 59944, 60309, 60675, 61040, 61405, 61770, 62135, 62501, 62866, 63231, 63596]
	year_ref    = [1993,1994,1995,1996,1997,1998,1999,2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026, 2027, 2028, 2029, 2030, 2031, 2032, 2033]

	
	for i in range(len(year_ref)):
		year = year_ref[i]
		mjd  = mjd_ref[i]
		
		if (yy == year):
			mjd += dd - 1
			
			break
			
	return mjd
	

def afficher(met,mjd,calen,year,doy):
	print("%20s %f" % ("Date MET (s) :",met))
	print("%20s %f" % ("Date MJD :",mjd))
	print("%20s %s" % ("Date calendrier :",calen))
	print("%20s %f" % ("Annee :",year))
	print("%20s %s" % ("DOY :",doy))
	
	return

