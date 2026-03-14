# source /data/upanda/.bashrc

from priwo.sigproc.hdr import readhdr, writehdr		
import sys
from astropy.time import Time


'''
file_name = sys.argv[1]	# name of the output .fil file
gmrt_hdr_ = sys.argv[2]	# name of gmrt_hdr file to extractr information from
mjd_


with open(gmrt_hdr_, 'r') as f:
	dict_info ={}
	for line_ in f:
		if not ':' in line_:
			continue
		key, value = line_.strip().split(':')[0], ''.join.(line_.strip().split(':')[1:])
		dict_info[key.strip()] = value.strip()
'''
RA = sys.argv[1]
DEC = sys.argv[2]
time_stmp_header = sys.argv[3]
file_name = sys.argv[4]


with open(str(time_stmp_header), 'r') as f:
	
	for line_ in f:
		if 'IST Time' in line_:
			#hrs, minute, sec  = line_.strip().split(':')[1:]
			ist_time =':'.join(line_.strip().split(':')[1:]).strip()
		if 'Date'  in line_:
			d, m ,y  = line_.strip().split(':')[1:]
			date = y.strip() + '-' + m.strip() + '-' + d.strip()


MJD_UTC = Time(date +' '+ ist_time,  format='iso').mjd - (5.5/24)
fil_header = {}

fil_header['rawdatafile'] = file_name
fil_header['source_name'] = file_name.split('_')[0]
fil_header['machine_id'] = 14
fil_header['telescope_id'] = 7
fil_header['src_raj'] = float(RA) # float(dict_info['Coordinates'].split(',')[0])
fil_header['src_dej'] = float(DEC) # float(dict_info['Coordinates'].split(',')[1])
fil_header['az_start'] = 0.0
fil_header['za_start'] = 0.0
fil_header['data_type'] = 1
fil_header['fch1'] = 500.00 #float(dict_info['Frequency Ch.1'])
fil_header['foff'] = -0.195315 #float(dict_info['Channel width'])
fil_header['nchans'] = 1024 #int(dict_info['Num Channels'])
fil_header['nbeams'] = 1
fil_header['ibeam'] = 1
fil_header['nbits'] = 16
fil_header['tstart'] = MJD_UTC
fil_header['tsamp'] = 4.096e-05
fil_header['nifs'] = 1
fil_header['size'] = 380



modified_header = writehdr(fil_header, str(file_name))


