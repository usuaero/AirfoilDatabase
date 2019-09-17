import subprocess
import math
import re
from collections import OrderedDict, Iterable
import json
import os


def generate_machup_from_xfoil(airfoils, alpha_min, alpha_max, delta_alpha, reynolds_number, grid_size, niter, ncrit, xtrt, xtrb, flap=0.0, xf=1.0, yf=0.0):
	if not isinstance(airfoils, Iterable): airfoils = [airfoils]
	for airfoil in airfoils:
		drag_polar_filename = execute(airfoil, 0.0, alpha_min, delta_alpha, reynolds_number, grid_size, niter, ncrit, xtrt, xtrb, flap, xf, yf)
		results1 = read_results_from_DP_file(drag_polar_filename)
		results1.remove(results1[-1])

		drag_polar_filename = execute(airfoil, 0.0, alpha_max, delta_alpha, reynolds_number, grid_size, niter, ncrit, xtrt, xtrb, flap, xf, yf)
		results2 = read_results_from_DP_file(drag_polar_filename)

		airfoil_name = 'NACA{}_Re={:.2E}_N={}_xtrt={}_xtrb={}'.format(airfoil, reynolds_number, ncrit, xtrt, xtrb)
		nskip = int(round(1.0 / delta_alpha))
		write_machup_airfoil(airfoil_name, results1, results2)
		for i in results2:
			results1.append(i)
	return results1

def execute(airfoil, alpha_start, target_alpha, delta_alpha, reynolds_number, grid_size, n_iter, ncrit, xtrt, xtrb, flap, xf, yf):
	case_name = '{}_Re={:.1E}_N={:3.1f}_xtrt={:3.1f}_xtrb={:3.1f}_{:0<5.1f}deg'.format(
		airfoil, reynolds_number, ncrit, xtrt, xtrb, target_alpha)
	drag_polar_filename = 'XFOIL_{0}_DP.txt'.format(case_name)
	if os.path.exists(drag_polar_filename): os.remove(drag_polar_filename)
	
	# Start execution of XFOIL
	with subprocess.Popen(['C:\\Users\\Zach\\Documents\\xfoil\\xfoil.exe'],
						  stdin = subprocess.PIPE,
						  stdout = subprocess.PIPE
						 ) as ps:
		commands = []
		isfile = False
		try:
			float(airfoil)
		except(ValueError):
			isfile = True
		if((len(airfoil)==4) and not isfile):
			commands = ['NACA {}'.format(airfoil).encode('utf-8')]
		else:
			commands = ['LOAD {}'.format(airfoil).encode('utf-8'),
						b'new']
		# Set up commands
		commands += [b'GDES',
					b'FLAP',
					'{}'.format(xf).encode('utf-8'),
					'{}'.format(yf).encode('utf-8'),
					'{}'.format(flap).encode('utf-8'),
					# b'GSET',
					b'EXEC',
					b'',
					b'PPAR',
					'N {}'.format(grid_size).encode('utf-8'),
					b'T 1',
					b'',
					b'',
					b'PANE',
					b'PPAR',
					b'',
					b'OPER',
					b'HINC',
					'VISC {}'.format(reynolds_number).encode('utf-8'),
					'ITER {}'.format(n_iter).encode('utf-8'),
					b'VPAR',
					'N {}'.format(ncrit).encode('utf-8'),
					'XTR {} {}'.format(xtrt, xtrb).encode('utf-8'),
					b'',
					b'PACC',
					drag_polar_filename.encode('utf-8'),
					b''
					]

		# Iterate to the desired alpha using the specified delta
		step = int(math.copysign(1, target_alpha))
		nstep = int(round(abs((target_alpha - alpha_start) / delta_alpha))) + 1
		alphas = [i * delta_alpha + alpha_start for i in range(0, nstep, step)]
		for i in range(0, int(target_alpha / delta_alpha) + step, step):
			commands.append('ALFA {}'.format(i * delta_alpha).encode('utf-8'))

		# Finish commands
		commands += [b'PACC',
					 b'',
					 b'QUIT'
					]
		# print('commands:')
		# print(commands)

		# Execute the commands and record the output
		response = ps.communicate(b'\n'.join(commands))[0].decode('utf-8')

	# Write the output to a file
	response = response.replace('\r\n', '\n')
	response = response.replace('a =-', 'a = -')
	# print('response:')
	print(response)
	output_filename = 'XFOIL_{}_dump.txt'.format(case_name)
	with open(output_filename, 'w') as file:
		file.write(response)

	return drag_polar_filename


def write_machup_airfoil(airfoil_name, results1, results2):
	airfoil_filename = airfoil_name + '.txt'
	with open(airfoil_filename, 'w') as file:
		fmt = '{:>7}\t {:>7}\t {:>7}\t {:>7}\n'
		file.write(fmt.format('alpha', 'CL', 'CD', 'Cm'))
		for i in range(0, len(results1)):
			file.write(fmt.format(results1[i][0], results1[i][1], results1[i][2], results1[i][3]))
		for i in range(0, len(results2)):
			file.write(fmt.format(results2[i][0], results2[i][1], results2[i][2], results2[i][3]))

		#always write out the last data point
#        l = len(results2) - 1
#        file.write(fmt.format(results2[l][0], results2[l][1], results2[l][2], results2[l][3]))

	json_properties = OrderedDict()
	json_properties['type'] = "datafile"
	json_properties['filename'] = airfoil_filename
	json_properties['Comments'] = "Tabular reference data from XFOIL"

	json_airfoil = OrderedDict()
	json_airfoil['properties'] = json_properties

	airfoil_name = airfoil_name.replace('.', 'p')
	json_data = OrderedDict()
	json_data[airfoil_name] = json_airfoil

	with open(airfoil_name + '.json', 'w') as json_file:
		json.dump(json_data, json_file, indent = 4)


def read_results_from_DP_file(results_filename):
	# Parse the output and extract alpha, CL, CD, and Cm
	results = []

	with open(results_filename, 'r') as results_file:
		results_text = results_file.readlines()

	header_regex = re.compile('.*alpha.*CL.*CD.*CDp.*CM.*Chinge.*Top_Xtr.*Bot_Xtr.*')
	header_iter = (item for item in results_text if re.match(header_regex, item))
	header = next(header_iter, None)
	if header is not None:
		ind = results_text.index(header)
		for i in range(ind + 2, len(results_text)):
			data = results_text[i].split()
			alpha = float(data[0])
			cL = float(data[1])
			cD = float(data[2])
			cm = float(data[4])
			cH = float(data[5])
			results.append([alpha, cL, cD, cm, cH])

	return sorted(results, key = lambda result: result[0])


def read_results_from_dump_file(dump_filename):
	# Parse the output and extract alpha, CL, CD, and Cm
	results = []

	with open(dump_filename, 'r') as dump_file:
		dump = dump_file.readlines()

	oper_regex = re.compile('\.OPERva   c>.*')
	alpha_regex = re.compile('.* a = .* CL = .*')
	cm_regex = re.compile('.* Cm = .* CD = .*')
	iter = (item for item in dump if re.match(oper_regex, item))
	ind = 0
	item = next(iter, None)
	while item is not None:
		ind = dump.index(item, ind + 1)
		alpha_line = dump[ind - 6]
		cm_line = dump[ind - 5]
		if re.match(alpha_regex, alpha_line) and re.match(cm_regex, cm_line):
			alpha_line = alpha_line.split()
			alpha = float(alpha_line[2])
			cL = float(alpha_line[5])

			cm_line = cm_line.split()
			cm = float(cm_line[2])
			cD = float(cm_line[5])

			results.append([alpha, cL, cD, cm])

		item = next(iter, None)

	return results
