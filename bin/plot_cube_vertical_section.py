#!/usr/bin/env python

import os
import argparse
import glob
import math
import numpy as np
import matplotlib.pyplot as plt
from elastic_stresses_py import PyCoulomb

DEBUG=True

def map_vert_profile(stress_cube, horiz_profile, depths, outfile):
	"""Display a vertical profile of stresses. Default colors for now."""

	print('horiz_profile.shape', horiz_profile.shape)
	lon2d = np.reshape(horiz_profile.lon1d, horiz_profile.shape);
	print('lon2d.shape', lon2d.shape)
	lat2d = np.reshape(horiz_profile.lat1d, horiz_profile.shape);
	print('lat2d.shape', lat2d.shape)
	
	max_nlon, max_nlat, max_ndepth = stress_cube.shape
	print(max_nlon, max_nlat, max_ndepth)
	
	# test vertical E-W slice near lat/2 (midLat)
	ilat = int(max_nlat / 2 - 1)
	#ilat = int(max_nlat / 2) # section at exactly lat/2 may give singularity
	#ilat += 6 # test vertical E-W slice at lat/2 + 3km (sourcePatch out/cube/20231208A/EastDipping/)
	coulomb_stress = stress_cube[:, ilat, :]
	coulomb_stress = coulomb_stress.transpose()
	coulomb_stress /= 1000.0   # convert to MPa
	
	lon_values = lon2d[ilat, :]
	X = np.tile(lon_values, (max_ndepth, 1))
	Y = np.tile(depths, (max_nlon, 1)).transpose()
	# convert to km
	X = X - horiz_profile.centerlon
	lon2km = 111.320 * math.cos(horiz_profile.centerlat * math.pi / 180)
	X = X * lon2km

	print('X.shape', X.shape)
	print('X', X)
	print('Y.shape', Y.shape)
	print('Y', Y)
	print('coulomb_stress.shape', coulomb_stress.shape)
	print('coulomb_stress', coulomb_stress)

	# Figure of stresses.
	plt.rcParams['figure.dpi'] = 300;
	plt.rcParams['savefig.dpi'] = 300;
	fig = plt.figure(figsize=(17, 8), dpi=300);
	#dislay_map = plt.contourf(X, Y, coulomb_stress, cmap='bwr');
	median_stress = np.mean(coulomb_stress)
	if median_stress >= 0.0:
		clevels = np.arange(-10.0*median_stress, 11.0*median_stress, median_stress)
	else:
		clevels = np.arange(10.0*median_stress, -11.0*median_stress, -median_stress)
	print("plot_cube_vertical_section.py: median_stress, clevels", median_stress, clevels)
	dislay_map = plt.contourf(X, Y, coulomb_stress, levels=clevels, extend='both', cmap='coolwarm');
	plt.title('Coulomb stresses on vertical profile, fixed reciever strike/dip/rake/depth of '+str(horiz_profile.strike)+', ' +
			  str(horiz_profile.dip)+', '+str(horiz_profile.rake)+', '+str(horiz_profile.depth_km));

	cb = fig.colorbar(dislay_map, ax=plt.gca(), location='right');
	cb.set_label('Coulomb Stress Change (MPa)', fontsize=22);
	for label in cb.ax.yaxis.get_ticklabels():
		label.set_size(18)
	plt.xlabel('Distance', fontsize=18)
	plt.ylabel('Depth (km)', fontsize=18)
	plt.gca().tick_params(labelsize=16)
	plt.gca().invert_yaxis()
	plt.axis('equal')
	plt.ioff()
	print("Saving figure %s " % outfile);
	plt.savefig(outfile, facecolor="w");
	return;


def welcome_and_parse_runstring():
	print("\n\nWelcome to a simple tool for plotting vertical section through a cube of stress results. ");
	parser = argparse.ArgumentParser(description='Run elastic stress models in Python',
									 epilog='\U0001f600 \U0001f600 \U0001f600 ');
	parser.add_argument('config', type=str, help='name of a config file for cube calculation. Required.')
	args = parser.parse_args()
	return args;

if __name__ == "__main__":   # CONFIGURE, INPUT, COMPUTE, OUTPUT

	args = welcome_and_parse_runstring()
	exp_config = PyCoulomb.configure_calc.configure_stress_calculation(args.config);
	if DEBUG:
		print('DEBUG: exp_config:')
		for key in exp_config.__dict__.keys():
			print('\t', key, exp_config.__dict__[key])
	
	[inputs, obs_disp_points, _] = PyCoulomb.input_values.read_inputs(exp_config);
	if DEBUG:
		print('DEBUG: inputs:')
		for key in inputs.__dict__.keys():
			print('\t', key, inputs.__dict__[key])
	if DEBUG:
		print('DEBUG: obs_disp_points:', obs_disp_points)
	#out_object = PyCoulomb.run_dc3d.do_stress_computation(exp_config, inputs, obs_disp_points, ());
	#PyCoulomb.output_manager.produce_outputs(exp_config, inputs, obs_disp_points, (), out_object);

	if DEBUG:
		print('DEBUG: receiver_horiz_profile:')
		for key in inputs.receiver_horiz_profile.__dict__.keys():
			print('\t', key, inputs.receiver_horiz_profile.__dict__[key])
		
	# construct cube of stress results for all horizontal profiles
	out_path_root = os.path.split(exp_config.__dict__["outdir"][ :-1])[0]
	if DEBUG:
		print('DEBUG: out_path_root:', out_path_root)
	horiz_profile_files = sorted(glob.glob(out_path_root + "/*/stresses_horiz_profile.txt"))
	if DEBUG:
		print('DEBUG: horiz_profile_files:', horiz_profile_files)
		
	# read each horizontal profile and populate 3D cube list of Coulomb stresses
# # lon lat depth_km normal_kPa shear_kPa coulomb_kPa
# # strike 6.000000, dip 28.000000, rake -93.000000
# -150.564602 61.129730 45.000000 -0.140972 1.829841 1.773452
# -150.562350 61.129730 45.000000 -0.138325 1.850203 1.794873
# -150.560097 61.129730 45.000000 -0.135666 1.870689 1.816423
# ...
	import csv
	max_nlon = inputs.receiver_horiz_profile.__dict__["shape"][1]
	max_nlat = inputs.receiver_horiz_profile.__dict__["shape"][0]
	max_ndepth = len(horiz_profile_files)
	stress_cube = [[[0 for k in range(max_ndepth)] for j in range(max_nlat)] for i in range(max_nlon)]
	depths = []
	ndepth = 0;
	for horiz_profile_file in horiz_profile_files:
		print('Reading:', horiz_profile_file)
		with open(horiz_profile_file, 'r') as infile:
			reader = csv.reader(infile, delimiter=' ')
			next(reader, None)  # skip the headers
			next(reader, None)  # skip the headers
			nlon = -1
			nlat = 0
			for row in reader:
				nlon += 1
				if (nlon == max_nlon):
					nlon = 0
					nlat += 1
				if (nlat == max_nlat):
					print('ERROR: reading horiz_profile_file: too many latitudes', nlat)
				# put stress value in list
				#print(nlon, "/", max_nlon, nlat, "/", max_nlat, ndepth, "/", max_ndepth)
				stress_cube[nlon][nlat][ndepth] = float(row[5])
		print('   nlon:', nlon, "/", max_nlon)
		print('   nlat:', nlat, "/", max_nlat)
		cdepth = os.path.basename(os.path.dirname(horiz_profile_file))
		print('	Depth:', cdepth, str(float(cdepth)))
		depths.append(float(cdepth))
		ndepth += 1

	stress_cube = np.asarray(stress_cube)
	depths = np.asarray(depths)
	print('depths:', depths)
		
	map_vert_profile(stress_cube, inputs.receiver_horiz_profile, depths, 
		os.path.join(out_path_root, 'vertical_profile_stresses.png'))









