# 20240712 - Anthony Lomax, ALomax Scientific

# BASH script to construct a 3D point-source deltaCFS cube using tools from https://github.com/kmaterna/elastic_stresses_py


RUN_NAME=20240303A_EastWest
#RUN_NAME=20240303A_EastEast

PROJECT=Anchorage_2018

# path to output
OUTPATH=out/cube/${RUN_NAME}
# MUST BE ABSOLUTE PATH!
OUTPATH="$(pwd)/${OUTPATH}"
mkdir -p ${OUTPATH}


# select point-source specification
if [ YES == YES ]; then 

	SOURCE_PATCH=PointSource
	# PointSource
	# https://earthquake.usgs.gov/earthquakes/eventpage/ak018fcnsk91/moment-tensor
	# Plane	Strike	Dip	Rake
	# NP1	6°	28°	-93°
	# NP2	189°	62°	-88°
	FAULT_PATCHES_FILE=${OUTPATH}/PointSource.txt
	cat > ${FAULT_PATCHES_FILE}  << END
# Source_FM: strike rake dip lon lat depth_km magnitude
# IMPORTANT! depth best centered between DEPTH_MIN and DEPTH_MAX, Source_FM centered in cube
#
# same fault slip as mainshock **East dipping**
Source_FM: 6.0 -93.0 28.0  -149.954 61.424  47.5  7.1  # https://earthquake.usgs.gov/earthquakes/eventpage/ak018fcnsk91/moment-tensor
#
# same fault slip as mainshock **West dipping** (FOR A PIONT SOURCE NECESSARILY GIVES IDENTICAL RESULTS TO ABOVE)
##Source_FM: 189.0 -88.0 62.0  -149.954 61.424  47.5  7.1  # https://earthquake.usgs.gov/earthquakes/eventpage/ak018fcnsk91/moment-tensor
END
	cat >> ${FAULT_PATCHES_FILE}  << END
# true centroid of fault in stress cube (not needed when Source_FM centered in cube)
# FaultSlipObject_Geographic_Center: lon lat depth
FaultSlipObject_Geographic_Center: -149.954 61.424  47.5
END

fi


OUTPUT_DIR=${OUTPATH}/${SOURCE_PATCH}
mkdir -p ${OUTPUT_DIR}
# save this file
cp -p $0 ${OUTPUT_DIR}



# construct point-source deltaCFS cube as stack of elastic_stresses_py Receiver_Horizontal_Profile's

# depth in meter units
declare -i DEPTH_STEP=500	# 0.5 km
declare -i DEPTH_MIN=25000	# 35 km
declare -i DEPTH_MAX=70000	# 60 km

declare -i DEPTH=DEPTH_MIN
declare -i INDEX=0
while [ $((${DEPTH} <= ${DEPTH_MAX})) == 1 ]; do

	DEPTH_USE=$(echo "scale=3; ${DEPTH}/1000" | bc)
	CDEPTH_USE='000000'${DEPTH_USE}      # get number, pack with zeros
	CDEPTH_USE=${CDEPTH_USE:(-7)}       # the last seven characters
	if [ ${CDEPTH_USE} == 0000000 ]; then 
		CDEPTH_USE=000.000
	fi	
	echo "Running: Depth=${CDEPTH_USE} ==================================="
	mkdir -p ${OUTPUT_DIR}/${CDEPTH_USE}/config
	CONFIG_FILE=${OUTPUT_DIR}/${CDEPTH_USE}/config/config_${CDEPTH_USE}.txt
	INPUT_FILE=${OUTPUT_DIR}/${CDEPTH_USE}/config/inputs_${CDEPTH_USE}.intxt

	# work in tmp sub-directories to prevent file access conflicts during parallel runs of elastic_stresses_driver.py
	rm -rf tmp_work_${CDEPTH_USE}/*
	mkdir tmp_work_${CDEPTH_USE}
	cd tmp_work_${CDEPTH_USE}

	cat << END > ${CONFIG_FILE}
[io-config]
exp_name = ${CDEPTH_USE}
input_file = ${INPUT_FILE}
output_dir = ${OUTPUT_DIR}
plot_stress = 1
plot_grd_disp = 0
gps_disp_points = 
aftershocks = 
strain_file = 

[compute-config]
strike_num_receivers = 1
dip_num_receivers = 1
# mu lambda from MuLambda_ak135-f.xlsx for depth ~45 km
mu = 72e9
lame1 = 87e9
B = 0
fixed_rake = 
END

# select receiver specification
	cat << END > ${INPUT_FILE}
# General: poissons_ratio friction_coef lon_min lon_max lon_zero lat_min lat_max lat_zero
General: 0.250 0.400 -150.4 -149.6 -150.0  61.2 61.6 61.4 

# Receiver_Horizontal_Profile: depth_km strike dip rake centerlon centerlat length_km width_km inc_km
#
# EastEast
# same fault slip as mainshock **East dipping**
##Receiver_Horizontal_Profile: ${DEPTH_USE}  6.0 28.0 -93.0  -149.954 61.424  25.0 25.0 0.5
#
# EastWest
# same fault slip as mainshock **West dipping**
Receiver_Horizontal_Profile: ${DEPTH_USE}  189.0 62.0 -88.0  -149.954 61.424  25.0 25.0 0.5
END

	cat << END >> ${INPUT_FILE}
# ${SOURCE_PATCH}
END
	cat < ${FAULT_PATCHES_FILE} >>  ${INPUT_FILE}
	
	#elastic_stresses_driver ${CONFIG_FILE}		# run sequentially
	nice elastic_stresses_driver ${CONFIG_FILE} &		# run in parallel
	PIDS[${INDEX}]=$!
	INDEX=INDEX+1

	DEPTH=DEPTH+DEPTH_STEP

	cd ..

done

# wait for all PIDS
for PID in ${PIDS[*]}; do
	wait $PID
	status=$?
	echo "Finished: PID=${PID} status=${status} ================================="
done

# save representative config file for use by plot_cube_vertical_section.py
cp -p ${CONFIG_FILE} ${OUTPUT_DIR}/config.txt

rm -r tmp_work_*


echo "Output written to: ${OUTPATH}"

