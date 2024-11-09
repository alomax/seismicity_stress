
# activate conda environment for https://github.com/kmaterna/elastic_stresses_py
conda activate elastic_py

export SEIS_STRESS_PYPATH=../../bin

# set out path
mv out out_ORIG
#
ln -s <my_out_path> out
# or
mkdir out

# construct point-source deltaCFS cube as stack of elastic_stresses_py Receiver_Horizontal_Profile's
# NOTE: if you do not have enough memory, use the "# run sequentially" command line and not "# run in parallel"
run_elastic_stresses_cube.bash

# plot vertical section of cube
python ${SEIS_STRESS_PYPATH}/plot_cube_vertical_section.py out/cube/20240212A/PointSource/config.txt
# visualize output *.png with your favorite image viewer

# perform seismicity-stress inversion
python ${SEIS_STRESS_PYPATH}/inv_aftershocks_stress_ffault.py --config out/cube/20240212A/PointSource/config.txt --aftershocks NLL_csv_catalogs/20231201A_coherence_Parkfield_2022_4hours.csv --max_aftershock_count 1 --diam_mask_stress_center 1.0 --decay_mask_stress_center 10.0 --stress_depth_top 0.0 --plot_scale 2.0 --save_cubes --out_path out/cube/20240212A/PointSource/FF_4hours

# visualize results
# requires that Java 8 or higher is installed and on your PATH 
java -jar SeismicityViewer50.jar net.alomax.seismicity.Seismicity -event.filetype CSV_COMPACT NLL_csv_catalogs/20231201A_coherence_Parkfield_2022_4hours.csv &




# synthetic fault slip on rectangular fault patches
# 20240614A/PatchSynth synthetic
#
# construct point-source deltaCFS cube as stack of elastic_stresses_py Receiver_Horizontal_Profile's
# IMPORTANT: in run_elastic_stresses_cube.bash, uncomment: "RUN_NAME=20240614A ; SOURCE_PATCH=PatchSynth"
run_elastic_stresses_cube.bash
#
# plot vertical section of cube
python ${SEIS_STRESS_PYPATH}/plot_cube_vertical_section.py out/cube/20240614A/PatchSynth/config.txt
#
# generate synthetic aftershocks
python ${SEIS_STRESS_PYPATH}/inv_aftershocks_stress_ffault.py --config out/cube/20240614A/PatchSynth/config.txt --aftershocks NLL_csv_catalogs/20231201A_coherence_Parkfield_2022_1week.csv --max_aftershock_count 1 --diam_mask_stress_center -1.0 --decay_mask_stress_center 10.0 --stress_depth_top 0.0 --plot_scale 4.0 --aftershock_write_rate 0.01 --out_path out/cube/20240614A/PatchSynth/Synth_AS
#
# visualize synthetic aftershocks
# IMPORTANT: in seismicitydefaults, change all: "out/cube/20240212A/PointSource/FF_4hours" -> "out/cube/20240614A/PatchSynth/Synth_AS"
#    and uncomment "seismicity.polygons.0 = maps/FiniteFault_PatchSynth.sdr;SOLID_BORDER;NO_FILL;255,255,255,255;2.0;PatchSynth"
sv -event.filetype BASIC_CSV out/cube/20240614A/PatchSynth/Synth_AS/masked_stress_aftershocks.csv &
#
# perform seismicity-stress inversion with invert synthetic aftershocks
python ${SEIS_STRESS_PYPATH}/inv_aftershocks_stress_ffault.py --config out/cube/20240212A/PointSource/config.txt --aftershocks out/cube/20240614A/PatchSynth/Synth_AS/masked_stress_aftershocks.csv --max_aftershock_count 1 --diam_mask_stress_center 1.0 --decay_mask_stress_center 10.0 --stress_depth_top 0.0 --plot_scale 2.0 --save_cubes --out_path out/cube/20240614A/PatchSynth/Synth_FF
#
# visualize results
# IMPORTANT: in seismicitydefaults, change all: "out/cube/20240614A/PatchSynth/Synth_AS" -> "out/cube/20240614A/PatchSynth/Synth_FF"
java -jar SeismicityViewer50.jar net.alomax.seismicity.Seismicity -event.filetype BASIC_CSV out/cube/20240614A/PatchSynth/Synth_AS/masked_stress_aftershocks.csv &
