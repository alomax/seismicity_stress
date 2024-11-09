
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
run_elastic_stresses_cube.bash 		# run with RUN_NAME=20240303A_EastWest

# plot vertical section of cube
python ${SEIS_STRESS_PYPATH}/plot_cube_vertical_section.py out/cube/20240303A_EastWest/PointSource/config.txt
# visualize output *.png with your favorite image viewer

# perform seismicity-stress inversion
# best (west-dipping receivers)
python ${SEIS_STRESS_PYPATH}/inv_aftershocks_stress_ffault.py --config out/cube/20240303A_EastWest/PointSource/config.txt --aftershocks NLL_csv_catalogs/20231125A_coherence_Anchorage_2018.1day.csv --max_aftershock_count 1 --diam_mask_stress_center 1.0 --decay_mask_stress_center 10.0  --plot_scale 2.0 --save_cubes --out_path out/cube/20240303A_EastWest/PointSource/FF_1day

# visualize results
# requires that Java 8 or higher is installed and on your PATH 
java -jar SeismicityViewer50.jar net.alomax.seismicity.Seismicity -event.filetype CSV_COMPACT NLL_csv_catalogs/20231125A_coherence_Anchorage_2018.1day.csv &
