
# activate conda environment for https://github.com/kmaterna/elastic_stresses_py
conda activate elastic_py

export SEIS_STRESS_PYPATH=../../bin

# set out path
mv out out_ORIG
ln -s <my_out_path> out
# or
mkdir out

# construct point-source deltaCFS cube as stack of elastic_stresses_py Receiver_Horizontal_Profile's
run_elastic_stresses_cube.bash

# plot vertical section of cube
python ${SEIS_STRESS_PYPATH}/plot_cube_vertical_section.py out/cube/20240710A_MSMS/PointSource/config.txt
# visualize output *.png with your favorite image viewer

# perform seismicity-stress inversion
python ${SEIS_STRESS_PYPATH}/inv_aftershocks_stress_ffault.py --config out/cube/20240710A_MSMS/PointSource/config.txt --aftershocks NLL_csv_catalogs/20231223A_coherence_4days.csv --max_aftershock_count 1 --diam_mask_stress_center 1.0 --decay_mask_stress_center 10.0 --stress_depth_top 0.0 --plot_scale 2.0 --save_cubes --out_path out/cube/20240710A_MSMS/PointSource/FF_4days

# visualize results
# requires that Java 8 or higher is installed and on your PATH 
java -jar SeismicityViewer50.jar net.alomax.seismicity.Seismicity -event.filetype CSV_COMPACT NLL_csv_catalogs/20231223A_coherence_4days.csv &

