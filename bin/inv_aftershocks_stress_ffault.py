#!/usr/bin/env python

import os
import pickle
import sys
import argparse
import glob
import copy
import math
import random
import numpy as np
import scipy
import matplotlib.pyplot as plt
from elastic_stresses_py import PyCoulomb
from Tectonic_Utils.geodesy import fault_vector_functions as fvf

DEBUG = True


# https://stackoverflow.com/questions/55917328/numpy-trim-zeros-in-2d-or-3d
def trim_zeros(arr):
    """Returns a trimmed view of an n-D array excluding any outer
    regions which contain only zeros.
    """
    slices = tuple(slice(idx.min(), idx.max() + 1) for idx in np.nonzero(arr))
    return arr[slices]


def load_aftershocks_to_cube(aftershocks, aftershock_cube, max_aftershock_count,
                             decay_mask_stress_center):
    """Load aftershocks into cube with pseudo-density as function of number of aftershocks in each cell."""

    # aftershock cells centered on stress, ffault grid points
    x_min = -aftershock_cube["xlen"] - aftershock_cube["xy_inc"] / 2.0
    x_max = aftershock_cube["xlen"] + aftershock_cube["xy_inc"] / 2.0
    y_min = -aftershock_cube["ylen"] - aftershock_cube["xy_inc"] / 2.0
    y_max = aftershock_cube["ylen"] + aftershock_cube["xy_inc"] / 2.0
    depth_min = aftershock_cube["centerdepth"] - aftershock_cube["zlen"] - aftershock_cube["z_inc"] / 2.0
    depth_max = aftershock_cube["centerdepth"] + aftershock_cube["zlen"] + aftershock_cube["z_inc"] / 2.0
    corner_ll = fvf.xy2lonlat_single(x_min, y_min, aftershock_cube["centerlon"], aftershock_cube["centerlat"])
    corner_ul = fvf.xy2lonlat_single(x_min, y_max, aftershock_cube["centerlon"], aftershock_cube["centerlat"])
    corner_lr = fvf.xy2lonlat_single(x_max, y_min, aftershock_cube["centerlon"], aftershock_cube["centerlat"])
    corner_ur = fvf.xy2lonlat_single(x_max, y_max, aftershock_cube["centerlon"], aftershock_cube["centerlat"])
    lon_min = min(corner_ll[0], corner_ul[0])
    lon_max = max(corner_lr[0], corner_ur[0])
    lat_min = min(corner_ll[1], corner_lr[1])
    lat_max = max(corner_ul[1], corner_ur[1])
    if False and DEBUG:
        print("lon:", lon_min, lon_max)
        print("lat:", lat_min, lat_max)
        print("depth:", depth_min, depth_max)
        print("x:", x_min, x_max)
        print("y:", y_min, y_max)

    cube = aftershock_cube["cube"]
    numx = cube.shape[0]
    numy = cube.shape[1]
    numz = cube.shape[2]

    # initalize cube
    cube.fill(0)

    aftershocks_in_cube = []

    for lon, lat, depth in aftershocks:
        # check hypocenter in bounding box of aftershocks cube
        if (
                lon >= lon_min and lon <= lon_max and lat >= lat_min and lat <= lat_max and depth > depth_min and depth < depth_max):
            # get exact x/y coords of aftershock
            x, y = fvf.latlon2xy_single(lon, lat, aftershock_cube["centerlon"], aftershock_cube["centerlat"])
            # check hypocenter is inside x/y limits of aftershock cube
            if (x > x_min and x < x_max and y > y_min and y < y_max):
                # get aftershock_cube cell indices
                nx = int((x - x_min) / aftershock_cube["xy_inc"])
                if nx >= 0 and nx < numx:
                    ny = int((y - y_min) / aftershock_cube["xy_inc"])
                    if ny >= 0 and ny < numy:
                        nz = int((depth - depth_min) / aftershock_cube["z_inc"])
                        if nz >= 0 and nz < numz:
                            # increment cell count if cell count less than max count
                            if max_aftershock_count > 0:
                                if (cube[nx][ny][nz] < max_aftershock_count):
                                    if (cube[nx][ny][nz] < 1):  # no previous aftershock in cell
                                        cube[nx][ny][nz] = 1
                                    else:
                                        cube[nx][ny][nz] += 1
                            else:
                                cube[nx][ny][nz] += 1
                            # add aftershock to list of aftershocks in cube
                            aftershocks_in_cube.append([lon, lat, depth])

    # fill empty cells with negative value that balances sum of negative and positive cells
    # DO NOT USE, DOES NOT WORK!
    if False:
        # trimmed view of aftershock_cube excluding any outer regions which contain only zeros
        if DEBUG:
            print("")
            print("aftershock_cube:", aftershock_cube["cube"].shape)
            for k, v in aftershock_cube.items():
                if k != "cube":
                    print(k, v, " ", end="")
            print("");
            print("")
        print("aftershock_cube size:", np.size(cube), cube.shape)
        trimmed_cube = trim_zeros(cube)
        trimmed_aftershock_cube_size = np.size(trimmed_cube)
        print("trimmed_aftershock_cube size:", trimmed_aftershock_cube_size, trimmed_cube.shape)

        trimmed_aftershock_cube_size = pow(decay_mask_stress_center / aftershock_cube["xy_inc"], 3)
        print("trimmed_aftershock_cube size:", trimmed_aftershock_cube_size)

        # print(cube)
        sum_pos = np.sum(cube)
        n_pos = np.count_nonzero(cube)
        print("DEBUG: sum_pos, n_pos", sum_pos, n_pos)
        if sum_pos > 0.0 and n_pos > 0:
            # neg_fill_value = -sum_pos / float(cube.size - n_pos)	# no difference
            # neg_fill_value = -1.0 / sum_pos						# very little difference
            # neg_fill_value = -1.0									# very large spreading of ffault prob
            ##neg_fill_value = -sum_pos / pow(np.size(cube), 2.0 / 3.0)			#
            ##neg_fill_value = -sum_pos / trimmed_aftershock_cube_size			#
            neg_fill_value = -0.02  #
            print("DEBUG: sum_pos, n_pos, neg_fill_value", sum_pos, n_pos, neg_fill_value)
            cube[cube < 0.5] = neg_fill_value
    # print(cube)

    return aftershock_cube, aftershocks_in_cube


def mask_cube_aftershocks(ffault_cube, aftershocks, no_mask_depth):
    if DEBUG:
        print("Masking finite-fault cube to convex hull of aftershocks inside aftershock cube...")
        print("")

    masked_cube = copy.deepcopy(ffault_cube)
    mask_cube = masked_cube["cube"]
    centerlon_masked_cube = masked_cube["centerlon"]
    centerlat_masked_cube = masked_cube["centerlat"]

    # print("aftershocks:", aftershocks)
    aftershock_array = np.array(aftershocks)  # lon, lat, depth
    # print("aftershock_array:", aftershock_array)

    # https://stackoverflow.com/questions/16750618/whats-an-efficient-way-to-find-if-a-point-lies-in-the-convex-hull-of-a-point-cl
    from scipy.spatial import Delaunay
    aftershock_array_hull = aftershock_array
    if no_mask_depth:
        aftershock_array_hull = np.delete(aftershock_array, 2, 1)
    hull = Delaunay(aftershock_array_hull)
    # print("hull:", hull)

    # mask if outside convex hull of aftershocks
    xpos = -masked_cube["xlen"]
    nx = 0
    for array1 in mask_cube:
        ypos = -masked_cube["ylen"]
        ny = 0
        for array2 in array1:
            depth = -masked_cube["zlen"] + masked_cube["centerdepth"]
            nz = 0
            for value in array2:
                lon, lat = fvf.xy2lonlat_single(xpos, ypos, centerlon_masked_cube, centerlat_masked_cube)
                # in_hull(points, x) check if inside
                if no_mask_depth:
                    point = np.array([lon, lat])
                else:
                    point = np.array([lon, lat, depth])
                if hull.find_simplex(point) < 0:
                    mask_cube[nx][ny][nz] = 0.0
                # else:
                #	print("lon, lat, depth, point, inside:", lon, lat, depth, point, hull.find_simplex(point))
                depth += masked_cube["z_inc"]
                nz += 1
            ypos += masked_cube["xy_inc"]
            ny += 1
        xpos += masked_cube["xy_inc"]
        nx += 1

    return masked_cube


def make_same_shape(bigger_cube, smaller_cube):
    current_shape = bigger_cube.shape
    ncenter = tuple(int(ti / 2) for ti in current_shape)

    sshape = smaller_cube.shape
    nlen = tuple(int(ti / 2) for ti in sshape)

    print("make_same_shape: bigger_cube.shape, smaller_cube.shape", bigger_cube.shape, sshape)
    print("make_same_shape: ncenter, nlen", ncenter, nlen)

    shrunk_cube = bigger_cube[
                  ncenter[0] - nlen[0]: ncenter[0] + nlen[0] + sshape[0] % 2,
                  ncenter[1] - nlen[1]: ncenter[1] + nlen[1] + sshape[1] % 2,
                  ncenter[2] - nlen[2]: ncenter[2] + nlen[2] + sshape[2] % 2
                  ]
    print("make_same_shape: shrunk_cube.shape", shrunk_cube.shape)

    return shrunk_cube


# https://stackoverflow.com/questions/17473917/is-there-a-equivalent-of-scipy-signal-deconvolve-for-2d-arrays
from scipy import fftpack


def convolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft * psf_fft)))


def deconvolve(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft / psf_fft)))


def pre_process_stress_cube(stress_cube):
    stress_cube_process_sign = np.sign(stress_cube)
    stress_cube_process_abs = np.abs(stress_cube)
    stress_cube_process = stress_cube_process_abs

    if False:
        # log stress cube
        print("np.percentile(stress_cube_process_abs, 1)", np.percentile(stress_cube_process_abs, 1))
        print("np.percentile(stress_cube_process_abs, 2)", np.percentile(stress_cube_process_abs, 2))
        print("np.percentile(stress_cube_process_abs, 5)", np.percentile(stress_cube_process_abs, 5))
        print("np.percentile(stress_cube_process_abs, 10)", np.percentile(stress_cube_process_abs, 10))
        print("np.percentile(stress_cube_process_abs, 50)", np.percentile(stress_cube_process_abs, 50))
        print("np.percentile(stress_cube_process_abs, 60)", np.percentile(stress_cube_process_abs, 60))
        print("np.percentile(stress_cube_process_abs, 70)", np.percentile(stress_cube_process_abs, 70))
        print("np.percentile(stress_cube_process_abs, 80)", np.percentile(stress_cube_process_abs, 80))
        print("np.percentile(stress_cube_process_abs, 90)", np.percentile(stress_cube_process_abs, 90))
        print("np.percentile(stress_cube_process_abs, 95)", np.percentile(stress_cube_process_abs, 95))
        print("np.percentile(stress_cube_process_abs, 99)", np.percentile(stress_cube_process_abs, 99))
        print("np.percentile(stress_cube_process_abs, 100)", np.percentile(stress_cube_process_abs, 100))
        waterline = np.percentile(stress_cube_process_abs, 90)
        stress_cube_waterline = np.maximum(stress_cube_process_abs, waterline)
        stress_cube_process = np.log(stress_cube_waterline)

    if False:
        # root stress cube
        stress_cube_process = np.cbrt(stress_cube_process_abs)

    stress_cube_process = np.multiply(stress_cube_process_sign, stress_cube_process)
    stress_cube_process = np.nan_to_num(stress_cube_process, copy=False, nan=0.0, posinf=None, neginf=None)

    return stress_cube_process


def correlate_stress_aftershock_cubes(stress_cube, aftershock_cube):
    # aftershock_cube contains stress_cube, so aftershock_cube is like input, stress_cube like weights

    # output is same size as aftershock_cube (?)

    do_correlation = True
    do_convolve = False
    do_deconvolve = False

    stress_cube_process = stress_cube["cube"]
    # log or root stress cube
    # stress_cube_process = pre_process_stress_cube(stress_cube["cube"])

    if do_correlation:
        # correlate
        ffault_cube_array = scipy.signal.correlate(aftershock_cube["cube"], stress_cube_process, mode='valid',
                                                   method='auto')
    elif do_convolve:
        ffault_cube_array = scipy.signal.convolve(aftershock_cube["cube"], stress_cube_process, mode='valid',
                                                  method='auto')
    elif do_deconvolve:
        # deconvolve
        # make sure cubes have same dimensions
        reduced_aftershock_cube = make_same_shape(aftershock_cube["cube"], stress_cube_process)
        normalized_stress_cube_process = stress_cube_process / np.max(
            np.abs(stress_cube_process))  # normalize to abs max
        waterlevel = 0.01
        normalized_stress_cube_process += waterlevel
        ffault_cube_array = deconvolve(reduced_aftershock_cube, normalized_stress_cube_process)
        ffault_cube_array = np.real(ffault_cube_array)

    print("stress_cube.shape", stress_cube_process.shape)
    print("aftershock_cube.shape", aftershock_cube["cube"].shape)
    print("ffault_cube_array.shape", ffault_cube_array.shape)

    ffault_cube = {
        "cube": ffault_cube_array,
        "centerlon": aftershock_cube["centerlon"],
        "centerlat": aftershock_cube["centerlat"],
        "centerdepth": aftershock_cube["centerdepth"],
        "ylen": (ffault_cube_array.shape[1] - 1) * aftershock_cube["xy_inc"] / 2.0,
        "xlen": (ffault_cube_array.shape[0] - 1) * aftershock_cube["xy_inc"] / 2.0,
        "xy_inc": aftershock_cube["xy_inc"],
        "zlen": (ffault_cube_array.shape[2] - 1) * aftershock_cube["z_inc"] / 2.0,
        "z_inc": aftershock_cube["z_inc"]
    }
    if DEBUG:
        print("")
        print("ffault_cube:", ffault_cube["cube"].shape)
        for k, v in ffault_cube.items():
            if k != "cube":
                print(k, v, " ", end="")
        print("");
        print("")

    return ffault_cube


def inv_aftershocks_stress_ffault(stress_cube, aftershocks, max_aftershock_count,
                                  decay_mask_stress_center):
    """Invert for finite-fault response in cube from mapping of Coulomb stress cube over aftershock distribution."""

    # assemble aftershock density cube
    # aftershock cube must include volume defined by stress cube when its center is at any edge point of finte faulting cube
    xy_inc_aftershock = stress_cube["xy_inc"]
    aftershock_cube_xlen = 2.0 * stress_cube["xlen"]
    aftershock_cube_ylen = 2.0 * stress_cube["ylen"]
    num_x_aftershock = int((2.0 * aftershock_cube_xlen + 0.001) / xy_inc_aftershock) + 1
    num_y_aftershock = int((2.0 * aftershock_cube_ylen + 0.001) / xy_inc_aftershock) + 1
    z_inc_aftershock = stress_cube["z_inc"]
    aftershock_cube_zlen = 2.0 * stress_cube["zlen"]
    num_z_aftershock = int((2.0 * aftershock_cube_zlen + 0.001) / z_inc_aftershock) + 1
    # initialize aftershock cube array with -1, indicates no aftershock in cell; positive values indicate aftershocks present; zero never used
    aftershock_cube_array = np.full([num_x_aftershock, num_y_aftershock, num_z_aftershock], -1, dtype=float)
    aftershock_cube = {
        "cube": aftershock_cube_array,
        "centerlon": stress_cube["centerlon"],
        "centerlat": stress_cube["centerlat"],
        "centerdepth": stress_cube["centerdepth"],
        "ylen": aftershock_cube_ylen,
        "xlen": aftershock_cube_xlen,
        "xy_inc": xy_inc_aftershock,
        "zlen": aftershock_cube_zlen,
        "z_inc": z_inc_aftershock
    }

    # load aftershock density to cube
    aftershock_cube, aftershocks_in_cube = \
        load_aftershocks_to_cube(aftershocks, aftershock_cube, max_aftershock_count,
                                 decay_mask_stress_center)

    # print("DEBUG: len(aftershocks), len(aftershocks_in_cube), ", len(aftershocks), len(aftershocks_in_cube))

    if DEBUG:
        print("")
        print("aftershock_cube:", aftershock_cube["cube"].shape)
        for k, v in aftershock_cube.items():
            if k != "cube":
                print(k, v, " ", end="")
        print("");
        print("")
        hist, bin_edges = np.histogram(aftershock_cube["cube"], bins=[-1.5, -0.5, 0.5, 1.5, 2.5, 3.5])  # -1,0,1,2,3
        print("bin_edges", bin_edges)
        print("           -1       0       1       2       3")
        print("hist", hist)

    # correlate stress cube and aftershock cube over volume defined by aftershock_cube, contains original finite-fault cube
    if DEBUG:
        print("")
        print("Correlating stress and aftershock cubes...")
    ffault_cube = correlate_stress_aftershock_cubes(stress_cube, aftershock_cube)

    return ffault_cube, aftershock_cube, aftershocks_in_cube


def gaussian(x, sigma):
    """ Return the unnormalized Gaussian (max amplitude == 1) with standard deviation sigma. """
    return np.exp(-0.5 * (x / sigma) ** 2)


def mask_cube_center_CONSTANT_TEST(stress_cube, diam_mask_stress_center, decay_mask_stress_center,
                                   fault_geographic_centerlon, fault_geographic_centerlat,
                                   fault_geographic_centerdepth):
    masked_stress_cube = copy.deepcopy(stress_cube)

    cube = stress_cube["cube"]
    masked_cube = masked_stress_cube["cube"]

    # adjust center for fault center (default is center at center of cube)
    xcent = 0.0
    ycent = 0.0
    zcent = 0.0
    if fault_geographic_centerlon != None and fault_geographic_centerlat != None and fault_geographic_centerdepth != None:
        xcent, ycent = fvf.latlon2xy_single(fault_geographic_centerlon, fault_geographic_centerlat,
                                            stress_cube["centerlon"], stress_cube["centerlat"])
        zcent = fault_geographic_centerdepth - stress_cube["centerdepth"]

    # find value closest to diam_mask_stress_center
    mask_value_pos = 0.0
    mask_value_neg = 0.0
    x = -stress_cube["xlen"] - xcent
    nx = 0
    for array1 in cube:
        y = -stress_cube["ylen"] - ycent
        ny = 0
        for array2 in array1:
            depth_offset = -stress_cube["zlen"] - zcent
            nz = 0
            for value in array2:
                distance = math.sqrt(x * x + y * y + depth_offset * depth_offset)
                if distance > diam_mask_stress_center / 2.0:
                    mask_value_test = value
                    if mask_value_test > mask_value_pos:
                        mask_value_pos = mask_value_test
                    elif mask_value_test < mask_value_neg:
                        mask_value_neg = mask_value_test
                depth_offset += stress_cube["z_inc"]
                nz += 1
            y += stress_cube["xy_inc"]
            ny += 1
        x += stress_cube["xy_inc"]
        nx += 1
    print("mask_cube_center: mask_value_neg", mask_value_neg, "mask_value_pos:", mask_value_pos)
    print("")

    # mask central values
    x = -stress_cube["xlen"] - xcent
    nx = 0
    for array1 in cube:
        y = -stress_cube["ylen"] - ycent
        ny = 0
        for array2 in array1:
            depth_offset = -stress_cube["zlen"] - zcent
            nz = 0
            for value in array2:
                distance = math.sqrt(x * x + y * y + depth_offset * depth_offset)
                if distance < diam_mask_stress_center / 2.0:
                    if value > mask_value_pos:
                        masked_cube[nx, ny, nz] = mask_value_pos
                    elif value < mask_value_neg:
                        masked_cube[nx, ny, nz] = mask_value_neg
                    else:
                        masked_cube[nx, ny, nz] = value
                else:
                    masked_cube[nx, ny, nz] = value
                depth_offset += stress_cube["z_inc"]
                nz += 1
            y += stress_cube["xy_inc"]
            ny += 1
        x += stress_cube["xy_inc"]
        nx += 1

    return masked_stress_cube


def mask_cube_center_BINARY_TEST(stress_cube, diam_mask_stress_center, decay_mask_stress_center,
                                 fault_geographic_centerlon, fault_geographic_centerlat, fault_geographic_centerdepth):
    masked_stress_cube = copy.deepcopy(stress_cube)

    cube = stress_cube["cube"]
    masked_cube = masked_stress_cube["cube"]

    # adjust center for fault center (default is center at center of cube)
    xcent = 0.0
    ycent = 0.0
    zcent = 0.0
    if fault_geographic_centerlon != None and fault_geographic_centerlat != None and fault_geographic_centerdepth != None:
        xcent, ycent = fvf.latlon2xy_single(fault_geographic_centerlon, fault_geographic_centerlat,
                                            stress_cube["centerlon"], stress_cube["centerlat"])
        zcent = fault_geographic_centerdepth - stress_cube["centerdepth"]

    # mask central values
    x = -stress_cube["xlen"] - xcent
    nx = 0
    for array1 in cube:
        y = -stress_cube["ylen"] - ycent
        ny = 0
        for array2 in array1:
            depth_offset = -stress_cube["zlen"] - zcent
            nz = 0
            for value in array2:
                distance = math.sqrt(x * x + y * y + depth_offset * depth_offset)
                mask = 1.0
                radius_mask_stress_center = diam_mask_stress_center / 2.0
                if distance <= radius_mask_stress_center:
                    mask = np.sign(value) * 0.0
                else:
                    # mask = np.sign(value) * (radius_mask_stress_center * radius_mask_stress_center) / (distance * distance)
                    mask = value * (radius_mask_stress_center * radius_mask_stress_center) / (distance * distance)
                masked_cube[nx, ny, nz] = mask
                depth_offset += stress_cube["z_inc"]
                nz += 1
            y += stress_cube["xy_inc"]
            ny += 1
        x += stress_cube["xy_inc"]
        nx += 1

    return masked_stress_cube


def mask_cube_center(stress_cube, diam_mask_stress_center, decay_mask_stress_center,
                     fault_geographic_centerlon, fault_geographic_centerlat, fault_geographic_centerdepth):
    masked_stress_cube = copy.deepcopy(stress_cube)

    cube = stress_cube["cube"]
    masked_cube = masked_stress_cube["cube"]

    # adjust center for fault center (default is center at center of cube)
    xcent = 0.0
    ycent = 0.0
    zcent = 0.0
    if fault_geographic_centerlon != None and fault_geographic_centerlat != None and fault_geographic_centerdepth != None:
        xcent, ycent = fvf.latlon2xy_single(fault_geographic_centerlon, fault_geographic_centerlat,
                                            stress_cube["centerlon"], stress_cube["centerlat"])
        zcent = fault_geographic_centerdepth - stress_cube["centerdepth"]

    # mask central values
    radius_mask_stress_center = diam_mask_stress_center / 2.0
    x = -stress_cube["xlen"] - xcent
    nx = 0
    for array1 in cube:
        y = -stress_cube["ylen"] - ycent
        ny = 0
        for array2 in array1:
            depth_offset = -stress_cube["zlen"] - zcent
            nz = 0
            for value in array2:
                distance = math.sqrt(x * x + y * y + depth_offset * depth_offset)
                mask = 1.0
                if distance < radius_mask_stress_center:
                    mask = 0.0
                else:
                    if decay_mask_stress_center > 1.0e-6:
                        mask = 1.0 - gaussian(distance - radius_mask_stress_center, decay_mask_stress_center)
                masked_cube[nx, ny, nz] = value * mask
                depth_offset += stress_cube["z_inc"]
                nz += 1
            y += stress_cube["xy_inc"]
            ny += 1
        x += stress_cube["xy_inc"]
        nx += 1

    return masked_stress_cube


def mask_cube_quantile(stress_cube, clip_quantile):
    masked_stress_cube = copy.deepcopy(stress_cube)

    masked_cube = masked_stress_cube["cube"]

    # get quantile independently from negative and positive values in cube
    flat_cube = np.reshape(masked_cube, (np.size(masked_cube),))
    quantile_high = np.quantile(flat_cube[flat_cube > 0.0], clip_quantile)
    quantile_low = np.quantile(flat_cube[flat_cube < 0.0], 1.0 - clip_quantile)

    if quantile_high > 0.0 and quantile_low < 0.0:
        clip = min(quantile_high, -quantile_low)
    elif quantile_high > 0.0:
        clip = quantile_high
    elif quantile_low < 0.0:
        clip = quantile_low

    print("clip, quantile_low, quantile_high", clip, quantile_low, quantile_high)
    masked_cube = np.clip(masked_cube, -clip, clip)

    masked_stress_cube["cube"] = masked_cube

    return masked_stress_cube


def mask_cube_binary(stress_cube, clip_quantile):
    masked_stress_cube = copy.deepcopy(stress_cube)

    masked_cube = masked_stress_cube["cube"]

    # get quantile independently from negative and positive values in cube
    flat_cube = np.reshape(masked_cube, (np.size(masked_cube),))
    quantile_high = np.quantile(flat_cube[flat_cube > 0.0], clip_quantile)
    quantile_low = np.quantile(flat_cube[flat_cube < 0.0], 1.0 - clip_quantile)

    if quantile_high > 0.0 and quantile_low < 0.0:
        clip = min(quantile_high, -quantile_low)
    elif quantile_high > 0.0:
        clip = quantile_high
    elif quantile_low < 0.0:
        clip = quantile_low

    print("clip, quantile_low, quantile_high", clip, quantile_low, quantile_high)
    masked_cube = np.clip(masked_cube, -clip, clip)
    masked_cube_int = np.zeros_like(masked_cube)
    masked_cube_int[masked_cube >= 0.99999 * clip] = 1
    masked_cube_int[masked_cube <= -(0.99999 * clip)] = -1

    masked_stress_cube["cube"] = masked_cube_int

    return masked_stress_cube


def read_aftershock_table_NLL_csv_compact(infile):
    """
    NLL compact CSV catalog format:
    date-time, latitude, longitude, depth, RMS, Nphs, Gap, Dist, errH, errZ, Mamp, Mdur, expect_lat, expect_lon, expect_z, EllipsoidAz1, EllipsoidDip1, EllipsoidLen1, EllipsoidAz2, EllipsoidDip2, EllipsoidLen2, EllipsoidLen3, EllipseAz1, EllipseLen1, EllipseAz2, EllipseLen2, publicId
2014-01-01T23:45:56.2806,61.4544,-150.0162,46.3281,0.17,14,111,38.87,2.4794,4.9719,1.20,-9.90,61.4552,-150.0180,46.7870, 7.94,-1.11,1.92,-82.14,-4.52,3.04,6.18,261.62,1.55,171.62,2.48, quakeml:earthquake.alaska.edu/event/01421i2zj
    ...

    modified from Elastic_stresses_py/PyCoulomb/io_additionals.py
    """

    print("Reading aftershocks from file %s " % infile);

    aftershocks = []

    nexcept = 0
    with open(infile) as ifile:
        for line in ifile:
            time = "ERROR";
            lat = -9999;
            lon = -9999;
            depth = -9999;
            magnitude = -9999;
            tokens = line.split(',');
            if tokens[0][0] == '#':
                continue;
            else:
                try:
                    time = tokens[0];
                    lat = float(tokens[1]);
                    lon = float(tokens[2]);
                    depth = float(tokens[3]);
                    magnitude = float(tokens[10]);
                    # aftershocks.append((lon, lat, depth, magnitude, time))
                    aftershocks.append((lon, lat, depth))
                except:
                    if not "latitude" in line:  # header line OK
                        print("Error reading line:", line)
                        nexcept += 1
                    continue

    print("Number aftershocks read:", len(aftershocks))

    if nexcept > 0:
        print("Warning:", nexcept, "errors reading aftershock table:", infile)

    return aftershocks;


def read_aftershock_table(infile):
    """
    Simple catalog format: time, lon, lat, depth, magnitude
    modified from Elastic_stresses_py/PyCoulomb/io_additionals.py
    """
    print("Reading aftershocks from file %s " % infile);
    aftershocks = []

    ifile = open(infile);
    for line in ifile:
        temp = line.split();
        if temp[0][0] == '#':
            continue;
        else:
            time = temp[0];
            lon = float(temp[1]);
            lat = float(temp[2]);
            depth = float(temp[3]);
            magnitude = float(temp[4]);
            # aftershocks.append((lon, lat, depth, magnitude, time))
            aftershocks.append((lon, lat, depth))
    ifile.close();

    return aftershocks;


def write_output(a_ffault_cube, a_masked_ffault_cube, plot_ffault_norm, a_stress_cube, a_masked_stress_cube,
                 aftershock_write_rate, out_path_root, no_output_neg, plot_scale, random_seed):

    write_gmt_points_cube(a_ffault_cube, out_path_root, "relative_fault_slip",
                          True, plot_ffault_norm, no_output_neg,
                          plot_scale, False, random_seed)

    if a_masked_ffault_cube != None:
        write_gmt_points_cube(a_masked_ffault_cube, out_path_root, "relative_fault_slip_masked",
                              True, plot_ffault_norm,
                              no_output_neg, plot_scale, False, random_seed)

    write_gmt_points_cube(a_stress_cube, out_path_root, "stress",
                          False, -1.0, no_output_neg, plot_scale, False, random_seed)
    # stress_cube_process = copy.deepcopy(a_stress_cube)
    # stress_cube_process["cube"] = pre_process_stress_cube(stress_cube_process["cube"])
    # write_gmt_points_cube(stress_cube_process, out_path_root, "stress_processsed", False, -1.0, no_output_neg, plot_scale, False, random_seed)
    write_gmt_points_cube(a_masked_stress_cube, out_path_root, "masked_stress",
                          False, -1.0, no_output_neg, plot_scale, False, random_seed)
    # stress_cube_process = copy.deepcopy(a_masked_stress_cube)
    # stress_cube_process["cube"] = pre_process_stress_cube(stress_cube_process["cube"])
    # write_gmt_points_cube(stress_cube_process, out_path_root, "masked_stress_processsed", False, -1.0, no_output_neg, plot_scale, False, random_seed)

    write_synthetic_aftershocks(a_masked_stress_cube, aftershock_write_rate, out_path_root, "masked_stress")

def write_synthetic_aftershocks(a_cube, write_rate, out_path_root, cube_name):

    cube = a_cube["cube"]

    value_max = np.amax(cube)  # will be positive if positive values present
    print("value_max", value_max)
    value_max = np.quantile(cube, 0.99)  # will be positive if positive values present
    print("value_max_99", value_max)

    centerlon = a_cube["centerlon"]
    centerlat = a_cube["centerlat"]

    rand_range_xy = a_cube["xy_inc"] / 2.0
    rand_range_z = a_cube["z_inc"] / 2.0
    #random.seed(54321)
    random.seed(12345)

    # write aftershocks in SeismicityViewer BASIC_CSV format
    asfile = os.path.join(out_path_root, cube_name + "_aftershocks.csv")
    print("Writing:", asfile)
    n_written = 0
    with open(asfile, 'w') as ofile:
        ofile.write("date-time, latitude, longitude, depth, RMS, Nphs, Gap, Dist, errH, errZ, Mamp\n")
        ofile.write("# " + out_path_root + "\n")
        x = -a_cube["xlen"]
        for array1 in cube:
            y = -a_cube["ylen"]
            for array2 in array1:
                depth = a_cube["centerdepth"] - a_cube["zlen"]
                for value in array2:
                    if random.random() < write_rate:  # reduce number saved
                        if value > sys.float_info.min:
                            value = value / value_max  # 0.0 -> 1.0
                            if value > 1.0 or value > random.random():  # save in proportion to value
                                x_rand = x + random.uniform(-rand_range_xy, rand_range_xy)
                                y_rand = y + random.uniform(-rand_range_xy, rand_range_xy)
                                lon, lat = fvf.xy2lonlat_single(x_rand, y_rand, centerlon, centerlat)
                                depth_rand = depth + random.uniform(-rand_range_z, rand_range_z)
                                ofile.write("2000-01-01T00:00:00.00Z, " + \
                                            str(lat) + ", " + str(lon) + ", " + str(
                                    depth_rand) + ", -, -, - , -, -, - ," + str(value) + "\n")
                                n_written += 1
                    depth += a_cube["z_inc"]
                y += a_cube["xy_inc"]
            x += a_cube["xy_inc"]

    print("Number synthetic aftershocks written:", n_written)


def write_gmt_points_cube(a_cube, out_path_root, cube_name, output_50percent, value_norm,
                          no_output_neg, plot_scale, square_amplitudes, random_seed):
    cube = a_cube["cube"]

    print("write_gmt_points_cube: value_norm", value_norm)
    value_min = np.amin(cube)  # will be negative if negative values present
    value_max = np.amax(cube)  # will be positive if positive values present
    print("write_gmt_points_cube: value_min, value_max", value_min, value_max)
    if value_norm <= 0.0:
        # value_norm = max(abs(value_min), abs(value_max))
        value_norm = abs(value_max)  # normalize to positive max
    print("write_gmt_points_cube: value_norm", value_norm)

    with open(os.path.join(out_path_root, cube_name + "_value_norm.log"), "w") as ofile_log:
        ofile_log.write(str(value_min) + ", " + str(value_max) + ", " + str(value_norm) + "\n")

    centerlon = a_cube["centerlon"]
    centerlat = a_cube["centerlat"]

    rand_range_xy = a_cube["xy_inc"] / 2.0
    rand_range_z = a_cube["z_inc"] / 2.0
    random.seed(random_seed)
    symbol_scale = 4.0 * min(a_cube["xy_inc"], a_cube["z_inc"])  # size of largest symbol in km

    plot_scale_std = plot_scale[0]
    if len(plot_scale) > 1:
        plot_scale_50 = plot_scale[1]
    else:
        plot_scale_50 = plot_scale[0]

    # write ffault correlation values in GMT_LATLONDEPTH format
    filepos = os.path.join(out_path_root, cube_name + "_pos.xyz")
    print("Writing:", filepos)
    fileneg = os.path.join(out_path_root, cube_name + "_neg.xyz")
    if no_output_neg:
        fileneg = os.path.join(out_path_root, "no_output_neg.xyz")
    else:
        print("Writing:", fileneg)
    # filepos_cutoff = sys.float_info.min
    filepos_cutoff = 0.1 * value_norm
    if output_50percent:
        filepos_50 = os.path.join(out_path_root, cube_name + "_pos_50.xyz")
        ##flat_cube = np.reshape(cube, (np.size(cube),))
        ##filepos_50_cutoff = np.quantile(flat_cube[flat_cube > sys.float_info.min], 0.99)
        ##filepos_50_cutoff = 0.5 * value_norm
        # filepos_50_cutoff = gaussian(1.0, 1.0) * value_norm
        # 50% makes sense since "negative" aftershocks missing, e.g. for cosine, sum is 50% instead of 0 when shifted by 1/2 cycle
        filepos_50_cutoff = 0.5 * value_norm
        print("write_gmt_points_cube: filepos_50_cutoff=", filepos_50_cutoff)
        ofile_pos_50 = open(filepos_50, 'w')
        print("Writing:", filepos_50)
        ofile_pos_50.write("> GMT_LATLONDEPTH\n")
        ofile_pos_50.write("> " + out_path_root + "\n")
    # max stats
    value_max_lon = 0.0
    value_max_lat = 0.0
    value_max_z = 0.0
    value_write_max = 0.0
    with (open(filepos, 'w') as ofile_pos, \
            open(fileneg, 'w') as ofile_neg):
        ofile_pos.write("> GMT_LATLONDEPTH\n")
        ofile_pos.write("> " + out_path_root + "\n")
        if not no_output_neg:
            ofile_neg.write("> GMT_LATLONDEPTH\n")
            ofile_neg.write("> " + out_path_root + "\n")
        x = -a_cube["xlen"]
        for array1 in cube:
            y = -a_cube["ylen"]
            for array2 in array1:
                depth = a_cube["centerdepth"] - a_cube["zlen"]
                for value in array2:
                    x_rand = x + random.uniform(-rand_range_xy, rand_range_xy)
                    y_rand = y + random.uniform(-rand_range_xy, rand_range_xy)
                    lon, lat = fvf.xy2lonlat_single(x_rand, y_rand, centerlon, centerlat)
                    depth_rand = depth + random.uniform(-rand_range_z, rand_range_z)
                    if value > filepos_cutoff:
                        value_write = value / value_norm  # 0.0 -> 1.0   # normalize to positive max
                        if square_amplitudes:
                            value_write *= value_write
                        value_write *= symbol_scale
                        ofile_pos.write(
                            str(lat) + " " + str(lon) + " " + str(depth_rand) + " " + str(value_write * plot_scale_std) + "\n")
                        if value_write > value_write_max:
                            value_write_max = value_write
                            value_max_lon = lon
                            value_max_lat = lat
                            value_max_z = depth
                        if output_50percent and value > filepos_50_cutoff:
                            # value_50_write = value / value_norm  # 0.0 -> 1.0   # normalize to positive max
                            value_50_write = (value - filepos_50_cutoff) / (value_norm - filepos_50_cutoff)  # 0.0 -> 1.0
                            ofile_pos_50.write(
                                str(lat) + " " + str(lon) + " " + str(depth_rand) + " " + str(value_50_write * plot_scale_50) + "\n")
                    elif not no_output_neg and value < -filepos_cutoff:
                        value_write = value / value_norm  # 0.0 -> 1.0   # normalize to positive max
                        if square_amplitudes:
                            value_write *= value_write
                        value_write *= symbol_scale
                        ofile_neg.write(
                            str(lat) + " " + str(lon) + " " + str(depth_rand) + " " + str(value_write * plot_scale_std) + "\n")
                    depth += a_cube["z_inc"]
                y += a_cube["xy_inc"]
            x += a_cube["xy_inc"]

    if output_50percent:
        ofile_pos_50.close()

    file_max = os.path.join(out_path_root, cube_name + "_max.xyz")
    with open(file_max, 'w') as ofile_cent:
        ofile_cent.write("> GMT_LATLONDEPTH\n")
        ofile_cent.write(str(value_max_lat) + " " + str(value_max_lon) + " " + str(value_max_z) + " " + str(2.0 * plot_scale_std) + "\n")

def write_cubes_output(a_cube, cube_name, out_path_root):
    picklefile = os.path.join(out_path_root, cube_name + ".pkl")
    with open(picklefile, 'wb') as fileout:
        print('Writing: ', picklefile)
        fileout.write(pickle.dumps(a_cube))


def welcome_and_parse_runstring():
    print("\n\nTool to invert aftershock seismicity and point-source Coulomb stress for finite-faulting. ");
    print("   Requires:");
    print(
        "	  - a set of Receiver_Horizontal_Profile stress files for a centered point source (see run_elastic_stresses_cube.bash)");
    print("		 used to construct a Coulomb stress cube with point source at center. ");
    print("	  - a list of aftershocks");
    print("	  - definition of target finite-faulting volume (cube)");
    print("   Outputs GMT format 3D point plot data for stress and finite-faulting results. ");
    print("   Optional output of 3D cubes as pickled structures: aftershock, stress and finite-faulting cubes. ");
    parser = argparse.ArgumentParser(description="Invert aftershocks and stress for finite-faulting",
                                     epilog="\U0001f600 \U0001f600 \U0001f600 ");
    parser.add_argument("--config", type=str, required=True,
                        help="name of a config.txt file used for a Receiver_Horizontal_Profile. Required.")
    parser.add_argument("--out_path", type=str, default=None,
                        help="Path for output. Default is path to config.txt file specified in --config.")
    parser.add_argument("--aftershocks", type=str, required=True,
                        help="file with aftershock hypocenters in simple format: date lon lat depth mag. Required.")
    parser.add_argument("--stress_measure", type=str, default="coulomb",
                        help="Stress measure to use to construct point source kernel: coulomb shear normal.")
    parser.add_argument("--max_aftershock_count", type=int, default=1,
                        help="Maximum count of afteshocks in single cube cell. < 0 for no maximum.")
    parser.add_argument('--no_mask_to_aftershocks', action='store_true', default=False,
                        help='Do not mask finite-faulting map to hull of aftershocks.')
    parser.add_argument('--no_mask_depth_aftershocks', action='store_true', default=False,
                        help='Do not mask finite-faulting map to depth of aftershocks, only mask lat/lon.')
    parser.add_argument("--stress_depth_top", type=float, default=None,
                        help="Target depth of top of stress cube depth; allows calcuation of stress cube with deep nominal depth to avoid free surface effects on point-source stress field. Default is input nominal depth of stress cube from config.")
    parser.add_argument("--diam_mask_stress_center", type=float, default=-1.0,
                        help="Diamater in of mask of center of stress cube (km). < 0.0 disables masking.")
    parser.add_argument("--decay_mask_stress_center", type=float, default=-1.0,
                        help="Gaussian std_dev decay outside mask of center of stress cube (km).")
    parser.add_argument("--clip_stres_quantile", type=float, default=-1.0,
                        help="Quantile beyond which to clip +/- stress cube.")
    parser.add_argument("--plot_ffault_norm", type=float, default=-1.0,
                        help="Normalization value for plotting ffault cube points.")
    parser.add_argument("--plot_scale", nargs='+', type=float, default=[1.0,1.0],
                        help="Scale factor for output of cube points; second value, if present, applied to 50 percent output.")
    parser.add_argument('--no_output_neg', action='store_true', default=False,
                        help='Do not output (possibly large) negative stress fields.')
    parser.add_argument('--save_cubes', action='store_true', default=False,
                        help='Write pickled stress and finite faulting cubes to output.')
    parser.add_argument('--save_aftershock_cube', action='store_true', default=False,
                        help='Write pickled aftershock cube to output.')
    parser.add_argument("--aftershock_write_rate", type=float, default=0.01,
                        help="0.0-1.0 rate for output of synthetic aftershocks distributed following Coulomb stress cube.")
    parser.add_argument("--random_seed", type=int, default=12345,
                        help="Random seed for generating samples of cube values for plot output.")

    args = parser.parse_args()

    return args;


if __name__ == "__main__":  # CONFIGURE, INPUT, COMPUTE, OUTPUT

    args = welcome_and_parse_runstring()

    exp_config = PyCoulomb.configure_calc.configure_stress_calculation(args.config);
    if DEBUG:
        print("")
        print("DEBUG: exp_config:")
        for k, v in exp_config.__dict__.items():
            if k != "cube":
                print(k, v, " ", end="")

    [inputs, obs_disp_points, _] = PyCoulomb.input_values.read_inputs(exp_config);
    if DEBUG:
        print("")
        print("DEBUG: inputs:")
        for k, v in inputs.__dict__.items():
            if k != "cube":
                print(k, v, " ", end="")
        print("");
        print("")
    if DEBUG:
        print("DEBUG: obs_disp_points:", obs_disp_points)
    # out_object = PyCoulomb.run_dc3d.do_stress_computation(exp_config, inputs, obs_disp_points, ());
    # PyCoulomb.output_manager.produce_outputs(exp_config, inputs, obs_disp_points, (), out_object);

    receiver_horiz_profile_params = inputs.receiver_horiz_profile.__dict__
    if DEBUG:
        print("")
        print("DEBUG: receiver_horiz_profile:")
        for k, v in receiver_horiz_profile_params.items():
            if k != "cube":
                print(k, v, " ", end="")
        print("lon1d shape:", receiver_horiz_profile_params["lon1d"].shape)
        print("lat1d shape:", receiver_horiz_profile_params["lat1d"].shape)

    horiz_profile_path_root = os.path.split(exp_config.__dict__["outdir"][:-1])[0]
    horiz_profile_files = sorted(glob.glob(horiz_profile_path_root + "/*/stresses_horiz_profile.txt"))
    # if DEBUG:
    #	print("DEBUG: horiz_profile_files:", horiz_profile_files)

    out_path_root = horiz_profile_path_root
    if args.out_path != None:
        out_path_root = args.out_path
    os.makedirs(out_path_root, exist_ok=True)
    if DEBUG:
        print("DEBUG: out_path_root:", out_path_root)

    # construct stress and finite-faulting cubes

    # read each horizontal profile and populate 3D cube list of Coulomb stresses
    # # lon lat depth_km normal_kPa shear_kPa coulomb_kPa
    # # strike 6.000000, dip 28.000000, rake -93.000000
    # -150.564602 61.129730 45.000000 -0.140972 1.829841 1.773452
    # -150.562350 61.129730 45.000000 -0.138325 1.850203 1.794873
    # -150.560097 61.129730 45.000000 -0.135666 1.870689 1.816423
    # ...
    import csv

    measure_field = 5
    if args.stress_measure == "coulomb":
        measure_field = 5
    elif args.stress_measure == "shear":
        measure_field = 4
    elif args.stress_measure == "normal":
        measure_field = 3
    else:
        print("Error unrecognized stress_measure:", args.stress_measure)
        exit()

    num_x = receiver_horiz_profile_params["shape"][1]
    num_y = receiver_horiz_profile_params["shape"][0]
    num_z = len(horiz_profile_files)
    stress_cube_array = np.zeros([num_x, num_y, num_z], dtype=float)
    depths = []
    ndepth = 0;
    for horiz_profile_file in horiz_profile_files:
        print("Reading:", horiz_profile_file, end="")
        with open(horiz_profile_file, "r") as infile:
            reader = csv.reader(infile, delimiter=" ")
            next(reader, None)  # skip the headers
            next(reader, None)  # skip the headers
            nlon = -1
            nlat = 0
            for row in reader:
                nlon += 1
                if (nlon == num_x):
                    nlon = 0
                    nlat += 1
                if (nlat == num_y):
                    print("Error: reading horiz_profile_file: too many latitudes", nlat)
                # put stress value in list
                # print(nlon, "/", num_x, nlat, "/", num_y, ndepth, "/", num_z)
                stress_cube_array[nlon][nlat][ndepth] = float(row[measure_field])
        cdepth = os.path.basename(os.path.dirname(horiz_profile_file))
        print("   nlon:", nlon, "/", num_x, "  nlat:", nlat, "/", num_y, "  Depth:", cdepth, str(float(cdepth)), "\r",
              end="")
        depths.append(float(cdepth))
        ndepth += 1
    print("")
    stress_depth_shift = 0.0
    if args.stress_depth_top != None:
        stress_depth_shift = depths[0] - args.stress_depth_top
    stress_cube = {
        "cube": stress_cube_array,
        "centerlon": receiver_horiz_profile_params["centerlon"],
        # not absolute, since stress cube will be convolved over finite-faulting cube
        "centerlat": receiver_horiz_profile_params["centerlat"],
        # not absolute, since stress cube will be convolved over finite-faulting cube
        "centerdepth": (depths[0] + depths[-1]) / 2.0 - stress_depth_shift,
        # not absolute, since stress cube will be convolved over finite-faulting cube
        "ylen": receiver_horiz_profile_params["width"],
        "xlen": receiver_horiz_profile_params["length"],
        "xy_inc": receiver_horiz_profile_params["inc"],
        "zlen": (depths[-1] - depths[0]) / 2.0,
        "z_inc": depths[1] - depths[0]
    }
    if DEBUG:
        print("")
        print("stress_cube:", stress_cube["cube"].shape)
        for k, v in stress_cube.items():
            if k != "cube":
                print(k, v, " ", end="")
        print("");
        print("")

    # save FaultSlipObject_Geographic_Center if available
    fault_geographic_centerlon = None
    fault_geographic_centerlat = None
    fault_geographic_centerdepth = None
    with open(exp_config.__dict__["input_file"]) as ifile:
        for line in ifile:
            if line[0] == "#":
                continue;
            if "FaultSlipObject_Geographic_Center" in line:
                temp = line.split();
                fault_geographic_centerlon = float(temp[1])
                fault_geographic_centerlat = float(temp[2])
                fault_geographic_centerdepth = float(temp[3]) - stress_depth_shift

    masked_stress_cube = stress_cube

    # mask center (around source) of stress cube
    if args.diam_mask_stress_center > 0.0:
        masked_stress_cube = mask_cube_center(masked_stress_cube, args.diam_mask_stress_center,
                                              args.decay_mask_stress_center,
                                              fault_geographic_centerlon, fault_geographic_centerlat,
                                              fault_geographic_centerdepth)

    # mask cube binary (seems to spread ffault much too much)
    # masked_stress_cube = mask_cube_binary(masked_stress_cube, 0.999)

    # clip stress cube
    if args.clip_stres_quantile > 0.0:
        masked_stress_cube = mask_cube_quantile(masked_stress_cube, args.clip_stres_quantile)

    aftershocks = read_aftershock_table_NLL_csv_compact(args.aftershocks)

    # invert (correlate) for finite-fault response in cube from mapping of Coulomb stress cube over aftershock distribution
    ffault_cube, aftershock_cube, aftershocks_in_cube = inv_aftershocks_stress_ffault(masked_stress_cube,
                                                                                      aftershocks,
                                                                                      args.max_aftershock_count,
                                                                                      args.decay_mask_stress_center)

    # mask finite-fault cube based on aftershock density
    masked_ffault_cube = None
    ##mask_cells_extend = 1  # mask cells within this number of cells from cells that contain aftershocks
    ##min_num_ashock_count = 1  # mask around aftershock cells containing at least this number of aftershocks
    ##masked_ffault_cube = mask_cube_aftershocks(ffault_cube, aftershock_cube, mask_cells_extend, min_num_ashock_count)
    if not args.no_mask_to_aftershocks:
        masked_ffault_cube = mask_cube_aftershocks(ffault_cube, aftershocks_in_cube, args.no_mask_depth_aftershocks)

    # adjust center of finite-fault cube if FaultSlipObject_Geographic_Center available in input_file
    if fault_geographic_centerlon != None and fault_geographic_centerlat != None and fault_geographic_centerdepth != None:
        ffault_cube["centerlon"] = fault_geographic_centerlon
        ffault_cube["centerlat"] = fault_geographic_centerlat
        ffault_cube["centerdepth"] = fault_geographic_centerdepth

    # write output
    aftershock_write_rate = args.aftershock_write_rate
    #masked_stress_cube_pos = masked_stress_cube["cube"][masked_stress_cube["cube"] > 0.0]
    #aftershock_write_rate = len(aftershocks) / ((np.sum(masked_stress_cube_pos) / np.amax(masked_stress_cube_pos)))
    # print("DEBUG: aftershock_write_rate, len(aftershocks), np.sum(masked_stress_cube_pos), np.amax(masked_stress_cube_pos)", aftershock_write_rate, len(aftershocks), np.sum(masked_stress_cube_pos), np.amax(masked_stress_cube_pos))
    write_output(ffault_cube, masked_ffault_cube, args.plot_ffault_norm, stress_cube, masked_stress_cube,
                 aftershock_write_rate, out_path_root, args.no_output_neg, args.plot_scale, args.random_seed)
    if args.save_cubes:
        write_cubes_output(ffault_cube, "ffault_cube", out_path_root)
        write_cubes_output(masked_ffault_cube, "masked_ffault_cube", out_path_root)
        write_cubes_output(stress_cube, "stress_cube", out_path_root)
        write_cubes_output(masked_stress_cube, "masked_stress_cube", out_path_root)
    if args.save_aftershock_cube:
        write_cubes_output(aftershock_cube, "aftershock_cube", out_path_root)
