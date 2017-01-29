# -*- coding: utf-8 -*-

# Functions for analysing the distribution of mitochondria

import os
import errno
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import bioformats
import pyome as ome

from matplotlib import gridspec
from mito import plot
from scipy import ndimage
from skimage import filters
from skimage.color import label2rgb
from skimage.exposure import histogram
from skimage.feature import peak_local_max
from skimage.morphology import (closing, square, watershed,
                                remove_small_objects)
from skimage.segmentation import clear_border
from skimage.measure import regionprops


def segment_dapi(dapi):
    dapi_g3 = filters.gaussian(dapi, sigma=3)
    dapi_thresh = filters.threshold_otsu(dapi_g3)
    return ndimage.binary_fill_holes(closing(dapi_g3 > dapi_thresh, square(5)))


def segment_mito(mt):
    mt_g1 = filters.gaussian(mt, sigma=1)
    mt_thresh = filters.threshold_li(mt_g1)
    return closing(mt_g1 > mt_thresh, square(5))


def segment_gammatub(gt):
    gt_g2 = filters.gaussian(gt, sigma=2)
    gt_thresh = np.mean(gt_g2)
    gt_bin = gt_g2 > gt_thresh
    return remove_small_objects(gt_bin, 50)


def localise_mtoc(gt, mask, bbox):
    gt_crop = filters.gaussian(gt[bbox[0]:bbox[2],
                                  bbox[1]:bbox[3]] * mask,
                               sigma=2)
    y, x = list(peak_local_max(gt_crop, indices=True, num_peaks=1)[0])

    return x + bbox[1], y + bbox[0]


def split_nuclei(dapi, mask, min_distance=100, min_size=1000):
    dapi_dm = ndimage.distance_transform_edt(mask)
    dapi_peaks = peak_local_max(filters.gaussian(dapi_dm, sigma=2),
                                min_distance=min_distance,
                                exclude_border=False,
                                indices=False)
    dapi_markers = ndimage.label(dapi_peaks)[0]
    dapi_labels = remove_small_objects(
        watershed(-dapi_dm, dapi_markers, mask=mask),
        min_size)

    return dapi_labels, dapi_markers * (dapi_labels > 0), dapi_dm


def get_polar_angle_map(xc, yc, offset_x, offset_y, mask):
    height, width = mask.shape
    dpolar = np.zeros(mask.shape)
    for y in range(height):
        for x in range(width):
            dpolar[y, x] = np.arctan2(y-(yc-offset_y), x-(xc-offset_x))

    return np.ma.array(dpolar, mask=~mask)


def get_nucleus_distmap(img, centroid):
    xc, yc = centroid
    cent_bin = np.zeros_like(img)
    cent_bin[np.floor(xc), np.floor(yc)] = 1

    return ndimage.distance_transform_edt(~cent_bin)


def calculate_r90(centroid, bounding_box, cell_mt_mask):
    yc, xc = centroid
    miny, minx, _, _ = bounding_box
    axc = int(np.floor(xc) - minx)
    ayc = int(np.floor(yc) - miny)

    cent_bin = np.zeros_like(cell_mt_mask)
    cent_bin[ayc, axc] = 1

    distmap = ndimage.distance_transform_edt(~cent_bin) * cell_mt_mask
    hist, bins = histogram(distmap)
    cdf = np.cumsum(hist[1:].astype(np.float)/np.float(np.sum(hist[1:])))

    return bins[np.argmax(cdf > 0.9)]


# DataFrame functions
def cell_status_classifier(rsv, threshold):
    if rsv > threshold:
        return "Infected"
    else:
        return "Normal"


# Clean angle data
def frequency90(angles):
    bet = angles[np.logical_and(angles > -45, angles <= 45)]
    return float(bet.size)/float(angles.size)


def process_images(name, dapi, mito, min_nucleus_distance=100,
                   min_nucleus_size=1000, min_mt_size=10,
                   output_folder=None):
    # Segment and split nuclei using the DAPI channel
    dapi_bin = segment_dapi(dapi)
    dapi_labels, dapi_markers, dapi_dm = split_nuclei(
        dapi, dapi_bin,
        min_distance=min_nucleus_distance,
        min_size=min_nucleus_size)

    # Segment mitochondrial stain. Exclude anything that touches image border
    mt_bin = segment_mito(mito)
    mt_dapi_bin = ndimage.binary_fill_holes(mt_bin+dapi_bin)

    # Label mitochondria per nuclei
    mt_dapi_labels = clear_border(
        watershed(-dapi_dm, dapi_markers, mask=mt_dapi_bin))
    mt_labels = remove_small_objects(mt_dapi_labels * mt_bin, min_mt_size)

    # Filter nuclei where no mitochodria were segemented
    # (e.g. where mito touch border of image)
    included = [mtl.label for mtl in regionprops(mt_labels)
                if mtl.area > min_mt_size or mtl.label != 0]
    mt_mask = np.in1d(mt_dapi_labels, included).reshape(mt_dapi_labels.shape)
    dapi_filt = dapi_labels * mt_mask
    mt_dapi_labels = mt_dapi_labels * mt_mask

    if output_folder:
        if not os.path.exists(output_folder):
            try:
                os.makedirs(output_folder)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(output_folder):
                    pass
                else:
                    raise

        fig, (ax0, ax1) = plt.subplots(nrows=2)
        ax0.imshow(label2rgb(dapi_filt, image=dapi, bg_label=0))
        ax0.axis('off')
        ax1.imshow(label2rgb(mt_labels, image=mito, bg_label=0))
        ax1.axis('off')

        fig.set_size_inches(5, 9)
        fig.savefig(os.path.join(output_folder,
                                 "{!s}.pdf".format(name)),
                    dpi=150)
        plt.close(fig)

    return [(dapi_filt, mt_labels)]


def mitochondria_angles(name, treatment, time, rsv, tub, dapi_labels,
                        mt_labels, infection_threshold, output_folder=None):
    # Get properties for each segmented nucleus
    dapi_rp = regionprops(dapi_labels)

    # Get properties for mitochondria surrounding each segmented nucleus.
    mt_rp = regionprops(mt_labels, intensity_image=rsv)

    mtocs = [localise_mtoc(tub, rp.image, rp.bbox) for rp in mt_rp]

    if len(mtocs) == 0:
        print("No MTOCs detected in {!s}".format(name))
        return
    elif len(mtocs) > len(dapi_rp):
        print("More than one MTOC (i.e., gamma-tubulin)"
              " per nucleus detected in {!s}".format(name))
        return

    irsv1 = [mt.mean_intensity for mt in mt_rp]

    # Calculate MTOC angles for each nuclei
    mtoc_angs = [np.arctan2(y-dp.centroid[0], x-dp.centroid[1])
                 for dp, (x, y) in zip(dapi_rp, mtocs)]

    # Generate Polar angle maps for mitochondria around each nuclei
    polar_maps = [get_polar_angle_map(dp.centroid[1], dp.centroid[0],
                                      mp.bbox[1], mp.bbox[0], mp.image)
                  for dp, mp in zip(dapi_rp, mt_rp)]

    def normalise_angle(pm, mtoc):
        npm = pm - mtoc
        if npm > np.pi:
            npm = npm - 2*np.pi
        elif npm < -np.pi:
            npm = npm + 2*np.pi

        return npm

    vnorm_angle = np.vectorize(normalise_angle)

    polar_maps = [vnorm_angle(pm, mtoca)
                  for pm, mtoca in zip(polar_maps, mtoc_angs)]

    # Function to create a dataframe from mtoc and mitoc angles.
    def mtoc_mt_df(n, mtoca, mt_map, irsv1):
        mt = mt_map.compressed()
        return pd.DataFrame({"Cell ID": np.repeat(n, len(mt)),
                             "MTOC angle": np.repeat(mtoca, len(mt)),
                             "RSV1 intensity": np.repeat(irsv1, len(mt)),
                             "Mitochondria angle": mt})

    results = pd.concat([mtoc_mt_df(i, mtoca, mt_map, irsv)
                         for i, (mtoca, mt_map, irsv)
                         in enumerate(zip(mtoc_angs, polar_maps, irsv1))])

    results["Name"] = name
    results["Pathogen"] = treatment
    results["Treatment"] = time
    results["Infected"] = results["RSV1 intensity"] > infection_threshold

    if output_folder:
        if not os.path.exists(output_folder):
            try:
                os.makedirs(output_folder)
            except OSError as exc:
                if exc.errno == errno.EEXIST and os.path.isdir(output_folder):
                    pass
                else:
                    raise

        fig = plt.figure(figsize=(8, 6))
        gs = gridspec.GridSpec(len(polar_maps), 2, width_ratios=[3, 1])
        ax = plt.subplot(gs[:, 0])
        ax.imshow(tub)
        ax.scatter([x for x, y in mtocs], [y for x, y in mtocs], marker=u'+')
        ax.axis('off')
        if len(polar_maps) > 1:
            for n, pm, mtoc in zip(range(len(polar_maps)), polar_maps, mtocs):
                axn = plt.subplot(gs[n, 1])
                axn.imshow(pm, cmap="RdBu")
                axn.axis('off')
        else:
            axs = plt.subplot(gs[0, 1])
            axs.imshow(polar_maps[0], cmap="RdBu")
            axs.axis('off')
        fig.savefig(os.path.join(output_folder, "{!s}_angle.pdf".format(name)),
                    dpi=150)
        plt.close(fig)

    return results


def mitochondria_params(name, treatment, time, dapi, mito, rsv, dapi_labels,
                        mt_labels, meta, infection_threshold,
                        output_folder=None):
    """ Return mitochondrial parameters

    Parameters
    ----------
    dapi: ndarray
        DAPI channel 2 x 2 array.
    rsv1: ndarray
        RSV1 channel 2 x 2 array.
    mt: ndarray
        MitoTracker channel 2 x 2 array.

    Returns
    -------
    data: DataFrame
        Dataframe containing mitochondrial params per cells
    """
    # Get properties for each segmented nucleus
    dapi_rp = regionprops(dapi_labels)

    # Get properties for mitochondria surrounding each segmented nucleus.
    mt_rp = regionprops(mt_labels, intensity_image=rsv)

    # calibration = np.sqrt(meta.voxel_size_x**2 + meta.voxel_size_y**2)
    calibration = meta.voxel_size_x

    area_r90 = pd.DataFrame(
        np.array([[rp.area * (meta.voxel_size_x * meta.voxel_size_y),
                   rp.mean_intensity]
                  for rp in mt_rp]), columns=['Area', 'RSV1 Intensity'])
    centroids = [rp.centroid for rp in dapi_rp]
    area_r90["yc"], area_r90["xc"] = zip(*centroids)
    area_r90["R90%"] = np.array(
        [calculate_r90(dp.centroid, mp.bbox, mp.image) * calibration
         for dp, mp in zip(dapi_rp, mt_rp)])
    area_r90["Name"] = name
    area_r90["Pathogen"] = treatment
    area_r90["Treatment"] = time
    area_r90["Infected"] = area_r90["RSV1 Intensity"] > infection_threshold

    if output_folder:
        output_path = os.path.join(output_folder, "{!s}_r90.pdf".format(name))

        plot.output_r90_img(np.array([dapi, mito, rsv]), area_r90,
                            calibration, output_path, hues=[0.6, 0, 0.4])

    return area_r90


def filter_data(name, dapi_labels, mt_labels):
    # Get properties for each segmented nucleus
    dapi_rp = regionprops(dapi_labels)
    if len(dapi_rp) == 0:
        print("Unable to segment any entire cells in {!s}".format(name))
        return False

    # Get properties for mitochondria surrounding each segmented nucleus.
    mt_rp = regionprops(mt_labels)
    if len(mt_rp) == 0:
        print("Mitochondria not successfully segmented in {!s}".format(name))
        return False
    elif len(mt_rp) != len(dapi_rp):
        print("Nuclei labels and mito labels are not equal lengths "
              "in {!s}".format(name))
        return False

    return True


def read_channels(img_path, series):
    dapi = bioformats.load_image(img_path, c=0, series=series, rescale=False)
    mt = bioformats.load_image(img_path, c=1, series=series, rescale=False)
    rsv1 = bioformats.load_image(img_path, c=2, series=series, rescale=False)
    tub = bioformats.load_image(img_path, c=3, series=series, rescale=False)

    return dapi, mt, rsv1, tub


def read_series(img_path, img_meta):
    with bioformats.ImageReader(img_path) as rdr:
        for i, m in enumerate(img_meta):
            img = np.array([rdr.read(c=c, series=i, rescale=False)
                            for c in range(m.sizec)])
            yield img, m


def process_file(img_path, treatment, time, output_dir, dapi_index=0,
                 mito_index=1, rsv_index=2, tub_index=3,
                 min_nucleus_distance=10, min_nucleus_area=100,
                 min_mito_area=10, infection_threshold=5):
    img_meta = ome.read(img_path)
    img_name = os.path.basename(img_path)[:-4]

    # Get mito area, r90 and angles
    return ((mitochondria_params("{0}_{1}".format(img_name, meta.name),
                                 treatment, time, img[dapi_index],
                                 img[mito_index], img[rsv_index],
                                 dapi_labels, mt_labels, meta,
                                 infection_threshold,
                                 output_folder=output_dir),
             mitochondria_angles("{0}_{1}".format(img_name, meta.name),
                                 treatment, time, img[rsv_index],
                                 img[tub_index], dapi_labels,
                                 mt_labels, infection_threshold,
                                 output_folder=output_dir))
            for img, meta in read_series(img_path, img_meta)
            for dapi_labels, mt_labels in
            process_images("{0}_{1}".format(img_name, meta.name),
                           img[dapi_index], img[mito_index],
                           min_nucleus_distance/np.sqrt(meta.voxel_size_x**2 + meta.voxel_size_y**2),
                           min_nucleus_area/meta.voxel_size_x * meta.voxel_size_y,
                           min_mito_area/meta.voxel_size_x * meta.voxel_size_y,
                           output_folder=output_dir)
            if filter_data("{0}_{1}".format(img_name, meta.name),
                           dapi_labels, mt_labels))


def process_r90(img_path, treatment, time, output_dir, dapi_index=0,
                mito_index=1, rsv_index=2, min_nucleus_distance=10,
                min_nucleus_area=100, min_mito_area=10, infection_threshold=5):
    img_meta = ome.read(img_path)
    img_name = os.path.basename(img_path)[:-4]

    # Get mito area, r90 and angles
    return (mitochondria_params("{0}_{1}".format(img_name, meta.name),
                                treatment, time, img[dapi_index],
                                img[mito_index], img[rsv_index],
                                dapi_labels, mt_labels, meta,
                                infection_threshold, output_folder=output_dir)
            for img, meta in read_series(img_path, img_meta)
            for dapi_labels, mt_labels in
            process_images("{0}_{1}".format(img_name, meta.name),
                           img[dapi_index], img[mito_index],
                           min_nucleus_distance/np.sqrt(meta.voxel_size_x**2 + meta.voxel_size_y**2),
                           min_nucleus_area/meta.voxel_size_x * meta.voxel_size_y,
                           min_mito_area/meta.voxel_size_x * meta.voxel_size_y,
                           output_folder=output_dir)
            if filter_data("{0}_{1}".format(img_name, meta.name),
                           dapi_labels, mt_labels))
