# Description -------------------------------------------------------------
## This file contrains functions to quantify stat1 within islets

# Load libraries ----------------------------------------------------------
import os
from glob import glob
import sys
from pyhere import here
import gc

import numpy as np
import pandas as pd

from microfilm import microplot
from microfilm import microanim
import matplotlib
import matplotlib.pyplot as plt

import nd2
from nd2reader import ND2Reader

import skimage as ski
from scipy import ndimage as ndi

# Helpers -----------------------------------------------------------------
def clear_matplotlib_cache():
    """Clearmatplotlib caches to free memory"""
    plt.close('all')
    # Clear the figure manager
    matplotlib.pyplot.rcParams.clear()
    # Force garbage collection
    gc.collect()
  
# Load data ----------------------------------------------------------------
def load_image(path):
  """
    Load a multichannel ND2 microscopy image and extract specific fluorescence channels.

    Parameters
    ----------
    path : str or list of str
        Path to the ND2 image file. If a list is given, the first path in the
        list is used and the rest are ignored.

    Returns
    -------
    insulin : numpy.ndarray
        Image data for the FITCprc20x (insulin) channel.
    stat1 : numpy.ndarray
        Image data for the TRITCprv20x (STAT1) channel.
    dapi : numpy.ndarray
        Image data for the DAPIprv20x1 (DAPI/nuclear) channel.
"""
    if isinstance(path, list):
        path = path[0]

    img = nd2.imread(path, xarray=True)
    insulin = img.sel(C="FITCprc20x")
    stat1   = img.sel(C="TRITCprv20x")
    dapi    = img.sel(C="DAPIprv20x1")
    insulin = insulin.values
    stat1   = stat1.values
    dapi    = dapi.values

    return insulin, stat1, dapi
  
def get_pixel_size(path, default=1.0):
    """
    Retrieve the pixel size (in microns) from an ND2 image's metadata.

    Parameters
    ----------
    path : str
        Path to the ND2 file.
    default : float, optional
        Value to return if pixel size metadata is missing (default 1.0).

    Returns
    -------
    float
        Pixel size in microns per pixel.
    """

    img = ND2Reader(path)
    px_size = img.metadata.get("pixel_microns")
    img.close()

    if px_size is None:
        px_size = default

    return px_size

# Remove background ------------------------------------------------------
def iterative_average(img, n_iter=10, window=15, epsilon=0.5):
  """
  Create a baseline image for background correction using iterative smoothing.
  
  Parameters
  ----------
  img : numpy.ndarray
    Input image data.
  n_iter : int, optional
      Maximum number of smoothing iterations to perform (default 10).
  window : int, optional
      Size of the uniform filter window used for smoothing (default 15).
  epsilon : float, optional
      Convergence threshold: if the mean absolute pixel change between
      the current and previous baseline is at or below this value,
      iteration stops early (default 0.5).
      
      
    Returns
    -------
    numpy.ndarray
        Estimated baseline (background) image, same shape as `img`.
  """
    baseline = img.copy()
    
    for i in range(n_iter):
        baseline_prev = baseline.copy()
        
        smooth = ndi.uniform_filter(baseline, size=window)
        baseline = np.minimum(smooth, img)
        
        change = np.sum(np.abs(baseline - baseline_prev)) / baseline.size
        
        if change <= epsilon:
            print(f'Change is {epsilon}, iterations {i}')
            break
    
    return baseline

# Segmentation ----------------------------------------------------------
def segment_insulin(insulin):
  """
  Segment insulin-positive regions from a fluorescence image.

  Parameters
  ----------
  insulin : numpy.ndarray
      2D grayscale image of the insulin fluorescence channel.

  Returns
  -------
  mask : numpy.ndarray (bool)
      Binary mask of segmented insulin-positive regions after cleanup.
  mask_label : numpy.ndarray (int)
      Label image where each connected region in `mask` has a unique
      integer ID.
  th : float
      Otsu threshold value used to generate the initial binary mask.
  """

    th = ski.filters.threshold_otsu(insulin)

    mask = insulin > th
    mask = ski.morphology.closing(mask, ski.morphology.disk(10))
    mask = ndi.binary_fill_holes(mask)
    mask = ski.morphology.remove_small_objects(mask, max_size=10000)
    mask_label = ski.measure.label(mask)

    return mask, mask_label, th

def segment_dapi(dapi):
  """
  Segment individual nuclei from a DAPI fluorescence image.

  Parameters
  ----------
  dapi : numpy.ndarray
      2D grayscale image of the DAPI (nuclear) channel.

  Returns
  -------
  mask : numpy.ndarray (bool)
      Binary mask of thresholded DAPI-positive regions.
  labels : numpy.ndarray (int)
      Label image of individually segmented nuclei after watershed
      splitting and small-object removal.
  expanded : numpy.ndarray (int)
      Label image with each nucleus label expanded outward by 8 pixels.
  th : float
      Otsu threshold value used for the initial binary mask.
  """
    th = ski.filters.threshold_otsu(dapi)
    mask = dapi > th

    distance = ndi.distance_transform_edt(mask)

    coords = ski.feature.peak_local_max(distance, min_distance=7, labels=mask)

    markers = np.zeros_like(distance, dtype=int)
    markers[tuple(coords.T)] = np.arange(1, len(coords) + 1)

    labels = ski.segmentation.watershed(-distance, markers, mask=mask)

    labels = ski.morphology.remove_small_objects(labels, max_size=30)

    expanded = ski.segmentation.expand_labels(labels, distance=8)

    return mask, labels, expanded, th

# Plotting ----------------------------------------------------------
def background_removal(insulin, stat1, dapi, insulin_base, stat1_base, dapi_base, insulin_bg, stat1_bg, dapi_bg, pixel_size, save_path):
    # Raw
    p1r = microplot.microshow(images=insulin, cmaps = ['grey'], label_text='Insulin raw',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p1r.fig)
    p2r = microplot.microshow(images=stat1, cmaps = ['grey'], label_text='STAT1 raw',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p2r.fig)
    p3r = microplot.microshow(images=dapi, cmaps = ['grey'], label_text='DAPI raw',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p3r.fig)
    
    # Base
    p1b = microplot.microshow(images=insulin_base, cmaps = ['grey'], label_text='Insulin background',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p1b.fig)
    p2b = microplot.microshow(images=stat1_base, cmaps = ['grey'], label_text='STAT1 background',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p2b.fig)
    p3b = microplot.microshow(images=dapi_base, cmaps = ['grey'], label_text='DAPI background',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p3b.fig)
    
    # Corrected
    p1c = microplot.microshow(images=insulin_bg, cmaps = ['grey'], label_text='Insulin corrected',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p1c.fig)
    p2c = microplot.microshow(images=stat1_bg, cmaps = ['grey'], label_text='STAT1 corrected',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p2c.fig)
    p3c = microplot.microshow(images=dapi_bg, cmaps = ['grey'], label_text='DAPI corrected',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p3c.fig)
    
    micropanel = microplot.Micropanel(rows=4, cols=3, figscaling=3)
    micropanel.add_element(pos=[0,0], microim=p1r)
    micropanel.add_element(pos=[0,1], microim=p2r)
    micropanel.add_element(pos=[0,2], microim=p3r)
    
    micropanel.add_element(pos=[1,0], microim=p1b)
    micropanel.add_element(pos=[1,1], microim=p2b)
    micropanel.add_element(pos=[1,2], microim=p3b)
    
    micropanel.add_element(pos=[2,0], microim=p1c)
    micropanel.add_element(pos=[2,1], microim=p2c)
    micropanel.add_element(pos=[2,2], microim=p3c)
    
    # Choose the middle row to display intensities
    row_idx = int(insulin.shape[0]/2)
    for col, (raw, base, corr, label) in enumerate([
        (insulin, insulin_base, insulin_bg, 'Insulin'),
        (stat1, stat1_base, stat1_bg, 'STAT1'),
        (dapi, dapi_base, dapi_bg, 'DAPI')
    ]):
        ax = micropanel.ax[3, col]
        ax.plot(raw[row_idx, :], color="blue", label="Raw", linewidth = 1)
        ax.plot(base[row_idx, :], color="red", label="Background", linewidth = 1)
        ax.plot(corr[row_idx, :], color="lightblue", label="Corrected", linewidth = 1)
        ax.set_title(label)
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        ax.legend(loc="upper right")
    
    # Add horizontal line to corrected images
    for microim in [p1r, p2r, p3r, p1b, p2b, p3b, p1c, p2c, p3c]:
        microim.ax.axhline(y=row_idx, color='yellow', linewidth=1, linestyle='--', alpha=0.7)
        
    plt.close(micropanel.fig)
    micropanel.savefig(save_path, bbox_inches = 'tight', pad_inches = 0, dpi=100)
    clear_matplotlib_cache


def multi_plot(insulin_bg, stat1_bg, dapi_bg,  pixel_size, save_path):
    # Plot multi channel image
    p1 = microplot.microshow(images=insulin_bg, cmaps = ['yellow'], label_text='Insulin',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p1.fig)
    p2 = microplot.microshow(images=stat1_bg, cmaps = ['red'], label_text='STAT1',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p2.fig)
    p3 = microplot.microshow(images=dapi_bg, cmaps = ['blue'], label_text='DAPI',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p3.fig)
    p4 = microplot.microshow(images=[dapi_bg, insulin_bg, stat1_bg], cmaps = ['blue', 'yellow', 'red'], proj_type='sum', 
                             label_text='Combined', 
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit = 'um')
    plt.close(p4.fig)
    
    micropanel = microplot.Micropanel(rows=1, cols=4, figscaling=3)
    micropanel.add_element(pos=[0,0], microim=p1)
    micropanel.add_element(pos=[0,1], microim=p2)
    micropanel.add_element(pos=[0,2], microim=p3)
    micropanel.add_element(pos=[0,3], microim=p4)
    plt.close(micropanel.fig)
    micropanel.savefig(save_path, bbox_inches = 'tight', pad_inches = 0, dpi=300)
    clear_matplotlib_cache

def segmentation_plot(ins_mask, dapi_mask, ins_labels, dapi_labels, insulin_bg, dapi_bg, ins_th, dapi_th,  pixel_size, save_path):
    p1 = microplot.microshow(images=ins_mask, 
                         cmaps=['red'], 
                         label_text='Islet mask', 
                         scalebar_unit_per_pix=pixel_size, 
                         scalebar_size_in_units=100, unit='um', fig_scaling = 3)
    plt.close(p1.fig)
    
    p2 = microplot.microshow(images=dapi_mask, 
                             cmaps=['red'], 
                             label_text='Nuclei mask', 
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit='um', fig_scaling = 3)
    plt.close(p2.fig)
    
    p3 = microplot.microshow(images=ins_labels, 
                             cmaps=['grey'], 
                             label_text='Islet segmentation',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit='um', fig_scaling=3)
    
    plt.close(p3.fig)
    
    p4 = microplot.microshow(images=dapi_labels, 
                             cmaps=['grey'], 
                             label_text='Nuclei segmentation',
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit='um', fig_scaling=3)
    plt.close(p4.fig)
    
    micropanel = microplot.Micropanel(rows=4, cols=2, figscaling=3)  
    micropanel.add_element(pos=[0,0], microim=p1)
    micropanel.add_element(pos=[0,1], microim=p2)
    micropanel.add_element(pos=[1,0], microim=p3)
    micropanel.add_element(pos=[1,1], microim=p4)
    
    # Add histograms to row 1
    for col, (corr, th, label) in enumerate([
        (insulin_bg, ins_th, 'Insulin'),
        (dapi_bg, dapi_th, 'DAPI')
    ]):
        ax = micropanel.ax[2, col]
        ax.hist(corr.flatten(), bins=100, color = "red")
        ax.axvline(th, color='black', linewidth=1, label=f'Threshold: {th}')  
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        ax.legend(loc="upper right")
    
    for col, (corr, th, label) in enumerate([
        (insulin_bg, ins_th, 'Insulin'),
        (dapi_bg, dapi_th, 'DAPI')
    ]):
        ax = micropanel.ax[3, col]
        ax.hist(corr.flatten(), bins=100, color = "red")
        ax.axvline(th, color='black', linewidth=1, label=f'Threshold: {th}') 
        ax.set_yscale('log') 
        ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
        ax.legend(loc="upper right")
    plt.close(micropanel.fig)
    micropanel.savefig(save_path, bbox_inches = 'tight', pad_inches = 0, dpi=100)
    clear_matplotlib_cache

def measure_plot(dapi_bg, stat1_bg, ins_mask, nuclei_in_islet, nuclei_removed, cells_in_islet, outside_mask,  pixel_size, save_path):
    p1 = microplot.microshow(images=[dapi_bg, stat1_bg], 
                         cmaps=['blue', 'red'], 
                         label_text='Islet mask', 
                         scalebar_unit_per_pix=pixel_size, 
                         scalebar_size_in_units=100, unit='um')
    plt.close(p1.fig)
    
    p2 = microplot.microshow(images=[dapi_bg, stat1_bg], 
                             cmaps=['blue', 'red'], 
                             label_text='Nuclei in islets', 
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit='um')
    plt.close(p2.fig)
    
    p3 = microplot.microshow(images=[dapi_bg, stat1_bg], 
                             cmaps=['blue', 'red'], 
                             label_text='Perinuclear area', 
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit='um')
    plt.close(p3.fig)
    
    p4 = microplot.microshow(images=[dapi_bg, stat1_bg], 
                             cmaps=['blue', 'red'], 
                             label_text='Background area', 
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit='um')
    plt.close(p4.fig)
    
    p5 = microplot.microshow(images=[dapi_bg, stat1_bg], 
                             cmaps=['blue', 'red'], 
                             label_text='Combined', 
                             scalebar_unit_per_pix=pixel_size, 
                             scalebar_size_in_units=100, unit='um')
    plt.close(p5.fig)
    
    micropanel = microplot.Micropanel(rows=1, cols=5, figscaling=5)
    micropanel.add_element(pos=[0,0], microim=p1)
    micropanel.add_element(pos=[0,1], microim=p2)
    micropanel.add_element(pos=[0,2], microim=p3)
    micropanel.add_element(pos=[0,3], microim=p4)
    micropanel.add_element(pos=[0,4], microim=p5)
    
    # Add contours to the micropanel axes
    bg_overlay = np.zeros((*outside_mask.shape, 4))
    bg_overlay[outside_mask] = [0, 1, 0, 0.25]  # Green, semi-transparent
    
    micropanel.ax[0, 0].contour(ins_mask, colors='red', linewidths=1)
    micropanel.ax[0, 1].contour(ins_mask, colors='red', linewidths=1)
    micropanel.ax[0, 1].contour(nuclei_in_islet, colors='white', linewidths=1)
    micropanel.ax[0, 1].contour(nuclei_removed, colors='yellow', linewidths=1)
    micropanel.ax[0, 2].contour(cells_in_islet, colors='cornflowerblue', linewidths=1)
    micropanel.ax[0, 3].imshow(bg_overlay)
    micropanel.ax[0, 3].contour(ins_mask, colors='red', linewidths=1)
    micropanel.ax[0, 4].contour(cells_in_islet, colors='cornflowerblue', linewidths=1)
    micropanel.ax[0, 4].contour(ins_mask, colors='red', linewidths=1)
    micropanel.ax[0, 4].contour(nuclei_in_islet, colors='white', linewidths=1)
    plt.close(micropanel.fig)
    micropanel.savefig(save_path, bbox_inches = 'tight', pad_inches = 0, dpi=100)
    clear_matplotlib_cache


# Master_function ----------------------------------------------------------
def extract_features(path, folder_name, image_name):

    img_path = str(here(f'{path}/{folder_name}/{image_name}.nd2'))
    bg_path = os.path.join(str(here()), f'data/background/{folder_name}')
    seg_path = os.path.join(str(here()), f'data/segmentation/{folder_name}')
    comb_path = os.path.join(str(here()), f'data/combined/{folder_name}')
    measure_path = os.path.join(str(here()), f'data/measure/{folder_name}')
    islet_path = os.path.join(str(here()), f'data/measure_islet/{folder_name}')
    cell_path = os.path.join(str(here()), f'data/measure_cell/{folder_name}')
    file_path_islet = os.path.join(str(here()), f'data/files/{folder_name}')
    file_path_cell = os.path.join(str(here()), f'data/files/{folder_name}')

    # Create directories
    os.makedirs(bg_path, exist_ok=True) 
    os.makedirs(seg_path, exist_ok=True) 
    os.makedirs(comb_path, exist_ok=True) 
    os.makedirs(measure_path, exist_ok=True) 
    os.makedirs(islet_path, exist_ok=True)
    os.makedirs(cell_path, exist_ok=True)
    os.makedirs(file_path_islet, exist_ok=True)
    os.makedirs(file_path_cell, exist_ok=True)
    
    # Load image
    insulin, stat1, dapi = load_image(img_path)
    
    # Remove background
    insulin_base = iterative_average(insulin, n_iter=300, window=15, epsilon=0.05)
    stat1_base = iterative_average(stat1, n_iter=300, window=15, epsilon=0.05)
    dapi_base = iterative_average(dapi, n_iter=300, window=15, epsilon=0.05)
    insulin_bg = insulin - insulin_base
    stat1_bg = stat1 - stat1_base
    dapi_bg = dapi - dapi_base
    
    # Pixel size
    pixel_size = get_pixel_size(img_path)

    # Segmentation
    ins_mask, ins_labels, ins_th = segment_insulin(insulin_bg)
    dapi_mask, dapi_labels, expanded, dapi_th = segment_dapi(dapi_bg)

    # Label regions
    ins_props = ski.measure.regionprops(ins_labels)
    cell_props = ski.measure.regionprops(expanded)

    # Skip images with no islets
    if len(ins_props) == 0:
        print("Skipping: no islets detected")
        return None

    # only keep cells that overlap 80 % with the insulin positive area
    overlap_threshold = 0.8
    filtered_cell_labels = []
    
    for cell in cell_props:
        cell_mask = expanded == cell.label
        cell_mask_sum = cell_mask.sum()
        if cell_mask_sum == 0:
            continue
        overlap_fraction = (cell_mask & ins_mask).sum() / cell_mask_sum
        if overlap_fraction >= overlap_threshold:
            filtered_cell_labels.append(cell.label)
    
    # Create filtered masks
    filtered_cells_mask = np.isin(expanded, filtered_cell_labels)
    nuclei_in_islet = (dapi_labels > 0) & ins_mask & filtered_cells_mask
    cells_in_islet = filtered_cells_mask & ins_mask
    outside_mask = ins_mask & ~filtered_cells_mask
    
    # Nuclei that was removed
    nuclei_all = (dapi_labels > 0) & ins_mask
    nuclei_removed = nuclei_all & ~nuclei_in_islet
          
    # Plot background
    background_removal(insulin, stat1, dapi, insulin_base, stat1_base, dapi_base, insulin_bg, stat1_bg, dapi_bg,  pixel_size, os.path.join(bg_path, f'{image_name}.png'))

    # Multichannel plot
    multi_plot(insulin_bg, stat1_bg, dapi_bg,  pixel_size, os.path.join(comb_path, f'{image_name}.png'))

    # Plot segmentation
    segmentation_plot(ins_mask, dapi_mask, ins_labels, dapi_labels, insulin_bg, dapi_bg, ins_th, dapi_th,  pixel_size, os.path.join(seg_path, f'{image_name}.png'))

    # Multichannel and measurements
    measure_plot(dapi_bg, stat1_bg, ins_mask, nuclei_in_islet, nuclei_removed, cells_in_islet, outside_mask,  pixel_size, os.path.join(measure_path, f'{image_name}.png'))

    # Feature extraction
    islet_df = []
    cell_df = []
    
    for islet in ins_props:
        # Islet id
        islet_id = islet.label
        # Islet area
        ins_area = ins_labels == islet_id
        # Nuclei in islets
        nuclei = ins_area & nuclei_in_islet
        # Perinular in islets
        peri = ins_area & cells_in_islet & ~ nuclei
        # Define background in islets
        outside = ins_area & outside_mask
    
        # Plot areas we will measure in
        p1 = microplot.microshow(images=ins_area,
                                 cmaps=['red'],
                                 label_text=f'Islet: {islet_id}')
        plt.close(p1.fig)
        p2 = microplot.microshow(images=nuclei,
                                 cmaps=['red'],
                                 label_text='Nuclei in islet')
        plt.close(p2.fig)
        p3 = microplot.microshow(images=peri,
                                 cmaps=['red'],
                                 label_text='Perinuclear in islet')
        plt.close(p3.fig)
        p4 = microplot.microshow(images=outside,
                                 cmaps=['red'],
                                 label_text='Background in islet')
        plt.close(p4.fig)
    
        micropanel = microplot.Micropanel(rows=1, cols=4, figscaling=3)
        micropanel.add_element(pos=[0,0], microim=p1)
        micropanel.add_element(pos=[0,1], microim=p2)
        micropanel.add_element(pos=[0,2], microim=p3)
        micropanel.add_element(pos=[0,3], microim=p4)
        plt.close(micropanel.fig)
        micropanel.savefig(os.path.join(islet_path, f'{image_name}_{islet_id}.png'), bbox_inches = 'tight', pad_inches = 0, dpi=100)

        del p1, p2, p3, p4, micropanel
        gc.collect()
    
        # Pixel size squared
        pixel_area = pixel_size ** 2
        
        # Measurements per islet ----------------------------------
        islet_df.append({
        "slide_name": folder_name,
        "image_name": image_name,
        "islet_id": islet_id,
        "pixel_size": pixel_size,
        "pixel_area_um2": pixel_area,
        # Insulin
        # Bakground
        "insulin_bg_mean": insulin_bg[outside].mean(),
        "insulin_bg_std": insulin_bg[outside].std(),
        "insulin_bg_area_um2": outside.sum() * pixel_area,
        # Peri nuclear
        "insulin_peri_mean": insulin_bg[peri].mean(),
        "insulin_peri_std": insulin_bg[peri].std(),
        "insulin_peri_area_um2": peri.sum() * pixel_area,
        # Nuclei
        "insulin_nuclei_mean": insulin_bg[nuclei].mean(),
        "insulin_nuclei_std": insulin_bg[nuclei].std(),
        "insulin_nuclei_area_um2": nuclei.sum() * pixel_area,
        # Islet
        "insulin_mean": insulin_bg[ins_area].mean(),
        "insulin_std": insulin_bg[ins_area].std(),
        "insulin_area_um2": ins_area.sum() * pixel_area,
        # STAT1
        "stat1_bg_mean": stat1_bg[outside].mean(),
        "stat1_bg_std": stat1_bg[outside].std(),
        "stat1_bg_area_um2": outside.sum() * pixel_area,
        "stat1_peri_mean": stat1_bg[peri].mean(),
        "stat1_peri_std": stat1_bg[peri].std(),
        "stat1_peri_area_um2": peri.sum() * pixel_area,
        "stat1_nuclei_mean": stat1_bg[nuclei].mean(),
        "stat1_nuclei_std": stat1_bg[nuclei].std(),
        "stat1_nuclei_area_um2": nuclei.sum() * pixel_area,
        "stat1_mean": stat1_bg[ins_area].mean(),
        "stat1_std": stat1_bg[ins_area].std(),
        "stat1_area_um2": ins_area.sum() * pixel_area,
        # DAPI
        "dapi_bg_mean": dapi_bg[outside].mean(),
        "dapi_bg_std": dapi_bg[outside].std(),
        "dapi_bg_area_um2": outside.sum() * pixel_area,
        "dapi_peri_mean": dapi_bg[peri].mean(),
        "dapi_peri_std": dapi_bg[peri].std(),
        "dapi_peri_area_um2": peri.sum() * pixel_area,
        "dapi_nuclei_mean": dapi_bg[nuclei].mean(),
        "dapi_nuclei_std": dapi_bg[nuclei].std(),
        "dapi_nuclei_area_um2": nuclei.sum() * pixel_area,
        "dapi_mean": dapi_bg[ins_area].mean(),
        "dapi_std": dapi_bg[ins_area].std(),
        "dapi_area_um2": ins_area.sum() * pixel_area})
    
    
        # Measurements per cell ----------------------------------
        # For the first nucleus, make a plot
        # Plot areas we will measure in
        for cell_id in filtered_cell_labels[:1]:
            cell_mask = expanded == cell_id
            cell_in_islet = cell_mask & ins_area
            
            if cell_in_islet.sum() == 0:
                continue
            
            nucleus = cell_in_islet & (dapi_labels > 0)
            peri = cell_in_islet & ~nucleus
                
            p1 = microplot.microshow(images=[ins_area, nucleus],
                                     cmaps=['red', 'yellow'],
                                     label_text=f'Islet: {islet_id}')
            p1.ax.contour(peri, colors='blue', linewidths=1)
            plt.close(p1.fig)
            p1.savefig(os.path.join(cell_path, f'{image_name}.png') , bbox_inches = 'tight', pad_inches = 0, dpi=100)
            del p1
    
        for cell_id in filtered_cell_labels:
            cell_mask = expanded == cell_id
            cell_in_islet = cell_mask & ins_area
            
            if cell_in_islet.sum() == 0:
                continue
            
            nucleus = cell_in_islet & (dapi_labels > 0)
            peri = cell_in_islet & ~nucleus
            
            cell_df.append({
                "slide_name": folder_name,
                "image_name": image_name,
                "islet_id": islet_id,
                "cell_id": cell_id,
                # Insulin
                "insulin_nucleus_mean": insulin_bg[nucleus].mean() if nucleus.sum() > 0 else np.nan,
                "insulin_nucleus_std": insulin_bg[nucleus].std() if nucleus.sum() > 0 else np.nan,
                "insulin_peri_mean": insulin_bg[peri].mean() if peri.sum() > 0 else np.nan,
                "insulin_peri_std": insulin_bg[peri].std() if peri.sum() > 0 else np.nan,
                # STAT1
                "stat1_nucleus_mean": stat1_bg[nucleus].mean() if nucleus.sum() > 0 else np.nan,
                "stat1_nucleus_std": stat1_bg[nucleus].std() if nucleus.sum() > 0 else np.nan,
                "stat1_peri_mean": stat1_bg[peri].mean() if peri.sum() > 0 else np.nan,
                "stat1_peri_std": stat1_bg[peri].std() if peri.sum() > 0 else np.nan,
                # DAPI
                "dapi_nucleus_mean": dapi_bg[nucleus].mean() if nucleus.sum() > 0 else np.nan,
                "dapi_nucleus_std": dapi_bg[nucleus].std() if nucleus.sum() > 0 else np.nan,
                "dapi_peri_mean": dapi_bg[peri].mean() if peri.sum() > 0 else np.nan,
                "dapi_peri_std": dapi_bg[peri].std() if peri.sum() > 0 else np.nan,
            })

            del cell_mask, cell_in_islet, nucleus, peri
    
    clear_matplotlib_cache
        
    pd.DataFrame(islet_df).to_csv(os.path.join(file_path_islet, f'{image_name}_islet.csv'), index = False)
    pd.DataFrame(cell_df).to_csv(os.path.join(file_path_cell, f'{image_name}_cell.csv'), index = False)


def extract_features_spatial(path, folder_name, image_name):

    img_path = str(here(f'{path}/{folder_name}/{image_name}.nd2'))
    bg_path = os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/background/{folder_name}')
    spatial_path = os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/spatial/{folder_name}')
    seg_path = os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/segmentation/{folder_name}')
    comb_path = os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/combined/{folder_name}')
    measure_path = os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/measure/{folder_name}')
    islet_path = os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/measure_islet/{folder_name}')
    cell_path = os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/measure_cell/{folder_name}')
    file_path_islet = os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/files/{folder_name}')
    file_path_cell = os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/files/{folder_name}')

    # Create directories
    os.makedirs(os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/'), exist_ok=True)
    os.makedirs(bg_path, exist_ok=True) 
    os.makedirs(spatial_path, exist_ok=True) 
    os.makedirs(seg_path, exist_ok=True) 
    os.makedirs(comb_path, exist_ok=True) 
    os.makedirs(measure_path, exist_ok=True) 
    os.makedirs(islet_path, exist_ok=True)
    os.makedirs(cell_path, exist_ok=True)
    os.makedirs(file_path_islet, exist_ok=True)
    os.makedirs(file_path_cell, exist_ok=True)
    
    # Load image
    insulin, stat1, dapi = load_image(img_path)
    
    # Remove background
    insulin_base = iterative_average(insulin, n_iter=300, window=15, epsilon=0.05)
    stat1_base = iterative_average(stat1, n_iter=300, window=15, epsilon=0.05)
    dapi_base = iterative_average(dapi, n_iter=300, window=15, epsilon=0.05)
    insulin_bg = insulin - insulin_base
    stat1_bg = stat1 - stat1_base
    dapi_bg = dapi - dapi_base
    
    # Pixel size
    pixel_size = get_pixel_size(img_path)

    # Segmentation
    ins_mask, ins_labels, ins_th = segment_insulin(insulin_bg)
    dapi_mask, dapi_labels, expanded, dapi_th = segment_dapi(dapi_bg)

    # Label regions
    ins_props = ski.measure.regionprops(ins_labels)
    cell_props = ski.measure.regionprops(expanded)

    # Make dictoinary of cell proporties
    cell_props_dict = {prop.label: prop for prop in cell_props}
    
    # Skip images with no islets
    if len(ins_props) == 0:
        print("Skipping: no islets detected")
        return None

    # only keep cells that overlap 80 % with the insulin positive area
    overlap_threshold = 0.8
    filtered_cell_labels = []
    
    for cell in cell_props:
        cell_mask = expanded == cell.label
        cell_mask_sum = cell_mask.sum()
        if cell_mask_sum == 0:
            continue
        overlap_fraction = (cell_mask & ins_mask).sum() / cell_mask_sum
        if overlap_fraction >= overlap_threshold:
            filtered_cell_labels.append(cell.label)
    
    # Create filtered masks
    filtered_cells_mask = np.isin(expanded, filtered_cell_labels)
    nuclei_in_islet = (dapi_labels > 0) & ins_mask & filtered_cells_mask
    cells_in_islet = filtered_cells_mask & ins_mask
    outside_mask = ins_mask & ~filtered_cells_mask
    
    # Nuclei that was removed
    nuclei_all = (dapi_labels > 0) & ins_mask
    nuclei_removed = nuclei_all & ~nuclei_in_islet
          
    # Plot background
    background_removal(insulin, stat1, dapi, insulin_base, stat1_base, dapi_base, insulin_bg, stat1_bg, dapi_bg,  pixel_size, os.path.join(bg_path, f'{image_name}.png'))

    # Multichannel plot
    multi_plot(insulin_bg, stat1_bg, dapi_bg,  pixel_size, os.path.join(comb_path, f'{image_name}.png'))

    # Plot segmentation
    segmentation_plot(ins_mask, dapi_mask, ins_labels, dapi_labels, insulin_bg, dapi_bg, ins_th, dapi_th,  pixel_size, os.path.join(seg_path, f'{image_name}.png'))

    # Multichannel and measurements
    measure_plot(dapi_bg, stat1_bg, ins_mask, nuclei_in_islet, nuclei_removed, cells_in_islet, outside_mask,  pixel_size, os.path.join(measure_path, f'{image_name}.png'))

    # Feature extraction
    islet_df = []
    cell_df = []
    
    for islet in ins_props:
        # Islet id
        islet_id = islet.label
        # Islet area
        ins_area = ins_labels == islet_id
        # Nuclei in islets
        nuclei = ins_area & nuclei_in_islet
        # Perinular in islets
        peri = ins_area & cells_in_islet & ~ nuclei
        # Define background in islets
        outside = ins_area & outside_mask
        
        # Compute distance from islet edge once per islet
        islet_distance = ndi.distance_transform_edt(ins_area)

        p1 = microplot.microshow(images=islet_distance, cmaps = ['RdYlBu_r'], label_text='Distance from edge',
                                 scalebar_unit_per_pix=pixel_size, 
                                 show_colorbar= True,
                                 scalebar_size_in_units=100, unit = 'um')
        p1.ax.contour(ins_area, colors = "red")
        plt.close(p1.fig)
        p1.savefig(os.path.join(spatial_path, f'{image_name}_{islet_id}.png'), bbox_inches = 'tight', pad_inches = 0, dpi=100)
    
    
        # Plot areas we will measure in
        p1 = microplot.microshow(images=ins_area,
                                 cmaps=['red'],
                                 label_text=f'Islet: {islet_id}')
        plt.close(p1.fig)
        p2 = microplot.microshow(images=nuclei,
                                 cmaps=['red'],
                                 label_text='Nuclei in islet')
        plt.close(p2.fig)
        p3 = microplot.microshow(images=peri,
                                 cmaps=['red'],
                                 label_text='Perinuclear in islet')
        plt.close(p3.fig)
        p4 = microplot.microshow(images=outside,
                                 cmaps=['red'],
                                 label_text='Background in islet')
        plt.close(p4.fig)
    
        micropanel = microplot.Micropanel(rows=1, cols=4, figscaling=3)
        micropanel.add_element(pos=[0,0], microim=p1)
        micropanel.add_element(pos=[0,1], microim=p2)
        micropanel.add_element(pos=[0,2], microim=p3)
        micropanel.add_element(pos=[0,3], microim=p4)
        plt.close(micropanel.fig)
        micropanel.savefig(os.path.join(islet_path, f'{image_name}_{islet_id}.png'), bbox_inches = 'tight', pad_inches = 0, dpi=100)

        del p1, p2, p3, p4, micropanel
        gc.collect()
    
        # Pixel size squared
        pixel_area = pixel_size ** 2
        
        # Measurements per islet ----------------------------------
        islet_df.append({
        "slide_name": folder_name,
        "image_name": image_name,
        "islet_id": islet_id,
        "pixel_size": pixel_size,
        "pixel_area_um2": pixel_area,
        # Insulin
        # Bakground
        "insulin_bg_mean": insulin_bg[outside].mean(),
        "insulin_bg_median": np.median(insulin_bg[outside]),
        "insulin_bg_std": insulin_bg[outside].std(),
        "insulin_bg_area_um2": outside.sum() * pixel_area,
        # Peri nuclear
        "insulin_peri_mean": insulin_bg[peri].mean(),
        "insulin_peri_median": np.median(insulin_bg[peri]),
        "insulin_peri_std": insulin_bg[peri].std(),
        "insulin_peri_area_um2": peri.sum() * pixel_area,
        # Nuclei
        "insulin_nuclei_mean": insulin_bg[nuclei].mean(),
        "insulin_nuclei_median": np.median(insulin_bg[nuclei]),
        "insulin_nuclei_std": insulin_bg[nuclei].std(),
        "insulin_nuclei_area_um2": nuclei.sum() * pixel_area,
        # Islet
        "insulin_mean": insulin_bg[ins_area].mean(),
        "insulin_median": np.median(insulin_bg[ins_area]),
        "insulin_std": insulin_bg[ins_area].std(),
        "insulin_area_um2": ins_area.sum() * pixel_area,
        # STAT1
        "stat1_bg_mean": stat1_bg[outside].mean(),
        "stat1_bg_median": np.median(stat1_bg[outside]),
        "stat1_bg_std": stat1_bg[outside].std(),
        "stat1_bg_area_um2": outside.sum() * pixel_area,
        "stat1_peri_mean": stat1_bg[peri].mean(),
        "stat1_peri_median": np.median(stat1_bg[peri]),
        "stat1_peri_std": stat1_bg[peri].std(),
        "stat1_peri_area_um2": peri.sum() * pixel_area,
        "stat1_nuclei_mean": stat1_bg[nuclei].mean(),
        "stat1_nuclei_median": np.median(stat1_bg[nuclei]),
        "stat1_nuclei_std": stat1_bg[nuclei].std(),
        "stat1_nuclei_area_um2": nuclei.sum() * pixel_area,
        "stat1_mean": stat1_bg[ins_area].mean(),
        "stat1_median": np.median(stat1_bg[ins_area]),
        "stat1_std": stat1_bg[ins_area].std(),
        "stat1_area_um2": ins_area.sum() * pixel_area,
        # DAPI
        "dapi_bg_mean": dapi_bg[outside].mean(),
        "dapi_bg_median": np.median(dapi_bg[outside]),
        "dapi_bg_std": dapi_bg[outside].std(),
        "dapi_bg_area_um2": outside.sum() * pixel_area,
        "dapi_peri_mean": dapi_bg[peri].mean(),
        "dapi_peri_median": np.median(dapi_bg[peri]),
        "dapi_peri_std": dapi_bg[peri].std(),
        "dapi_peri_area_um2": peri.sum() * pixel_area,
        "dapi_nuclei_mean": dapi_bg[nuclei].mean(),
        "dapi_nuclei_median": np.median(dapi_bg[nuclei]),
        "dapi_nuclei_std": dapi_bg[nuclei].std(),
        "dapi_nuclei_area_um2": nuclei.sum() * pixel_area,
        "dapi_mean": dapi_bg[ins_area].mean(),
        "dapi_median": np.median(dapi_bg[ins_area]),
        "dapi_std": dapi_bg[ins_area].std(),
        "dapi_area_um2": ins_area.sum() * pixel_area})
    
    
        # Measurements per cell ----------------------------------
        # For the first nucleus, make a plot
        # Plot areas we will measure in
        for cell_id in filtered_cell_labels[:1]:
            cell_mask = expanded == cell_id
            cell_in_islet = cell_mask & ins_area
            
            if cell_in_islet.sum() == 0:
                continue
            
            nucleus = cell_in_islet & (dapi_labels > 0)
            peri = cell_in_islet & ~nucleus
                
            p1 = microplot.microshow(images=[ins_area, nucleus],
                                     cmaps=['red', 'yellow'],
                                     label_text=f'Islet: {islet_id}')
            p1.ax.contour(peri, colors='blue', linewidths=1)
            plt.close(p1.fig)
            p1.savefig(os.path.join(cell_path, f'{image_name}.png') , bbox_inches = 'tight', pad_inches = 0, dpi=100)
            del p1
    
        for cell_id in filtered_cell_labels:
            cell_mask = expanded == cell_id
            cell_in_islet = cell_mask & ins_area
            
            if cell_in_islet.sum() == 0:
                continue
            
            nucleus = cell_in_islet & (dapi_labels > 0)
            peri = cell_in_islet & ~nucleus

            # Get spatial features
            cell_prop = cell_props_dict[cell_id]
            centroid_y, centroid_x = cell_prop.centroid
            centroid_y_int = int(np.round(centroid_y))
            centroid_x_int = int(np.round(centroid_x))
            
            # Get distance from edge
            if (0 <= centroid_y_int < islet_distance.shape[0] and 
                0 <= centroid_x_int < islet_distance.shape[1]):
                distance_from_edge_px = islet_distance[centroid_y_int, centroid_x_int]
            else:
                distance_from_edge_px = np.nan

            cell_df.append({
                "slide_name": folder_name,
                "image_name": image_name,
                "islet_id": islet_id,
                "cell_id": cell_id,
                # SPATIAL FEATURES
                "centroid_x_px": centroid_x,
                "centroid_y_px": centroid_y,
                "centroid_x_um": centroid_x * pixel_size,
                "centroid_y_um": centroid_y * pixel_size,
                "distance_from_edge_px": distance_from_edge_px,
                "distance_from_edge_um": distance_from_edge_px * pixel_size,
                "cell_area_px": cell_prop.area,
                "cell_area_um2": cell_prop.area * pixel_area,
                # Insulin
                "insulin_nucleus_mean": insulin_bg[nucleus].mean() if nucleus.sum() > 0 else np.nan,
                "insulin_nucleus_median": np.median(insulin_bg[nucleus]) if nucleus.sum() > 0 else np.nan,
                "insulin_nucleus_std": insulin_bg[nucleus].std() if nucleus.sum() > 0 else np.nan,
                "insulin_peri_mean": insulin_bg[peri].mean() if peri.sum() > 0 else np.nan,
                "insulin_peri_median": np.median(insulin_bg[peri]) if peri.sum() > 0 else np.nan,
                "insulin_peri_std": insulin_bg[peri].std() if peri.sum() > 0 else np.nan,
                # STAT1
                "stat1_nucleus_mean": stat1_bg[nucleus].mean() if nucleus.sum() > 0 else np.nan,
                "stat1_nucleus_median": np.median(stat1_bg[nucleus]) if nucleus.sum() > 0 else np.nan,
                "stat1_nucleus_std": stat1_bg[nucleus].std() if nucleus.sum() > 0 else np.nan,
                "stat1_peri_mean": stat1_bg[peri].mean() if peri.sum() > 0 else np.nan,
                "stat1_peri_median": np.median(stat1_bg[peri]) if peri.sum() > 0 else np.nan,
                "stat1_peri_std": stat1_bg[peri].std() if peri.sum() > 0 else np.nan,
                # DAPI
                "dapi_nucleus_mean": dapi_bg[nucleus].mean() if nucleus.sum() > 0 else np.nan,
                "dapi_nucleus_median": np.median(dapi_bg[nucleus]) if nucleus.sum() > 0 else np.nan,
                "dapi_nucleus_std": dapi_bg[nucleus].std() if nucleus.sum() > 0 else np.nan,
                "dapi_peri_mean": dapi_bg[peri].mean() if peri.sum() > 0 else np.nan,
                "dapi_peri_median": np.median(dapi_bg[peri]) if peri.sum() > 0 else np.nan,
                "dapi_peri_std": dapi_bg[peri].std() if peri.sum() > 0 else np.nan
            })
            
            del cell_mask, cell_in_islet, nucleus, peri
    
    clear_matplotlib_cache
        
    pd.DataFrame(islet_df).to_csv(os.path.join(file_path_islet, f'{image_name}_islet.csv'), index = False)
    pd.DataFrame(cell_df).to_csv(os.path.join(file_path_cell, f'{image_name}_cell.csv'), index = False)

import os
import json
import numpy as np
from scipy import ndimage as ndi
import skimage as ski
from pyhere import here

def extract_insulin_masks(path, folder_name, image_name):
    """
    Extract insulin masks and convert to polygons for each islet.
    Saves as JSON for easy loading in R.
    """
    
    img_path = str(here(f'{path}/{folder_name}/{image_name}.nd2'))
    masks_path = os.path.join(str(here()), f'data/revisions/spatial_lfd_hfd/insulin_masks/{folder_name}')
    
    # Create directory
    os.makedirs(masks_path, exist_ok=True)
    
    # Load image
    insulin, _, _ = load_image(img_path)
    
    # Get pixel size
    pixel_size = get_pixel_size(img_path)
    
    # Remove background
    insulin_base = iterative_average(insulin, n_iter=300, window=15, epsilon=0.05)
    insulin_bg = insulin - insulin_base
    
    # Segment insulin
    ins_mask, ins_labels, _ = segment_insulin(insulin_bg)
    
    # Get islet properties
    ins_props = ski.measure.regionprops(ins_labels)
    
    # Convert masks to polygons
    masks_dict = {}
    
    for islet in ins_props:
        islet_id = int(islet.label)
        
        # Get binary mask for this islet
        islet_mask = (ins_labels == islet_id)
        
        # Find contours (polygon boundary)
        contours = ski.measure.find_contours(islet_mask, 0.5)
        
        if len(contours) > 0:
            # Use the longest contour (main boundary)
            contour = max(contours, key=len)
            
            # Convert to x, y coordinates (flip row/col to x/y)
            polygon = {
                'islet_id': islet_id,
                'x': contour[:, 1].tolist(),  # column = x
                'y': contour[:, 0].tolist(),  # row = y
                'area_px': int(islet.area),
                'area_um2': float(islet.area * (pixel_size ** 2))
            }
            
            masks_dict[islet_id] = polygon
    
    # Save as JSON (easy to load in R)
    save_file = os.path.join(masks_path, f'{image_name}_insulin_masks.json')
    
    with open(save_file, 'w') as f:
        json.dump(masks_dict, f)
    
    return len(masks_dict)
        
