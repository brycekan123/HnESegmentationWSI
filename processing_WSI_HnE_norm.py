# https://youtu.be/tNfcvgPKgyU
"""
Here, we use openslide to read a whole slide image. 
We will then extract a lower reolution version of the image to normalize it
and then to extract H and E signals separately. 

We will also perform the exact operation on the entire whole slide image by 
extracting tilee, processing them, and saving processed images separately. 

Please note that this code will not cover putting tiles back into a 
whole slide image (image pyramid). You can explore pyvips or similar package
to put together tiles into an image pyramid. 

For an introduction to openslide, please watch video 266: https://youtu.be/QntLBvUZR5c
    
For details about H&E normalization, please watch my video 122: https://youtu.be/yUrwEYgZUsA
    
Useful references:
A method for normalizing histology slides for quantitative analysis. M. Macenko et al., ISBI 2009
http://wwwx.cs.unc.edu/~mn/sites/default/files/macenko2009.pdf

Efficient nucleus detector in histopathology images. J.P. Vink et al., J Microscopy, 2013

Other useful references:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5226799/
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169875

"""
import os
from openslide import open_slide
from PIL import Image
import numpy as np
from matplotlib import pyplot as plt
import tifffile as tiff
import math
import cv2

# Load the slide file (svs) into an object
slide = open_slide("images/image1.svs")
slide_props = slide.properties
slide_dimensions = slide.dimensions
print(slide_dimensions)

from normalize_HnE import norm_HnE

def count_he_pixels(norm_img, H_img, E_img, process, threshold=0.1):
    """
    Count pixels for pink, purple, background from H&E normalized images
    """
    if process == False:
        total_pixels = norm_img.shape[0] * norm_img.shape[1]
        return {
            'purple_pixels': 0,
            'pink_pixels': 0,  
            'background_pixels': total_pixels,
            'total_pixels': total_pixels
        }
    else:
        # Convert to grayscale for thresholding
        H_gray = cv2.cvtColor((H_img * 255).astype(np.uint8), cv2.COLOR_RGB2GRAY)
        E_gray = cv2.cvtColor((E_img * 255).astype(np.uint8), cv2.COLOR_RGB2GRAY)
        norm_gray = cv2.cvtColor((norm_img * 255).astype(np.uint8), cv2.COLOR_RGB2GRAY)
        
        # Create masks
        mask_H = H_gray > (threshold * 255)  # Purple/nuclei regions
        mask_E = E_gray > (threshold * 255)  # Pink/cytoplasm regions  
        mask_background = norm_gray < (threshold * 255)  # Background (white/light areas)
        
        # Count pixels
        total_pixels = norm_img.shape[0] * norm_img.shape[1]
        num_purple_pixels = cv2.countNonZero(mask_H.astype(np.uint8))
        num_pink_pixels = cv2.countNonZero(mask_E.astype(np.uint8)) - num_purple_pixels  # E minus overlap
        num_background_pixels = cv2.countNonZero(mask_background.astype(np.uint8))
        
        return {
            'purple_pixels': num_purple_pixels,
            'pink_pixels': max(0, num_pink_pixels),  # Ensure non-negative
            'background_pixels': num_background_pixels,
            'total_pixels': total_pixels
        }
#FOR SELECT PART OF SLIDE
slide_width, slide_height = slide_dimensions
tile_size = 1024

# Start from middle of slide
start_x = slide_width // 2  # Middle x coordinate
start_y = slide_height // 2  # Middle y coordinate

# Convert to tile coordinates
start_tile_x = start_x // tile_size
start_tile_y = start_y // tile_size

# Process only 50 tiles (7x8 grid) starting from middle
tiles_x = 7
tiles_y = 8

##FOR ENTIRE SLIDE
# slide_width, slide_height = slide_dimensions
# tile_size = 1024

# # Start from middle of slide
# start_x = 0
# start_y = 0

# # Convert to tile coordinates
# start_tile_x = start_x // tile_size
# start_tile_y = start_y // tile_size

# tiles_x = math.ceil(slide_width / tile_size)
# tiles_y = math.ceil(slide_height / tile_size)

# print(f"Slide dimensions: {slide_width} x {slide_height} pixels")
# print(f"Starting from middle: pixel ({start_x}, {start_y}) = tile ({start_tile_x}, {start_tile_y})")
# print(f"Processing {tiles_x}x{tiles_y} grid, stopping at 50 tiles")

# Create output directory if it doesn't exist
os.makedirs("output_images", exist_ok=True)

# Initialize pixel counters for entire slide
total_purple_pixels = 0
total_pink_pixels = 0 
total_background_pixels = 0
total_slide_pixels = 0  # Fixed variable name

# Process each tile across the entire slide
tile_count = 0
for tile_x in range(start_tile_x, start_tile_x + tiles_x):
    for tile_y in range(start_tile_y, start_tile_y + tiles_y):
        tile_count += 1
        
        # Calculate coordinates
        x = tile_x * tile_size
        y = tile_y * tile_size
        
        # Adjust tile size for edge tiles
        current_tile_width = min(tile_size, slide_width - x)
        current_tile_height = min(tile_size, slide_height - y)
        
        print(f"Processing tile {tile_count}/{tiles_x * tiles_y} at ({x}, {y})")

        # Extract tile - process ALL tiles (including blank ones)
        smaller_region = slide.read_region((x, y), 0, (current_tile_width, current_tile_height))
        smaller_region_RGB = smaller_region.convert('RGB')
        smaller_region_np = np.array(smaller_region_RGB)
        
        # Call norm_HnE function - this handles blank tiles internally
        norm_img, H_img, E_img, process = norm_HnE(smaller_region_np, Io=240, alpha=1, beta=0.15, process=True)
        
        # Count pixels for this tile
        pixel_counts = count_he_pixels(norm_img, H_img, E_img, process)
        
        # Add pixel counts to totals
        total_purple_pixels += pixel_counts['purple_pixels']
        total_pink_pixels += pixel_counts['pink_pixels']
        total_background_pixels += pixel_counts['background_pixels']
        total_slide_pixels += pixel_counts['total_pixels']  # Fixed variable name
        
        print(f"Tile {tile_count} at ({x}, {y}): {pixel_counts['purple_pixels']:,} purple pixels")
        print(f"Tile {tile_count} at ({x}, {y}): {pixel_counts['pink_pixels']:,} pink pixels ")
        print(f"Tile {tile_count} at ({x}, {y}): {pixel_counts['background_pixels']:,} background pixels")

        # Save the 3 images for each tile (needed for stitching)
        plt.imsave(f"output_images/normalized_{x}_{y}.png", norm_img)
        plt.imsave(f"output_images/H_{x}_{y}.png", H_img)
        plt.imsave(f"output_images/E_{x}_{y}.png", E_img)

# Calculate true percentages based on total pixel counts
true_purple_percent = (total_purple_pixels / total_slide_pixels) * 100
true_pink_percent = (total_pink_pixels / total_slide_pixels) * 100  
true_background_percent = (total_background_pixels / total_slide_pixels) * 100

summary_text = f"""TOTAL SLIDE ANALYSIS ({total_slide_pixels:,} pixels):
  Purple: {true_purple_percent:.1f}% ({total_purple_pixels:,} pixels)
  Pink: {true_pink_percent:.1f}% ({total_pink_pixels:,} pixels)
  Background: {true_background_percent:.1f}% ({total_background_pixels:,} pixels)
"""

with open("summary.txt", "w") as f:
    f.write(summary_text)

print(f"Analysis complete! Results saved to summary.txt")

