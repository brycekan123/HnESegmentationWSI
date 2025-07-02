

import os
import glob
from PIL import Image
import numpy as np

def stitchback(imagedir):
    
    
    os.makedirs("stitched", exist_ok=True)
    
    # Get all image files for each type
    image_types = ['normalized', 'H', 'E']
    
    for img_type in image_types:
        print(f"Stitching {img_type} images...")
        
        # Find all files of this type
        pattern = os.path.join(imagedir, f"{img_type}_*.png")
        files = glob.glob(pattern)
        
        if not files:
            print(f"No {img_type} files found!")
            continue
        
        # Extract coordinates from filenames and store file info
        tile_info = []
        for file in files:
            basename = os.path.basename(file)
            coords_part = basename.replace(f"{img_type}_", "").replace(".png", "")
            x, y = map(int, coords_part.split("_"))
            tile_info.append((x, y, file))
            
        #Does not affect Functionality but sorts x and y coordinates so stiching can happen in unform order.
        tile_info.sort(key=lambda item: (item[1], item[0]))
        # Get ALL unique x and y coordinates
        x_coords = sorted(list(set([item[0] for item in tile_info])))
        y_coords = sorted(list(set([item[1] for item in tile_info])))
        
        print(f"Found {len(tile_info)} tiles")
        print(f"X coordinates: {x_coords}")
        print(f"Y coordinates: {y_coords}")
        
        first_img = Image.open(tile_info[0][2])
        tile_width, tile_height = first_img.size
    
        final_width = len(x_coords) * tile_width
        final_height = len(y_coords) * tile_height
        
        print(f"Creating {img_type} image: {final_width}x{final_height}")
        #StitchedImage
        final_image = Image.new('RGB', (final_width, final_height))
        
        # Place each tile in correct position
        for x_coord, y_coord, filepath in tile_info:
            # Calculate position in final image onto the main tile
            x_index = x_coords.index(x_coord)
            y_index = y_coords.index(y_coord)
            
            paste_x = x_index * tile_width
            paste_y = y_index * tile_height
            
            #Pasting tile to stitched image!
            tile_img = Image.open(filepath)
            final_image.paste(tile_img, (paste_x, paste_y))
            
            print(f"Placed tile ({x_coord}, {y_coord}) at position ({paste_x}, {paste_y})")
        
        # Save stitched image
        output_path = os.path.join("stitched", f"{img_type}_stitched.png")
        final_image.save(output_path)
        print(f"Saved {output_path}")
        
    print("Stitching complete!")
    
    
stitchback("output_images")
