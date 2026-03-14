import re
import subprocess
import glob
from PIL import Image
import os

def extract_bounding_box(ps_file):
    # Read the file to find the BoundingBox at the end
    with open(ps_file, 'r') as file:
        lines = file.readlines()

    bounding_box = None
    for line in reversed(lines):
        match = re.match(r'%%BoundingBox: (\d+) (\d+) (\d+)', line)
        if match:
            bounding_box = line.strip()
            break

    if not bounding_box:
        raise ValueError("BoundingBox not found in the PostScript file.")
    return bounding_box

def create_modified_ps(ps_file, new_ps_file, bounding_box):
    # Create a new file with the BoundingBox at the beginning
    with open(ps_file, 'r') as file:
        lines = file.readlines()

    with open(new_ps_file, 'w') as file:
        for line in lines:
            if '%%BoundingBox:' in line:
                continue  # Skip the original BoundingBox line
            file.write(line)
        file.write(f'{bounding_box}\n')

def ps_to_png(ps_file, output_file):
    args = [
        "gs",
        "-dSAFER",
        "-dBATCH",
        "-dNOPAUSE",
        "-sDEVICE=png16m",
        "-r300",  # Set resolution to 300 DPI
        f"-sOutputFile={output_file}",
        ps_file
    ]
    
    result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    
    if result.returncode != 0:
        print("Error occurred:")
        print(result.stderr)
    else:
        print("Conversion successful")
        print(result.stdout)

def rotate_png(image_file, output_file, angle=90):
    with Image.open(image_file) as img:
        rotated_img = img.rotate(-angle, expand=True)
        rotated_img.save(output_file)
        print(f"Image rotated by {angle} degrees and saved as {output_file}")
# Process all .ps files in the current directory
ps_files = glob.glob("*.ps")

for ps_file in ps_files:
    base_name = os.path.splitext(ps_file)[0]
    modified_ps_file = f"{base_name}_modified.ps"
    output_png_file = f"{base_name}.png"
    rotated_png_file = f"{base_name}_rotated.png"
    
    try:
        bounding_box = extract_bounding_box(ps_file)
        create_modified_ps(ps_file, modified_ps_file, bounding_box)
        ps_to_png(modified_ps_file, output_png_file)
        rotate_png(output_png_file, rotated_png_file, angle=90)
        # Remove the original and intermediate files
        #os.remove(ps_file)
        os.remove(modified_ps_file)
        os.remove(output_png_file)
    except Exception as e:
        print(f"An error occurred with file {ps_file}: {e}")

