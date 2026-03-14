import re

def convert_to_ra_dec(name):
    # Regular expression to capture RA and Dec from the name
    match = re.match(r'HR_(\d{2})(\d{2})-(\d{2})', name)
    
    if match:
        ra_hours = match.group(1)
        ra_minutes = match.group(2)
        dec_degrees = match.group(3)
        
        # Construct RA and Dec strings
        ra = f"{ra_hours}h{ra_minutes}m0.00s"
        dec = f"-{dec_degrees}d00'0.00\""
        return ra, dec
    else:
        return None

def convert_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            name = line.strip()
            result = convert_to_ra_dec(name)
            if result:
                ra, dec = result
                outfile.write(f"{name}  {ra} {dec} 2000.\n")

# Example usage
input_file = 'names.txt'
output_file = 'converted_ra_dec.txt'
convert_file(input_file, output_file)

print(f"Conversion completed! Output written to {output_file}.")

