import argparse
import pysam

def split_bam_by_tag(input_bam, tag, output_prefix):
  """
  Splits a BAM file based on the values of a specific tag.

  Args:
    input_bam: Path to the input BAM file (must be indexed).
    tag: Name of the tag to use for splitting.
    output_prefix: Prefix for the output BAM files.
  """
  # Open input and create a dictionary for output files
  samfile = pysam.AlignmentFile(input_bam, "rb")
  output_files = {}
  
  # Loop through reads in the BAM file
  for read in samfile:
    # Get the tag value
    tag_value = read.get_tag(tag)
    
    # Check if output file exists for the tag value
    if tag_value not in output_files:
      output_filename = f"{output_prefix}_{tag_value}.bam"
      output_files[tag_value] = pysam.AlignmentFile(output_filename, "wb", template=samfile)
    
    # Write the read to the corresponding output file
    output_files[tag_value].write(read)
  
  # Close input and output files
  samfile.close()
  for filename in output_files.values():
    filename.close()

if __name__ == "__main__":
  # Parse command-line arguments
  parser = argparse.ArgumentParser(description="Split BAM file based on a tag")
  parser.add_argument("input_bam", help="Path to the input BAM file")
  parser.add_argument("tag", help="Name of the tag to use for splitting")
  parser.add_argument("output_prefix", help="Prefix for the output BAM files")
  args = parser.parse_args()

  # Call the split function with parsed arguments
  split_bam_by_tag(args.input_bam, args.tag, args.output_prefix)

  print(f"Split BAM file completed. Output files prefixed with: {args.output_prefix}")
