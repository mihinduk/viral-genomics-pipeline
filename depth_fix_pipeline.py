import re

# Read the pipeline script
with open("scripts/run_viral_pipeline.py", "r") as f:
    content = f.read()

# Find the old depth command
old_depth_block = """            depth_cmd = f\"\"\"\/home\/mihindu\/miniforge3\/bin\/mamba run -n viral_genomics \\\\
                samtools depth {bam_file} > {output_dir}\/{args.output}_depth.txt\"\"\"
            run_command(depth_cmd, \"Generating depth file\")"""

# Create the new improved depth command with header and proper format
new_depth_block = """            # Generate depth file with proper header and 3-column format
            depth_cmd = f\"\"\"\/home\/mihindu\/miniforge3\/bin\/mamba run -n viral_genomics bash -c 
