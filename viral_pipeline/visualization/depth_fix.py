# Fix for the depth visualization - add proper size constraints
import re

# Read the original file
with open("visualize_depth.py", "r") as f:
    content = f.read()

# Add size constraint function after imports
size_constraint_code = """
def constrain_figure_size(width, height, max_pixels=30000):
    \"\"\"Constrain figure size to prevent matplotlib errors\"\"\"
    # Maximum pixels per dimension (well below 2^16 = 65536)
    max_width = max_pixels
    max_height = max_pixels
    
    # Constrain dimensions
    if width > max_width:
        print(f"⚠️  Constraining width from {width} to {max_width}")
        width = max_width
    if height > max_height:
        print(f"⚠️  Constraining height from {height} to {max_height}")
        height = max_height
    
    # Ensure reasonable minimums
    width = max(4, min(width, max_width))
    height = max(3, min(height, max_height))
    
    return width, height

"""

# Find the imports section and add our function
import_end = content.find("warnings.filterwarnings")
import_end = content.find("\n", import_end) + 1

# Insert size constraint function
content = content[:import_end] + size_constraint_code + content[import_end:]

# Fix the figure creation line to use constraints
old_fig_line = "fig = plt.figure(figsize=figsize)"
new_fig_line = """# Apply size constraints to prevent matplotlib errors
    constrained_width, constrained_height = constrain_figure_size(figsize[0], figsize[1])
    print(f"📐 Figure size: {constrained_width}x{constrained_height}")
    fig = plt.figure(figsize=(constrained_width, constrained_height))"""

content = content.replace(old_fig_line, new_fig_line)

# Also reduce DPI to be safe
content = content.replace("dpi=300", "dpi=150")

# Write the fixed version
with open("visualize_depth_fixed.py", "w") as f:
    f.write(content)

print("✅ Created fixed depth visualization script")
