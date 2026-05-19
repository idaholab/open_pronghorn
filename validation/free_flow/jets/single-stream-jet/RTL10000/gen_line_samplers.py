# gen_line_samplers.py
L = 0.3
nz = 150
y_end = 0.1
variables = "vel_z"

dz = L / nz

with open("radial.i", "w") as f:
    for i in range(nz):
        z = (i + 0.5) * dz
        name = f"radial_{i:03d}"

        f.write(f"  [{name}]\n")
        f.write("    type = LineValueSampler\n")
        f.write(f"    start_point = '${{eps}} ${{eps}} {z:.6f}'\n")
        f.write(f"    end_point   = '${{eps}} {y_end:.6f} {z:.6f}'\n")
        f.write("    num_points  = ${ny}\n")
        f.write(f"    variable    = '{variables}'\n")
        f.write("    sort_by     = 'y'\n")
        f.write("    execute_on  = 'FINAL'\n")
        f.write("  []\n\n")
