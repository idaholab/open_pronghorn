# gen_line_samplers.py
n = 151
z_start = 0.0
z_end   = 0.3
y_end   = 0.1
num_points = 41
variables = "mu_eff mu_t pressure TKE TKED vel_x vel_y vel_z yplus"

with open("radial.i", "w") as f:  # "w" = overwrite
    for i in range(n):
        z = z_start + i * (z_end - z_start) / (n - 1)
        name = f"radial_{i:03d}"

        f.write(f"  [{name}]\n")
        f.write("    type = LineValueSampler\n")
        f.write(f"    start_point = '0 0.0 {z:.6f}'\n")
        f.write(f"    end_point   = '0 {y_end:.6f} {z:.6f}'\n")
        f.write("    num_points  = ${ny}\n")  # <-- here
        f.write(f"    variable    = '{variables}'\n")
        f.write("    sort_by     = 'y'\n")
        f.write("    execute_on  = 'timestep_end'\n")
        f.write("  []\n\n")
