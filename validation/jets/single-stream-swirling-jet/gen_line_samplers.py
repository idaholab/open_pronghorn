# gen_line_samplers.py
n = 101
z_start = 0.0
z_end   = 0.3
y_end   = 0.1
num_points = 101
variables = "mu_eff mu_t pressure TKE TKED vel_x vel_y vel_z yplus"

for i in range(n):
    z = z_start + i * (z_end - z_start) / (n - 1)
    name = f"radial_{i:03d}"
    print(f"  [{name}]")
    print("    type = LineValueSampler")
    print(f"    start_point = '0 0.0 {z:.6f}'")
    print(f"    end_point   = '0 {y_end:.6f} {z:.6f}'")
    print(f"    num_points  = {num_points}")
    print(f"    variable    = '{variables}'")
    print("    sort_by     = 'y'")
    print("    execute_on  = 'timestep_end'")
    print("  []\n")
