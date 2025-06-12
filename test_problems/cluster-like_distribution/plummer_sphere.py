import numpy as np
import csv
import os

file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
dir_path = os.path.dirname(dir_path)
path = dir_path + "/generated_particles/"

G = 1.0  # Gravitational constant

def plummer_sphere_3d(
    n_particles,
    scale_radius=1.0,
    total_mass=1.0,
    mass_std=None,
    filename="particles",
):
    """
    生成 Plummer sphere 分布的 3D 粒子資料，含自洽速度分布，輸出 binary + CSV 檔

    每顆粒子資料格式為：
    [mass, x, y, z, vx, vy, vz]

    Parameters:
        n_particles (int): 粒子數
        scale_radius (float): Plummer 尺度半徑 a
        total_mass (float): 總質量
        mass_std (float): 若給定則質量為 N(mass, mass_std)
        filename (str): 檔名（不含副檔名）
    """

    # ======= Sample positions ========
    def sample_radius(n):
        u = np.random.uniform(0.0, 1.0, n)
        return scale_radius / np.sqrt(u ** (-2.0 / 3.0) - 1.0)

    r = sample_radius(n_particles)
    theta = np.arccos(np.random.uniform(-1, 1, n_particles))
    phi = np.random.uniform(0, 2 * np.pi, n_particles)

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    # ======= Sample velocities ========
    def plummer_escape_velocity(r):
        return np.sqrt(2 * G * total_mass / np.sqrt(r**2 + scale_radius**2))

    def sample_velocity(r):
        v_esc = plummer_escape_velocity(r)
        v_ratio = np.sqrt(2.0/9.0)
        f_max =  v_ratio**2 * (1 - v_ratio**2)**(7.0/2.0) * v_esc**2
        while True:
            v = np.random.uniform(0, v_esc)
            f_v = v**2 * (1 - (v / v_esc) ** 2) ** 3.5
            if np.random.uniform(0, 1) < f_v / f_max:
                return v

    vx = np.zeros(n_particles)
    vy = np.zeros(n_particles)
    vz = np.zeros(n_particles)

    for i in range(n_particles):
        v = sample_velocity(r[i])

        # Random direction on a sphere
        costheta_v = np.random.uniform(-1, 1)
        sintheta_v = np.sqrt(1 - costheta_v**2)
        phi_v = np.random.uniform(0, 2 * np.pi)

        vx[i] = v * sintheta_v * np.cos(phi_v)
        vy[i] = v * sintheta_v * np.sin(phi_v)
        vz[i] = v * costheta_v

    # ======= Mass distribution ========
    if mass_std is not None:
        m = np.random.normal(total_mass / n_particles, mass_std, n_particles).astype(np.float64)
    else:
        m = np.full(n_particles, total_mass / n_particles, dtype=np.float64)

    # ======= Output files ========
    particles = np.stack([m, x, y, z, vx, vy, vz], axis=1)

    with open(filename + ".bin", "wb") as f:
        particles.astype(np.float64).tofile(f)

    with open(filename + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(particles)

    print(f"Saved {n_particles} particles to '{filename}.bin' and '{filename}.csv'.")

def double_plummer_sphere_3d(
    n_particles_total,
    scale_radius1=50.0,
    scale_radius2=50.0,
    total_mass1=None,
    total_mass2=None,
    offset_pos=(200.0, 0.0, 0.0),
    offset_vel=(-0.5, 0.0, 0.0),
    filename="particles_double",
):
    """
    生成 Double Plummer 分布的粒子群，合併輸出。

    offset_pos: 第二個 Plummer 球的相對位置
    offset_vel: 第二個 Plummer 球的相對速度
    """

    if total_mass1 is None:
        total_mass1 = n_particles_total / 2.0
    if total_mass2 is None:
        total_mass2 = n_particles_total / 2.0

    n1 = n_particles_total // 2
    n2 = n_particles_total - n1

    # 第一個 Plummer 球（在原點）
    plummer_sphere_3d(
        n_particles=n1,
        scale_radius=scale_radius1,
        total_mass=total_mass1,
        filename="plummer1_temp"
    )
    p1 = np.fromfile("plummer1_temp.bin", dtype=np.float64).reshape(-1, 7)
    p1[:, 1:4] -= np.array(offset_pos) / 2
    p1[:, 4:7] -= np.array(offset_vel) / 2

    # 第二個 Plummer 球（偏移位置和速度）
    plummer_sphere_3d(
        n_particles=n2,
        scale_radius=scale_radius2,
        total_mass=total_mass2,
        filename="plummer2_temp"
    )
    p2 = np.fromfile("plummer2_temp.bin", dtype=np.float64).reshape(-1, 7)
    p2[:, 1:4] += np.array(offset_pos) / 2
    p2[:, 4:7] += np.array(offset_vel) / 2

    # 合併
    particles = np.vstack([p1, p2])

    # 寫入檔案
    with open(filename + ".bin", "wb") as f:
        particles.astype(np.float64).tofile(f)

    with open(filename + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(particles)

    # 清理暫存檔
    os.remove("plummer1_temp.bin")
    os.remove("plummer2_temp.bin")
    os.remove("plummer1_temp.csv")
    os.remove("plummer2_temp.csv")

    print(f"Saved double Plummer model to '{filename}.bin' and '{filename}.csv'.")



# # 範例呼叫
# if __name__ == "__main__":
#     plummer_sphere_3d(
#     n_particles=int(1e5),
#     scale_radius=50.0,
#     total_mass=1.0*int(1e5),
#     mass_std=None,
#     filename=path+"plummer_vel_3d_1e5"
#     )

if __name__ == "__main__":
    double_plummer_sphere_3d(
        n_particles_total=int(2e5),
        scale_radius1=50.0,
        scale_radius2=50.0,
        offset_pos=(200.0, 0.0, 0.0),
        offset_vel=(-0.5, 0.0, 0.0),
        filename=path + "double_plummer_2e5"
    )