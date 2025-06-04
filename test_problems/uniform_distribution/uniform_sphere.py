import numpy as np
import csv
import os

file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
dir_path = os.path.dirname(dir_path)
path = dir_path + "/generated_particles/"

def uniform_sphere_3d(
    n_particles,
    radius=1.0,
    mass=1.0,
    mass_std=None,
    velocity=0.0,
    velocity_std=None,
    filename="particles",
):
    """
    生成均勻分布於球體內的 3D 粒子資料，輸出 binary + CSV 檔

    每顆粒子資料格式為：
    [mass, x, y, z, vx, vy, vz]

    Parameters:
        n_particles (int): 粒子數
        radius (float): 球體半徑（粒子均勻分布於半徑為 radius 的球內）
        mass (float): 質量均值（固定或高斯）
        mass_std (float): 若給定則質量為 N(mass, mass_std)
        velocity (float or tuple): 速度均值 (vx, vy)
        velocity_std (float): 若給定則速度為 N(v, velocity_std)
        filename (str): 檔名（不含副檔名）
    """

    # 均勻球體內部隨機取點的方法（球坐標轉直角座標）
    u = np.random.uniform(0, 1, n_particles)
    costheta = np.random.uniform(-1, 1, n_particles)
    phi = np.random.uniform(0, 2 * np.pi, n_particles)
    theta = np.arccos(costheta)
    r = radius * u**(1/3)

    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)

    # 質量分布
    if mass_std is not None:
        m = np.random.normal(mass, mass_std, n_particles).astype(np.float64)
    else:
        m = np.full(n_particles, mass, dtype=np.float64)

    # 速度分布
    if isinstance(velocity, tuple):
        v_mean_x, v_mean_y = velocity
    else:
        v_mean_x = v_mean_y = velocity

    if velocity_std is not None:
        vx = np.random.normal(v_mean_x, velocity_std, n_particles).astype(np.float64)
        vy = np.random.normal(v_mean_y, velocity_std, n_particles).astype(np.float64)
    else:
        vx = np.full(n_particles, v_mean_x, dtype=np.float64)
        vy = np.full(n_particles, v_mean_y, dtype=np.float64)

    vz = np.zeros(n_particles, dtype=np.float64)

    # 合併成 (n, 7)
    particles = np.stack([m, x, y, z, vx, vy, vz], axis=1)

    # 儲存 binary 檔
    with open(filename + ".bin", "wb") as f:
        particles.astype(np.float64).tofile(f)

    # 儲存 CSV 檔
    with open(filename + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(particles)

    print(f"Saved {n_particles} particles to '{filename}.bin' and '{filename}.csv'.")


def uniform_disk_2d(
    n_particles,
    radius=1.0,
    mass=1.0,
    mass_std=None,
    velocity=0.0,
    velocity_std=None,
    filename="particles",
):
    """
    生成均勻分布的 2D 粒子資料（在圓盤內，以 3D 格式儲存），輸出 binary + CSV 檔

    每顆粒子資料格式為：
    [mass, x, y, z=0, vx, vy, vz=0]

    Parameters:
        n_particles (int): 粒子數
        radius (float): 圓盤半徑
        mass (float): 質量均值（固定或高斯）
        mass_std (float): 若給定則質量為 N(mass, mass_std)
        velocity (float or tuple): 速度均值 (vx, vy)
        velocity_std (float): 若給定則速度為 N(v, velocity_std)
        filename (str): 檔名
    """

    # 均勻分布在圓盤內的方法：r = R * sqrt(u)
    u = np.random.uniform(0, 1, n_particles)
    theta = np.random.uniform(0, 2 * np.pi, n_particles)
    r = radius * np.sqrt(u)

    x = r * np.cos(theta).astype(np.float32)
    y = r * np.sin(theta).astype(np.float32)
    z = np.zeros(n_particles, dtype=np.float32)

    # 質量分布
    if mass_std is not None:
        m = np.random.normal(mass, mass_std, n_particles).astype(np.float32)
    else:
        m = np.full(n_particles, mass, dtype=np.float32)

    # 速度分布
    if isinstance(velocity, tuple):
        v_mean_x, v_mean_y = velocity
    else:
        v_mean_x = v_mean_y = velocity

    if velocity_std is not None:
        vx = np.random.normal(v_mean_x, velocity_std, n_particles).astype(np.float32)
        vy = np.random.normal(v_mean_y, velocity_std, n_particles).astype(np.float32)
    else:
        vx = np.full(n_particles, v_mean_x, dtype=np.float32)
        vy = np.full(n_particles, v_mean_y, dtype=np.float32)

    vz = np.zeros(n_particles, dtype=np.float32)

    # 合併成 (n, 7)
    particles = np.stack([m, x, y, z, vx, vy, vz], axis=1)

    # 儲存 binary 檔
    with open(filename + ".bin", "wb") as f:
        particles.astype(np.float32).tofile(f)

    # 儲存 CSV 檔
    with open(filename + ".csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(particles)

    print(f"Saved {n_particles} particles to '{filename}.bin' and '{filename}.csv'.")

# 範例呼叫
if __name__ == "__main__":
    uniform_sphere_3d(
        n_particles=int(1e6),
        radius=50.0,
        mass=1.0,
        mass_std=0.1,
        velocity=(0.0, 0.0),
        velocity_std=0.5,
        filename=path+"uniform_sphere_3d_1e6"
    )
