import numpy as np
from numba import cuda
import time
import math # Import math for sqrt
import os
import pandas as pd

alpha = 1.25 # v_cut = alpha * v_esc
method = "fmm" # fmm, tree, naive

@cuda.jit
def potential_energy_kernel(masses, positions, potential_energies):
    """
    Numba CUDA 核心，用於計算每個粒子的重力位能。
    每個執行緒計算一個粒子 (i) 的位能貢獻。
    """

    i = cuda.grid(1)
    if i >= positions.shape[0]:
        return

    # G = 6.67430e-11  # 萬有引力常數
    G = 1.0
    
    local_energy = 0.0
    m_i = masses[i]
    pos_i = positions[i]

    # 計算粒子 i 與所有粒子 j (j > i) 之間的位能
    for j in range(i + 1, positions.shape[0]):
        m_j = masses[j]
        pos_j = positions[j]

        dx = pos_i[0] - pos_j[0]
        dy = pos_i[1] - pos_j[1]
        dz = pos_i[2] - pos_j[2]
        
        # 避免計算與自身的距離 (雖然迴圈已避免，但這是個好習慣)
        # 並加上一個很小的數 (epsilon) 來避免 r_sq == 0 的情況
        r_sq = dx*dx + dy*dy + dz*dz
        if r_sq < 1e-12: # epsilon^2
            r_sq = 1e-12
        inv_r = 1.0 / math.sqrt(r_sq)
        local_energy -= G * m_i * m_j * inv_r

    # 使用原子操作將局部能量加到總能量中
    cuda.atomic.add(potential_energies, 0, local_energy)

@cuda.jit
def kinetic_energy_kernel(masses, velocities, kinetic_energies, escape_num):
    """
    Numba CUDA 核心，用於計算每個粒子的動能。
    每個執行緒計算一個粒子 (i) 的動能貢獻。
    """
    i = cuda.grid(1)
    if i >= velocities.shape[0]:
        return
    
    local_escape = 0
    m_i = masses[i]
    vec_i = velocities[i]

    local_energy = 0.5 * m_i * (vec_i[0] * vec_i[0]
                              + vec_i[1] * vec_i[1]
                              + vec_i[2] * vec_i[2])
    
    # 去除太大的動能(逃逸動能=2e4)
    if local_energy > 2e4 * alpha:
        local_energy = 0
        local_escape = 1

    # 使用原子操作將局部能量加到總能量中
    cuda.atomic.add(kinetic_energies, 0, local_energy)
    cuda.atomic.add(escape_num, 0, local_escape)

@cuda.jit
def momentum_x_kernel(masses, velocities, momentum):
    """
    Numba CUDA 核心，用於計算每個粒子的動能。
    每個執行緒計算一個粒子 (i) 的動能貢獻。
    """
    i = cuda.grid(1)
    if i >= velocities.shape[0]:
        return

    m_i = masses[i]
    vec_i = velocities[i]

    local_momentum = m_i * vec_i[0]

    # 使用原子操作將局部能量加到總能量中
    cuda.atomic.add(momentum, 0, local_momentum)

@cuda.jit
def momentum_y_kernel(masses, velocities, momentum):
    """
    Numba CUDA 核心，用於計算每個粒子的動能。
    每個執行緒計算一個粒子 (i) 的動能貢獻。
    """
    i = cuda.grid(1)
    if i >= velocities.shape[0]:
        return

    m_i = masses[i]
    vec_i = velocities[i]

    local_momentum = m_i * vec_i[1]

    # 使用原子操作將局部能量加到總能量中
    cuda.atomic.add(momentum, 0, local_momentum)

@cuda.jit
def momentum_z_kernel(masses, velocities, momentum):
    """
    Numba CUDA 核心，用於計算每個粒子的動能。
    每個執行緒計算一個粒子 (i) 的動能貢獻。
    """
    i = cuda.grid(1)
    if i >= velocities.shape[0]:
        return

    m_i = masses[i]
    vec_i = velocities[i]

    local_momentum = m_i * vec_i[2]

    # 使用原子操作將局部能量加到總能量中
    cuda.atomic.add(momentum, 0, local_momentum)

def calculate_potential_energy_gpu(filename="uniform_sphere_3d_1e6.bin", n_particles=1_000_000):
    """
    使用 Numba 在 GPU 上計算總重力位能。
    """
    # 每個粒子的數據結構：[mass, x, y, z, vx, vy, vz] (7個 float64)
    
    # 讀取二進位檔案

    def load_particles_bin(filename, n_particles):
        with open(filename, "rb") as f:
            data = np.fromfile(f, dtype=np.float64).reshape(n_particles, 7)
        return data
    particles = load_particles_bin(filename, n_particles)

    num_particles = len(particles)
    print(f"成功讀取 {num_particles} 個粒子。")

    # 提取質量和位置
    masses = particles[:, 0]
    masses = np.ascontiguousarray(masses)
    positions = np.ascontiguousarray(np.vstack((particles[:, 1], particles[:, 2], particles[:, 3])).T.copy())
    velocities = np.ascontiguousarray(np.vstack((particles[:, 4], particles[:, 5], particles[:, 6])).T.copy())

    start_time = time.time()

    # 將數據傳輸到 GPU
    d_masses = cuda.to_device(masses)
    d_positions = cuda.to_device(positions)
    d_velocities = cuda.to_device(velocities)
    d_potential_energies = cuda.to_device(np.zeros(1, dtype=np.float64))
    d_kinetic_energies = cuda.to_device(np.zeros(1, dtype=np.float64))
    d_escape_num = cuda.to_device(np.zeros(1, dtype=np.float64))
    d_momentum_x = cuda.to_device(np.zeros(1, dtype=np.float64))
    d_momentum_y = cuda.to_device(np.zeros(1, dtype=np.float64))
    d_momentum_z = cuda.to_device(np.zeros(1, dtype=np.float64))

    # 設定 CUDA 核心的執行配置
    threads_per_block = 256
    blocks_per_grid = (num_particles + (threads_per_block - 1)) // threads_per_block

    # 執行 CUDA 核心
    potential_energy_kernel[blocks_per_grid, threads_per_block](
        d_masses, d_positions, d_potential_energies
    )
    kinetic_energy_kernel[blocks_per_grid, threads_per_block](
        d_masses, d_velocities, d_kinetic_energies, d_escape_num
    )
    momentum_x_kernel[blocks_per_grid, threads_per_block](
        d_masses, d_velocities, d_momentum_x
    )
    momentum_y_kernel[blocks_per_grid, threads_per_block](
        d_masses, d_velocities, d_momentum_y
    )
    momentum_z_kernel[blocks_per_grid, threads_per_block](
        d_masses, d_velocities, d_momentum_z
    )
    
    # 將結果從 GPU 傳回 CPU
    total_potential_energy = d_potential_energies.copy_to_host()
    total_kinetic_energy = d_kinetic_energies.copy_to_host()
    total_escape_num = d_escape_num.copy_to_host()
    total_momentum_x = d_momentum_x.copy_to_host()
    total_momentum_y = d_momentum_y.copy_to_host()
    total_momentum_z = d_momentum_z.copy_to_host()

    end_time = time.time()

    print(f"計算完成。")
    print(f"總重力位能: {total_potential_energy[0]:.6e} J")
    print(f"總動能: {total_kinetic_energy[0]:.6e} J，逃逸粒子數: {total_escape_num[0]:.6e}")
    print(f"總動量: {total_momentum_x[0]:.6e} {total_momentum_y[0]:.6e} {total_momentum_z[0]:.6e}")
    print(f"GPU 計算耗時: {end_time - start_time:.4f} 秒")

    return total_potential_energy[0], total_kinetic_energy[0], total_momentum_x[0], total_momentum_y[0], total_momentum_z[0], total_escape_num[0]

if __name__ == "__main__":

    summary = []

    alpha = 1.25


    method = 'fmm'
    for i in range(201):
        print("執行次數:", i, "方法", method)
        datafile = f'particles_{method}/plummer_vel_3d_1e5_step{i}_var.bin'
        U, K, Px, Py, Pz, Ne =calculate_potential_energy_gpu(filename=datafile, n_particles=100_000)
        summary.append({
            "i": i,
            "U": U,
            "K": K,
            "Px": Px,
            "Py": Py,
            "Pz": Pz,
            "Ne": Ne,
         })
    
    df_summary = pd.DataFrame(summary)
    df_summary.to_csv(f"conservation/conservation_results_alpha{alpha}_{method}_var.csv", index=False)
    print(df_summary)