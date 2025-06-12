import numpy as np
import pandas as pd
import subprocess
import os

exe_name = "./tree2fmm_parallel_high_acc"

def generate_input(datafile, outfile):
    with open("tree.in", "w") as f:
        f.write(f"""# G
1.0
# THETA
0.3
# MAX_PARTICLES_PER_CELL
10
# datafile
{datafile}
# outfile
{outfile}
# NUM_OF_THREADS
24
""")

def run_KDK_step(dt=0.01, infile='particles/plummer_vel_3d_1e5_step0.bin', outfile='particles/plummer_vel_3d_1e5_step1.bin'):
    """
    執行一步 KDK 軌道積分。
    """
    # --- 讀取資料 ---
    print("\n--- 開始 KDK 積分步驟 ---")

    # 1. 讀取粒子初始狀態 (mass, pos, vel)
    # 使用 np.fromfile 高效讀取二進位檔
    try:
        particle_data = np.fromfile(infile, dtype=np.float64)
        # 將一維陣列重塑為 (N, 7) 的形狀
        particles = particle_data.reshape(-1, 7)
    except FileNotFoundError:
        print(f"錯誤: {infile} 不存在。")
        return

    mass = particles[:, 0:1]  # 維持 (N, 1) 形狀以便廣播
    pos = particles[:, 1:4]   # 位置 (x, y, z)
    vel = particles[:, 4:7]   # 速度 (vx, vy, vz)
    print(f"成功讀取 {len(particles)} 顆粒子。")

    if os.path.exists('force/force_tree2fmm_K1.csv'):
        os.remove('force/force_tree2fmm_K1.csv')

    # 1.5. run ./tree2fmm_parallel
    generate_input(infile, 'force/force_tree2fmm_K1.csv')
    subprocess.run([exe_name])

    # 2. 讀取 tree2fmm_parallel 計算出的力，並計算加速度 a = F/m
    try:
        force_df = pd.read_csv('force/force_tree2fmm_K1.csv', header=None)
        force = force_df.to_numpy()
        # 加速度 a = F/m。NumPy的廣播機制會自動處理維度。
        accel = force / mass
    except FileNotFoundError:
        print("錯誤: force_tree2fmm_K1.csv 不存在。請確保 tree2fmm_parallel 已成功執行。")
        return
    print("成功讀取力，並計算出加速度。")


    # --- KDK (Kick-Drift-Kick) 演算法 ---
    # KDK 分為三步: Kick (半步) -> Drift (全步) -> Kick (半步)
    # v(t + dt/2) = v(t) + a(t) * dt/2
    # x(t + dt)   = x(t) + v(t + dt/2) * dt
    # v(t + dt)   = v(t + dt/2) + a(t + dt) * dt/2

    # 3. Kick (第一步): 使用目前的加速度 a(t) 更新速度半個時間步
    vel_half = vel + accel * (dt / 2.0)

    # 4. Drift: 使用更新後的速度 vel_half 更新位置一整個時間步
    pos_new = pos + vel_half * dt

    # 5. 更新粒子狀態並儲存，為下一次計算力做準備
    if os.path.exists('particles/half_step.bin'):
        os.remove('particles/half_step.bin')
    updated_particles = np.hstack((mass, pos_new, vel_half))
    updated_particles.tofile('particles/half_step.bin')

    # 5.5 run ./tree2fmm_parallel
    if os.path.exists('force/force_tree2fmm_K2.csv'):
        os.remove('force/force_tree2fmm_K2.csv')
    generate_input('particles/half_step.bin', 'force/force_tree2fmm_K2.csv')
    subprocess.run([exe_name])

    # 6. 讀取 tree2fmm_parallel 計算出的力，並計算加速度 a = F/m
    try:
        force_df = pd.read_csv('force/force_tree2fmm_K2.csv', header=None)
        force = force_df.to_numpy()
        # 加速度 a = F/m。NumPy的廣播機制會自動處理維度。
        accel = force / mass
    except FileNotFoundError:
        print("錯誤: force_tree2fmm_K2.csv 不存在。請確保 tree2fmm_parallel 已成功執行。")
        return
    print("成功讀取力，並計算出加速度。")

    # 7. Kick (第二步): 使用目前的加速度 a(t) 更新速度半個時間步
    vel_final = vel_half + accel * (dt / 2.0)

    # 8. 將更新後的粒子資料寫回檔案，以便進行下一個循環
    updated_particles = np.hstack((mass, pos_new, vel_final))
    updated_particles.tofile(outfile)   

    # 完整的模擬循環是：
    #   a. 計算 a(t) (透過 tree2fmm)
    #   b. 更新速度半步、位置一步 (如上)
    #   c. 用新的位置 pos_new 再次執行 tree2fmm 得到 a(t+dt)
    #   d. 用 a(t+dt) 完成速度的後半步更新
    
    print("\n--- 積分計算完成 ---")
    print(f"時間步長 (dt): {dt}")
    print("\n前3顆粒子的位置變化:")
    for i in range(3):
        print(f"  粒子 {i}:")
        print(f"  舊位置: {pos[i]}")
        print(f"  新位置: {pos_new[i]}")

if __name__ == '__main__':
    # 執行一次積分步驟，時間步長為 0.1
    for i in range(200):
        infile=f"particles/plummer_vel_3d_1e5_step{i}.bin"
        outfile=f"particles/plummer_vel_3d_1e5_step{i+1}.bin"
        run_KDK_step(dt=0.1, infile=infile, outfile=outfile)