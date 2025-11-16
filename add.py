import numpy as np
from scipy.optimize import fsolve

# --- CO2 Wagner 方程参数 ---
T_c = 304.1282  # 临界温度
P_c = 7.3773    # 临界压力

n1 = -7.0602087
n2 = 1.9391218
n3 = -1.6412611
n4 = 3.8966446
n5 = -2.7490843

# --- 核心计算函数 ---

def temperature_from_pressure(P_MPa):
    """
    主函数：根据Wagner方程，从压力(MPa)迭代求解饱和温度(K)
    """
    # 目标函数：找到T，使得 Wagner(T) - P = 0
    objective = lambda T: _wagner_pressure_from_temp(T) - P_MPa
    
    # 提供一个合理的初始猜测值，使用简化的Antoine方程估算
    T_guess_K = 1301.679 / (6.81228 - np.log10(P_MPa * 10)) - 3.494 + 273.15
    
    # 使用fsolve求解
    solution_K = fsolve(objective, T_guess_K)[0]
    return solution_K

# --- 使用示例 ---
if __name__ == "__main__":
    print("--- CO2 已知压力求温度计算器 ---")
    
    # 您之前的问题中的压力列表
    pressures_MPa = [4.612, 4.837, 5.021, 5.339, 5.532, 5.841, 6.093, 6.307, 6.548, 6.885, 7.136]
    
    print("压力(MPa)\t温度(°C)")
    print("-" * 30)
    
    for p in pressures_MPa:
        temp_K = temperature_from_pressure(p)
        temp_C = temp_K - 273.15
        print(f"{p:.3f}\t\t{temp_C:.2f}")

