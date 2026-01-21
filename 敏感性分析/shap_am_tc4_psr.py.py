import numpy as np
import pandas as pd
import shap
import matplotlib.pyplot as plt

# ================== 1. 读取参数文件 ==================
filename = 'AM-TC4-GRO.xlsx'
df = pd.read_excel(filename)
param_names = df.columns.tolist()
X_full = df.values    # N × 4
d = X_full.shape[1]

# ================== 2. PSR 裂纹扩展模型（log10 输出） ==================
def crack_model_PSR_log(X, deltaK=6.0, R_ratio=0.8):
    log_theta1 = X[:, 0]
    theta2 = X[:, 1]
    theta3 = X[:, 2]
    k2 = X[:, 3]

    theta1 = 10 ** log_theta1
    Kmax = deltaK / (1 - R_ratio)

    bracket_term = np.maximum(Kmax / k2 - 1.0, 0.0)
    da_dN = theta1 * (deltaK ** theta2) * (bracket_term ** theta3)

    return np.log10(da_dN + 1e-12)

# ================== 3. 样本数量收敛性设置 ==================
sample_list = np.round(np.linspace(500, 10000, 10)).astype(int)
shap_convergence = np.zeros((d, len(sample_list)))

# ================== 4. SHAP 收敛性分析 ==================
for k, N in enumerate(sample_list):
    print(f'Computing SHAP for N = {N}')
    X = X_full[:N, :]

    explainer = shap.KernelExplainer(
        lambda x: crack_model_PSR_log(x),
        X
    )

    shap_values = explainer.shap_values(
        X,
        nsamples=1000,
        l1_reg="num_features(4)"
    )

    # Mean |SHAP|
    shap_convergence[:, k] = np.mean(np.abs(shap_values), axis=0)

# ================== 5. 归一化 SHAP（便于与 Sobol 对比） ==================
shap_convergence_norm = shap_convergence / np.sum(shap_convergence, axis=0)

# ================== 6. SHAP 收敛性可视化 ==================
colors = plt.cm.tab10.colors
markers = ['o', 's', '^', 'd']

plt.figure(figsize=(11, 4))

# --- 子图 1：Mean |SHAP| ---
plt.subplot(1, 2, 1)
for i in range(d):
    plt.plot(sample_list, shap_convergence[i, :],
             '-', color=colors[i], linewidth=1.5)
    plt.scatter(sample_list, shap_convergence[i, :],
                color=colors[i], marker=markers[i], s=40,
                label=param_names[i])

plt.xlabel('Sample Number')
plt.ylabel('Mean |SHAP|')
plt.title('SHAP Importance Convergence')
plt.grid(True)
plt.legend()

# --- 子图 2：归一化 SHAP ---
plt.subplot(1, 2, 2)
for i in range(d):
    plt.plot(sample_list, shap_convergence_norm[i, :],
             '-', color=colors[i], linewidth=1.5)
    plt.scatter(sample_list, shap_convergence_norm[i, :],
                color=colors[i], marker=markers[i], s=40)

plt.xlabel('Sample Number')
plt.ylabel('Normalized Mean |SHAP|')
plt.title('Normalized SHAP Convergence')
plt.grid(True)

plt.suptitle('SHAP Convergence Analysis for PSR Crack Growth Model')
plt.tight_layout()
plt.show()

# ================== 7. 最大样本数下 SHAP 结果表 ==================
final_shap = shap_convergence[:, -1]
final_norm = final_shap / np.sum(final_shap)

shap_table = pd.DataFrame({
    'Parameter': param_names,
    'MeanAbsSHAP': final_shap,
    'NormalizedSHAP': final_norm
})

print('\n--- Final SHAP Importance (Max Sample Size) ---')
print(shap_table)

# ================== 8. SHAP Summary Plot（最大样本） ==================
X_max = X_full[:sample_list[-1], :]
explainer = shap.KernelExplainer(
    lambda x: crack_model_PSR_log(x),
    X_max
)

shap_values = explainer.shap_values(X_max, nsamples=1000)

shap.summary_plot(
    shap_values,
    X_max,
    feature_names=param_names,
    show=True
)
