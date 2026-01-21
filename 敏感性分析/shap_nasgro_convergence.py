import numpy as np
import pandas as pd
import shap
import matplotlib.pyplot as plt

# ================== 1. 读取参数文件 ==================
filename = 'NASGRO(H-S).xlsx'
df = pd.read_excel(filename)
param_names = df.columns.tolist()
X_full = df.values  # N x 4

# ================== 2. NASGRO 模型定义 ==================
def nasgro_model(X, deltaK=6, R_ratio=0.5):
    log10_D = X[:, 0]
    p = X[:, 1]
    K_thr = X[:, 2]
    A = X[:, 3]

    D = 10 ** log10_D
    Kmax = deltaK / (1 - R_ratio)

    effective_dK = np.maximum(deltaK - K_thr, 1e-6)
    denominator = np.maximum(1 - Kmax / A, 1e-6)

    da_dN = D * (effective_dK) ** p / (denominator) ** (p / 2)
    return np.log10(da_dN + 1e-12)

# ================== 3. SHAP 收敛性分析设置 ==================
sample_list = np.round(np.linspace(200, len(X_full), 8)).astype(int)
d = X_full.shape[1]
shap_convergence = np.zeros((d, len(sample_list)))

# ================== 4. 循环不同样本数计算 SHAP ==================
for k, N in enumerate(sample_list):
    print(f"Running SHAP for N = {N}")

    X = X_full[:N, :]

    explainer = shap.KernelExplainer(
        lambda x: nasgro_model(x),
        X
    )

    shap_values = explainer.shap_values(
        X,
        nsamples=1000,
        l1_reg="num_features(4)"
    )

    # 平均绝对 SHAP
    shap_convergence[:, k] = np.mean(np.abs(shap_values), axis=0)

# ================== 5. 归一化 SHAP（可选，但推荐） ==================
shap_convergence_norm = shap_convergence / np.sum(shap_convergence, axis=0)

# ================== 6. 收敛性可视化 ==================
colors = plt.cm.tab10.colors
markers = ['o', 's', '^', 'd']

plt.figure(figsize=(10, 4))

# ---- 子图 1：原始 Mean |SHAP| ----
plt.subplot(1, 2, 1)
for i in range(d):
    plt.plot(
        sample_list,
        shap_convergence[i, :],
        '-',
        color=colors[i],
        linewidth=1.5
    )
    plt.scatter(
        sample_list,
        shap_convergence[i, :],
        color=colors[i],
        marker=markers[i],
        s=40,
        label=param_names[i]
    )

plt.xlabel('Sample Number')
plt.ylabel('Mean |SHAP|')
plt.title('SHAP Importance Convergence')
plt.grid(True)
plt.legend()

# ---- 子图 2：归一化 SHAP ----
plt.subplot(1, 2, 2)
for i in range(d):
    plt.plot(
        sample_list,
        shap_convergence_norm[i, :],
        '-',
        color=colors[i],
        linewidth=1.5
    )
    plt.scatter(
        sample_list,
        shap_convergence_norm[i, :],
        color=colors[i],
        marker=markers[i],
        s=40
    )

plt.xlabel('Sample Number')
plt.ylabel('Normalized Mean |SHAP|')
plt.title('Normalized SHAP Convergence')
plt.grid(True)

plt.suptitle('SHAP Convergence Analysis for NASGRO Model')
plt.tight_layout()
plt.show()

# ================== 7. 最大样本数下 SHAP 汇总 ==================
final_shap = shap_convergence[:, -1]
final_norm = final_shap / np.sum(final_shap)

shap_table = pd.DataFrame({
    'Parameter': param_names,
    'MeanAbsSHAP': final_shap,
    'NormalizedSHAP': final_norm
})

print("\n--- Final SHAP Importance (Max Sample Size) ---")
print(shap_table)

# ================== 8. SHAP Summary Plot（最大样本） ==================
X_max = X_full[:sample_list[-1], :]
explainer = shap.KernelExplainer(nasgro_model, X_max)
shap_values = explainer.shap_values(X_max, nsamples=1000)

shap.summary_plot(shap_values, X_max, feature_names=param_names)
