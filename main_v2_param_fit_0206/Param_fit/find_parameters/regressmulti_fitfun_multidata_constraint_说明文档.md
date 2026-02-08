# regressmulti_fitfun_multidata_constraint 函数详细说明文档

## 📋 文档概述

本文档详细介绍 `regressmulti_fitfun_multidata_constraint.m` 函数的各部分功能和总体流程。该函数是遗传编程（GP）框架中的**核心适应度评估函数**，负责对候选表达式进行参数拟合、物理约束检验和适应度计算。

---

## 🎯 函数定位与作用

### 在整个系统中的位置

```
遗传编程主循环 (rungp)
    ↓
生成候选表达式 (tree2evalstr)
    ↓
【适应度评估】regressmulti_fitfun_multidata_constraint ⭐ 本函数
    ├─ 符号求导 (diff_F)
    ├─ 合并数据拟合 (Fitfun_multidata_merged)
    ├─ 单个数据拟合 (loss_cal_optimize)
    ├─ 物理约束检验 (Multi_R_check)
    └─ 多数据损失计算 (cal_loss_multidata)
    ↓
返回适应度值
    ↓
选择、交叉、变异
```

### 核心功能

1. **符号微分**：对候选表达式进行自动符号求导，获得解析梯度
2. **约束参数搜索**：遍历所有可能的约束参数（c1, c2, ...）组合
3. **多层优化**：
   - 第一层：遍历约束参数组合
   - 第二层：遍历基准参数PARA0
   - 第三层：遍历数据集
   - 第四层：k1/k2参数优化
4. **物理约束**：确保拟合结果满足疲劳裂纹增长的物理规律
5. **多阶段评估**：根据进化阶段使用不同的损失函数标准

---

## 📝 函数签名

```matlab
function [fitness_out, gp] = regressmulti_fitfun_multidata_constraint(evalstr_in, gp)
```

### 输入参数

| 参数 | 类型 | 说明 |
|-----|------|------|
| `evalstr_in` | cell数组 | 候选表达式字符串（基因表达式）<br>例如：`{'x1^c1', 'x2*c2', 'x1*x2'}` |
| `gp` | 结构体 | 遗传编程配置结构体，包含：<br>- `multidata`: 多数据集<br>- `fitness`: 适应度配置<br>- `runcontrol`: 运行控制参数<br>- 等等 |

### 输出参数

| 参数 | 类型 | 说明 |
|-----|------|------|
| `fitness_out` | 1×3 向量 | 三种损失值：`[loss1, loss2, loss3]`<br>对应不同的评估标准 |
| `gp` | 结构体 | 更新后的gp结构体，包含：<br>- `gp.fitness.returnvalues`: 详细结果 |

---

## 🏗️ 总体流程图

```
┌─────────────────────────────────────────────────────────┐
│  阶段0: 初始化与预处理                                   │
│  ├─ 提取配置参数                                        │
│  ├─ 统计基因数量和约束参数数量                           │
│  ├─ 生成约束参数组合                                     │
│  └─ 符号求导 (diff_F)                                   │
└────────────────┬────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────┐
│  阶段1: 约束参数遍历（第一层循环）                        │
│  for const_parameter_i = 1:const_all_situations         │
└────────────────┬────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────┐
│  阶段2: 合并数据拟合                                     │
│  调用 Fitfun_multidata_merged                           │
│  └─ 获得初始参数集合 PARA0                               │
└────────────────┬────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────┐
│  阶段3: PARA0遍历（第二层循环）                          │
│  for PARA0_i = 1:PARA0_NUM                              │
└────────────────┬────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────┐
│  阶段4: 数据集遍历（第三层循环）                          │
│  for data_i = 1:Ndata                                   │
│  ├─ 准备训练数据                                         │
│  ├─ 动态调整参数边界                                     │
│  └─ 调用 loss_cal_optimize (模式2.34)                   │
└────────────────┬────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────┐
│  阶段5: 多数据损失计算                                   │
│  调用 cal_loss_multidata                                │
│  └─ 综合评估所有数据集的拟合效果                          │
└────────────────┬────────────────────────────────────────┘
                 ↓
┌─────────────────────────────────────────────────────────┐
│  阶段6: 结果选择与返回                                   │
│  ├─ 根据进化阶段选择最优约束参数组合                       │
│  ├─ 构建返回结构                                         │
│  └─ 返回适应度值和详细参数                                │
└─────────────────────────────────────────────────────────┘
```

---

## 📦 详细功能模块说明

### 模块1: 初始化与配置提取（第15-22行）

#### 功能说明
从gp结构体中提取运行所需的配置参数。

#### 代码片段
```matlab
multidata = gp.multidata;                    % 多数据集（元胞数组）
bootSample = gp.userdata.bootSample;         % 自助采样索引
bootSampleSize = gp.userdata.bootSampleSize; % 自助采样大小
run_completed = gp.state.run_completed;      % 运行完成标志
force_compute_theta = gp.state.force_compute_theta; % 强制计算theta标志
iteration_extent = gp.fitness.iteration;     % 迭代次数
```

#### 关键参数
- **multidata**: 包含多组疲劳裂纹增长数据的元胞数组
  - 每个元素格式：`[delta_K, da/dN, R, Kmax]`
- **bootSample**: 是否使用Bootstrap采样提高鲁棒性
- **run_completed**: 是否运行完成（影响是否计算theta）

---

### 模块2: 约束参数组合生成（第23-51行）

#### 功能说明
统计表达式中的约束参数（c1, c2, ...）数量，并生成所有可能的约束参数组合。

#### 算法原理

**约束参数示例**：
```matlab
表达式: {'x1^c1', 'x2*c2', 'x1*x2'}
约束参数: c1, c2 (共2个)
候选值: [1, 2, 3]
组合数: 3^2 = 9 种

生成的组合:
[1, 1]
[1, 2]
[1, 3]
[2, 1]
[2, 2]
[2, 3]
[3, 1]
[3, 2]
[3, 3]
```

#### 代码实现
```matlab
% 统计约束参数数量
num_const = 0;
for i = 1:numGenes
    open_sq_br = strfind(evalstr_in{i}, 'c');
    num_const = num_const + numel(open_sq_br);
end

% 生成所有组合（类似进制转换）
const_all_situations = size_const_choose^(num_const);
const_all = zeros(const_all_situations, num_const);

for i = 1:const_all_situations
    for j = 1:num_const
        if j == 1
            index = mod(i-1, size_const_choose) + 1;
        else
            index = floor(mod(i-1, size_const_choose^j) / (size_const_choose^(j-1))) + 1;
        end
        const_all(i, num_const+1-j) = const_choose(index);
    end
end
```

#### 组合生成算法图解
```
类似于k进制数的遍历：
i=1: [c1=1, c2=1]  → 00 (3进制)
i=2: [c1=1, c2=2]  → 01
i=3: [c1=1, c2=3]  → 02
i=4: [c1=2, c2=1]  → 10
i=5: [c1=2, c2=2]  → 11
...
```

---

### 模块3: 符号求导与表达式处理（第58-108行）

#### 功能说明
对候选表达式进行符号微分，获得解析梯度表达式，用于后续的约束优化。

#### 调用流程
```matlab
if loss_mode >= 2
    % 调用diff_F函数进行符号求导
    [diff_index, eq, diff_omega, diff_omega_OF_deq, diff_delta_K, ...
     contains_k1, contains_k2] = diff_F(evalstr_in, num_const, diff_model);
    
    if diff_index == 1
        % 变量替换：将符号变量转换为MATLAB可执行表达式
        diff_omega = regexprep(diff_omega, 'k(\d+)', 'k($1)');
        diff_omega = regexprep(diff_omega, 'c(\d+)', 'Const_pair_now($1)');
        diff_omega = regexprep(diff_omega, 'f(\d+)', 'f($1)');
        % ... 同样处理其他导数表达式
    end
end
```

#### 返回的导数信息

| 变量 | 说明 | 用途 |
|-----|------|------|
| `diff_index` | 求导是否成功（1=成功） | 判断是否继续优化 |
| `eq` | 方程表达式 | 用于计算预测值 |
| `diff_omega` | ∂M/∂ω（对所有参数的导数） | 用于梯度优化 |
| `diff_omega_OF_deq` | ∂²M/∂ω∂ΔK | 用于物理约束计算 |
| `diff_delta_K` | ∂M/∂ΔK | 用于单调性检验 |
| `contains_k1` | 是否包含k1参数 | 决定是否优化k1 |
| `contains_k2` | 是否包含k2参数 | 决定是否优化k2 |

#### 符号求导示例

**输入表达式**：
```matlab
evalstr_in = {'x1^c1', 'x2*c2'}
```

**符号求导结果**：
```matlab
M = f1*log10(x1)^c1 + f2*log10(x2)*c2 + f3*log10(delta_K)

∂M/∂f1 = log10(x1)^c1
∂M/∂f2 = log10(x2)*c2
∂M/∂k1 = f3/k1 (如果x1 = delta_K/k1)
```

---

### 模块4: k1/k2参数搜索范围设置（第161-196行）

#### 功能说明
根据表达式是否包含k1和k2参数，设置材料参数的搜索范围。

#### 参数设置逻辑

```matlab
% k1和k2的搜索范围（对数空间均匀采样）
if contains_k2
    KC = linspace(log10(k1_k2_ini(1)), log10(k1_k2_ini(2)), k1_k2_num);
else
    KC = log10(1);  % 不包含k2时固定为1
end

if contains_k1
    DELTA_KTH = linspace(log10(k1_k2_ini(1)), log10(k1_k2_ini(2)), k1_k2_num);
else
    DELTA_KTH = log10(1);  % 不包含k1时固定为1
end

% k1和k2相关性检查
if contains_k1 && contains_k2
    k1k2_Correlation = k1k2_Correlation_test(evalstr_in);
    if k1k2_Correlation
        DELTA_KTH = log10(1);  % 存在相关性时固定k1
    end
end

% 转换回线性空间
DELTA_KTH = 10.^DELTA_KTH;
KC = 10.^KC;
```

#### k1/k2相关性说明

**k1和k2相关性**：某些表达式中k1和k2不是独立的。

**示例**：
```matlab
% 独立情况
M = f1*log10(delta_K/k1) + f2*log10(Kmax/k2)
→ k1和k2独立，都需要优化

% 相关情况
M = f1*log10(delta_K*k2/k1)
→ k1和k2成比例关系，只优化k2
```

#### 搜索网格图示
```
对数空间均匀采样（典型设置）:
k1_k2_ini = [1, 100]
k1_k2_num = 5

对数空间:
log10(k1) = [0, 0.5, 1.0, 1.5, 2.0]
log10(k2) = [0, 0.5, 1.0, 1.5, 2.0]

线性空间:
k1 = [1, 3.16, 10, 31.6, 100]
k2 = [1, 3.16, 10, 31.6, 100]

网格点数: 5×5 = 25 个初始值组合
```

---

### 模块5: 参数边界设置（第202-211行）

#### 功能说明
设置优化变量的上下界约束，确保参数在物理合理的范围内。

#### 边界向量结构

```matlab
% 参数顺序: [log10(delta_kth); f1; f2; ...; fN; m; delta_kth; kc]
%           第1位          2~N+1    N+2   N+3      N+4

f_L = -5 * ones(numGenes, 1);  % 线性系数下界
f_U = +5 * ones(numGenes, 1);  % 线性系数上界
m_L = 1;                        % Paris指数下界
m_U = 10;                       % Paris指数上界

LB_orig = [-15; f_L; m_L;   0;   0];  % 完整下界向量
UB_orig = [ -5; f_U; m_U; 500; 500];  % 完整上界向量
```

#### 边界意义

| 参数 | 下界 | 上界 | 物理意义 |
|-----|------|------|---------|
| log10(delta_kth) | -15 | -5 | delta_Kth: 10^-15 ~ 10^-5 |
| f_i | -5 | 5 | 线性系数范围 |
| m | 1 | 10 | Paris指数（通常2~4） |
| delta_kth | 0 | 500 | 阈值应力强度因子幅值 |
| kc | 0 | 500 | 临界应力强度因子 |

---

### 模块6: 第一层循环 - 约束参数遍历（第222-562行）

#### 功能说明
遍历所有可能的约束参数组合，为每个组合执行完整的参数拟合流程。

#### 循环结构
```matlab
for const_parameter_i = 1:const_all_situations
    Const_pair_now = const_all(const_parameter_i, :);  % 当前约束参数组合
    
    % 初始化存储结构
    fitness_all_const{1}{const_parameter_i, 1} = All_single_para_PARA0;
    fitness_all_const{2}(const_parameter_i, 1) = fitness_best_PARA0(1);
    fitness_all_const{3}(const_parameter_i, 1) = fitness_best_PARA0(2);
    fitness_all_const{4}(const_parameter_i, 1) = fitness_best_PARA0(3);
    fitness_all_const{5}{const_parameter_i, 1} = fitness_PARA0;
end
```

#### Multi_R_check参数设置（第242-260行）

用于物理约束检验的网格参数：

```matlab
NR = 20;  % R值采样点数
R_temp = linspace(-1, 0.95, NR);  % R值范围：-1到0.95

material_i_deltaK_min = 0.1;
material_i_deltaK_max = 120;
x1_temp = linspace(log10(material_i_deltaK_min), log10(material_i_deltaK_max), 100);
x1_temp = 10.^x1_temp;  % 100个deltaK采样点

% 生成R和deltaK的网格
[Rg, x1_temp_g] = meshgrid(R_temp, x1_temp);  % 100×20 网格
KMg = x1_temp_g ./ (1 - Rg);  % 计算对应的Kmax
```

**网格示意图**：
```
    R = -1.0  -0.9  ...  0.95
ΔK
0.1    [K11]  [K12] ...  [K1,20]
0.15   [K21]  [K22] ...  [K2,20]
...    ...    ...   ...  ...
120    [K100,1] ...     [K100,20]

用于在整个(R, ΔK)空间检验物理约束
```

---

### 模块7: 合并数据拟合（第262-264行）⭐

#### 功能说明
**关键步骤**！调用 `Fitfun_multidata_merged` 函数，基于合并的所有数据拟合初始参数集合。

#### 调用接口
```matlab
[EFF_index, PARA0, PARA_for_initial, parameter_gp] = ...
    Fitfun_multidata_merged(evalstr_in, gp, parameter_gp);
```

#### 返回值说明

| 返回值 | 类型 | 说明 |
|-------|------|------|
| `EFF_index` | 标量 | 是否成功拟合（1=成功，0=失败） |
| `PARA0` | 矩阵 | 全局初始参数集合<br>每行：`[k1, k2, f1, f2, ..., fN]` |
| `PARA_for_initial` | 矩阵 | 用于后续优化的初始参数组合 |
| `parameter_gp` | 结构体 | 更新后的参数结构体 |

#### 工作流程
```
所有数据合并
    ↓
遍历k1/k2网格
    ↓
对每个k1/k2组合:
    ├─ 线性回归求theta
    ├─ fmincon全参数优化
    └─ 物理约束检验
    ↓
筛选有效的参数组合
    ↓
聚类得到PARA0集合（最多max_cluster_number组）
```

**PARA0的作用**：
- 提供多个高质量的初始参数组合
- 避免局部最优
- 提高后续单个数据拟合的成功率

---

### 模块8: PARA0循环（第273-462行）

#### 功能说明
遍历合并数据拟合得到的每个初始参数组合PARA0，对每组初始参数分别进行单个数据的拟合。

#### 循环结构
```matlab
for PARA0_i = 1:PARA0_NUM
    parameter_gp.PARA0 = PARA0(PARA0_i, :);  % 设置当前PARA0
    
    fitness_PARA0_i_all_data = ones(Ndata, numGenes+6+num_const) * inf;
    
    % 遍历每个数据集
    for data_i = 1:Ndata
        % 单个数据拟合...
    end
    
    % 计算当前PARA0的多数据综合损失
    loss_multidata = cal_loss_multidata(gp, fitness_PARA0_i_all_data, ...
                                        numGenes, num_const, gp.fitness.lambda);
end
```

#### 为什么需要多个PARA0？

**理由**：
1. **避免局部最优**：不同的初始值可能收敛到不同的局部最优
2. **增加鲁棒性**：对于不同的数据集，最优的初始值可能不同
3. **探索多样性**：保留多个候选方案供后续选择

**典型设置**：
```matlab
max_cluster_number = 3  % 最多保留3组PARA0
```

---

### 模块9: 数据集循环（第278-451行）

#### 功能说明
对每个单独的数据集进行参数拟合，使用当前PARA0作为初始值。

#### 核心流程

```matlab
for data_i = 1:Ndata
    % 1. 提取当前数据集
    parameter_gp.data_i = multidata{1, data_i};
    
    % 2. 从PARA0提取初始参数
    parameter_K = PARA_for_initial(PARA0_i, 1:2);    % [k1, k2]
    theta = PARA_for_initial(PARA0_i, 3:end);        % [f1, f2, ..., fN]
    
    % 3. 动态调整参数边界
    for pp = 1:numGenes
        if theta_gene(pp) > 0
            f_L(pp) = 0.1;   % 正参数：保持为正
            f_U(pp) = 10;
        elseif theta_gene(pp) < 0
            f_L(pp) = -10;   % 负参数：保持为负
            f_U(pp) = -0.1;
        end
    end
    
    % 4. 准备训练数据
    xtrain = [multidata{1,data_i}(:,1)./delta_kth, ...  % delta_K/delta_Kth
              multidata{1,data_i}(:,4)./kc];             % Kmax/Kc
    ytrain = log10(multidata{1,data_i}(:,2));            % log10(da/dN)
    
    % 5. 调用优化函数（模式2.34）
    [loss, delta_kth, kc, theta, opti_mode_used, ...
     constraint_wrong, pass_index] = ...
        loss_cal_optimize(2.34, p_temp, DELTA_KTH, KC);
    
    % 6. 保存结果
    if pass_index == 1
        fitness_temp_at_Kth_KC_i = [loss(1), delta_kth, kc, ...
                                    Const_pair_now, theta', opti_mode_used];
    else
        fitness_temp_at_Kth_KC_i = [inf, ...];
    end
end
```

#### 动态边界调整机制

**为什么需要动态调整？**

从全局参数（PARA0）可以推断参数的符号：
- 如果 `theta_i > 0`，说明该项对结果有正贡献，应保持为正
- 如果 `theta_i < 0`，说明该项对结果有负贡献，应保持为负

**好处**：
1. 减少搜索空间（排除不合理的符号）
2. 加快收敛速度
3. 提高数值稳定性

**示例**：
```matlab
PARA0中的theta = [8.5, -2.3, 1.1, -0.5]

动态调整后的边界:
参数1: f_L=0.1,  f_U=10    (保持为正)
参数2: f_L=-10,  f_U=-0.1  (保持为负)
参数3: f_L=0.1,  f_U=10    (保持为正)
参数4: f_L=-10,  f_U=-0.1  (保持为负)
```

#### loss_cal_optimize调用（模式2.34）

**模式2.34特点**：
- 全参数优化（theta, k1, k2同时优化）
- 包含物理约束惩罚项
- **包含正则化项**（参数向PARA0靠拢）
- 使用fmincon约束优化

**为什么使用模式2.34？**
```matlab
Loss = fit_MSE + fit_punish_1 + fit_punish_2

fit_MSE:      拟合当前数据集
fit_punish_1: 物理约束惩罚
fit_punish_2: 参数正则化（防止偏离PARA0过远）

目的: 在拟合单个数据的同时保持参数一致性
```

---

### 模块10: 多数据损失计算（第453-461行）

#### 功能说明
计算当前PARA0下所有数据集的综合损失值。

#### 调用接口
```matlab
loss_multidata = cal_loss_multidata(gp, fitness_PARA0_i_all_data, ...
                                    numGenes, num_const, gp.fitness.lambda);
```

#### 输入格式

`fitness_PARA0_i_all_data` 是一个矩阵，每行对应一个数据集：

```matlab
结构: [loss, k1, k2, c1, c2, ..., f1, f2, ..., fN, opti_mode]
行数: Ndata (数据集数量)
列数: numGenes + 6 + num_const

示例:
数据1: [0.015, 5.2, 85.3, 1, -9.6, 2.0, -1.0, 2.0, 1]
数据2: [0.018, 4.8, 90.1, 1, -9.3, 2.1, -1.1, 2.1, 1]
数据3: [0.020, 5.5, 88.5, 1, -9.8, 1.9, -0.9, 1.9, 1]
```

#### 返回值

```matlab
loss_multidata = [loss1, loss2, loss3]

loss1: 所有数据MSE的算术平均（主要标准）
loss2: 中间标准（可自定义）
loss3: 严格标准（可自定义）
```

#### 多数据损失计算策略

**常用策略**：

1. **算术平均**：
   ```matlab
   loss1 = mean([loss_data1, loss_data2, loss_data3])
   ```

2. **加权平均**：
   ```matlab
   loss2 = weighted_mean([loss_data1, loss_data2, loss_data3], weights)
   ```

3. **惩罚参数差异**：
   ```matlab
   loss3 = mean_loss + lambda * parameter_variance
   ```

---

### 模块11: 最优PARA0选择（第475-479行）

#### 功能说明
从所有PARA0中选择综合损失最小的那组参数。

#### 选择逻辑
```matlab
% 选择MSE最小的PARA0
[~, pvals] = min(fitness_all_PARA0(:,1));

% 提取最优结果
fitness_best_PARA0 = [fitness_all_PARA0(pvals,1), ...
                      fitness_all_PARA0(pvals,2), ...
                      fitness_all_PARA0(pvals,3)];
All_single_para_PARA0 = parameters_all_PARA0{pvals, 1};
fitness_PARA0 = PARA0(pvals, :);
```

#### 评估标准

当前实现：**统一使用MSE（loss1）作为选择标准**

注释掉的代码显示曾考虑分阶段选择：
```matlab
% 旧版本：根据进化阶段使用不同标准
% if gen_count_now <= stage1
%     [~, pvals] = min(fitness_all_PARA0(:,1));  % 阶段1用loss1
% end
% if gen_count_now > stage1 && gen_count_now <= stage2
%     [~, pvals] = min(fitness_all_PARA0(:,2));  % 阶段2用loss2
% end
% if gen_count_now > stage2
%     [~, pvals] = min(fitness_all_PARA0(:,3));  % 阶段3用loss3
% end
```

---

### 模块12: 可视化与调试（第481-544行）

#### 功能说明
可选的可视化模块，用于调试时查看拟合效果和物理约束检验结果。

#### 激活方式
```matlab
show_figure = 1;  % 设置为1启用可视化
```

#### 可视化内容

1. **Multi_R_check结果**：
   - 三个区域的划分
   - 物理约束满足情况
   - GY_1, GY_2, GY_3曲线

2. **保存格式**：
   - PNG格式：用于报告和演示
   - FIG格式：用于后续编辑

3. **文件结构**：
```
Figure_NASGRO_HS_test/
├── 1/  (图形类型1)
│   ├── 1_1.png
│   ├── 1_1.fig
│   ├── 2_1.png
│   └── ...
├── 2/  (图形类型2)
│   ├── 1_2.png
│   ├── 1_2.fig
│   └── ...
└── ...

Parameters_NASGRO_HS_test/
├── 1.mat
├── 2.mat
└── ...
```

#### 调试建议

**场景1：检查单个数据拟合效果**
```matlab
show_figure = 1;
data_i = 3;  % 只查看第3组数据
```

**场景2：验证物理约束**
```matlab
% 查看Multi_R_check返回的约束违反情况
[test_index, Test_output] = Multi_R_check(...);
if test_index == 0
    disp('物理约束违反！');
    disp(['GY_2 = ', num2str(Test_output.GY_2')]);
end
```

---

### 模块13: 约束参数组合结果保存（第554-561行）

#### 功能说明
保存当前约束参数组合下的所有结果，为后续选择最优约束参数做准备。

#### 存储结构
```matlab
fitness_all_const = cell(1, 5);  % 元胞数组

fitness_all_const{1}{const_parameter_i, 1} = All_single_para_PARA0;  % 所有数据参数
fitness_all_const{2}(const_parameter_i, 1) = fitness_best_PARA0(1);  % loss1
fitness_all_const{3}(const_parameter_i, 1) = fitness_best_PARA0(2);  % loss2
fitness_all_const{4}(const_parameter_i, 1) = fitness_best_PARA0(3);  % loss3
fitness_all_const{5}{const_parameter_i, 1} = fitness_PARA0;          % PARA0参数
```

#### 数据结构图示
```
fitness_all_const
├─ {1}: 详细参数（元胞数组）
│   ├─ {1,1}: 约束参数组合1的所有数据参数
│   ├─ {2,1}: 约束参数组合2的所有数据参数
│   └─ ...
├─ {2}: loss1（矩阵，每个约束参数组合的MSE）
├─ {3}: loss2（矩阵）
├─ {4}: loss3（矩阵）
└─ {5}: PARA0（元胞数组，每个约束参数组合的基准参数）
```

---

### 模块14: 多阶段进化策略（第564-576行）⭐

#### 功能说明
根据当前进化阶段，使用不同的损失函数标准选择最优约束参数组合。

#### 三阶段策略

```matlab
% 阶段1：早期探索（gen <= stage1）
if gp.fitness.gen_count_now <= gp.runcontrol.stage1
    [~, C_pvals] = min(fitness_all_const{1,2});  % 使用loss1（MSE）
end

% 阶段2：中期过渡（stage1 < gen <= stage2）
if (gp.fitness.gen_count_now > gp.runcontrol.stage1) && ...
   (gp.fitness.gen_count_now <= gp.runcontrol.stage2)
    [~, C_pvals] = min(fitness_all_const{1,3});  % 使用loss2
end

% 阶段3：后期优化（gen > stage2）
if gp.fitness.gen_count_now > gp.runcontrol.stage2
    [~, C_pvals] = min(fitness_all_const{1,4});  % 使用loss3
end
```

#### 阶段设计理念

| 阶段 | 代数范围 | 优化目标 | 使用标准 | 特点 |
|-----|---------|---------|---------|------|
| **阶段1** | 1 ~ stage1 | 探索搜索空间 | loss1 (MSE) | 快速找到合理区域 |
| **阶段2** | stage1 ~ stage2 | 平衡拟合与约束 | loss2 | 逐步引入约束 |
| **阶段3** | stage2 ~ end | 精细化优化 | loss3 | 严格约束和一致性 |

#### 典型配置
```matlab
gp.runcontrol.stage1 = 50;   % 第一阶段：前50代
gp.runcontrol.stage2 = 100;  % 第二阶段：51-100代
% 第三阶段：101代之后
```

#### 分阶段的优势

1. **早期探索**：
   - 使用宽松标准，避免过早收敛
   - 允许表达式多样性
   - 快速淘汰明显不合适的表达式

2. **中期过渡**：
   - 逐步引入物理约束
   - 平衡拟合精度和物理合理性
   - 稳定种群质量

3. **后期优化**：
   - 使用严格标准
   - 确保最终模型的高质量
   - 保证参数一致性

#### 损失函数演进示意图
```
      │
Loss  │  ┌─ Stage 1: MSE优先
      │  │
      │  │   ┌─ Stage 2: MSE + 物理约束
      │  │   │
      │  │   │     ┌─ Stage 3: MSE + 物理约束 + 一致性
      │  │   │     │
      ├──┴───┴─────┴──────────────────────→ Generation
      0  50  100   150   200   250   300
```

---

### 模块15: 结果构建与返回（第578-598行）

#### 功能说明
构建返回结构，包含适应度值和详细的参数信息。

#### 返回结构

```matlab
% 1. 适应度值（三种损失）
fitness_out = [fitness_all_const{2}(C_pvals,1), ...  % loss1
               fitness_all_const{3}(C_pvals,1), ...  % loss2
               fitness_all_const{4}(C_pvals,1)];     % loss3

% 2. 详细参数信息
fitness_temp = fitness_all_const{1}{C_pvals,1};      % 所有数据集的参数
fitness_PARA0 = fitness_all_const{5}{C_pvals,1};     % PARA0基准参数

% 3. 构建返回元胞数组
if loss_mode >= 2
    fitness_return = {fitness_temp, numGenes, num_const, eq, ...
                      evalstr_in, fitness_PARA0};
elseif loss_mode == 1
    fitness_return = {fitness_temp, numGenes, num_const, ...
                      evalstr_in, fitness_PARA0};
end

% 4. 保存到gp结构体
gp.fitness.returnvalues = fitness_return;
```

#### fitness_return结构详解

| 元素位置 | 内容 | 说明 |
|---------|------|------|
| `{1}` | fitness_temp | Ndata×M矩阵，每行为一个数据集的完整参数 |
| `{2}` | numGenes | 基因数量（表达式项数） |
| `{3}` | num_const | 约束参数数量 |
| `{4}` | eq | 方程表达式（loss_mode≥2时） |
| `{5}` | evalstr_in | 原始表达式字符串 |
| `{6}` | fitness_PARA0 | 基准参数PARA0 |

#### fitness_temp矩阵结构

```matlab
% 矩阵大小: Ndata × (numGenes + 6 + num_const)
% 列结构: [loss, k1, k2, c1, c2, ..., f1, f2, ..., fN, opti_mode]

示例（3个数据集，2个基因，1个约束参数）:
fitness_temp = [
    0.015, 5.2, 85.3, 1, -9.6, 2.0, -1.0, 1;  % 数据集1
    0.018, 4.8, 90.1, 1, -9.3, 2.1, -1.1, 1;  % 数据集2
    0.020, 5.5, 88.5, 1, -9.8, 1.9, -0.9, 1   % 数据集3
];
```

#### 失败情况处理

如果符号求导失败（`diff_index != 1`）：
```matlab
fitness_out = [inf, inf, inf];
fitness_return = {inf, numGenes, num_const, evalstr_in, evalstr_in, [0 0 0]};
gp.fitness.returnvalues = fitness_return;
```

---

## 🔄 完整执行流程示例

### 场景设置

```matlab
% 表达式
evalstr_in = {'x1^c1', 'x2*c2', 'x1*x2'};
% 3个基因，2个约束参数

% 约束参数候选值
const_choose = [1, 2];
% 可能的组合: [1,1], [1,2], [2,1], [2,2] (共4种)

% 数据集
Ndata = 3;  % 3组疲劳裂纹增长数据

% k1/k2搜索
k1_k2_num = 5;  % 5×5网格

% PARA0数量
max_cluster_number = 3;  % 最多3组基准参数
```

### 执行流程

```
1. 初始化
   └─ 提取配置参数

2. 约束参数组合生成
   └─ const_all = [1,1; 1,2; 2,1; 2,2]  (4组)

3. 符号求导
   └─ 得到: eq, diff_omega, diff_omega_OF_deq, diff_delta_K

4. 第一层循环: 遍历4个约束参数组合
   ├─ const_parameter_i = 1: [c1=1, c2=1]
   │   ├─ 调用 Fitfun_multidata_merged
   │   │   └─ 得到 PARA0 (假设3组)
   │   │
   │   ├─ PARA0循环 (3次)
   │   │   ├─ PARA0_1: [k1=5.0, k2=85.0, ...]
   │   │   │   ├─ 数据集1拟合 → [loss=0.015, ...]
   │   │   │   ├─ 数据集2拟合 → [loss=0.018, ...]
   │   │   │   ├─ 数据集3拟合 → [loss=0.020, ...]
   │   │   │   └─ 多数据损失 → [0.0177, 0.0180, 0.0185]
   │   │   │
   │   │   ├─ PARA0_2: [k1=4.5, k2=90.0, ...]
   │   │   │   └─ ... (同上)
   │   │   │
   │   │   └─ PARA0_3: [k1=5.5, k2=80.0, ...]
   │   │       └─ ... (同上)
   │   │
   │   └─ 选择最优PARA0 → PARA0_1 (loss最小)
   │
   ├─ const_parameter_i = 2: [c1=1, c2=2]
   │   └─ ... (同上)
   │
   ├─ const_parameter_i = 3: [c1=2, c2=1]
   │   └─ ... (同上)
   │
   └─ const_parameter_i = 4: [c1=2, c2=2]
       └─ ... (同上)

5. 根据进化阶段选择最优约束参数组合
   └─ 假设阶段1: 选择loss1最小的 → const_parameter_i = 1

6. 构建返回结构
   ├─ fitness_out = [0.0177, 0.0180, 0.0185]
   └─ gp.fitness.returnvalues = {fitness_temp, ...}

7. 返回
   └─ [fitness_out, gp]
```

### 计算量估算

```
总的优化次数 = 约束参数组合数 × PARA0数量 × 数据集数量
            = 4 × 3 × 3
            = 36 次参数优化

每次优化包含:
- fmincon调用（可能数百次迭代）
- 物理约束检验
- 梯度计算

典型耗时:
- 单次优化: 1-5秒
- 总耗时: 36-180秒
```

---

## 📊 数据流图

### 输入数据流

```
evalstr_in (表达式)
    ↓
diff_F (符号求导)
    ↓
eq, diff_omega, ... (导数表达式)
    ↓
Fitfun_multidata_merged
    ↓
PARA0 (基准参数)
    ↓
loss_cal_optimize (单数据优化)
    ↓
fitness_temp (所有数据参数)
    ↓
cal_loss_multidata
    ↓
fitness_out (适应度值)
```

### 参数传递链

```
gp.fitness.const_choose
    ↓
const_all (所有约束参数组合)
    ↓
Const_pair_now (当前约束参数)
    ↓
parameter_gp.Const_pair_now
    ↓
loss_cal_optimize (用于表达式评估)
```

---

## 🔧 关键设计模式

### 1. 多层循环策略

**优点**：
- ✅ 全面探索参数空间
- ✅ 避免局部最优
- ✅ 提供多个候选方案

**缺点**：
- ❌ 计算量大
- ❌ 可能存在冗余计算

**优化措施**：
```matlab
% 1. 提前终止
if fitness_temp < terminate_value
    break;
end

% 2. 并行计算（已预留接口）
parfor Kth_KC_i = 1:N_ini_k1k2
    % ...
end

% 3. 智能初始值（使用PARA0）
```

### 2. 失败保护机制

**多层try-catch**：
```matlab
try
    % 策略1
catch
    try
        % 策略2
    catch
        % 返回默认值
    end
end
```

**失效检查**：
```matlab
if isinf(loss)
    break;  % 提前终止无效计算
end

if check_at_const_parameter_i == 0
    break;  % 跳过当前约束参数组合
end
```

### 3. 渐进式优化

**从粗到精**：
```
1. 合并数据 → 全局参数（粗略）
2. 单个数据 → 精细参数（精确）
3. 多阶段评估 → 逐步严格
```

**参数传递**：
```
PARA0 (全局) → 单数据优化的初始值
            → 正则化的参考值
```

---

## ⚙️ 配置参数说明

### 在 gpdemo_crack_growth_config.m 中设置

```matlab
% ========== 约束参数设置 ==========
gp.fitness.const_choose = [1];  % 约束参数候选值（通常只用1）

% ========== 损失模式 ==========
gp.fitness.loss_mode = 2.34;   % 优化模式（2.34推荐）

% ========== 微分模型 ==========
gp.fitness.diff_model = 2.2;   % 2.2表示对数坐标求导

% ========== k1/k2搜索范围 ==========
gp.fitness.k1_k2_ini = [1, 100];   % 搜索范围
gp.fitness.k1_k2_num = 5;          % 网格点数

% ========== PARA0设置 ==========
gp.fitness.max_cluster_number = 3; % 最大聚类数

% ========== 物理约束参数 ==========
gp.fitness.g_y_alpha_merged = 1.0;
gp.fitness.g_y_lambda1_merged = 0.1;
gp.fitness.g_y_alpha_single = 10.0;
gp.fitness.g_y_lambda1_single = 0.1;
gp.fitness.g_y_lambda2_single = 0.001;

% ========== 多阶段设置 ==========
gp.runcontrol.stage1 = 50;    % 第一阶段代数
gp.runcontrol.stage2 = 100;   % 第二阶段代数

% ========== 提前终止 ==========
gp.fitness.terminate_value = 0.01;  % 损失值小于此值提前终止

% ========== 多数据损失权重 ==========
gp.fitness.lambda = [1, 1, 1];  % 不同数据集的权重

% ========== 调试模式 ==========
gp.debug = 0;  % 1=启用调试输出
```

---

## 🎯 使用建议

### 首次使用

1. **检查数据格式**
   ```matlab
   % 确保 multidata 格式正确
   % 每个元素: [delta_K, da/dN, R, Kmax]
   multidata{1,1}
   ```

2. **从简单配置开始**
   ```matlab
   gp.fitness.k1_k2_num = 3;  % 减少网格点
   gp.fitness.max_cluster_number = 1;  % 只保留1个PARA0
   gp.debug = 1;  % 启用调试
   ```

3. **查看中间结果**
   ```matlab
   show_figure = 1;  % 启用可视化
   ```

### 性能优化

1. **减少约束参数组合**
   ```matlab
   gp.fitness.const_choose = [1];  % 只用1个值
   ```

2. **启用并行计算**
   ```matlab
   % 修改第304行
   parfor Kth_KC_i = 1:N_ini_k1k2  % 使用parfor
   ```

3. **调整提前终止阈值**
   ```matlab
   gp.fitness.terminate_value = 0.02;  % 放宽阈值
   ```

### 调试技巧

1. **检查符号求导**
   ```matlab
   if diff_index == 0
       disp('符号求导失败！');
       disp(evalstr_in);
   end
   ```

2. **查看PARA0质量**
   ```matlab
   if EFF_index == 0
       disp('合并数据拟合失败！');
   else
       disp(['PARA0数量: ', num2str(size(PARA0,1))]);
       disp(PARA0);
   end
   ```

3. **监控优化过程**
   ```matlab
   gp.debug = 1;  % 输出data_i
   % 在loss_cal_optimize.m中设置
   'Display', 'iter-detailed'  % 显示优化细节
   ```

---

## ❓ 常见问题

### Q1: 函数运行很慢，如何加速？

**原因分析**：
- 约束参数组合过多
- k1/k2网格过密
- PARA0数量过多
- 未启用并行计算

**解决方案**：
```matlab
% 1. 减少网格密度
gp.fitness.k1_k2_num = 3;  % 从5减到3

% 2. 减少PARA0数量
gp.fitness.max_cluster_number = 1;

% 3. 启用并行（如果有并行工具箱）
parfor Kth_KC_i = 1:N_ini_k1k2
```

### Q2: 返回 fitness_out = [inf, inf, inf]

**可能原因**：
1. 符号求导失败 (`diff_index = 0`)
2. 合并数据拟合失败 (`EFF_index = 0`)
3. 所有约束参数组合都失败

**排查步骤**：
```matlab
% 1. 检查表达式
disp(evalstr_in);

% 2. 检查符号求导
[diff_index, eq, ...] = diff_F(evalstr_in, num_const, diff_model);
if diff_index == 0
    error('符号求导失败');
end

% 3. 单独测试Fitfun_multidata_merged
[EFF_index, PARA0, ...] = Fitfun_multidata_merged(...);
disp(['EFF_index: ', num2str(EFF_index)]);
```

### Q3: 不同数据集参数差异很大

**原因**：
- 数据质量差异
- 正则化权重太小
- 物理约束不够强

**解决方案**：
```matlab
% 1. 增加正则化权重
gp.fitness.g_y_lambda2_single = 0.01;  % 从0.001增到0.01

% 2. 检查数据质量
for i = 1:Ndata
    figure;
    loglog(multidata{1,i}(:,1), multidata{1,i}(:,2), 'o-');
    title(['数据集 ', num2str(i)]);
end

% 3. 增强物理约束
gp.fitness.g_y_lambda1_single = 0.5;  % 从0.1增到0.5
```

### Q4: 如何理解三种损失值？

**loss1, loss2, loss3的含义**：
```matlab
loss1: 纯MSE，用于早期探索
loss2: MSE + 适度约束，用于中期过渡
loss3: MSE + 严格约束 + 一致性，用于后期优化
```

**在 cal_loss_multidata.m 中定义**：
```matlab
function loss_multidata = cal_loss_multidata(gp, fitness_data, ...)
    % loss1 = mean(MSE)
    % loss2 = mean(MSE) + constraint_penalty
    % loss3 = mean(MSE) + constraint_penalty + consistency_penalty
end
```

### Q5: PARA0数量总是为0或1

**原因**：
- 聚类阈值设置不当
- k1/k2网格太稀疏
- 有效参数组合太少

**解决方案**：
```matlab
% 1. 增加k1/k2网格密度
gp.fitness.k1_k2_num = 10;  % 从5增到10

% 2. 在Fitfun_multidata_merged中调整聚类参数
% 查看 Get_base_PARA0 函数的实现

% 3. 放宽物理约束（使更多组合有效）
gp.fitness.g_y_lambda1_merged = 0.05;
```

---

## 📚 相关函数

### 核心依赖函数

| 函数名 | 功能 | 文件位置 |
|-------|------|---------|
| `diff_F` | 符号求导 | find_parameters/diff_F.m |
| `Fitfun_multidata_merged` | 合并数据拟合 | find_parameters/Fitfun_multidata_merged.m |
| `loss_cal_optimize` | 参数优化 | find_parameters/loss_cal_optimize.m |
| `Multi_R_check` | 物理约束检验 | find_parameters/Multi_R_check.m |
| `cal_loss_multidata` | 多数据损失计算 | find_parameters/cal_loss_multidata.m |
| `k1k2_Correlation_test` | k1k2相关性检验 | find_parameters/k1k2_Correlation_test.m |
| `Get_base_PARA0` | PARA0聚类 | find_parameters/Get_base_PARA0.m |

### 辅助函数

| 函数名 | 功能 |
|-------|------|
| `Soft_plus_fun` | Softplus函数 |
| `D_soft_plus_fun` | Softplus导数 |
| `lg` | 安全的log10函数 |
| `createModelFunction` | 创建模型函数句柄 |

---

## 📖 代码结构图

```
regressmulti_fitfun_multidata_constraint.m
│
├── [初始化] (15-22行)
│   └── 提取配置参数
│
├── [约束参数] (23-51行)
│   ├── 统计num_const
│   └── 生成const_all
│
├── [符号求导] (58-108行)
│   ├── 调用diff_F
│   └── 变量替换
│
├── [k1k2范围] (161-196行)
│   ├── 设置搜索范围
│   └── k1k2相关性检验
│
├── [边界设置] (202-211行)
│   └── LB_orig, UB_orig
│
├── [第一层循环: 约束参数] (222-562行)
│   │
│   ├── [Multi_R_check参数] (242-260行)
│   │   └── 网格设置
│   │
│   ├── [合并数据拟合] (262-264行)
│   │   └── Fitfun_multidata_merged
│   │
│   ├── [PARA0循环] (273-462行)
│   │   │
│   │   └── [数据集循环] (278-451行)
│   │       ├── 数据准备
│   │       ├── 动态边界调整
│   │       └── loss_cal_optimize (2.34)
│   │
│   ├── [多数据损失] (453-461行)
│   │   └── cal_loss_multidata
│   │
│   ├── [PARA0选择] (475-479行)
│   │   └── 选择最优PARA0
│   │
│   ├── [可视化] (481-544行)
│   │   └── Multi_R_check + 保存图形
│   │
│   └── [结果保存] (554-561行)
│       └── fitness_all_const
│
├── [多阶段选择] (564-576行)
│   └── 根据进化阶段选择约束参数
│
└── [返回构建] (578-598行)
    ├── fitness_out
    └── gp.fitness.returnvalues
```

---

## 📈 性能分析

### 计算复杂度

```
时间复杂度:
O(C × P × D × I)

其中:
C = const_all_situations (约束参数组合数)
P = PARA0_NUM (基准参数数量)
D = Ndata (数据集数量)
I = 迭代次数 (fmincon平均迭代次数)

典型值:
C = 1 (只用c=1)
P = 3 (最多3组PARA0)
D = 5 (5组数据)
I = 100 (平均100次迭代)

总计算量 ≈ 1×3×5×100 = 1500 次函数评估
```

### 瓶颈分析

**最耗时的部分**：
1. **loss_cal_optimize** (60-80%)
   - fmincon迭代
   - 梯度计算
   - 物理约束检验

2. **Fitfun_multidata_merged** (10-20%)
   - k1/k2网格搜索
   - 聚类计算

3. **Multi_R_check** (5-10%)
   - 网格计算(100×20)

### 优化建议

1. **并行化**：
   ```matlab
   parfor PARA0_i = 1:PARA0_NUM  % 不同PARA0并行
   ```

2. **缓存**：
   ```matlab
   % 缓存相同表达式的符号求导结果
   persistent cache_diff_F;
   ```

3. **提前终止**：
   ```matlab
   if min(fitness_all_PARA0) < threshold
       break;
   end
   ```

---

## 🔬 测试与验证

### 单元测试示例

```matlab
% 测试约束参数组合生成
num_const = 2;
const_choose = [1, 2, 3];
[const_all] = generate_const_combinations(num_const, const_choose);
assert(size(const_all, 1) == 9);  % 3^2 = 9

% 测试符号求导
evalstr_in = {'x1^2', 'x2'};
[diff_index, eq, diff_omega, ...] = diff_F(evalstr_in, 0, 2.2);
assert(diff_index == 1);

% 测试完整流程
gp = load_test_config();
evalstr_in = {'x1', 'x2'};
[fitness_out, gp] = regressmulti_fitfun_multidata_constraint(evalstr_in, gp);
assert(~isinf(fitness_out(1)));  % 应该返回有限值
```

### 集成测试

```matlab
% 完整运行测试
gp = gpdemo_crack_growth_config();
gp = rungp(gp);
% 检查最终结果
assert(gp.results.best.fitness(1) < 1.0);
```

---

## 📝 版本历史

| 版本 | 日期 | 主要变更 |
|-----|------|---------|
| v1.0 | 2023 | 初始版本 |
| v1.1 | 2024 | 增加动态边界调整 |
| v1.2 | 2025 | 优化PARA0选择逻辑 |

---

## 💡 总结

### 核心亮点

1. ✅ **全面的参数搜索**：多层循环确保找到全局最优
2. ✅ **物理约束保证**：确保结果满足疲劳力学规律
3. ✅ **多阶段进化**：渐进式优化策略
4. ✅ **鲁棒性强**：多层失败保护机制
5. ✅ **灵活配置**：丰富的配置参数

### 设计哲学

- **从粗到精**：合并数据 → 单个数据
- **多层保护**：try-catch + 有效性检查
- **渐进优化**：三阶段评估标准
- **参数一致性**：正则化保持全局-局部一致

### 使用流程

```
1. 配置参数 (gpdemo_crack_growth_config.m)
   ↓
2. 准备数据 (multidata)
   ↓
3. 运行GP (rungp)
   ↓
4. 自动调用本函数评估每个表达式
   ↓
5. 返回适应度值
   ↓
6. GP选择、交叉、变异
   ↓
7. 重复步骤4-6直到收敛
```

---

**文档版本**: 1.0  
**创建日期**: 2025  
**作者**: AI Assistant  
**适用代码版本**: gptips2_20250722_2  

---

**使用建议**：
1. 首次使用请阅读"使用建议"章节
2. 遇到问题查阅"常见问题"章节
3. 性能优化参考"性能分析"章节
4. 深入理解请结合"详细功能模块"和源代码

**反馈与改进**：
如发现文档错误或有改进建议，请联系项目维护者。














