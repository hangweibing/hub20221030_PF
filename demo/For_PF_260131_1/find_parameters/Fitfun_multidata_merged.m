function [EFF_index,PARA_for_merged_data,PARA_for_initial,parameter_gp]=Fitfun_multidata_merged(evalstr_in,gp,parameter_gp)
% FITFUN_MULTIDATA_MERGED 多数据合并拟合函数
% 
% 功能：对遗传编程生成的表达式进行参数拟合和优化
% 
% 输入参数：
%   evalstr_in - 基因表达式字符串（元胞数组）
%   gp - GP结构体，包含配置信息
%   parameter_gp - 参数结构体，包含已计算的导数等信息
% 
% 输出参数：
%   EFF_index - 有效性索引（1表示有效，0表示无效）
%   PARA_for_merged_data - 合并数据的最优参数
%   PARA_for_initial - 初始参数集合（用于后续优化）
%   parameter_gp - 更新后的参数结构体


%% ==================== 第一部分：提取并行计算和配置参数 ====================
Data_merged=gp.multidata_merged;              % 合并后的训练数据
bootSample=gp.userdata.bootSample;              % 是否使用Bootstrap采样
bootSampleSize=gp.userdata.bootSampleSize;      % Bootstrap采样大小
run_completed=gp.state.run_completed;           % 运行是否完成标志
force_compute_theta=gp.state.force_compute_theta;  % 是否强制计算theta
iteration_extent=gp.fitness.iteration;         % 迭代次数设置
loss_mode=gp.fitness.loss_mode;                % 损失函数模式

%% ==================== 第二部分：统计常数数量和生成常数组合 ====================
numGenes = numel(evalstr_in);                   % 基因数量（表达式数量）
num_const=0;                                    % 初始化常数计数器

% 统计所有基因中常数c的数量（c1, c2, c3...）
for i=1:numGenes
    open_sq_br = strfind(evalstr_in{i},'c');    % 查找所有'c'字符位置
    num_const = num_const+numel(open_sq_br);    % 累加常数数量
end

% 获取常数候选值列表 const_choose为[1]
const_choose=gp.fitness.const_choose;
size_const_choose=size(const_choose,2);         % 候选常数数量

% 计算所有可能的常数组合情况数（笛卡尔积）
const_all_situations=size_const_choose^(num_const);
const_all=zeros(const_all_situations,num_const);  % 初始化常数组合矩阵

% 生成所有可能的常数组合（全排列）
for i=1:const_all_situations
    for j=1:num_const
        if j==1
            % 第一个常数的索引
            index=mod(i-1,size_const_choose)+1;      
            const_all(i,num_const+1-j)=const_choose(index);
        else
            % 后续常数的索引（按位计算）
            index=floor(mod(i-1,size_const_choose^j)/(size_const_choose^(j-1)))+1;
            const_all(i,num_const+1-j)=const_choose(index);
        end
    end
end

%% ==================== 第三部分：准备材料参数（k1, k2）的初始值网格 ====================
Ndata=size(Data_merged,2);                     % 数据组数
k1_k2_num=gp.fitness.k1_k2_num;               % k1和k2的网格点数

% 根据是否包含k1和k2，生成参数网格
if parameter_gp.contains_k2
    % 如果包含k2，生成KC的对数空间网格
    KC=linspace(log10(gp.fitness.k1_k2_ini(1)),log10(gp.fitness.k1_k2_ini(2)),k1_k2_num);  
else
    KC=log10(1);                                % 不包含k2时设为1
end

if parameter_gp.contains_k1
    % 如果包含k1，生成DELTA_KTH的对数空间网格
    DELTA_KTH=linspace(log10(gp.fitness.k1_k2_ini(1)),log10(gp.fitness.k1_k2_ini(2)),k1_k2_num); 
else
    DELTA_KTH=log10(1);                         % 不包含k1时设为1
end

% 如果k1和k2相关，则DELTA_KTH设为1
if parameter_gp.contains_k1  &&  parameter_gp.contains_k2
    if parameter_gp.k1k2_Correlation
        DELTA_KTH=log10(1);
    end    
end                      

% 将对数空间转换回线性空间
DELTA_KTH=10.^DELTA_KTH;
KC=10.^KC;

%% ==================== 第四部分：构建频繁调用的函数句柄 ====================
% 准备测试用的表达式字符串
evalstr_test=evalstr_in;     
pat11 = 'x(\d+)';                              % 匹配输入变量x1, x2等的正则表达式
pat22 = 'c(\d+)';                              % 匹配常数c1, c2等的正则表达式

% 将表达式中的c1, c2替换为Const_pair_now(1), Const_pair_now(2)等
evalstr_test = regexprep(evalstr_test,pat22,'Const_pair_now($1)');
% 将表达式中的x1, x2替换为xtest{1}, xtest{2}等（用于测试）
evalstr_test = regexprep(evalstr_test,pat11,'xtest{$1}'); 

vars_test = {'xtest','Const_pair_now'};        % 测试函数的变量列表
parameter_gp.evalstr_test_fun=cell(1, numel(evalstr_test));

% 为每个基因表达式创建函数句柄
for i = 1:numel(evalstr_test)
    parameter_gp.evalstr_test_fun{i} = createModelFunction(evalstr_test{i}, vars_test);
end

% 构建其他相关函数的句柄
vars = {'delta_K','Const_pair_now','f','k','R'};  % 变量列表
parameter_gp.diff_delta_K_fun = createModelFunction(parameter_gp.diff_delta_K, vars);

vars1 = {'delta_K','Const_pair_now','f','k','Kmax'};  % 变量列表1
vars2 = {'delta_K_test','Const_pair_now','f','k','R'}; % 变量列表2

parameter_gp.eq_fun = createModelFunction(parameter_gp.eq, vars1);

% 构建omega导数的函数句柄
parameter_gp.diff_omega_fun=cell(1, numel(parameter_gp.diff_omega));
for i = 1:numel(parameter_gp.diff_omega)
    parameter_gp.diff_omega_fun{i} = createModelFunction(parameter_gp.diff_omega{i}, vars1);
end

% 处理用于曲面测试的omega导数（使用delta_K_test而非delta_K）
diff_omega_OF_deq=parameter_gp.diff_omega_OF_deq;
pat_test = 'delta_K';
diff_omega_OF_deq = regexprep(diff_omega_OF_deq,pat_test,'delta_K_test');  % 替换为测试变量

parameter_gp.diff_omega_OF_deq_fun=cell(1, numel(diff_omega_OF_deq));
for i = 1:numel(diff_omega_OF_deq)
    parameter_gp.diff_omega_OF_deq_fun{i} = createModelFunction(diff_omega_OF_deq{i}, vars2);
end

%% ==================== 第五部分：设置参数优化范围 ====================
f_L=-5*ones(numGenes,1);                       % f参数的下界（除了Paris部分的两个参数）
f_U= 5*ones(numGenes,1);                       % f参数的上界
m_L=1;                                          % m参数的下界
m_U=10;                                         % m参数的上界

% 组合所有参数的下界和上界：[log10(f1); f2...fN; m; k1; k2]
LB_orig = [-15; f_L; m_L;   0;   0];           % 下界向量
UB_orig = [ -5; f_U; m_U; 500; 500];          % 上界向量
parameter_gp.LB_orig=LB_orig;
parameter_gp.UB_orig=UB_orig;

% 传递正则化参数
parameter_gp.g_y_alpha_merged=gp.fitness.g_y_alpha_merged;
parameter_gp.g_y_lambda1_merged=gp.fitness.g_y_lambda1_merged;  

parameter_gp.g_y_alpha_single=gp.fitness.g_y_alpha_single;
parameter_gp.g_y_lambda1_single=gp.fitness.g_y_lambda1_single;      
parameter_gp.g_y_lambda2_single=gp.fitness.g_y_lambda2_single;
   
%% ==================== 第六部分：第一层循环 - 遍历所有常数组合 ====================
for const_parameter_i=1:const_all_situations    % 遍历所有常数组合情况
    check_at_const_parameter_i=1;              % 有效性检查标志（0表示发现loss=inf，放弃该常数组合）
    Const_pair_now=const_all(const_parameter_i,:);  % 当前常数组合
    parameter_gp.Const_pair_now=Const_pair_now;     % 存储到参数结构体
    fitness_best_constC_all_data=zeros(Ndata,numGenes+6+num_const)*inf;  % 初始化最优适应度矩阵

%% ==================== 第七部分：第二层循环 - 遍历数据组 ====================
    for data_i=1:1                             % 第二层循环（当前只处理一组数据）
        if  check_at_const_parameter_i==0
            break;                              % 如果常数组合无效，跳出循环
        end

        % 设置当前数据组
        parameter_gp.data_i=Data_merged{1,data_i};
        [num_of_data_i,~]=size(Data_merged{1,data_i}(:,1));
        parameter_gp.num_of_data_i=num_of_data_i;        
    
        % 生成k1和k2的参数网格组合
        [~,n1]=size(DELTA_KTH);                % k1的网格点数
        [~,n2]=size(KC);                       % k2的网格点数
        N_ini_k1k2=n1*n2;                      % 总组合数
        Parameters_meterial=zeros(N_ini_k1k2,2);  % 初始化参数矩阵
        
        % 生成所有k1和k2的组合对
        for i=1:N_ini_k1k2
            index2=mod(i-1,n2)+1;              % k2的索引
            index1=floor(mod(i-1,n1*n2)/(n2))+1;  % k1的索引
            Parameters_meterial(i,1)= DELTA_KTH(index1);   % 存储k1值
            Parameters_meterial(i,2)= KC(index2);          % 存储k2值
        end

        % 随机打乱参数组合顺序（避免按固定顺序搜索）
        rng(2);                                 % 设置随机种子
        permIndex = randperm(N_ini_k1k2);      % 生成随机排列索引
        Parameters_meterial = Parameters_meterial(permIndex, :);  % 按随机顺序重排

        % 初始化存储数组
        loss_at_Kth_KC_all=zeros(N_ini_k1k2,1)*inf;  % 存储每个参数组合的损失值
        fitness_at_Kth_KC_all=cell(N_ini_k1k2,1);     % 存储每个参数组合的适应度
        fitness_at_Kth_KC_all_check=cell(N_ini_k1k2,1);  % 存储检查用的适应度
        geneOutputs=cell(N_ini_k1k2,1);               % 存储基因输出
        geneOutputs_real=cell(N_ini_k1k2,1);           % 存储实数基因输出

        %% ==================== 第八部分：准备参数和表达式 ====================
        parameter_gp.evalstr1=evalstr_in;      % 原始表达式
        parameter_gp.fmincon_option1=gp.fitness.fmincon_option1;  % 优化选项1
        parameter_gp.fmincon_option2=gp.fitness.fmincon_option2;  % 优化选项2
        parameter_gp.noise_on=gp.fitness.noise_on;      % 是否启用噪声
        parameter_gp.noise_level=gp.fitness.noise_level;  % 噪声水平
        parameter_gp.ridge_on=gp.fitness.ridge_on;     % 是否启用岭回归
        parameter_gp.ridge_k=gp.fitness.ridge_k;       % 岭回归参数
        parameter_gp.theta_end_limit=gp.fitness.theta_end_limit;  % theta的边界限制
        p_temp=parameter_gp;                   % 临时参数结构体
        evalstr=evalstr_in;                    % 临时表达式

        %% ==================== 第九部分：第三层循环 - 遍历k1和k2的初始值 ====================
        Kth_KC_i_EFF_NUM=0;                    % 有效参数组合计数器
       
        for Kth_KC_i=1:N_ini_k1k2              % 遍历所有k1和k2的组合
        %parfor Kth_KC_i=1:N_ini_k1k2          % 并行版本（已注释）
            theta=zeros(numGenes+2,1);         % 初始化theta（线性系数）
            fitness_temp_at_Kth_KC_i=zeros(1,5+num_const+numGenes);  % 初始化适应度向量
            check_fitness=1;                   % 适应度检查标志
            
            % 获取当前k1和k2的值
            parameter_K=Parameters_meterial(Kth_KC_i,:);  % 作为非线性优化的不同起点
            delta_kth=parameter_K(1);          % k1值（delta_K阈值）
            kc=parameter_K(2);                 % k2值（Kc值）
            
            % 准备训练数据
            xtrain=[];                         % 初始化输入矩阵
            ytrain=[];                         % 初始化输出向量
            
            % 构建输入特征：delta_K/delta_Kth, Kmax/Kc
            xtrain=[Data_merged{1,data_i}(:,1)./delta_kth,Data_merged{1,data_i}(:,4)./kc];
            % 构建输出：log10(da/dN)
            ytrain=log10(Data_merged{1,data_i}(:,2));
            
            % 处理表达式字符串，将c1替换为Const_pair_now(1)，x1替换为xtrain(:,1)等
            pat1 = 'x(\d+)';                   % 匹配x1, x2等的正则表达式
            pat2 = 'c(\d+)';                   % 匹配c1, c2等的正则表达式
            evalstr = regexprep(evalstr,pat2,'Const_pair_now($1)');  % 替换常数
            evalstr = regexprep(evalstr,pat1,'xtrain(:,$1)');        % 替换输入变量
            p_temp.evalstr2=evalstr;           % 存储处理后的表达式
            
            [numData,~] =size(ytrain);         % 数据点数量
            
            % 初始化基因输出矩阵（包含偏置列）
            geneOutputs{Kth_KC_i,1} = ones(numData,numGenes+2);
            
            %% ==================== 第十部分：计算基因输出并检查有效性 ====================
            % 初始化检验：如果起点parameter_K就存在基因中有无穷或复数，就放弃这个参数起点
            for i = 1:numGenes
                ind = i + 1;                   % 索引（第1列是偏置）
                
                % 尝试计算基因输出
                try
                    gene_temp=eval([evalstr{i} ';']);  % 执行表达式计算
                catch
                    disp('An error occurred.');  % 发生错误时提示
                end

                % 对基因输出取对数（lg函数）
                geneOutputs{Kth_KC_i,1}(:,ind)=lg(gene_temp);
            end
            
            % 添加最后一列：log10(delta_K)
            geneOutputs{Kth_KC_i,1}(:,numGenes+2)=lg(Data_merged{1,data_i}(:,1));

            % 查找实数点并检查有效性
            [realMask,check_fitness,Complex_index,DATA_OUT]=Find_real_points(geneOutputs{Kth_KC_i,1});
            y_real = ytrain(realMask, :);      % 实数点对应的输出
            
            if check_fitness==0
                % 如果无效，设置适应度为无穷大
                fitness_temp_at_Kth_KC_i=[Inf,delta_kth,kc,Const_pair_now,zeros(1,numGenes+2),0]; 
                fitness_temp_at_Kth_KC_i_check=[inf,inf,inf,delta_kth,kc,Const_pair_now,theta',0];
            elseif check_fitness==1
                geneOutputs_real{Kth_KC_i,1}=DATA_OUT;  % 存储实数输出
                loss_cal_optimize_in=1;         % 允许进行损失计算和优化
                
                %% ==================== 第十一部分：计算线性系数theta ====================
                % 只在实际运行或强制计算时计算权重系数
                if ~run_completed || force_compute_theta
            
                    % 准备最小二乘法的矩阵
                    if bootSample
                        % 使用Bootstrap采样
                        sampleInds = bootsample(geneOutputs_real{Kth_KC_i,1},bootSampleSize);
                        goptrans = geneOutputs_real{Kth_KC_i,1}(sampleInds,:)';
                        prj = goptrans * geneOutputs_real{Kth_KC_i,1}(sampleInds,:);
                        ysample = y_real(sampleInds);
                    else
                        % 使用全部数据
                        goptrans = geneOutputs_real{Kth_KC_i,1}';
                        prj = goptrans * geneOutputs_real{Kth_KC_i,1};
                    end
            
                    % 使用SVD基于最小二乘法计算树权重系数
                    try
                        if bootSample
                            theta = pinv(prj) * goptrans * ysample;  % 伪逆求解
                            if p_temp.ridge_on
                                %% 采用岭回归求解f（防止过拟合）
                                geneOutputs_real{Kth_KC_i,1}(:,1)=[];  % 移除偏置列
                                theta = ridge(ysample,geneOutputs_real{Kth_KC_i,1},p_temp.ridge_k,0);             
                            end                                 
                        else
                            theta = pinv(prj) * goptrans * y_real;
                            if p_temp.ridge_on
                                %% 采用岭回归求解f
                                geneOutputs_real{Kth_KC_i,1}(:,1)=[];  % 移除偏置列
                                theta = ridge(y_real,geneOutputs_real{Kth_KC_i,1},p_temp.ridge_k,0);             
                            end                                  
                        end                        
                    catch
                        % 计算失败时设置标志
                        loss_cal_optimize_in=0;
                        fitness_temp_at_Kth_KC_i=[Inf,delta_kth,kc,Const_pair_now,zeros(1,numGenes+2),0]; 
                        fitness_temp_at_Kth_KC_i_check=[inf,inf,inf,delta_kth,kc,Const_pair_now,theta',0];
                    end
            
                    % 如果系数中有NaN或Inf，标记为无效
                    if any(isinf(theta)) || any(isnan(theta))
                        loss_cal_optimize_in=0;
                        fitness_temp_at_Kth_KC_i=[Inf,delta_kth,kc,Const_pair_now,zeros(1,numGenes+2),0]; 
                        fitness_temp_at_Kth_KC_i_check=[inf,inf,inf,delta_kth,kc,Const_pair_now,theta',0];
                    end
            
                else % 如果是运行后，从返回值字段获取存储的系数
                    error('没有计算theta，需要重新给theta');
                end    
                
                %% ==================== 第十二部分：损失计算和参数优化 ====================
                if loss_cal_optimize_in==1  
       
                    % 设置优化所需的参数
                    p_temp.theta=theta;
                    p_temp.parameter_K=parameter_K;
                    p_temp.iteration_extent=iteration_extent;
                    p_temp.geneOutputs=geneOutputs{Kth_KC_i,1};
                    p_temp.geneOutputs_real=geneOutputs_real{Kth_KC_i,1};
                    p_temp.ytrain=ytrain;
                    
                    showfigure=0;              % 不显示图形

                    % 调用损失计算和优化函数
                    [loss,delta_kth,kc,theta,opti_mode_used,constraint_wrong,pass_index]=loss_cal_optimize(loss_mode,p_temp,DELTA_KTH,KC);
                
                    %% 对非线性优化后的theta与k1,k2进行初步泛化性能检验
                    if pass_index==1
                        % 如果通过检验，记录适应度
                        fitness_temp_at_Kth_KC_i=[loss(1),delta_kth,kc,Const_pair_now,theta',opti_mode_used];     % 格式：[loss, k1, k2, c1, c2..., f1, f2..., mode]
                        fitness_temp_at_Kth_KC_i_check=[loss,delta_kth,kc,Const_pair_now,theta',opti_mode_used]; 
                        Kth_KC_i_EFF_NUM=Kth_KC_i_EFF_NUM+1;  % 有效组合计数加1

                    else
                        % 如果未通过检验，设置适应度为无穷大
                        fitness_temp_at_Kth_KC_i=[inf,delta_kth,kc,Const_pair_now,theta',opti_mode_used];     
                        fitness_temp_at_Kth_KC_i_check=[inf,inf,inf,delta_kth,kc,Const_pair_now,theta',opti_mode_used]; 
                    end

                end
 
            end
            
            % 存储当前参数组合的结果
            loss_at_Kth_KC_all(Kth_KC_i,1)=fitness_temp_at_Kth_KC_i(1);
            fitness_at_Kth_KC_all{Kth_KC_i,1}=fitness_temp_at_Kth_KC_i; 
            fitness_at_Kth_KC_all_check{Kth_KC_i,1}=fitness_temp_at_Kth_KC_i_check; 

            % 如果损失值足够小，提前终止
            if fitness_temp_at_Kth_KC_i(1)<gp.fitness.terminate_value
                break;
            end
            % 如果有效组合数超过10个，提前终止
            if Kth_KC_i_EFF_NUM>10
                break;
            end
            
        end 

        %% ==================== 第十三部分：结果处理和筛选 ====================
        % 将元胞数组转换为矩阵
        fitness_at_Kth_KC_all_metrix=cell2mat(fitness_at_Kth_KC_all);              
        fitness_at_Kth_KC_all_CHECK_metrix=cell2mat(fitness_at_Kth_KC_all_check);

        % 筛选出有效的适应度值（非无穷大）
        mask_inf = isinf(fitness_at_Kth_KC_all_metrix(:,1));
        mask_valid = ~mask_inf;
        fitness_at_Kth_KC_all_EFF = fitness_at_Kth_KC_all_metrix(mask_valid, :);
        fitness_at_Kth_KC_all_EFF_sort=sortrows(fitness_at_Kth_KC_all_EFF, 1);  % 按损失值排序

        fitness_at_Kth_KC_all_CHECK_EFF = fitness_at_Kth_KC_all_CHECK_metrix(mask_valid, :);
        fitness_at_Kth_KC_all_CHECK_EFF_sort=sortrows(fitness_at_Kth_KC_all_CHECK_EFF, 1);
        
        % 确定是否有有效结果
        N_EFF=size(fitness_at_Kth_KC_all_CHECK_EFF,1);
        if N_EFF==0
            EFF_index=0;                       % 无有效结果
        else
            EFF_index=1;                       % 有有效结果
        end

        % 找出所有k1/k2起点中损失最小的那个
        [loss_min_at_Kth_KC_all,pvals]=min(loss_at_Kth_KC_all);      
        if  isinf(loss_min_at_Kth_KC_all)
            check_at_const_parameter_i=0;       % 如果最小损失也是无穷大，放弃这个常数组合
        end
        
        % 存储当前常数组合下的最优结果
        if check_at_const_parameter_i==1
            fitness_best_constC_all_data(data_i,1:numGenes+6+num_const)=fitness_at_Kth_KC_all{pvals,1};     
        end

        %% ==================== 第十四部分：生成初始参数集合 ====================
        if EFF_index==1        
            max_cluster_number=gp.fitness.max_cluster_number;  % 最大聚类数
            
            % 移除不需要的列（损失值、k1、k2、常数、模式标志）
            fitness_at_Kth_KC_all_CHECK_EFF_sort(:,2:3)=[];  % 移除k1和k2
            fitness_at_Kth_KC_all_CHECK_EFF_sort(:,4:3+num_const)=[];  % 移除常数
            fitness_at_Kth_KC_all_CHECK_EFF_sort(:,end)=[];  % 移除模式标志
            
            % 获取基础参数集合（用于后续优化）
            PARA_for_initial=Get_base_PARA0(fitness_at_Kth_KC_all_CHECK_EFF_sort,max_cluster_number);
            PARA_for_merged_data=PARA_for_initial;
        else
            PARA_for_initial=[];                % 无有效结果时返回空
            PARA_for_merged_data=[];
        end

    end 

end

% 函数结束
