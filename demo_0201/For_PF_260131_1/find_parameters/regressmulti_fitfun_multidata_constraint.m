function [fitness_out,gp]=regressmulti_fitfun_multidata_constraint(evalstr_in,gp)
% 多数据回归拟合函数（带约束条件）
% 输入参数:
%   evalstr_in: 备选方程的字符串表达式（基因表达式）
%   gp: 遗传规划结构体，包含所有配置参数和数据
% 输出参数:
%   fitness_out: 适应度输出值（损失函数值）
%   gp: 更新后的遗传规划结构体

% 默认值设置（用于提前退出的情况）
theta=[];ypredtrain=[];fitnessTest=[];ypredtest=[];
r2train=[];r2test=[];r2val=[];geneOutputs=[];geneOutputsTest=[];
geneOutputsVal=[];

%% 并行计算相关参数提取
multidata=gp.multidata;                    % 多数据集（元胞数组）
bootSample=gp.userdata.bootSample;         % 自助采样索引
bootSampleSize=gp.userdata.bootSampleSize; % 自助采样大小
run_completed=gp.state.run_completed;       % 运行完成标志
force_compute_theta=gp.state.force_compute_theta; % 强制计算theta标志
iteration_extent=gp.fitness.iteration;      % 迭代次数

%% 统计基因数量和约束参数数量
numGenes = numel(evalstr_in);              % 基因数量（方程项数）
num_const=0;                                % 初始化约束参数数量
% 遍历所有基因表达式，统计约束参数'c'的出现次数
for i=1:numGenes
    open_sq_br = strfind(evalstr_in{i},'c'); % 查找'c'字符位置
    num_const = num_const+numel(open_sq_br); % 累加约束参数数量
end

%% 生成所有可能的约束参数组合
const_choose=gp.fitness.const_choose;      % 约束参数的候选值
size_const_choose=size(const_choose,2);    % 候选值数量
const_all_situations=size_const_choose^(num_const); % 所有可能的组合数（排列组合）
const_all=zeros(const_all_situations,num_const);     % 初始化组合矩阵

% 生成所有约束参数组合（类似进制转换的方式）
for i=1:const_all_situations
    for j=1:num_const
        if j==1
            % 第一个约束参数：直接取模
            index=mod(i-1,size_const_choose)+1;      
            const_all(i,num_const+1-j)=const_choose(index);
        else
            % 其他约束参数：通过除法和取模计算索引
            index=floor(mod(i-1,size_const_choose^j)/(size_const_choose^(j-1)))+1;
            const_all(i,num_const+1-j)=const_choose(index);
        end
    end
end

Ndata=size(multidata,2);                    % 数据集数量
% gp.fitness.num_const=num_const;
% gp.fitness.numGenes=numGenes;
% gp.fitness.Ndata=Ndata;

%% 对备选方程F求导（用于约束条件计算）
loss_mode=gp.fitness.loss_mode;            % 损失函数模式
diff_model=gp.fitness.diff_model;          % 微分模型类型
parameter_gp.debug=gp.debug;               % 调试标志

% 当损失模式>=2时，需要对方程进行求导
if loss_mode>=2
    % 调用diff_F函数对备选方程进行符号求导
    % 返回: 微分索引、方程表达式、对omega的导数、对omega_OF_deq的导数、
    %      对delta_K的导数、是否包含k1、是否包含k2
    [diff_index,eq,diff_omega,diff_omega_OF_deq,diff_delta_K, contains_k1,contains_k2]=diff_F(evalstr_in,num_const,diff_model);

    % 如果微分索引有效（diff_index==1），进行变量替换
    if diff_index==1
        % 定义正则表达式模式
        pata = 'k(\d+)';  % 匹配k(数字)模式
        patb = 'c(\d+)';  % 匹配c(数字)模式
        patc = 'f(\d+)';  % 匹配f(数字)模式
        
        % 对diff_omega进行变量替换，将符号表达式转换为可执行的MATLAB表达式
        diff_omega = regexprep(diff_omega,pata,'k($1)');
        diff_omega = regexprep(diff_omega,patb,'Const_pair_now($1)'); % c替换为当前约束值
        diff_omega = regexprep(diff_omega,patc,'f($1)');

        % 对diff_omega_OF_deq进行同样的替换
        diff_omega_OF_deq = regexprep(diff_omega_OF_deq,pata,'k($1)');
        diff_omega_OF_deq = regexprep(diff_omega_OF_deq,patb,'Const_pair_now($1)');
        diff_omega_OF_deq = regexprep(diff_omega_OF_deq,patc,'f($1)');

        % 对diff_delta_K进行同样的替换
        diff_delta_K = regexprep(diff_delta_K,pata,'k($1)');
        diff_delta_K = regexprep(diff_delta_K,patb,'Const_pair_now($1)');
        diff_delta_K = regexprep(diff_delta_K,patc,'f($1)'); 

        % 对方程表达式eq进行同样的替换
        eq = regexprep(eq,pata,'k($1)');
        eq = regexprep(eq,patb,'Const_pair_now($1)');
        eq = regexprep(eq,patc,'f($1)');

        % 将处理后的导数表达式保存到parameter_gp结构体中
        parameter_gp.diff_omega=diff_omega;
        parameter_gp.diff_omega_OF_deq=diff_omega_OF_deq;
        parameter_gp.eq=eq;
        parameter_gp.diff_delta_K=diff_delta_K;
        parameter_gp.contains_k1=contains_k1;  % 标记是否包含k1参数
        parameter_gp.contains_k2=contains_k2;  % 标记是否包含k2参数
%     parameter_gp.dF_dk1=dF_dk1;
%     parameter_gp.dF_dk2=dF_dk2;

    end
end

%% 基于数据的差分（用于计算mode2下的loss，这部分代码暂未用到）        
if loss_mode>=10  % 这部分代码暂未用到
    % 模式2.1：基于线性坐标的差分
    if diff_model==2.1
        DELTA_y_DELTA_k1=cell(1,Ndata);  % y对k1的差分
        DELTA_y_DELTA_k2=cell(1,Ndata);  % y对k2的差分
        R=cell(1,Ndata);                 % 应力比R
        k1=cell(1,Ndata);                % k1的中间值
        k2=cell(1,Ndata);                % k2的中间值
        
        % 对每个数据集计算差分
        for data_j=1:Ndata
            D_i=multidata{1,data_j};
            [num_of_data_i,~]=size(D_i(:,1));
            R{data_j}=D_i(1,3);  % 提取应力比R
            
            % 计算y对k1的差分（线性坐标）
            DELTA_y_DELTA_k1{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   D_i(2:num_of_data_i,1)-D_i(1:num_of_data_i-1,1)   );
            % 计算y对k2的差分（线性坐标）
            DELTA_y_DELTA_k2{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   D_i(2:num_of_data_i,4)-D_i(1:num_of_data_i-1,4)   );
            % 计算k1和k2的中间值
            k1{data_j}=0.5*(  D_i(1:num_of_data_i-1,1) +D_i(2:num_of_data_i,1)  );
            k2{data_j}=0.5*(  D_i(1:num_of_data_i-1,4) +D_i(2:num_of_data_i,4)  );
        end

    % 模式2.2：基于对数坐标的差分
    elseif diff_model==2.2
        DELTA_y_DELTA_z1=cell(1,Ndata);  % y对z1的差分（对数坐标）
        DELTA_y_DELTA_z2=cell(1,Ndata);  % y对z2的差分（对数坐标）
        R=cell(1,Ndata);                 % 应力比R
        z1=cell(1,Ndata);                % z1的中间值（对数坐标）
        z2=cell(1,Ndata);                % z2的中间值（对数坐标）
        
        % 对每个数据集计算差分
        for data_j=1:Ndata
            D_i=multidata{1,data_j};
            [num_of_data_i,~]=size(D_i(:,1));
            R{data_j}=D_i(1,3);  % 提取应力比R
            
            % 计算y对z1的差分（对数坐标）
            DELTA_y_DELTA_z1{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   lg(D_i(2:num_of_data_i,1))-lg(D_i(1:num_of_data_i-1,1))   );
            % 计算y对z2的差分（对数坐标）
            DELTA_y_DELTA_z2{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   lg(D_i(2:num_of_data_i,4))-lg(D_i(1:num_of_data_i-1,4))   );
            % 计算z1和z2的中间值（对数坐标）
            z1{data_j}=0.5*(  lg(D_i(1:num_of_data_i-1,1)) +lg(D_i(2:num_of_data_i,1))  )  ;
            z2{data_j}=0.5*(  lg(D_i(1:num_of_data_i-1,4)) +lg(D_i(2:num_of_data_i,4))  )  ;
        end

    end
end

%% 参数拟合准备：设置k1和k2的搜索范围
if diff_index==1
    % 检查k1和k2是否存在相关性
    if contains_k1&&contains_k2
       % k1k2_Correlation表示k1和k2是否存在比例关系
       % 如果只有k1或k2就不用计算相关性
       k1k2_Correlation=k1k2_Correlation_test(evalstr_in);           
       parameter_gp.k1k2_Correlation=k1k2_Correlation;
    end

    %% 备选方程参数拟合：设置参数搜索范围
    % 如果包含k2参数，在对数空间均匀采样
    if contains_k2
        KC=linspace(log10(gp.fitness.k1_k2_ini(1)),log10(gp.fitness.k1_k2_ini(2)),gp.fitness.k1_k2_num);  
    else
        KC=log10(1);    % 如果不包含k2，设为1的对数
    end
    
    % 如果包含k1参数，在对数空间均匀采样
    if contains_k1
        DELTA_KTH=linspace(log10(gp.fitness.k1_k2_ini(1)),log10(gp.fitness.k1_k2_ini(2)),gp.fitness.k1_k2_num); 
    else
        DELTA_KTH=log10(1);  % 如果不包含k1，设为1的对数
    end
    
    % 如果k1和k2存在相关性，则DELTA_KTH固定为1
    if contains_k1&&contains_k2
        if k1k2_Correlation
            DELTA_KTH=log10(1);
        end    
    end
                           
    % 将对数空间的采样值转换回线性空间
    DELTA_KTH=10.^DELTA_KTH;
    KC=10.^KC;   
end


%% 主要拟合流程：多层循环结构
if diff_index==1
 
    %%% 参数范围设置（优化变量的上下界）
    f_L=-5*ones(numGenes,1);            % 除了paris部分的两个参数，其他的参数施加统一约束（下界）
    f_U= 5*ones(numGenes,1);            % 除了paris部分的两个参数，其他的参数施加统一约束（上界）
    m_L=1;                              % m参数的下界
    m_U=10;                             % m参数的上界
    % 参数边界向量：[log10(delta_kth); f参数; m参数; delta_kth; kc]
    LB_orig = [-15; f_L; m_L;   0;   0];    % 下界向量
    UB_orig = [ -5; f_U; m_U; 300; 300];    % 上界向量
    parameter_gp.LB_orig=LB_orig;
    parameter_gp.UB_orig=UB_orig;

    % 保存合并数据和单数据的正则化参数
    parameter_gp.g_y_alpha_merged=gp.fitness.g_y_alpha_merged;       % 合并数据的alpha参数
    parameter_gp.g_y_lambda1_merged=gp.fitness.g_y_lambda1_merged;  % 合并数据的lambda1参数

    parameter_gp.g_y_alpha_single=gp.fitness.g_y_alpha_single;      % 单数据的alpha参数
    parameter_gp.g_y_lambda1_single=gp.fitness.g_y_lambda1_single; % 单数据的lambda1参数
    parameter_gp.g_y_lambda2_single=gp.fitness.g_y_lambda2_single;  % 单数据的lambda2参数

    %% 第一层循环：遍历所有约束参数组合
    for const_parameter_i=1:const_all_situations 

        check_at_const_parameter_i=1;  % 若check_at_const_parameter_i=0（发现了loss=inf的情况），
                                       % 则抛弃方程在这个常数组下的探索
        Const_pair_now=const_all(const_parameter_i,:);  % 当前约束参数组合
        parameter_gp.Const_pair_now=Const_pair_now;     % 保存到参数结构体
        parameter_gp.evalstr1=evalstr_in;                % 保存原始表达式
        
        % 保存优化选项和正则化参数
        parameter_gp.fmincon_option1=gp.fitness.fmincon_option1;  % fmincon优化选项1
        parameter_gp.fmincon_option2=gp.fitness.fmincon_option2;  % fmincon优化选项2
        parameter_gp.noise_on=gp.fitness.noise_on;                % 是否启用噪声
        parameter_gp.noise_level=gp.fitness.noise_level;          % 噪声水平
        parameter_gp.ridge_on=gp.fitness.ridge_on;                % 是否启用岭回归
        parameter_gp.ridge_k=gp.fitness.ridge_k;                  % 岭回归系数
        parameter_gp.theta_end_limit=gp.fitness.theta_end_limit;   % theta的终止限制
        parameter_gp.num_const=num_const;                          % 约束参数数量
        parameter_gp.numGenes=numGenes;                            % 基因数量
        parameter_gp.c_log10=log(10);                              % 常用常数

        % Multi_R_check初始值设置（用于多R值检验）
        NR = 20;                                                    % R值采样点数
        R_temp =linspace(-1,0.95,NR);                              % R值范围：-1到0.95
        material_i_deltaK_min=0.1;                                 % deltaK最小值
        material_i_deltaK_max=120;                                 % deltaK最大值
        % 在对数空间均匀采样deltaK
        x1_temp = linspace (log10(material_i_deltaK_min), log10(material_i_deltaK_max),100); 
        x1_temp=10.^x1_temp;                                        % 转换回线性空间

        %—— 构造网格并向量化计算 y ——%
        [Rg, x1_temp_g] = meshgrid(R_temp, x1_temp);      % 生成R和deltaK的网格（大小 100×20）
        KMg       = x1_temp_g ./ (1 - Rg);                % 计算Kmax = deltaK/(1-R)

        % 保存网格参数
        parameter_gp.NR=NR;
        parameter_gp.R_temp=R_temp;
        parameter_gp.Rg=Rg;
        parameter_gp.x1_temp_g=x1_temp_g;
        parameter_gp.KMg=KMg;

        %% 计算模型对于所有合并数据的系数（获取初始参数）
        % 调用Fitfun_multidata_merged函数，基于合并数据拟合初始参数
        [EFF_index,PARA0,PARA_for_initial,parameter_gp]=Fitfun_multidata_merged(evalstr_in,gp,parameter_gp);    
        
        % 如果合并数据拟合成功（EFF_index==1），继续对每个数据集进行拟合
        if EFF_index==1
            PARA0_NUM=size(PARA0,1);                                    % 初始参数组合数量
            fitness_all_PARA0=ones(PARA0_NUM,1)*inf;                   % 初始化所有PARA0的适应度（设为无穷大）
            parameters_all_PARA0=cell(PARA0_NUM,1);                     % 存储所有PARA0对应的参数
            
            % 遍历每个初始参数组合PARA0
            for PARA0_i=1:PARA0_NUM  
                parameter_gp.PARA0=PARA0(PARA0_i,:);                    % 设置当前PARA0组合
                fitness_PARA0_i_all_data=ones(Ndata,numGenes+6+num_const)*inf;  % 初始化当前PARA0下所有数据的适应度
                
                %% 第二层循环：遍历每个数据集
                for data_i=1:Ndata
                    % 如果当前约束参数组合已经失效，跳出循环
                    if  check_at_const_parameter_i==0
                        break;
                    end

                    % 提取当前数据集
                    parameter_gp.data_i=multidata{1,data_i};
                    [num_of_data_i,~]=size(multidata{1,data_i}(:,1));
                    parameter_gp.num_of_data_i=num_of_data_i;          % 当前数据集的数据点数量
            
                    %%% 生成DELTA_KTH与KC的所有可能组合，在对数坐标系下均布采样
                    N_ini_k1k2=size(PARA_for_initial,1);               % 初始k1k2组合数量
                    %N_ini_k1k2=1;
            
                    % 初始化存储变量
                    loss_at_Kth_KC_all=zeros(N_ini_k1k2,1)*inf;        % 存储所有k1k2组合的损失值
                    fitness_at_Kth_KC_all=cell(N_ini_k1k2,1);          % 存储所有k1k2组合的适应度信息
                    geneOutputs=cell(N_ini_k1k2,1);                     % 存储基因输出
                    
                    %% 遍历delta_kth与kc（第三层循环的准备）
                    p_temp=parameter_gp;                                % 复制参数结构体
                    evalstr=evalstr_in;                                 % 复制表达式
                    evalstr_test=evalstr_in;                            % 复制测试表达式
                    
                    %% 第三层循环：遍历不同的k1k2初始值组合
                    for Kth_KC_i=1:1                 
                    %parfor Kth_KC_i=1:N_ini_k1k2      % 并行循环（已注释）
                        fitness_temp_at_Kth_KC_i=zeros(1,5+num_const+numGenes);  % 初始化当前组合的适应度
                        check_fitness=1;                                % 适应度检查标志
                        
                        % 从初始参数中提取k1和k2作为非线性优化的起点
                        parameter_K=PARA_for_initial(PARA0_i,1:2);      % [delta_kth, kc]
                        delta_kth=parameter_K(1);                       % 提取delta_kth
                        kc=parameter_K(2);                              % 提取kc
                        theta=PARA_for_initial(PARA0_i,3:end);          % 提取其他参数theta
                        xtrain=[];                                      % 初始化训练输入
                        ytrain=[];                                      % 初始化训练输出

                        %% 重新设定参数范围（根据当前theta值动态调整）
                        theta_gene=theta(2:end-1);                      % 提取基因参数（排除首尾）
                        % 初始化f_L和f_U
                        f_L = zeros(numGenes, 1);
                        f_U = zeros(numGenes, 1);
                        
                        % 根据theta_gene的符号设置f_L和f_U（保持符号一致性）
                        for pp = 1:numGenes
                            if theta_gene(pp) > 0
                                f_L(pp) = 0.1;   % 正元素：下界设为0.1
                                f_U(pp) = 10;   % 正元素：上界设为10
                            elseif theta_gene(pp) < 0
                                f_L(pp) = -10;  % 负元素：下界设为-10
                                f_U(pp) = -0.1;   % 负元素：上界设为-0.1
                            else
                                % 处理a(i)=0的情况（默认设为0）
                                f_L(pp) = 0;
                                f_U(pp) = 0;
                            end
                        end
                        % 设置其他参数边界
                        m_L = 1;                                        % m参数下界
                        m_U = 10;                                       % m参数上界
                        
                        % 构建最终边界向量：[log10(delta_kth); f参数; m参数; delta_kth; kc]
                        LB_orig = [-15; f_L; m_L;   0;   0];           % 下界向量
                        UB_orig = [ -5; f_U; m_U; 500; 500];           % 上界向量
                        parameter_gp.LB_orig=LB_orig;
                        parameter_gp.UB_orig=UB_orig;


            %           %%%delta_Kth/delta_K, delta_K/delta_Kth, Kmax/Kc, Kc/Kmax, R, log10(delta_K)
            %           xtrain=[delta_kth./multidata{1,data_i}(:,1),multidata{1,data_i}(:,1)./delta_kth,multidata{1,data_i}(:,4)./kc,kc./multidata{1,data_i}(:,4),multidata{1,data_i}(:,3)];
                        
                        % 准备训练数据：delta_K/delta_Kth, Kmax/Kc
                        xtrain=[multidata{1,data_i}(:,1)./delta_kth,multidata{1,data_i}(:,4)./kc];
                        ytrain=log10(multidata{1,data_i}(:,2));  % 目标值取对数
                        
                        % 处理表达式字符串，使用正则表达式替换以允许直接访问数据矩阵
                        pat1 = 'x(\d+)';  % 匹配x(数字)模式
                        pat2 = 'c(\d+)';  % 匹配c(数字)模式
                        evalstr = regexprep(evalstr,pat2,'Const_pair_now($1)');  % 将c替换为当前约束值
                        evalstr = regexprep(evalstr,pat1,'xtrain(:,$1)');        % 将x替换为训练数据列
                        p_temp.evalstr2=evalstr;  % 保存处理后的表达式
                        y = ytrain;
                        [numData,~] =size(ytrain);
                        % 设置矩阵存储树输出加上偏置列（全1列）
                        geneOutputs{Kth_KC_i,1} = ones(numData,numGenes+2);
                        
                        %% 初始化检验，如果起点parameter_K就存在基因中有无穷或者复数，就放弃这个参数起点
                        if check_fitness==1
            
                            % 注释掉的代码：初始值检验（已改为直接通过）
                            % if ~run_completed || force_compute_theta
                            %     showfigure=0;
                            %     [loss_cal_optimize_in,~] = Multi_R_check(p_temp,showfigure,Const_pair_now,theta,delta_kth,kc);           
                            % 
                            %     %初始值是否满足要求，满足要求则再进一步优化
                            %     if loss_cal_optimize_in==0
                            %         fitness_temp_at_Kth_KC_i=[Inf,delta_kth,kc,Const_pair_now,zeros(1,numGenes+2),0];                  
                            %     end                   
                            % end  
                            % 初始曲面对所有数据都是满足的，对单条数据肯定也是满足的
                            loss_cal_optimize_in=1;

                            % 如果初始值检验通过，进行参数优化
                            if loss_cal_optimize_in==1  
                   
                                % 设置优化所需的参数
                                p_temp.theta=theta;                      % 当前参数值
                                p_temp.parameter_K=parameter_K;         % k1和k2参数
                                p_temp.iteration_extent=iteration_extent; % 迭代次数
                                p_temp.geneOutputs=geneOutputs{Kth_KC_i,1};  % 基因输出矩阵
                                p_temp.ytrain=ytrain;                    % 训练目标值
          
                                % 调用损失计算和优化函数
                                % 返回：损失值、优化后的delta_kth、kc、theta、使用的优化模式、约束错误标志、通过标志
                                [loss,delta_kth,kc,theta,opti_mode_used,constraint_wrong,pass_index]=loss_cal_optimize(2.34,p_temp,DELTA_KTH,KC);
                                
                                %% 再次对优化后的参数进行检验
                                % showfigure=0;
                                % [pass_index,~] = Multi_R_check(p_temp,showfigure,Const_pair_now,theta,delta_kth,kc);
            
                                % 根据检验结果设置适应度值
                                if pass_index==1
                                    % 如果通过检验，保存损失值和参数
                                    fitness_temp_at_Kth_KC_i=[loss(1),delta_kth,kc,Const_pair_now,theta',opti_mode_used];     % 格式：3+num_const+numGenes+2+1
                                else
                                    % 如果未通过检验，损失值设为无穷大
                                    fitness_temp_at_Kth_KC_i=[inf,delta_kth,kc,Const_pair_now,theta',opti_mode_used];     % 格式：3+num_const+numGenes+2+1
                                end
            
                            end
             
                            %fitness_temp_at_Kth_KC_i=[loss,delta_kth,kc,Const_pair_now,theta',opti_mode_used];     %3+num_const+numGenes+2+1
                        end
                        
                        % 保存当前k1k2组合的结果
                        loss_at_Kth_KC_all(Kth_KC_i,1)=fitness_temp_at_Kth_KC_i(1);
                        fitness_at_Kth_KC_all{Kth_KC_i,1}=fitness_temp_at_Kth_KC_i;    
                        
                        % 如果损失值已经足够小，提前终止循环
                        if fitness_temp_at_Kth_KC_i(1)<gp.fitness.terminate_value
                            break
                        end
                        
                    end  % 第三层循环结束
                    
                    % 将所有k1k2组合的结果转换为矩阵
                    fitness_at_Kth_KC_all_check=cell2mat(fitness_at_Kth_KC_all);
                    % 挑出众多Kth_KC起点下最小的那个loss
                    % 返回回来的loss如果不是inf那么theta_end一定满足要求，这在loss_cal_optimize.m中已经实现了
                    [loss_min_at_Kth_KC_all,pvals]=min(loss_at_Kth_KC_all);      
                    
                    % 如果所有k1k2组合都无法拟合，跳出数据循环
                    if  isinf(loss_min_at_Kth_KC_all)
                        break;  % 当前PARA0的情况下，没法对第data_i条数据取得拟合值，后面的数据没有必要再试了
                    else
                        % 保存当前PARA0组合下，当前data，遍历delta_Kth、Kc后最优的delta_Kth、Kc、theta
                        fitness_PARA0_i_all_data(data_i,1:numGenes+6+num_const)=fitness_at_Kth_KC_all{pvals,1};     
                    end

                    % 
                    % if  isinf(loss_min_at_Kth_KC_all)
                    %     check_at_const_parameter_i=0;          %意味着这条材料没有办法实现theta_end满足要求，后面的材料没必要再试了，抛弃这个常数组
                    % end
                    % if check_at_const_parameter_i==1
                    %     fitness_best_constC_all_data(data_i,1:numGenes+6+num_const)=fitness_at_Kth_KC_all{pvals,1};     %当前const组合下,当前data，遍历delta_Kth，Kc后最优的delta_Kth，Kc，theta
                    % end
                    % 调试输出
                    if p_temp.debug==1
                        data_i 
                    end

                end  % 第二层循环结束（遍历所有数据集）
    
                % 计算当前PARA0下所有数据的综合损失
                loss_multidata=cal_loss_multidata(gp,fitness_PARA0_i_all_data,numGenes,num_const,gp.fitness.lambda);

                % 保存当前PARA0的结果
                %loss_PARA0_i(PARA0_i,1)=loss_multidata;
                parameters_all_PARA0{PARA0_i,1}=fitness_PARA0_i_all_data;  % 保存所有数据集的参数
                fitness_all_PARA0(PARA0_i,1)=loss_multidata(1);            % 保存损失值1（MSE）
                fitness_all_PARA0(PARA0_i,2)=loss_multidata(2);            % 保存损失值2
                fitness_all_PARA0(PARA0_i,3)=loss_multidata(3);            % 保存损失值3
            end  % PARA0循环结束
            
            % 注释掉的代码：根据不同阶段选择不同的损失标准（已改为统一使用MSE）
            % if gp.fitness.gen_count_now<=gp.runcontrol.stage1        
            %     [~,pvals]=min(fitness_all_PARA0(:,1));       %这里是备选方程在所有PARA0情况下最优的loss_multidata作为fitness_out
            % end
            % if (gp.fitness.gen_count_now>gp.runcontrol.stage1)&&(gp.fitness.gen_count_now<=gp.runcontrol.stage2)   
            %     [~,pvals]=min(fitness_all_PARA0(:,2)); 
            % end
            % if gp.fitness.gen_count_now>gp.runcontrol.stage2 
            %     [~,pvals]=min(fitness_all_PARA0(:,3));       %pvals代表最优的一组PARA0情况
            % end

            % 选择最优的PARA0组合（最优基准参数的选择就一个标准，MSE最小）
            [~,pvals]=min(fitness_all_PARA0(:,1));      
            fitness_best_PARA0=[fitness_all_PARA0(pvals,1),fitness_all_PARA0(pvals,2),fitness_all_PARA0(pvals,3)];  % 最优PARA0的三种损失值
            All_single_para_PARA0=parameters_all_PARA0{pvals,1};  % 最优PARA0对应的所有数据集参数
            fitness_PARA0=PARA0(pvals,:);                          % 最优PARA0的初始参数

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%可视化,debug用（调试时启用）
            show_figure=0;  % 是否显示和保存图形（0=不显示，1=显示并保存）
            if show_figure
                plot_temp=p_temp;
                % 对每个数据集进行可视化
                for data_i=1:Ndata   
                    % 提取当前数据集的参数
                    para_data_i=All_single_para_PARA0(data_i,2:end-1);
                    theta=para_data_i(2+num_const+1:end);  % 提取theta参数
                    delta_kth=para_data_i(1);                % 提取delta_kth
                    kc=para_data_i(2);                      % 提取kc
                    plot_temp.data_i=multidata{1,data_i};   % 设置当前数据集
                    
                    % 调用Multi_R_check进行多R值检验并绘图
                    [~,~] = Multi_R_check(plot_temp,show_figure,Const_pair_now,theta,delta_kth,kc);
        
                    % 设置保存文件夹路径
                    saveFolder1 = 'C:\Users\wzy59\Desktop\test250618\test5\Figure_NASGRO_HS_test';
                    if ~exist(saveFolder1, 'dir')
                        mkdir(saveFolder1);
                    end   
                    % 在主文件夹下创建 1、2 两个子文件夹
                    for s = 1:2
                        subFolder = fullfile(saveFolder1, num2str(s));
                        if ~exist(subFolder, 'dir')
                            mkdir(subFolder);
                        end
                    end
        
                    % 获取所有打开的图形句柄
                    figHandles = findall(0, 'Type', 'figure');
                    [~, order] = sort([figHandles.Number]);
                    figHandles = figHandles(order);
          
                    % 遍历并分别保存 .png 和 .fig 格式
                    for k = 1:numel(figHandles)
                        fh = figHandles(k);
                        
                        % 目标子文件夹
                        subFolder = fullfile(saveFolder1, num2str(k));
                        
                        % 构造基础文件名，例如 "5_1"
                        baseName = sprintf('%d_%d', data_i, k);            
                        % 完整路径：PNG
                        pngFile = fullfile(subFolder, [baseName, '.png']);
                        saveas(fh, pngFile);           % 保存为 PNG
                        
                        % 完整路径：FIG
                        figFile = fullfile(subFolder, [baseName, '.fig']);
                        savefig(fh, figFile);          % 保存为 FIG
                    end  
                    
                    % 保存参数数据
                    saveFolder2 = 'C:\Users\wzy59\Desktop\test250618\test5\Parameters_NASGRO_HS_test';
                    if ~exist(saveFolder2, 'dir')
                        mkdir(saveFolder2);
                    end 
                    fileName2 = fullfile(saveFolder2, sprintf('%d.mat', data_i));
                    save(fileName2, 'fitness_at_Kth_KC_all_check');
                    close all;  % 关闭所有图形

                end
            end
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%可视化,debug用

        else
            % 如果合并数据拟合失败，设置默认值
            fitness_best_PARA0=[inf,inf,inf];
            All_single_para_PARA0={inf,numGenes,num_const,evalstr_in};
            fitness_PARA0=[];

        end
        
        % 保存当前约束参数组合的结果
        %loss_multidata=cal_loss_multidata(gp,fitness_best_constC_all_data,numGenes,num_const,gp.fitness.lambda);
        fitness_all_const{1}{const_parameter_i,1}=All_single_para_PARA0;  % 保存所有数据集的参数
        fitness_all_const{2}(const_parameter_i,1)=fitness_best_PARA0(1);   % 保存损失值1
        fitness_all_const{3}(const_parameter_i,1)=fitness_best_PARA0(2);   % 保存损失值2
        fitness_all_const{4}(const_parameter_i,1)=fitness_best_PARA0(3);   % 保存损失值3
        fitness_all_const{5}{const_parameter_i,1}=fitness_PARA0;           % 保存PARA0参数
  
    end  % 第一层循环结束（遍历所有约束参数组合）
    
    %% 根据进化阶段选择最优的约束参数组合
    if gp.fitness.gen_count_now<=gp.runcontrol.stage1        
        % 阶段1：使用损失值1选择最优约束参数组合
        [~,C_pvals]=min(fitness_all_const{1,2});       % 这里是备选方程在所有给定C参数范围内最优的loss_multidata作为fitness_out
    end
    if (gp.fitness.gen_count_now>gp.runcontrol.stage1)&&(gp.fitness.gen_count_now<=gp.runcontrol.stage2)   
        % 阶段2：使用损失值2选择最优约束参数组合
        [~,C_pvals]=min(fitness_all_const{1,3}); 
    end
    if gp.fitness.gen_count_now>gp.runcontrol.stage2 
        % 阶段3：使用损失值3选择最优约束参数组合
        [~,C_pvals]=min(fitness_all_const{1,4});       % C_pvals代表最优的一组C常数，当C只取值1的时候，就只有一组参数，没有影响
    end
    
    % 提取最优约束参数组合的结果
    fitness_out=[fitness_all_const{2}(C_pvals,1),fitness_all_const{3}(C_pvals,1),fitness_all_const{4}(C_pvals,1)];  % 三种损失值
    fitness_temp=fitness_all_const{1}{C_pvals,1};      % 所有数据集的参数
    fitness_PARA0=fitness_all_const{5}{C_pvals,1};     % PARA0参数

    % 根据损失模式构建返回结构
    if loss_mode>=2
       % 损失模式>=2：包含方程表达式eq
       fitness_return={fitness_temp,numGenes,num_const,eq,evalstr_in,fitness_PARA0};
    elseif loss_mode==1
       % 损失模式=1：不包含方程表达式eq
       fitness_return={fitness_temp,numGenes,num_const,evalstr_in,fitness_PARA0};
    end
    gp.fitness.returnvalues = fitness_return;  % 保存返回结果到gp结构体

else
    % 如果diff_index无效，返回默认值
    fitness_out=[inf,inf,inf];
    fitness_return={inf,numGenes,num_const,evalstr_in,evalstr_in,[0 0 0]};
    gp.fitness.returnvalues = fitness_return;

end    %% diff_index有效才执行上述代码




end

%foldername=['W:\test231012' sprintf('%d',inpChoice)];



