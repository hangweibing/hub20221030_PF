

%%
    % 2026-1-27 改动日志：
    % (1)generate_stress_spectra用于生成不同应力比的载荷块谱AsteixSpectraData_fake.mat
    % (2)在gen_synthetic_data_spectrum.m中重新运行，生成新的检查数据，应力-载荷转化系数设置为60
    % (3）pred_a_N中的应力-载荷转化系数也需要做对应的修改
%%
clc
clear all
close all

% load('pop_for_debug_250722_test11.mat')
load('AM-TC4-GRO_260123_merged.mat')
%load('parameter_gp_AM_TC4_GRO.mat')
load('parameter_gp_AM-TC4-GRO_combine.mat')

data = load('AsteixSpectraData_fake.mat', 'spectra');
parameter_gp.spectra = data.spectra;

options_fmincon_1 = optimoptions('fmincon', ...
    'Algorithm', 'interior-point', ...       % 支持边界约束的算法
    'SpecifyObjectiveGradient',false, ...   % <-- 开启
    'OptimalityTolerance', 1e-8, ...
    'StepTolerance', 1e-10, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 2000, ...
    'MaxFunctionEvaluations', 2000, ...
    'Display', 'off');

% options_fmincon_1 = optimoptions('fmincon', ...
%     'Algorithm', 'interior-point', ...       % 支持边界约束的算法
%     'SpecifyObjectiveGradient',true, ...   % <-- 开启
%     'OptimalityTolerance', 1e-10, ...
%     'StepTolerance', 1e-12, ...
%     'FunctionTolerance', 1e-10, ...
%     'MaxIterations', 1000, ...
%     'MaxFunctionEvaluations', 1000, ...
%     'Display', 'off');



options_fmincon_2 = optimoptions('fmincon', ...
    'Algorithm', 'sqp', ...       % 支持边界约束的算法
    'SpecifyObjectiveGradient',false, ...
    'OptimalityTolerance', 1e-6, ...
    'StepTolerance', 1e-10, ...
    'FunctionTolerance', 1e-6, ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 1000, ...
    'Display', 'off');


parameter_gp.fmincon_option1=options_fmincon_1;
parameter_gp.fmincon_option2=options_fmincon_2;

t_check = [110.3872 220.7745 331.1617 441.5489 551.9362 662.3234 772.7106 883.0979 993.4851];
z = [10.9043 11.0747 12.2835 13.0287 14.2125 16.4101 17.4748 25.1182 30.2851];



parameter_gp.PARA0=pop_now{1,6};

% --- Added for k2 fix ---
parameter_gp.fix_k2 = true; % 固定 k2 开关，手动选择
if parameter_gp.fix_k2
    PARA0_default_temp = pop_now{1,6};
    pop_all_temp = pop_now{1,1};
    % 选择 PARA0_default(2) 和 pop_all(:, 3) 平均值中的较小者
    k2_fixed_val = min(PARA0_default_temp(2), mean(pop_all_temp(:, 3)));
    k2_fixed_val = 6.6750;
    fprintf('k2 Fixed Switch is ON. Fixed Value: %.4f\n', k2_fixed_val);
    parameter_gp.k2_fixed_val = k2_fixed_val;
end
% -----------------------

parameter_gp.data_a_N=[t_check',z'];
[num_of_data_a_N,~]=size(parameter_gp.data_a_N(:,1));
parameter_gp.num_of_data_a_N=num_of_data_a_N;

evalstr=parameter_gp.evalstr1;
evalstr_test=parameter_gp.evalstr1;


fitness_temp_at_Kth_KC_i=zeros(1,5+parameter_gp.num_const+parameter_gp.numGenes);
check_fitness=1;
parameter_K=parameter_gp.PARA0(1:2);                        %作为非线性优化不同的起点
delta_kth=parameter_K(1);
kc=parameter_K(2);
theta=parameter_gp.PARA0(3:end);
xtrain=[];
ytrain=[];


%% 重新设定参数范围
theta_gene=theta(2:end-1);
% 初始化f_L和f_U
f_L = zeros(parameter_gp.numGenes, 1);
f_U = zeros(parameter_gp.numGenes, 1);

% 根据a的符号设置f_L和f_U
for pp = 1:parameter_gp.numGenes
    if theta_gene(pp) > 0
        f_L(pp) = 0.1;   % 正元素：下界设为0.1
        f_U(pp) = 10;   % 正元素：上界设为5
    elseif theta_gene(pp) < 0
        f_L(pp) = -10;  % 负元素：下界设为-5
        f_U(pp) = -0.1;   % 负元素：上界设为-0.1
    else
        % 处理a(i)=0的情况（默认设为0）
        f_L(pp) = 0;
        f_U(pp) = 0;
    end
end
% 设置其他参数边界
m_L = 1;
m_U = 10;

% 构建最终边界向量
LB_orig = [-15; f_L; m_L;   0;   0];
UB_orig = [ -5; f_U; m_U; 500; 500];
parameter_gp.LB_orig=LB_orig;
parameter_gp.UB_orig=UB_orig;

xtrain=[15./delta_kth,20./kc];
ytrain=0.7;
% process evalstr with regex to allow direct access to data matrices
pat1 = 'x(\d+)';
pat2 = 'c(\d+)';
evalstr = regexprep(evalstr,pat2,'Const_pair_now($1)');
evalstr = regexprep(evalstr,pat1,'xtrain(:,$1)');
parameter_gp.evalstr2=evalstr;
y = ytrain;
[numData,~] =size(ytrain);
%set up a matrix to store the tree outputs plus a bias column of ones

parameter_gp.theta=theta;
parameter_gp.parameter_K=parameter_K;
parameter_gp.ytrain=ytrain;

% [loss,delta_kth,kc,theta,opti_mode_used,constraint_wrong,pass_index]=loss_cal_optimize_for_PF(parameter_gp);
kc
theta
N_starts = 10; % 手动设置起始点数量
[loss,delta_kth,kc,theta,opti_mode_used,constraint_wrong,pass_index]=loss_cal_optimize_for_PF_multistart(parameter_gp, pop_now, N_starts);


