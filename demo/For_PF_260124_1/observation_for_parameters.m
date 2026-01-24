
clc
clear all
close all
%% ===================================================================
%% 随机数种子设置
%% ===================================================================
SIM_SEED = 2023;
rng(SIM_SEED);

load('pop_for_debug_250722_test11.mat')
load('AM-TC4-GRO_260123.mat')
load('parameter_gp_AM_TC4_GRO.mat')

data = load('AsteixSpectraData.mat', 'spectra');
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

t_check = [202.5475  405.0949  607.6424  810.1899 1012.7374 1215.2848 1417.8323 1620.3798 1822.9273];
z = [10.5573 11.1964 11.9421 12.8326 13.9283 15.3444 17.3190 20.5362 30.0371];


parameter_gp.PARA0=pop_now{1,6};

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

N_starts = 10; % 手动设置起始点数量
[loss,delta_kth,kc,theta,opti_mode_used,constraint_wrong,pass_index]=loss_cal_optimize_for_PF_multistart(parameter_gp, pop_now, N_starts);




