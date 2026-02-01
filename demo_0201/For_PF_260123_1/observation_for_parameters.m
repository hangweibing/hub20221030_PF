
clc
clear all
close all

load('pop_for_debug_250722_test11.mat')
load('AM-TC4-GRO_260123.mat')
load('parameter_gp_AM_TC4_GRO.mat')

data = load('AsteixSpectraData.mat', 'spectra');
parameter_gp.spectra = data.spectra;

parameter_gp.fmincon_option1.SpecifyObjectiveGradient=false;
parameter_gp.fmincon_option2.SpecifyObjectiveGradient=false;

t_check = [202.5475  405.0949  607.6424  810.1899 1012.7374 1215.2848 1417.8323 1620.3798 1822.9273];
z = [10.5573 11.1964 11.9421 12.8326 13.9283 15.3444 17.3190 20.5362 30.0371];


parameter_gp.PARA0=pop_now{1,6};  
      
parameter_gp.data_a_N=[t_check',z'];     
[num_of_data_a_N,~]=size(parameter_gp.data_a_N(:,1));
parameter_gp.num_of_data_a_N=num_of_data_a_N;        


p_temp=parameter_gp;
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
p_temp.evalstr2=evalstr;
y = ytrain;
[numData,~] =size(ytrain);
%set up a matrix to store the tree outputs plus a bias column of ones

p_temp.theta=theta;
p_temp.parameter_K=parameter_K;
p_temp.ytrain=ytrain;

[loss,delta_kth,kc,theta,opti_mode_used,constraint_wrong,pass_index]=loss_cal_optimize_for_PF(p_temp);




