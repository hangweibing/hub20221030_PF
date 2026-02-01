function [fitness_out,gp]=regressmulti_fitfun_multidata_constraint(evalstr_in,gp)

%defaults in case of early exit
theta=[];ypredtrain=[];fitnessTest=[];ypredtest=[];
r2train=[];r2test=[];r2val=[];geneOutputs=[];geneOutputsTest=[];
geneOutputsVal=[];
%% for parallel
multidata=gp.multidata;
bootSample=gp.userdata.bootSample;
bootSampleSize=gp.userdata.bootSampleSize;
run_completed=gp.state.run_completed;
force_compute_theta=gp.state.force_compute_theta;
iteration_extent=gp.fitness.iteration;
%% 
numGenes = numel(evalstr_in); 
num_const=0;
for i=1:numGenes
    open_sq_br = strfind(evalstr_in{i},'c');
    num_const = num_const+numel(open_sq_br);
end
const_choose=gp.fitness.const_choose;
size_const_choose=size(const_choose,2);
const_all_situations=size_const_choose^(num_const);
const_all=zeros(const_all_situations,num_const);
for i=1:const_all_situations
    for j=1:num_const
        if j==1
            index=mod(i-1,size_const_choose)+1;      
            const_all(i,num_const+1-j)=const_choose(index);
        else
            index=floor(mod(i-1,size_const_choose^j)/(size_const_choose^(j-1)))+1;
            const_all(i,num_const+1-j)=const_choose(index);
        end
    end
end

Ndata=size(multidata,2);
% gp.fitness.num_const=num_const;
% gp.fitness.numGenes=numGenes;
% gp.fitness.Ndata=Ndata;
%% 对备选方程F求导
loss_mode=gp.fitness.loss_mode;
diff_model=gp.fitness.diff_model;
parameter_gp.debug=gp.debug;
if loss_mode>=2
    [diff_index,eq,diff_omega,diff_omega_OF_deq,diff_delta_K, contains_k1,contains_k2]=diff_F(evalstr_in,num_const,diff_model);

    if diff_index==1
    pata = 'k(\d+)';
    patb = 'c(\d+)';
    patc = 'f(\d+)';
    diff_omega = regexprep(diff_omega,pata,'k($1)');
    diff_omega = regexprep(diff_omega,patb,'Const_pair_now($1)');
    diff_omega = regexprep(diff_omega,patc,'f($1)');

    diff_omega_OF_deq = regexprep(diff_omega_OF_deq,pata,'k($1)');
    diff_omega_OF_deq = regexprep(diff_omega_OF_deq,patb,'Const_pair_now($1)');
    diff_omega_OF_deq = regexprep(diff_omega_OF_deq,patc,'f($1)');


    diff_delta_K = regexprep(diff_delta_K,pata,'k($1)');
    diff_delta_K = regexprep(diff_delta_K,patb,'Const_pair_now($1)');
    diff_delta_K = regexprep(diff_delta_K,patc,'f($1)'); 

    eq = regexprep(eq,pata,'k($1)');
    eq = regexprep(eq,patb,'Const_pair_now($1)');
    eq = regexprep(eq,patc,'f($1)');

    parameter_gp.diff_omega=diff_omega;
    parameter_gp.diff_omega_OF_deq=diff_omega_OF_deq;
    parameter_gp.eq=eq;
    parameter_gp.diff_delta_K=diff_delta_K;
    parameter_gp.contains_k1=contains_k1;
    parameter_gp.contains_k2=contains_k2;
%     parameter_gp.dF_dk1=dF_dk1;
%     parameter_gp.dF_dk2=dF_dk2;

    end
end

%% 基于数据的差分（用于算mode2下的loss）
if loss_mode>=10  %这部分代码暂未用到
    if diff_model==2.1
        DELTA_y_DELTA_k1=cell(1,Ndata);
        DELTA_y_DELTA_k2=cell(1,Ndata);
        R=cell(1,Ndata);
        k1=cell(1,Ndata);
        k2=cell(1,Ndata);
        for data_j=1:Ndata
            D_i=multidata{1,data_j};
            [num_of_data_i,~]=size(D_i(:,1));
            R{data_j}=D_i(1,3);
            DELTA_y_DELTA_k1{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   D_i(2:num_of_data_i,1)-D_i(1:num_of_data_i-1,1)   );
            DELTA_y_DELTA_k2{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   D_i(2:num_of_data_i,4)-D_i(1:num_of_data_i-1,4)   );
            k1{data_j}=0.5*(  D_i(1:num_of_data_i-1,1) +D_i(2:num_of_data_i,1)  );
            k2{data_j}=0.5*(  D_i(1:num_of_data_i-1,4) +D_i(2:num_of_data_i,4)  );
        end

    elseif diff_model==2.2
        DELTA_y_DELTA_z1=cell(1,Ndata);
        DELTA_y_DELTA_z2=cell(1,Ndata);
        R=cell(1,Ndata);
        z1=cell(1,Ndata);
        z2=cell(1,Ndata);
        for data_j=1:Ndata
            D_i=multidata{1,data_j};
            [num_of_data_i,~]=size(D_i(:,1));
            R{data_j}=D_i(1,3);
            DELTA_y_DELTA_z1{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   lg(D_i(2:num_of_data_i,1))-lg(D_i(1:num_of_data_i-1,1))   );
            DELTA_y_DELTA_z2{data_j}=(   lg(D_i(2:num_of_data_i,2))-lg(D_i(1:num_of_data_i-1,2))   )./(   lg(D_i(2:num_of_data_i,4))-lg(D_i(1:num_of_data_i-1,4))   );
            z1{data_j}=0.5*(  lg(D_i(1:num_of_data_i-1,1)) +lg(D_i(2:num_of_data_i,1))  )  ;
            z2{data_j}=0.5*(  lg(D_i(1:num_of_data_i-1,4)) +lg(D_i(2:num_of_data_i,4))  )  ;
        end

    end
end

if diff_index==1
    if contains_k1&&contains_k2
       k1k2_Correlation=k1k2_Correlation_test(evalstr_in);           %k1k2_Correlation意味着k1,k2是否存在比例关系，如果只有k1或k2就不用计算
       parameter_gp.k1k2_Correlation=k1k2_Correlation;
    end

    %% 备选方程参数拟合
    if contains_k2
        KC=linspace(log10(gp.fitness.k1_k2_ini(1)),log10(gp.fitness.k1_k2_ini(2)),gp.fitness.k1_k2_num);  
    else
        KC=log10(1);    
    end
    if contains_k1
        DELTA_KTH=linspace(log10(gp.fitness.k1_k2_ini(1)),log10(gp.fitness.k1_k2_ini(2)),gp.fitness.k1_k2_num); 
    else
        DELTA_KTH=log10(1);
    end
    if contains_k1&&contains_k2
        if k1k2_Correlation
            DELTA_KTH=log10(1);
        end    
    end
                           
    DELTA_KTH=10.^DELTA_KTH;
    KC=10.^KC;   
end


if diff_index==1
 
    %%%参数范围设置
    f_L=-5*ones(numGenes,1);            %除了paris部分的两个参数，其他的参数施加统一约束
    f_U= 5*ones(numGenes,1);            %除了paris部分的两个参数，其他的参数施加统一约束
    m_L=1;
    m_U=10;
    LB_orig = [-15; f_L; m_L;   0;   0];
    UB_orig = [ -5; f_U; m_U; 300; 300];   
    parameter_gp.LB_orig=LB_orig;
    parameter_gp.UB_orig=UB_orig;

    parameter_gp.g_y_alpha_merged=gp.fitness.g_y_alpha_merged;       
    parameter_gp.g_y_lambda1_merged=gp.fitness.g_y_lambda1_merged;

    parameter_gp.g_y_alpha_single=gp.fitness.g_y_alpha_single;       
    parameter_gp.g_y_lambda1_single=gp.fitness.g_y_lambda1_single;
    parameter_gp.g_y_lambda2_single=gp.fitness.g_y_lambda2_single;       


    %% 第一层循环
    for const_parameter_i=1:const_all_situations 

        check_at_const_parameter_i=1;                %若check_at_const_parameter_i=0（发现了loss=inf的情况），则抛弃方程在这个常数组下的探索
        Const_pair_now=const_all(const_parameter_i,:);
        parameter_gp.Const_pair_now=Const_pair_now;
        parameter_gp.evalstr1=evalstr_in;
        parameter_gp.fmincon_option1=gp.fitness.fmincon_option1;
        parameter_gp.fmincon_option2=gp.fitness.fmincon_option2;
        parameter_gp.noise_on=gp.fitness.noise_on;
        parameter_gp.noise_level=gp.fitness.noise_level; % 噪声水平
        parameter_gp.ridge_on=gp.fitness.ridge_on;        
        parameter_gp.ridge_k=gp.fitness.ridge_k;
        parameter_gp.theta_end_limit=gp.fitness.theta_end_limit;
        parameter_gp.num_const=num_const;
        parameter_gp.numGenes=numGenes;
        parameter_gp.c_log10=log(10);

        % Multi_R_check初始值
        NR = 20;
        R_temp =linspace(-1,0.95,NR);
        material_i_deltaK_min=0.1;
        material_i_deltaK_max=120; 
        x1_temp = linspace (log10(material_i_deltaK_min), log10(material_i_deltaK_max),100); 
        x1_temp=10.^x1_temp;

        %—— 构造网格并向量化计算 y ——%
        [Rg, x1_temp_g] = meshgrid(R_temp, x1_temp);      % 大小 400×40
        KMg       = x1_temp_g ./ (1 - Rg);

        parameter_gp.NR=NR;
        parameter_gp.R_temp=R_temp;
        parameter_gp.Rg=Rg;
        parameter_gp.x1_temp_g=x1_temp_g;
        parameter_gp.KMg=KMg;

        %% 计算模型对于所有合并数据的系数
        [EFF_index,PARA0,PARA_for_initial,parameter_gp]=Fitfun_multidata_merged(evalstr_in,gp,parameter_gp);    
        if EFF_index==1
            PARA0_NUM=size(PARA0,1);
            fitness_all_PARA0=ones(PARA0_NUM,1)*inf;
            parameters_all_PARA0=cell(PARA0_NUM,1);
            for PARA0_i=1:PARA0_NUM  
                parameter_gp.PARA0=PARA0(PARA0_i,:);  
                fitness_PARA0_i_all_data=ones(Ndata,numGenes+6+num_const)*inf;                
                %% 第二层循环    
                for data_i=1:Ndata                                             %第二层循环
                    if  check_at_const_parameter_i==0
                        break;
                    end

                    parameter_gp.data_i=multidata{1,data_i};
                    [num_of_data_i,~]=size(multidata{1,data_i}(:,1));
                    parameter_gp.num_of_data_i=num_of_data_i;        
                    %%%生成DELTA_KTH与KC的所有可能组合，在对数坐标系下均布采样
            
                    N_ini_k1k2=size(PARA_for_initial,1);
                    %N_ini_k1k2=1;
            
                    loss_at_Kth_KC_all=zeros(N_ini_k1k2,1)*inf;
                    fitness_at_Kth_KC_all=cell(N_ini_k1k2,1);
                    geneOutputs=cell(N_ini_k1k2,1);       
                    %% 遍历delta_kth与kc

                    p_temp=parameter_gp;
                    evalstr=evalstr_in;
                    evalstr_test=evalstr_in;     
                    %% 第三层循环          
                    for Kth_KC_i=1:1                 
                    %parfor Kth_KC_i=1:N_ini_k1k2      
                        fitness_temp_at_Kth_KC_i=zeros(1,5+num_const+numGenes);
                        check_fitness=1;
                        parameter_K=PARA_for_initial(PARA0_i,1:2);                        %作为非线性优化不同的起点
                        delta_kth=parameter_K(1);
                        kc=parameter_K(2);
                        theta=PARA_for_initial(PARA0_i,3:end);            
                        xtrain=[];
                        ytrain=[];


                        %% 重新设定参数范围
                        theta_gene=theta(2:end-1);
                        % 初始化f_L和f_U
                        f_L = zeros(numGenes, 1);
                        f_U = zeros(numGenes, 1);
                        
                        % 根据a的符号设置f_L和f_U
                        for pp = 1:numGenes
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


            %           %%%delta_Kth/delta_K, delta_K/delta_Kth, Kmax/Kc, Kc/Kmax, R, log10(delta_K)
            %           xtrain=[delta_kth./multidata{1,data_i}(:,1),multidata{1,data_i}(:,1)./delta_kth,multidata{1,data_i}(:,4)./kc,kc./multidata{1,data_i}(:,4),multidata{1,data_i}(:,3)];
                        %delta_K/delta_Kth, Kmax/Kc, log10(delta_K)
                        xtrain=[multidata{1,data_i}(:,1)./delta_kth,multidata{1,data_i}(:,4)./kc];
                        ytrain=log10(multidata{1,data_i}(:,2));
                        % process evalstr with regex to allow direct access to data matrices
                        pat1 = 'x(\d+)';
                        pat2 = 'c(\d+)';
                        evalstr = regexprep(evalstr,pat2,'Const_pair_now($1)');
                        evalstr = regexprep(evalstr,pat1,'xtrain(:,$1)');
                        p_temp.evalstr2=evalstr;
                        y = ytrain;
                        [numData,~] =size(ytrain);
                        %set up a matrix to store the tree outputs plus a bias column of ones
                        geneOutputs{Kth_KC_i,1} = ones(numData,numGenes+2);
                        %% 初始化检验，如果起点parameter_K就存在基因中有无穷或者复数，就放弃这个参数起点。
            
                        if check_fitness==1
            
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

                            if loss_cal_optimize_in==1  
                   
                                p_temp.theta=theta;
                                p_temp.parameter_K=parameter_K;
                                p_temp.iteration_extent=iteration_extent;
                                p_temp.geneOutputs=geneOutputs{Kth_KC_i,1};
                                p_temp.ytrain=ytrain;
          
                                [loss,delta_kth,kc,theta,opti_mode_used,constraint_wrong,pass_index]=loss_cal_optimize(2.34,p_temp,DELTA_KTH,KC);
                                %%再次对优化后的参数进行检验
                                % showfigure=0;
                                % [pass_index,~] = Multi_R_check(p_temp,showfigure,Const_pair_now,theta,delta_kth,kc);
            
                                if pass_index==1
                                    fitness_temp_at_Kth_KC_i=[loss(1),delta_kth,kc,Const_pair_now,theta',opti_mode_used];     %3+num_const+numGenes+2+1
                                else
                                    fitness_temp_at_Kth_KC_i=[inf,delta_kth,kc,Const_pair_now,theta',opti_mode_used];     %3+num_const+numGenes+2+1
                                end
            
                            end
             
                            %fitness_temp_at_Kth_KC_i=[loss,delta_kth,kc,Const_pair_now,theta',opti_mode_used];     %3+num_const+numGenes+2+1
                        end
                        loss_at_Kth_KC_all(Kth_KC_i,1)=fitness_temp_at_Kth_KC_i(1);
                        fitness_at_Kth_KC_all{Kth_KC_i,1}=fitness_temp_at_Kth_KC_i;    
                        if fitness_temp_at_Kth_KC_i(1)<gp.fitness.terminate_value
                            break
                        end
                        
                    end 
                    fitness_at_Kth_KC_all_check=cell2mat(fitness_at_Kth_KC_all);
                    [loss_min_at_Kth_KC_all,pvals]=min(loss_at_Kth_KC_all);      %挑出众多Kth_KC起点下最小的那个loss。返回回来的loss如果不是inf那么theta_end一定满足要求，这在loss_cal_optimize.m中已经实现了
                    if  isinf(loss_min_at_Kth_KC_all)
                        break;              %当前PARA0的情况下，没法对第data_i条数据取得拟合值 ,后面的数据没有必要再试了
                    else
                        fitness_PARA0_i_all_data(data_i,1:numGenes+6+num_const)=fitness_at_Kth_KC_all{pvals,1};     %当前PARA0组合下,当前data，遍历delta_Kth，Kc后最优的delta_Kth，Kc，theta
                    end

                    % 
                    % if  isinf(loss_min_at_Kth_KC_all)
                    %     check_at_const_parameter_i=0;          %意味着这条材料没有办法实现theta_end满足要求，后面的材料没必要再试了，抛弃这个常数组
                    % end
                    % if check_at_const_parameter_i==1
                    %     fitness_best_constC_all_data(data_i,1:numGenes+6+num_const)=fitness_at_Kth_KC_all{pvals,1};     %当前const组合下,当前data，遍历delta_Kth，Kc后最优的delta_Kth，Kc，theta
                    % end
                    if p_temp.debug==1
                        data_i 
                    end


                end 
    
                loss_multidata=cal_loss_multidata(gp,fitness_PARA0_i_all_data,numGenes,num_const,gp.fitness.lambda);

                %loss_PARA0_i(PARA0_i,1)=loss_multidata;
                parameters_all_PARA0{PARA0_i,1}=fitness_PARA0_i_all_data;
                fitness_all_PARA0(PARA0_i,1)=loss_multidata(1);   
                fitness_all_PARA0(PARA0_i,2)=loss_multidata(2);   
                fitness_all_PARA0(PARA0_i,3)=loss_multidata(3); 
            end 
            % if gp.fitness.gen_count_now<=gp.runcontrol.stage1        
            %     [~,pvals]=min(fitness_all_PARA0(:,1));       %这里是备选方程在所有PARA0情况下最优的loss_multidata作为fitness_out
            % end
            % if (gp.fitness.gen_count_now>gp.runcontrol.stage1)&&(gp.fitness.gen_count_now<=gp.runcontrol.stage2)   
            %     [~,pvals]=min(fitness_all_PARA0(:,2)); 
            % end
            % if gp.fitness.gen_count_now>gp.runcontrol.stage2 
            %     [~,pvals]=min(fitness_all_PARA0(:,3));       %pvals代表最优的一组PARA0情况
            % end

            [~,pvals]=min(fitness_all_PARA0(:,1));      %最优基准参数的选择就一个标准，MSE最小
            fitness_best_PARA0=[fitness_all_PARA0(pvals,1),fitness_all_PARA0(pvals,2),fitness_all_PARA0(pvals,3)];
            All_single_para_PARA0=parameters_all_PARA0{pvals,1}; 
            fitness_PARA0=PARA0(pvals,:);

            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%可视化,debug用  
            show_figure=0;
            if show_figure
                plot_temp=p_temp;
                for data_i=1:Ndata   
                    para_data_i=All_single_para_PARA0(data_i,2:end-1);
                    theta=para_data_i(2+num_const+1:end);
                    delta_kth=para_data_i(1);
                    kc=para_data_i(2);
                    plot_temp.data_i=multidata{1,data_i};
                    [~,~] = Multi_R_check(plot_temp,show_figure,Const_pair_now,theta,delta_kth,kc);
        
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
        
                    figHandles = findall(0, 'Type', 'figure');
                    [~, order] = sort([figHandles.Number]);
                    figHandles = figHandles(order);
          
                    % 遍历并分别保存 .png 和 .fig
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
                    saveFolder2 = 'C:\Users\wzy59\Desktop\test250618\test5\Parameters_NASGRO_HS_test';
                    if ~exist(saveFolder2, 'dir')
                        mkdir(saveFolder2);
                    end 
                    fileName2 = fullfile(saveFolder2, sprintf('%d.mat', data_i));
                    save(fileName2, 'fitness_at_Kth_KC_all_check');
                    close all; 

                end
            end
            %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%可视化,debug用

        else
            fitness_best_PARA0=[inf,inf,inf];
            All_single_para_PARA0={inf,numGenes,num_const,evalstr_in};
            fitness_PARA0=[];

        end
        %loss_multidata=cal_loss_multidata(gp,fitness_best_constC_all_data,numGenes,num_const,gp.fitness.lambda);
        fitness_all_const{1}{const_parameter_i,1}=All_single_para_PARA0;
        fitness_all_const{2}(const_parameter_i,1)=fitness_best_PARA0(1);   
        fitness_all_const{3}(const_parameter_i,1)=fitness_best_PARA0(2);   
        fitness_all_const{4}(const_parameter_i,1)=fitness_best_PARA0(3);   
        fitness_all_const{5}{const_parameter_i,1}=fitness_PARA0;
  
    end
    
    if gp.fitness.gen_count_now<=gp.runcontrol.stage1        
        [~,C_pvals]=min(fitness_all_const{1,2});       %这里是备选方程在所有给定C参数范围内最优的loss_multidata作为fitness_out
    end
    if (gp.fitness.gen_count_now>gp.runcontrol.stage1)&&(gp.fitness.gen_count_now<=gp.runcontrol.stage2)   
        [~,C_pvals]=min(fitness_all_const{1,3}); 
    end
    if gp.fitness.gen_count_now>gp.runcontrol.stage2 
        [~,C_pvals]=min(fitness_all_const{1,4});       %pvals代表最优的一组C常数，当C只取值1的时候，就只有一组参数，没有影响
    end
    
    fitness_out=[fitness_all_const{2}(C_pvals,1),fitness_all_const{3}(C_pvals,1),fitness_all_const{4}(C_pvals,1)];
    fitness_temp=fitness_all_const{1}{C_pvals,1};
    fitness_PARA0=fitness_all_const{5}{C_pvals,1};

    if loss_mode>=2
       fitness_return={fitness_temp,numGenes,num_const,eq,evalstr_in,fitness_PARA0};
    elseif loss_mode==1
       fitness_return={fitness_temp,numGenes,num_const,evalstr_in,fitness_PARA0};
    end
    gp.fitness.returnvalues = fitness_return;


else

    fitness_out=[inf,inf,inf];
    fitness_return={inf,numGenes,num_const,evalstr_in,evalstr_in,[0 0 0]};
    gp.fitness.returnvalues = fitness_return;

end    %% diff_index有效才执行上述代码




end

%foldername=['W:\test231012' sprintf('%d',inpChoice)];



