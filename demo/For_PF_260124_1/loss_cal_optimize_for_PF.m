function [loss_out,delta_kth_out,kc_out,theta_out,opti_mode_used,constraint_wrong,pass_index] = loss_cal_optimize_for_PF(p)

%% 对所有参数采用非线性优化求解，初始值均为1

    %numGenes = p.numGenes;
    Const_pair_now=p.Const_pair_now;
    f=p.theta;
    k=p.parameter_K;
    k0=k;
    f0=f;
    omega0=[f0,k0];                       %放缩系数
    [numData,~] =size(p.ytrain);
    delta_K=p.data_i(:,1);
    Kmax=p.data_i(:,2);
    ytrain=p.ytrain;
    num_p=size(p.diff_omega,2);
    num_f=num_p-2;
    p.num_p=num_p;
    p.num_f=num_f;
    p.numData=numData;
    
    z0=ones(1,num_p);
    
    PARA0=p.PARA0;
    K_merged=PARA0(1:2);
    %numConsts = numel(Const_pair_now);    
    f_merged=PARA0(3:end);
    omega_merged=[f_merged,K_merged];
    
    showfigure=0;
    theta=p.theta;
    
    [~,Test_output] = Multi_R_check(p,showfigure,Const_pair_now,theta,k(1),k(2));
    
    g_y_alpha=p.g_y_alpha_single;       
    g_y_lambda1=p.g_y_lambda1_single; 
    g_y_lambda2=p.g_y_lambda2_single;  
    
    rou_2=Soft_plus_fun(g_y_alpha,-1*Test_output.GY_2);
    R_physic_2=1/2/Test_output.NII*sum( rou_2.^2 );
    
    g_y=[Test_output.GY_1; Test_output.GY_3];
    GY_2_mean=Test_output.GY_2_mean*ones((Test_output.NI+Test_output.NIII),1);
    rou_1_3=Soft_plus_fun(g_y_alpha,(GY_2_mean-g_y));
    R_physic_1_3=1/2/Test_output.NI*sum( rou_1_3.^2 );
    fit_punish_1=g_y_lambda1*(R_physic_1_3+R_physic_2);
    
    omega0_punish_2 = [omega0(1:end-2), log10(omega0(end-1:end))];
    omega_merged_punish_2 = [omega_merged(1:end-2), log10(omega_merged(end-1:end))];
    fit_punish_2=g_y_lambda2*1/2/num_p*sum((omega0_punish_2-omega_merged_punish_2).^2);
    [fit_MSE, a_pred] = pred_a_N(p,f,k,Const_pair_now);
    %y_pre=p.eq_fun(delta_K, Const_pair_now, f, k, Kmax);                                %这个时候k,f是更新过的
    %fit_MSE=1/2/numData*sum((y_pre-ytrain).^2);
    loss_MSE_0=fit_MSE+fit_punish_1+fit_punish_2;
    
    fun = @(z) rosenbrockwithgrad_for_PF(z, p);
    
    % 2. 原始的上下界（在 x 空间里）
    % % 2. 原始的上下界（在 x 空间里）
    LB_orig = p.LB_orig;
    UB_orig = p.UB_orig;     
    
    % 3. 把它们映射到 z 空间：
    ratio   = [LB_orig./omega0',  UB_orig./omega0'];
    LB      = min(ratio,[],2);   % 下界要取二者之小
    UB      = max(ratio,[],2);   % 上界要取二者之大
    
    warning off;
    p.fmincon_option1.Display = 'iter';
    try
        % 这里是可能产生错误的代码
        [z_opt, fval, exitflag, output] = fmincon( ...
            fun,       ... % 目标函数
            z0,        ... % 初始点
            [], [],    ... % 线性不等式 A, b
            [], [],    ... % 线性等式 Aeq, beq
            LB, UB,    ... % 下界、上界
            [],        ... % 非线性约束（这里无）
            p.fmincon_option1);
        omega=z_opt.*omega0;
        loss_out_temp=fval;                 
    
        opti_mode_used=1;
    catch
        % 在发生错误时执行的代码
        try
            % 这里是可能产生错误的代码
            [z_opt, fval, exitflag, output] = fmincon( ...
                fun,       ... % 目标函数
                z0,        ... % 初始点
                [], [],    ... % 线性不等式 A, b
                [], [],    ... % 线性等式 Aeq, beq
                LB, UB,    ... % 下界、上界
                [],        ... % 非线性约束（这里无）
                p.fmincon_option2);
            omega=z_opt.*omega0;
            loss_out_temp=fval; 
    
            opti_mode_used=2;
        catch
            % 在发生错误时执行的代码
            z_opt=z0;
            fval=loss_MSE_0;
    
            omega=z_opt.*omega0;
            loss_out_temp=fval; 
    
            opti_mode_used=0;
        end
    end
    
    warning on;
    
    %omega=z_opt.*omega0
    %loss_out_temp=fval
    %omega=z_min.*omega0;
    
    showfigure=0;
    omega=omega(:);
    
    f=omega(1:num_f,1);
    k=omega(num_f+1:num_p,1);   
    
    % f=[-9.6345; 2.0; -1.0; 2.0];  k=[6.02; 87.3];  
    
    [test_index,Test_output] = Multi_R_check(p,showfigure,Const_pair_now,f,k(1),k(2));
    
    if test_index==1
        rou_2=Soft_plus_fun(g_y_alpha,-1*Test_output.GY_2);
        R_physic_2=1/2/Test_output.NII*sum( rou_2.^2 );
    
        g_y=[Test_output.GY_1; Test_output.GY_3];
        GY_2_mean=Test_output.GY_2_mean*ones((Test_output.NI+Test_output.NIII),1);
        rou_1_3=Soft_plus_fun(g_y_alpha,(GY_2_mean-g_y));
        R_physic_1_3=1/2/Test_output.NI*sum( rou_1_3.^2 );
    
        fit_punish_1=g_y_lambda1*(R_physic_1_3+R_physic_2);
    
        omega_punish_2 = [omega(1:end-2); log10(omega(end-1:end))]';
        omega_merged_punish_2 = [omega_merged(1:end-2), log10(omega_merged(end-1:end))];
        fit_punish_2=g_y_lambda2*1/2/num_p*sum((omega_punish_2-omega_merged_punish_2).^2);
    
        [fit_MSE, ~] = pred_a_N(p,f,k,Const_pair_now); % 使用优化后的参数重新计算MSE

        if isreal(fit_MSE) && isfinite(fit_MSE) && ~isnan(fit_MSE)
            pass_index = 1;
        else
            pass_index = 0;
        end
    
    else
        pass_index = 0;
    end
    
    
    % showfigure=1;
    % [~,Test_output] = Copy_of_Multi_R_check(p,showfigure,Const_pair_now,f,k(1),k(2));
    
    delta_kth_out=k(1);
    kc_out=k(2);            
    theta_out=f;
    
    %% theta_end检验与k1,k2检验
    %由于theta的最后一项对应delta_K，所以对其施加限制，必须落在区间
    constraint_wrong=0;
    if theta_out(end)<p.theta_end_limit(1) 
        pass_index=0;
        constraint_wrong=1;
    elseif theta_out(end)>p.theta_end_limit(2)
        pass_index=0;
        constraint_wrong=2;
    % elseif (delta_kth_out<0)||(kc_out<0)
    %     pass_index=0;
    %     constraint_wrong=3;
    % elseif delta_kth_out>kc_out
    %     pass_index=0;                        
    end
    
    if pass_index==1
        loss_out=[fit_MSE,fit_punish_1,fit_punish_2];
        %loss_out=loss_out_temp;
    
    else
        loss_out=[inf,inf,inf];
    end


end

%load('B:\Desktop\matlab.mat')
%loss_cal_optimize(loss_mode,p)