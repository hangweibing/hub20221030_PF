function [fit,G] = call_for_fminunc_3_1(z)
    %提取全局变量
    global global_p_temp
    global search_index
    global loss_MSE  
    global loss_min            
    global z_now
    global theta_now
    global z_min
    global theta_min

    numGenes = numel(global_p_temp.evalstr2);
    Const_pair_now=global_p_temp.Const_pair_now;
    f=global_p_temp.theta';
    k0=global_p_temp.parameter_K;
    %omega0=[f,k];                       %放缩系数
    %omega=omega0';
    [numData,~] =size(global_p_temp.ytrain);
    delta_K=global_p_temp.data_i(:,1);
    Kmax=global_p_temp.data_i(:,4);
    eq=global_p_temp.eq;
    ytrain=global_p_temp.ytrain;
    num_p=numGenes+4;
    num_f=num_p-2;
    %% 更新k
    k=z.*k0;
    if global_p_temp.contains_k2==0
        k(2)=k0(2);          %k2不更新
    end
    if global_p_temp.contains_k1==0
        k(1)=k0(1);          %k1不更新 
    end
    if global_p_temp.contains_k1&&global_p_temp.contains_k2
        if global_p_temp.k1k2_Correlation
            k(1)=k0(1);          %k1不更新 
        end    
    end

    
    [theta,y_pre]=get_f_lr(global_p_temp,k);
    %% 更新f
    f=theta';
    R=global_p_temp.data_i(:,3);
    g_y=eval(global_p_temp.diff_delta_K);
   
    g_y_alpha=100;
    g_y_lambda=1;
    % g_y=-0.02:0.01:0.5
    % ReLu_g_y=g_y_alpha*log(1+exp(-1*g_y/g_y_alpha));
    % figure
    % plot(g_y,ReLu_g_y)
    rou_g_y=1/g_y_alpha*log(1+exp(-1*g_y*g_y_alpha));     %代入的是负梯度
    fit_punish_1=g_y_lambda*sum(rou_g_y.^2);
    fit_MSE=1/2/numData*sum((y_pre-ytrain).^2);
    fit=fit_MSE+fit_punish_1;
    
    if nargout > 1 % gradient required
        for i=1:2
            %%%求loss对参数z1,z2的梯度，注意，迭代的是参数z，所有后面要乘以一项k0(i)
            dM_dki=eval(global_p_temp.diff_omega{num_f+i}).*k0(i);        
            ddM_dki=eval(global_p_temp.diff_omega_OF_deq{num_f+i}).*k0(i);
            G1=1/numData*sum( (y_pre-ytrain).*dM_dki );
            G2=-2*g_y_lambda*sum( rou_g_y./(1+exp(g_y*g_y_alpha)).*ddM_dki );
            G(i,1)=G1+G2;                %对参数i的批梯度
        end
    end
    search_index=search_index+1;
    loss_MSE(search_index)=fit;
    z_now=z;
    theta_now=theta;
    if fit<loss_min
        loss_min=fit;
        z_min=z;
        theta_min=theta;
    end


end

