function [fit,G] = call_for_fminunc_3_4(z)
%提取全局变量

    global global_p_temp
    global search_index
    global loss_MSE  
    global loss_min            
    global z_now
    global theta_now
    global z_min
    global theta_min
     
    PARA0=global_p_temp.PARA0;
    numGenes = numel(global_p_temp.evalstr2);
    Const_pair_now=global_p_temp.Const_pair_now;
    f0=global_p_temp.theta;
    k0=global_p_temp.parameter_K;
    omega0=[f0,k0];                       %放缩系数

    K_merged=PARA0(1:2);
    numConsts = numel(Const_pair_now);    
    f_merged=PARA0(2+numConsts+1:end);
    omega_merged=[f_merged,K_merged];


    [numData,~] =size(global_p_temp.ytrain);
    delta_K=global_p_temp.data_i(:,1);
    Kmax=global_p_temp.data_i(:,4);
    eq=global_p_temp.eq;
    ytrain=global_p_temp.ytrain;
    num_p=size(global_p_temp.diff_omega,2);
    num_f=num_p-2;

    %% 更新k,f
    omega=omega0'.*z';
    f=omega(1:num_f,1);
    k=omega(num_f+1:num_p,1);        
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

    %% 计算fit, G  
    theta=f;
    delta_kth=k(1);
    kc=k(2);    
    showfigure=0;
    %[~,delta_K_test,g_y,R] = Multi_R_check(global_p_temp,showfigure,Const_pair_now,theta,delta_kth,kc);
    [~,delta_K_test,g_y,R] = Multi_R_check(global_p_temp.data_i,global_p_temp.evalstr_test_fun,global_p_temp.diff_delta_K_fun,showfigure,Const_pair_now,theta,k(1),k(2));
    N_gy=size(g_y,1);
    g_y_alpha=global_p_temp.g_y_alpha;
    g_y_lambda1=global_p_temp.g_y_lambda1;
    g_y_lambda2=global_p_temp.g_y_lambda2;   
    % g_y=-10:0.01:5
    % ReLu_g_y=g_y_alpha*log(1+exp(-1*g_y/g_y_alpha));
    % figure
    % plot(g_y,ReLu_g_y)
    rou_g_y=1/g_y_alpha*log(1+exp(-1*g_y*g_y_alpha));     %代入的是负梯度
    fit_punish_1=g_y_lambda1*1/2/N_gy*sum(rou_g_y.^2);
    y_pre=eval(eq);                                %这个时候k,f是更新过的
    fit_MSE=1/2/numData*sum((y_pre-ytrain).^2);
    fit_punish_2=g_y_lambda2*1/2/num_p*sum((omega'-omega_merged).^2);

    fit=fit_MSE+fit_punish_1+fit_punish_2;
  
    if nargout > 1 % gradient required
        diff_omega=global_p_temp.diff_omega;
        diff_omega_OF_deq=global_p_temp.diff_omega_OF_deq;
        pat_test = 'delta_K';
        diff_omega_OF_deq = regexprep(diff_omega_OF_deq,pat_test,'delta_K_test');

        for i=1:num_p
            %%%求loss对所有参数的梯度，注意，迭代的是参数z，所有后面要乘以一项omega0(i)
            dM_dfi=eval(diff_omega{i}).*omega0(i); 
            ddM_dfi=eval(diff_omega_OF_deq{i})*log(10).*delta_K_test.*omega0(i);
            G1=1/numData*sum( (y_pre-ytrain).*dM_dfi );
            G2=-1/N_gy*g_y_lambda1*sum( rou_g_y./(1+exp(g_y*g_y_alpha)).*ddM_dfi );
            if i<=num_f
                G3=g_y_lambda2/num_p*(omega(i)-omega_merged(i)).*omega0(i);
            else
                G3=0;
            end

            G(i,1)=G1+G2+G3;              %对参数i的批梯度
        end
    end

    search_index=search_index+1;
    loss_MSE(search_index)=fit;
    z_now=z;
    theta_now=theta;
    if isreal(fit) && all(fit < loss_min) 
        loss_min=fit;
        z_min=z;
        theta_min=theta;
    end


end

