function [fit,G] = rosenbrockwithgrad_3_4_test(z, p)


    PARA0=p.PARA0;
    Const_pair_now=p.Const_pair_now;       
    omega0 = [p.theta(:); p.parameter_K(:)]; 

    K_merged=PARA0(1:2);
    %numConsts = numel(Const_pair_now);    
    f_merged=PARA0(3:end);
    omega_merged=[f_merged,K_merged];

    delta_K=p.data_i(:,1);
    Kmax=p.data_i(:,4);
    ytrain=p.ytrain;
    numData=p.numData;
    num_p = p.num_p;
    num_f = p.num_f;

    %% 参数提取（直接访问结构体字段）

    omega = omega0.*z';
    f = omega(1:num_f,1);
    k = omega(num_f+1:num_p,1); 

    if p.contains_k2==0
        k(2)=k0(2);          %k2不更新
    end
    if p.contains_k1==0
        k(1)=k0(1);          %k1不更新 
    end
    if p.contains_k1 && p.contains_k2
        if p.k1k2_Correlation
            k(1)=k0(1);          %k1不更新 
        end    
    end

    %% 计算fit, G  
    showfigure=0;
    g_y_alpha=p.g_y_alpha_single;
    g_y_lambda1=p.g_y_lambda1_single;
    g_y_lambda2=p.g_y_lambda2_single;

    [~,Test_output]=Multi_R_check(p,showfigure,Const_pair_now,f,k(1),k(2));        %这个时候k,f是更新过的

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

    y_pre=p.eq_fun(delta_K, Const_pair_now, f, k, Kmax); 
    fit_MSE = sum((y_pre - ytrain).^2) / (2*numData);

    fit=fit_MSE+fit_punish_1+fit_punish_2;
    % %% debug用
    % omega = omega0.*z'             
    % fit_MSE
    % fit_punish_1
    % fit_punish_2
    % %% debug用
  
    if nargout > 1 % gradient required
        G1 = zeros(num_p, 1);
        G2 = zeros(num_p, 1);
        G3 = zeros(num_p, 1);
        delta_K_test_1_3 = [Test_output.X1_1(:);Test_output.X1_3];  % 确保列向量
        delta_K_test_2 = Test_output.X1_2;  % 确保列向量
        R_1_3=[Test_output.R_1;Test_output.R_3];

        % 预计算公共部分
        y_diff = (y_pre - p.ytrain) / numData;
        D_rou_2=D_soft_plus_fun(g_y_alpha,-1*Test_output.GY_2);
        exp_term_2=-1*rou_2.*D_rou_2/Test_output.NII;

        D_rou_1_3=D_soft_plus_fun(g_y_alpha,(GY_2_mean-g_y));
        exp_term_1_3=rou_1_3.*D_rou_1_3/Test_output.NI;       

        for i=1:num_p
            %%%求loss对所有参数的梯度，注意，迭代的是参数z，所有后面要乘以一项omega0(i)
            %% loss第一项对参数i的梯度G1   
            dM_dki=p.diff_omega_fun{i}(delta_K, Const_pair_now, f, k, Kmax).*omega0(i);    %用的实验数据点
            G1(i,1)=sum(y_diff .* dM_dki);

            %% loss第二项对参数i的梯度G2
            ddM_dki_II=p.diff_omega_OF_deq_fun{i}(delta_K_test_2, Const_pair_now, f, k, Test_output.R_2)*(p.c_log10).*delta_K_test_2.*omega0(i);
            D_R_physic_2=sum(exp_term_2.*ddM_dki_II);

            ddM_dki_II_mean=mean(ddM_dki_II)*ones((Test_output.NI+Test_output.NIII),1);
            ddM_dki_I_III=p.diff_omega_OF_deq_fun{i}(delta_K_test_1_3, Const_pair_now, f, k, R_1_3)*(p.c_log10).*delta_K_test_1_3.*omega0(i);
            D_R_physic_1_3=sum(exp_term_1_3.*(ddM_dki_II_mean-ddM_dki_I_III));
            G2(i,1)=g_y_lambda1*(D_R_physic_2+D_R_physic_1_3);
            %% loss第三项对参数i的梯度G3
            if i<=num_f
                G3(i,1)=g_y_lambda2/num_p*(omega(i)-omega_merged(i)).*omega0(i);
            else
                G3(i,1)=g_y_lambda2/num_p*(omega(i)-omega_merged(i))/omega(i)/log(10).*omega0(i);
            end
        end
        G=G1+G2+G3;              %对参数i的批梯度
    end



end

