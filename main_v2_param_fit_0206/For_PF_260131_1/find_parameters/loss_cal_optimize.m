function [loss_out,delta_kth_out,kc_out,theta_out,opti_mode_used,constraint_wrong,pass_index] = loss_cal_optimize(loss_mode,p,DELTA_KTH,KC)

    if loss_mode==1
        %% 基于离散采样的线性回归
        [numData,~] =size(p.ytrain);
        ypredtrain = p.geneOutputs * p.theta;                          
        err = p.ytrain - ypredtrain;
        loss_out=sqrt(((err'*err)/numData));
        theta_out=p.theta';
        delta_kth_out=p.parameter_K(1);
        kc_out=p.parameter_K(2);
        delta_K=p.data_i(:,1);
        ytrain=p.ytrain;
    elseif loss_mode>=2   
        %% 非线性优化求解参数

        if loss_mode==2.1

            [loss_out,delta_kth_out,kc_out,theta_out] = user_defined_optimization_method_all_p(p);
       
        elseif loss_mode==2.2

            [loss_out,delta_kth_out,kc_out,theta_out] = user_defined_optimization_method_k1k2(p);

        elseif loss_mode==2.31
            %% 只对参数k1,k2采用非线性优化求解，初始值均为1
            numGenes = numel(p.evalstr2);
            Const_pair_now=p.Const_pair_now;
            f=p.theta';
            k=p.parameter_K;
            k0=k;
            f0=f;
            omega0=[f0,k0];                       %放缩系数
            [numData,~] =size(p.ytrain);
            delta_K=p.data_i(:,1);
            Kmax=p.data_i(:,4);
            num_p=size(p.diff_omega,2);
            num_f=num_p-2;
            loss_MSE_0=sum((eval(p.eq)-p.ytrain).^2)/2/numData;   
            z0 = [1 1];

            global loss_MSE
            loss_MSE=[];
            loss_MSE(1)=loss_MSE_0;      %记录loss的整个演变过程
            global loss_min
            loss_min=[];
            loss_min=loss_MSE_0; 

            global z_now           %now表示跟随当前优化状态的系数
            z_now=[];
            z_now=z0;         
            global theta_now
            theta_now=[];
            theta_now=p.theta;  
            global z_min           %min表示对应当前记录到的loss_min的系数
            z_min=[];
            z_min=z0;         
            global theta_min
            theta_min=[];
            theta_min=p.theta;  

            global search_index
            search_index=1;
            global global_p_temp
            global_p_temp=p;

            loss_now=loss_MSE_0;
            fun = @rosenbrockwithgrad_3_1;
            try
                % 这里是可能产生错误的代码
                [~,~,exitflag,output] = fminunc(fun,z0,p.fit_option1);
                opti_mode_used=1;
            catch
                % 在发生错误时执行的代码
                try
                    % 这里是可能产生错误的代码
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;  
                    [~,~,exitflag,output] = fminunc(fun,z0,p.fit_option2);        %经过测试，最后返回的z，会多一次梯度下降，导致下面的z和call_for_fminunc_5.m中最后的z相差一步
                    opti_mode_used=2;
                catch
                    % 在发生错误时执行的代码
                    %loss_min=loss_MSE(end);
                    %z=z0;
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;                      
                    opti_mode_used=0;
                end
            end

            loss_out_temp=loss_min;
            k=z_min.*k0;
            f=theta_min';
         
            theta_out=f';

            clear global_p_temp
            clear search_index
            clear loss_MSE  
            clear loss_min            
            clear z_now
            clear theta_now
            clear z_min
            clear theta_min
            %% theta_end检验与k1,k2检验
            %由于theta的最后一项对应delta_K，所以对其施加限制，必须落在区间
            pass_index=1;
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
                loss_out=loss_out_temp;
            else
                loss_out=inf;
            end

        elseif loss_mode==2.32
            %% 对所有参数采用非线性优化求解，初始值均为1
            % 将原全局变量封装为结构体
            params = struct(...
                'p_temp', p, ...
                'search_index', 0, ...
                'loss_MSE', [], ...
                'loss_min', inf, ...
                'z_min', [], ...
                'theta_min', []);

            numGenes = numel(p.evalstr2);
            Const_pair_now=p.Const_pair_now;
            f=p.theta';
            k=p.parameter_K;
            k0=k;
            f0=f;
            omega0=[f0,k0];                       %放缩系数
            [numData,~] =size(p.ytrain);
            delta_K=p.data_i(:,1);
            Kmax=p.data_i(:,4);
            ytrain=p.ytrain;
            num_p=size(p.diff_omega,2);
            num_f=num_p-2;
            z0=ones(1,num_p);

            showfigure=0;
            theta=p.theta;
            delta_kth=k(1);
            kc=k(2); 
            [~,delta_K_test,g_y,R] = Multi_R_check(p,showfigure,Const_pair_now,theta,delta_kth,kc);

            N_gy=size(g_y,1);

            g_y_alpha=p.g_y_alpha_merged;       
            g_y_lambda1=p.g_y_lambda1_merged;          
            rou_g_y=1/g_y_alpha*log(1+exp(-1*g_y*g_y_alpha));     %代入的是负梯度
            fit_punish_1=g_y_lambda1*1/2/N_gy*sum(rou_g_y.^2);

            vars1 = {'delta_K','Const_pair_now','f','k','Kmax'};
            vars2 = {'delta_K','Const_pair_now','f','k','R'};
            eq_fun = createModelFunction(p.eq, vars1);
            diff_omega_fun=createModelFunction(p.diff_omega, vars1);
            diff_omega_OF_deq_fun=createModelFunction(p.diff_omega_OF_deq, vars2);
            params.eq_fun=eq_fun;
            params.eq_fun=diff_omega_fun;
            params.eq_fun=diff_omega_OF_deq_fun;


            y_pre=eq_fun(delta_K, Const_pair_now, f, k, Kmax);
            %y_pre=eval(p.eq);                                %这个时候k,f是更新过的
            fit_MSE=1/2/numData*sum((y_pre-ytrain).^2);
            loss_MSE_0=fit_MSE+fit_punish_1;

            fun = @(z) rosenbrockwithgrad_3_2_test(z, params);


            global loss_MSE
            loss_MSE=[];
            loss_MSE(1)=loss_MSE_0;      %记录loss的整个演变过程
            global loss_min
            loss_min=[];
            loss_min=loss_MSE_0; 

            global z_now           %now表示跟随当前优化状态的系数
            z_now=[];
            z_now=z0;         
            global theta_now
            theta_now=[];
            theta_now=p.theta;  
            global z_min           %min表示对应当前记录到的loss_min的系数
            z_min=[];
            z_min=z0;         
            global theta_min
            theta_min=[];
            theta_min=p.theta;  

            global search_index
            search_index=1;
            global global_p_temp
            global_p_temp=p;


            options_test = optimoptions('fminunc',...           % 使用fminunc函数进行无约束最小化
                'Algorithm', 'quasi-newton',...            % 选择拟牛顿法作为算法
                'OptimalityTolerance', 1e-6, ...           % 最优性容忍度，较小值提高精度
                'StepTolerance', 1e-10,  ...               % 步长容忍度，非常小，以确保精细搜索
                'FunctionTolerance', 1e-6, ...             % 函数值容忍度，与最优性容忍度一致
                'MaxIterations', 1000, ...                  % 最大迭代次数
                'MaxFunctionEvaluations', 2000, ...        % 最大函数评估次数
                'Display', 'off');            % 显示每次迭代的详细信息


            fun = @rosenbrockwithgrad_3_2_old;
            try
                % 这里是可能产生错误的代码
                %[~,~,exitflag,output] = fminunc(fun,z0,p.fit_option1);
                [~,~,exitflag,output] = fminunc(fun,z0,options_test);
                opti_mode_used=1;
            catch
                % 在发生错误时执行的代码
                try
                    % 这里是可能产生错误的代码
                    loss_MSE=[];
                    z_now=[];
                    z_now=z0;  
                    theta_now=[]; 
                    loss_min=[];      
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    [~,~,exitflag,output] = fminunc(fun,z0,p.fit_option2);        %经过测试，最后返回的z，会多一次梯度下降，导致下面的z和call_for_fminunc_5.m中最后的z相差一步
                    opti_mode_used=2;
                catch
                    % 在发生错误时执行的代码
                    %loss_min=loss_MSE(end);
                    %z=z0;
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;                      
                    opti_mode_used=0;
                end
            end
            loss_out_temp=loss_min;
            omega=z_min.*omega0;
            omega=omega';

            clear global_p_temp  
            clear loss_MSE
            f=omega(1:num_f,1);
            k=omega(num_f+1:num_p,1);      
            delta_kth_out=k(1);
            kc_out=k(2);            
            theta_out=f;

            clear global_p_temp
            clear search_index
            clear loss_MSE  
            clear loss_min            
            clear z_now
            clear theta_now
            clear z_min
            clear theta_min

            %% theta_end检验与k1,k2检验
            %由于theta的最后一项对应delta_K，所以对其施加限制，必须落在区间
            pass_index=1;
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
                loss_out=loss_out_temp;
            else
                loss_out=inf;
            end
        elseif loss_mode==2.33
            %% 对所有参数采用非线性优化求解，初始值均为1

            Const_pair_now=p.Const_pair_now;
            f=p.theta';
            k=p.parameter_K;
            k0=k;
            f0=f;
            omega0=[f0,k0];                       %放缩系数
            [numData,~] =size(p.ytrain);
            delta_K=p.data_i(:,1);
            Kmax=p.data_i(:,4);
            ytrain=p.ytrain;
            num_p=size(p.diff_omega,2);
            num_f=num_p-2;
            p.num_p=num_p;
            p.num_f=num_f;
            p.numData=numData;


            z0=ones(1,num_p);

            showfigure=0;
            theta=p.theta;
            [~,Test_output] = Multi_R_check(p,showfigure,Const_pair_now,theta,k(1),k(2));

            g_y_alpha=p.g_y_alpha_merged;       
            g_y_lambda1=p.g_y_lambda1_merged; 

            rou_2=Soft_plus_fun(g_y_alpha,-1*Test_output.GY_2);
            R_physic_2=1/2/Test_output.NII*sum( rou_2.^2 );

            g_y=[Test_output.GY_1; Test_output.GY_3];
            GY_2_mean=Test_output.GY_2_mean*ones((Test_output.NI+Test_output.NIII),1);
            rou_1_3=Soft_plus_fun(g_y_alpha,(GY_2_mean-g_y));
            R_physic_1_3=1/2/Test_output.NI*sum( rou_1_3.^2 );

            fit_punish_1=g_y_lambda1*(R_physic_1_3+R_physic_2);

            y_pre=p.eq_fun(delta_K, Const_pair_now, f, k, Kmax);

            fit_MSE=1/2/numData*sum((y_pre-ytrain).^2);
            loss_MSE_0=fit_MSE+fit_punish_1;

            fun = @(z) rosenbrockwithgrad_3_2_test(z, p);


            % % 2. 原始的上下界（在 x 空间里）
            LB_orig = p.LB_orig;
            UB_orig = p.UB_orig;            
            % 3. 把它们映射到 z 空间：
            % 
            ratio   = [LB_orig./omega0',  UB_orig./omega0'];
            LB      = min(ratio,[],2);   % 下界要取二者之小
            UB      = max(ratio,[],2);   % 上界要取二者之大

            %fun = @rosenbrockwithgrad_3_2;

            warning off;

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
                %profile viewer
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
            % omega=z_opt.*omega0;
            % loss_out_temp=fval; 


            %% %%%%%%%%%figure for debug

            omega=omega(:);

            f=omega(1:num_f,1);
            k=omega(num_f+1:num_p,1);  


            showfigure=0; 

            [test_index,Test_output] = Multi_R_check(p,showfigure,Const_pair_now,f,k(1),k(2));
            if test_index==1
                rou_2=Soft_plus_fun(g_y_alpha,-1*Test_output.GY_2);
                R_physic_2=1/2/Test_output.NII*sum( rou_2.^2 );
    
                g_y=[Test_output.GY_1; Test_output.GY_3];
                GY_2_mean=Test_output.GY_2_mean*ones((Test_output.NI+Test_output.NIII),1);
                rou_1_3=Soft_plus_fun(g_y_alpha,(GY_2_mean-g_y));
                R_physic_1_3=1/2/Test_output.NI*sum( rou_1_3.^2 );
    
                fit_punish_1=g_y_lambda1*(R_physic_1_3+R_physic_2);
                y_pre=p.eq_fun(delta_K, Const_pair_now, f, k, Kmax);
                fit_MSE=1/2/numData*sum((y_pre-ytrain).^2);
    
                if isreal(fit_MSE) && isfinite(fit_MSE) && ~isnan(fit_MSE)
                    pass_index = 1;
                else
                    pass_index = 0;
                end

            else

                pass_index = 0;

            end


            %% %%%%%%%%%figure for debug
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
                %loss_out=[fit_MSE,fit_punish_1];
                loss_out=[loss_out_temp,fit_MSE,fit_punish_1];
                
            else
                loss_out=[inf,inf,inf];
            end


        elseif loss_mode==2.34
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
            Kmax=p.data_i(:,4);
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

            y_pre=p.eq_fun(delta_K, Const_pair_now, f, k, Kmax);                                %这个时候k,f是更新过的
            fit_MSE=1/2/numData*sum((y_pre-ytrain).^2);
            loss_MSE_0=fit_MSE+fit_punish_1+fit_punish_2;

            fun = @(z) rosenbrockwithgrad_3_4_test(z, p);

            % 2. 原始的上下界（在 x 空间里）
            % % 2. 原始的上下界（在 x 空间里）
            LB_orig = p.LB_orig;
            UB_orig = p.UB_orig;     
            
            % 3. 把它们映射到 z 空间：
            ratio   = [LB_orig./omega0',  UB_orig./omega0'];
            LB      = min(ratio,[],2);   % 下界要取二者之小
            UB      = max(ratio,[],2);   % 上界要取二者之大

            warning off;

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
    
                y_pre=p.eq_fun(delta_K, Const_pair_now, f, k, Kmax);                                %这个时候k,f是更新过的
                fit_MSE=1/2/numData*sum((y_pre-ytrain).^2);
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



        elseif loss_mode==2.4
            %% 对所有参数采用非线性优化求解，初始值均为1
            numGenes = numel(p.evalstr2);
            Const_pair_now=p.Const_pair_now;
            f=p.theta';
            k=p.parameter_K;
            k0=k;
            f0=f;
            omega0=[f0,k0];                       %放缩系数
            [numData,~] =size(p.ytrain);
            delta_K=p.data_i(:,1);
            Kmax=p.data_i(:,4);
            num_p=size(p.diff_omega,2);
            num_f=num_p-2;
            loss_MSE_0=sum((eval(p.eq)-p.ytrain).^2)/2/numData;   
            global loss_MSE
            loss_MSE(1)=loss_MSE_0;

            global global_p_temp
            global_p_temp=p;
            z0=ones(1,num_p);
            fun = @rosenbrockwithgrad_4;
            try
                % 这里是可能产生错误的代码
                [z,loss_now,exitflag,output] = fminunc(fun,z0,p.fit_option1);
            catch
                % 在发生错误时执行的代码
                try
                    % 这里是可能产生错误的代码
                    [z,loss_now,exitflag,output] = fminunc(fun,z0,p.fit_option2);
                catch
                    % 在发生错误时执行的代码
                    loss_now=loss_MSE(1);
                    z=z0;
                end
            end
            clear global_p_temp  
            clear loss_MSE

            if loss_now<loss_MSE(1)
                loss_out_temp=loss_now;
                omega=omega0'.*z';
                f=omega(1:num_f,1);
                k=omega(num_f+1:num_p,1);      

                delta_kth_out=k(1);
                kc_out=k(2);
                theta_out=f;        
            else
                loss_out_temp=loss_MSE(1);
                delta_kth_out=k0(1);
                kc_out=k0(2);
                theta_out=f0';               
            end
            %% theta_end检验
            %由于theta的最后一项对应delta_K，所以对其施加限制，必须落在区间
            if (theta_out(end)>theta_end_limit(1)) &&(theta_out(end)<theta_end_limit(2))
                loss_out=loss_out_temp;
            else
                loss_out=inf;                     %抛弃该方程
            end

        elseif loss_mode==2.5
            %% 只对参数k1,k2采用非线性优化求解，初始值均为1
            numGenes = numel(p.evalstr2);
            Const_pair_now=p.Const_pair_now;
            f=p.theta';
            k=p.parameter_K;
            k0=k;
            f0=f;
            omega0=[f0,k0];                       %放缩系数
            [numData,~] =size(p.ytrain);
            delta_K=p.data_i(:,1);
            Kmax=p.data_i(:,4);
            num_p=size(p.diff_omega,2);
            num_f=num_p-2;
            loss_MSE_0=sum((eval(p.eq)-p.ytrain).^2)/2/numData;   
            z0 = [1 1];

            global loss_MSE
            loss_MSE=[];
            loss_MSE(1)=loss_MSE_0;      %记录loss的整个演变过程
            global loss_min
            loss_min=[];
            loss_min=loss_MSE_0; 

            global z_now           %now表示跟随当前优化状态的系数
            z_now=[];
            z_now=z0;         
            global theta_now
            theta_now=[];
            theta_now=p.theta;  
            global z_min           %min表示对应当前记录到的loss_min的系数
            z_min=[];
            z_min=z0;         
            global theta_min
            theta_min=[];
            theta_min=p.theta;  

            global search_index
            search_index=1;
            global global_p_temp
            global_p_temp=p;



            loss_now=loss_MSE_0;
            fun = @rosenbrockwithgrad_5;
            try
                % 这里是可能产生错误的代码
                [~,~,exitflag,output] = fminunc(fun,z0,p.fit_option1);
                opti_mode_used=1;
            catch
                % 在发生错误时执行的代码
                try
                    % 这里是可能产生错误的代码
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;  
                    [~,~,exitflag,output] = fminunc(fun,z0,p.fit_option2);        %经过测试，最后返回的z，会多一次梯度下降，导致下面的z和call_for_fminunc_5.m中最后的z相差一步
                    opti_mode_used=2;
                catch
                    % 在发生错误时执行的代码
                    %loss_min=loss_MSE(end);
                    %z=z0;
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;                      
                    opti_mode_used=0;
                end
            end

            loss_out_temp=loss_min;
            k=z_min.*k0;
            f=theta_min';
            delta_kth_out=k(1);
            kc_out=k(2);            
            theta_out=f';

            clear global_p_temp
            clear search_index
            clear loss_MSE  
            clear loss_min            
            clear z_now
            clear theta_now
            clear z_min
            clear theta_min
            %% theta_end检验与k1,k2检验
            %由于theta的最后一项对应delta_K，所以对其施加限制，必须落在区间
            pass_index=1;
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
                loss_out=loss_out_temp;
            else
                loss_out=inf;
            end

        elseif loss_mode==2.6
            %% 只对参数k1,k2采用非线性优化求解，初始值均为1（此处为有约束的非线性优化）
            numGenes = numel(p.evalstr2);
            Const_pair_now=p.Const_pair_now;
            f=p.theta';
            k=p.parameter_K;
            k0=k;
            f0=f;
            omega0=[f0,k0];                       %放缩系数
            [numData,~] =size(p.ytrain);
            delta_K=p.data_i(:,1);
            Kmax=p.data_i(:,4);
            num_p=size(p.diff_omega,2);
            num_f=num_p-2;
            loss_MSE_0=sum((eval(p.eq)-p.ytrain).^2)/2/numData;   
            z0 = [1 1];
            lb = [min(DELTA_KTH), min(KC)];  % 下边界
            ub = [max(DELTA_KTH), max(KC)];  % 上边界


            global loss_MSE
            loss_MSE=[];
            loss_MSE(1)=loss_MSE_0;      %记录loss的整个演变过程
            global loss_min
            loss_min=[];
            loss_min=loss_MSE_0; 

            global z_now           %now表示跟随当前优化状态的系数
            z_now=[];
            z_now=z0;         
            global theta_now
            theta_now=[];
            theta_now=p.theta;  
            global z_min           %min表示对应当前记录到的loss_min的系数
            z_min=[];
            z_min=z0;         
            global theta_min
            theta_min=[];
            theta_min=p.theta;  

            global search_index
            search_index=1;
            global global_p_temp
            global_p_temp=p;

            loss_now=loss_MSE_0;
            fun = @rosenbrockwithgrad_5;
            try
                % 这里是可能产生错误的代码
                [~,~,exitflag,output]= fmincon(fun, z0, [], [], [], [], lb, ub, [], p.fit_option3);       
                opti_mode_used=3;
            catch
                % 在发生错误时执行的代码
                try
                    % 这里是可能产生错误的代码
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;  
                    [~,~,exitflag,output] = fminunc(fun,z0,p.fit_option1);        %经过测试，最后返回的z，会多一次梯度下降，导致下面的z和call_for_fminunc_5.m中最后的z相差一步
                    opti_mode_used=1;
                catch
                    % 在发生错误时执行的代码
                    %loss_min=loss_MSE(end);
                    %z=z0;
                    loss_MSE=[];
                    loss_MSE(1)=loss_MSE_0;     
                    z_now=[];
                    z_now=z0;  
                    theta_now=[];
                    theta_now=p.theta;     
                    loss_min=[];
                    loss_min=loss_MSE_0;         
                    z_min=[];
                    z_min=z0;   
                    theta_min=[];
                    theta_min=p.theta;                      
                    opti_mode_used=0;
                end
            end

            loss_out_temp=loss_min;
            k=z_min.*k0;
            f=theta_min';
            delta_kth_out=k(1);
            kc_out=k(2);            
            theta_out=f';

            clear global_p_temp
            clear search_index
            clear loss_MSE  
            clear loss_min            
            clear z_now
            clear theta_now
            clear z_min
            clear theta_min
            %% theta_end检验与k1,k2检验
            %由于theta的最后一项对应delta_K，所以对其施加限制，必须落在区间
            pass_index=1;
            % if theta_out(end)<p.theta_end_limit(1) 
            %     pass_index=0;
            % elseif theta_out(end)>p.theta_end_limit(2)
            %     pass_index=0;
            if (delta_kth_out<0)||(kc_out<0)
                pass_index=0;
            % elseif delta_kth_out>kc_out
            %     pass_index=0;                        
            end

            if pass_index==1
                loss_out=loss_out_temp;
            else
                loss_out=inf;
            end

        end

    else
        %%
        error('wrong loss_mode!')

    end
end

%load('B:\Desktop\matlab.mat')
%loss_cal_optimize(loss_mode,p)