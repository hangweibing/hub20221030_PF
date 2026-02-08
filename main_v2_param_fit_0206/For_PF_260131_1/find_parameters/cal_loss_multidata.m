function [loss_multidata] = cal_loss_multidata(gp,fitness_best_constC_all_data,numGenes,num_const,lambda)
    theta=fitness_best_constC_all_data(:,4+num_const:4+num_const+numGenes+1);           %theta的数量是numGenes+2
    lambda1=lambda(1);
    lambda2=lambda(2);
    lambda3=lambda(3);
    lambda4=lambda(4);

    multidata_loss=gp.fitness.multidata_loss_mode;

    if multidata_loss==1
%% 
        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=inf;                                   %如果对不同数据的拟合中存在异常状况
        else
            loss_multidata=max(fitness_best_constC_all_data(:,1));                   %取不同数据中拟合效果最差的loss作为fitness
        end


    elseif multidata_loss==2
%%         
        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=inf;
        else
            Ndata=size(theta,1);
            m=mean(theta);
            ss=std(theta);
            max_of_var=max(ss);
            if max_of_var>lambda1
                loss_multidata=inf;
            else
                loss_of_fit=lambda2*max(fitness_best_constC_all_data(:,1));
                loss_of_mean=lambda3/Ndata*(sum(m.^2));
                loss_of_var=lambda4*sum(ss);
                loss_multidata=loss_of_fit+loss_of_mean+loss_of_var;
            end        
    %         loss_of_mean=lambda1/Ndata*(sum(m.^2));
    %         loss_of_var=lambda2*(max(v));
    %         loss=loss_of_fit+loss_of_mean+loss_of_var;
           %fitness_at_constC_all_data(:,1)=loss;
        end
    
    elseif multidata_loss==3
        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=inf;
        else
            Ndata=size(theta,1);
            m=mean(theta);                        %theta的均值
            ss=std(theta);                            %各个系数的标准差
            max_of_var=max(ss);
            loss_of_fit=max(fitness_best_constC_all_data(:,1));
            %loss_of_var=sum(ss);
            loss_of_var=max_of_var;

            if gp.fitness.gen_count_now<gp.runcontrol.stage1                       %分阶段改变fitness考察要求
                loss_multidata=loss_of_fit;             
            else

                if loss_of_fit>gp.fitness.loss_MSE_limit          %到第二阶段的时候已经存在gp.fitness.loss_MSE_limit的值了，1e-10或上一阶段loss_of_fit最好的值再乘个系数
                    loss_multidata=inf;                                     %不允许loss_of_fit再有大幅上涨
                else
                    loss_multidata=loss_of_var;
                end
            end
        end

    elseif multidata_loss==4
        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=inf;
        else
            Ndata=size(theta,1);                               %有多少组数据
            loss_of_fit=max(fitness_best_constC_all_data(:,1));
            k_of_AIC=size(theta,2);

            if gp.fitness.gen_count_now<gp.runcontrol.stage1                       %分阶段改变fitness考察要求
                loss_multidata=loss_of_fit;             
            else

                if loss_of_fit>gp.fitness.loss_MSE_limit          %到第二阶段的时候已经存在gp.fitness.loss_MSE_limit的值了，1e-10或上一阶段loss_of_fit最好的值再乘个系数
                    loss_multidata=inf;                                     %不允许loss_of_fit再有大幅上涨
                else
                    for i=1:Ndata
                        n_of_AIC=size(gp.multidata{1,i},1);      %每组数据有多少个数据点
                        AIC(i)=2*k_of_AIC+n_of_AIC*log(fitness_best_constC_all_data(i,1));
                        AICc(i)=AIC(i)+2*k_of_AIC*(k_of_AIC+1)/(n_of_AIC-k_of_AIC-1);
                    end                    
                    loss_multidata=max(AICc);
                    if loss_multidata<0
                        loss_multidata=-1/loss_multidata;   
                    else
                        loss_multidata=inf;
                    end

                end

            end

        end
    elseif multidata_loss==5  

        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=[inf,inf,inf];
        else
            Ndata=size(theta,1);
            m=mean(theta);                        %theta的均值
            ss=std(theta);                            %各个系数的标准差
            max_of_var=max(ss);
            loss_of_fit=max(fitness_best_constC_all_data(:,1));
            %loss_of_var=sum(ss);
            loss_of_var=max_of_var;
            loss_of_parameter_num=numGenes+2;

            %第一阶段，选择loss_of_fit
            if gp.fitness.gen_count_now<=gp.runcontrol.stage1                       
                loss_multidata=[loss_of_fit,loss_of_var,loss_of_parameter_num];             
            end
            %第二阶段，选择loss_of_var，但是loss_of_fit不能超过gp.fitness.loss_MSE_limit
            if (gp.fitness.gen_count_now>gp.runcontrol.stage1)&&(gp.fitness.gen_count_now<=gp.runcontrol.stage2)   
                if loss_of_fit<gp.fitness.loss_MSE_limit          %到第二阶段的时候已经存在gp.fitness.loss_MSE_limit的值了，1e-10或上一阶段loss_of_fit最好的值再乘个系数
                    %loss_multidata=loss_of_var;
                    loss_multidata=[loss_of_fit,loss_of_var,loss_of_parameter_num];
                else
                    loss_multidata=[inf,inf,inf];      
                    %loss_multidata=inf;                                     %不允许loss_of_fit再有大幅上涨
                end   
            end
            %第三阶段，选择loss_of_parameter_num，但是loss_of_fit不能超过gp.fitness.loss_MSE_limit，且loss_of_var不能超过gp.fitness.loss_VAR_limit
            if gp.fitness.gen_count_now>gp.runcontrol.stage2   
                if (loss_of_fit<gp.fitness.loss_MSE_limit)&&(loss_of_var<gp.fitness.loss_VAR_limit)   
                    loss_multidata=[loss_of_fit,loss_of_var,loss_of_parameter_num];
                    %loss_multidata=loss_of_parameter_num;
                else
                    loss_multidata=[inf,inf,inf];
                    %loss_multidata=inf;                                     %不允许loss_of_fit和loss_of_var再有大幅上涨
                end   

            end
        end
%% 模式6：不筛选，直接三个值传回去
    elseif multidata_loss==6    
        if  any(~isfinite(fitness_best_constC_all_data(:,1))) || any(~isreal(theta))
            loss_multidata=[inf,inf,inf];
        else
            Ndata=size(theta,1);
            m=mean(theta);                        %theta的均值
            %ss=std(theta);                            %各个系数的标准差
            ss=var(theta);                            %各个系数的标准差
            if gp.fitness.cal_k_var
                k1=fitness_best_constC_all_data(:,2);
                k2=fitness_best_constC_all_data(:,3);
                %max_of_var=max([ss,std(k1),std(k2)]);
                max_of_var=max([ss,var(k1),var(k2)]);
            else
                max_of_var=max(ss);
            end
            loss_of_fit=max(fitness_best_constC_all_data(:,1));

            %%%判断是单材料情况还是多材料情况
            if Ndata==1
                loss_of_var=0.1;
            else
                %loss_of_var=sum(ss);
                loss_of_var=max_of_var;
            end

            loss_of_parameter_num=numGenes+2;
            loss_multidata=[loss_of_fit,loss_of_var,loss_of_parameter_num];
        end
    else

        error('wrong multidata_loss')

    end



end


%load('B:\Desktop\test240108.mat')
%cal_loss_multidata(gp,fitness_best_constC_all_data,numGenes,num_const,lambda)
