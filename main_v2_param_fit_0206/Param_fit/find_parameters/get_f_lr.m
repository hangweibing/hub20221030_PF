function [theta,y_pre] = get_f_lr(p,k)
%GET_F_LR  get f1,f2,...fn by linear regression

    numGenes = numel(p.evalstr2);
    Const_pair_now=p.Const_pair_now;
    [numData,~] =size(p.ytrain);
    delta_K=p.data_i(:,1);
    Kmax=p.data_i(:,4);
    xtrain=[delta_K./k(1),Kmax./k(2)];                 %这里输入的是新的k
    geneOutputs = ones(numData,numGenes+2);
    check_index=1;
    for i = 1:numGenes
        ind = i + 1;
        %gene_temp=eval_exe(evalstr{i},xtrain,Const_pair_now);
        try
            % 这里是可能产生错误的代码
            gene_temp=eval([p.evalstr2{i} ';']);
        catch
            % 在发生错误时执行的代码
            disp('出错：表达式中带入参数.');
        end                
        geneOutputs(:,ind)=lg(gene_temp);
        if  any(~isfinite(geneOutputs(:,ind))) || any(~isreal(geneOutputs(:,ind)))
            %disp('出错：表达式带入数据后出现了复数/无穷大.');
            check_index=0;
            break
        end
    end
    
    geneOutputs(:,numGenes+2)=lg(delta_K);          %补最后一列
    geneOutputs_temp=geneOutputs;
    if check_index==1
        if p.ridge_on
            %% 采用岭回归求解f
            geneOutputs(:,1)=[];                                                 %岭回归需要去掉第一列
            theta = ridge(p.ytrain,geneOutputs,p.ridge_k,0);     %输出的theta维数是一样的
        else
            goptrans = geneOutputs';
            prj = goptrans * geneOutputs;
            theta = pinv(prj) * goptrans * p.ytrain;                     %输出的theta维数是一样的
        end        
    else
        theta=zeros(numGenes+2,1);
    end
    y_pre=geneOutputs_temp*theta;

end

