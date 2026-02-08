function [diff_index, eq,diff_omega, diff_omega_OF_deq, diff_delta_K, contains_k1,contains_k2 ] = diff_F(evalstr_IN,numConst,diff_model)
%%此求导程序只针对delta_K和K_max两个独立变量
numGenes = numel(evalstr_IN); 
for i=1:numGenes+2
    syms(['f',num2str(i)]);
end
for i=1:numConst
    syms(['c',num2str(i)]);
end
syms delta_K Kmax k1 k2 R
%eq=['log10(f1)'];
eq=['f1'];
for i=1:numGenes
    x1_str=['(','delta_K/','k1',')'];
    x2_str=['(','Kmax/','k2',')'];
    evalstr_IN{i} = strrep(evalstr_IN{i},'x1',x1_str);
    evalstr_IN{i} = strrep(evalstr_IN{i},'x2',x2_str);
    %eq_part{i}=['f',num2str(i+1),'*',evalstr_IN{i}];
    eq_part{i}=['f',num2str(i+1),'*','log10(',evalstr_IN{i},')',];
    eq=[eq,'+',eq_part{i}];
end
eq=[eq,'+','f',num2str(numGenes+2),'*','log10(delta_K)'];

derivative_result=diff(eval(eq));
% 检测逻辑同上
if any(isnan(derivative_result(:)))
    diff_index = 0;
else
    diff_index = 1;
end

if diff_index==1

    if diff_model==2.1
    
        Kmax_str=['(','delta_K/','(','1-R',')',')'];
        eq_delta_K= strrep(eq,'Kmax',Kmax_str);
        variables_in_eq = symvar(eval(eq));
        % 检查 k1 是否在变量列表中
        contains_k1 = ismember('k1', variables_in_eq);
        contains_k2 = ismember('k2', variables_in_eq);
    %%%对f求导
        diff_omega=cell(1,numGenes+4);
        diff_omega_OF_deq=cell(1,numGenes+4);
        for i=1:numGenes+2
            diff_omega{i}=char(diff(eval(eq),['f',num2str(i)]));     
        end
        diff_omega{numGenes+3}=char(diff(eval(eq),k1));
        diff_omega{numGenes+4}=char(diff(eval(eq),k2));
    
        %转换以方便矩阵运算：
        diff_omega=strrep(diff_omega,'/','./');
        diff_omega=strrep(diff_omega,'^','.^');
        diff_omega=strrep(diff_omega,'*','.*');
    
        %%%eq对于delta_K求导，然后再对参数求导
        diff_delta_K = char(diff(eval(eq_delta_K),delta_K));
        for i=1:numGenes+2
            diff_omega_OF_deq{i}=char(diff(eval(diff_delta_K),['f',num2str(i)]));     
        end
        diff_omega_OF_deq{numGenes+3}=char(diff(eval(diff_delta_K),k1));
        diff_omega_OF_deq{numGenes+4}=char(diff(eval(diff_delta_K),k2));
        diff_delta_K=strrep(diff_delta_K,'/','./');
        diff_delta_K=strrep(diff_delta_K,'^','.^');
        diff_delta_K=strrep(diff_delta_K,'*','.*');  
        
        %转换以方便矩阵运算：
        diff_omega_OF_deq=strrep(diff_omega_OF_deq,'/','./');
        diff_omega_OF_deq=strrep(diff_omega_OF_deq,'^','.^');
        diff_omega_OF_deq=strrep(diff_omega_OF_deq,'*','.*');
      
        eq=strrep(eq,'/','./');
        eq=strrep(eq,'^','.^');
        eq=strrep(eq,'*','.*');
    
    elseif diff_model==2.2
        syms z1 z2 R pk1 pk2
        eq=['f',num2str(numGenes+2),'*','z1'];
        for i=1:numGenes
            x1_str=['(','pk1','/(10^z1)',')'];
            x2_str=['(','(10^z1)/','pk1',')'];
            x3_str=['(','(10^z2)/','pk2',')'];
            x4_str=['(','pk2','/(10^z2)',')'];
            x5_str='R';
            evalstr_IN{i} = strrep(evalstr_IN{i},'x1',x1_str);
            evalstr_IN{i} = strrep(evalstr_IN{i},'x2',x2_str);
            evalstr_IN{i} = strrep(evalstr_IN{i},'x3',x3_str);
            evalstr_IN{i} = strrep(evalstr_IN{i},'x4',x4_str);
            evalstr_IN{i} = strrep(evalstr_IN{i},'x5',x5_str);
            eq_part{i}=['f',num2str(i+1),'*',evalstr_IN{i}];
            eq=[eq,'+',eq_part{i}];
        end
        eq=[eq,'+','f1'];
        dF_dk1=char(diff(eval(eq),z1));                %dy/dk1
        dF_dk2=char(diff(eval(eq),z2));                %dy/dk2
        
        %转换以方便矩阵运算：
        dF_dk1=strrep(dF_dk1,'/','./');
        dF_dk2=strrep(dF_dk2,'/','./');
        dF_dk1=strrep(dF_dk1,'^','.^');
        dF_dk2=strrep(dF_dk2,'^','.^');
        dF_dk1=strrep(dF_dk1,'*','.*');
        dF_dk2=strrep(dF_dk2,'*','.*');    
    end

else

    diff_omega=[]; diff_omega_OF_deq=[]; diff_delta_K=[]; contains_k1=[]; contains_k2=[];


end

end

%load('B:\Desktop\test231018.mat')
%diff_F(evalstr_IN,numConst)

