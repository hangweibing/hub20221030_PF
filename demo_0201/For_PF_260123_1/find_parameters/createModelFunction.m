function [calcFun] = createModelFunction(expr_str, variables)
% 创建模型计算函数
% expr_str: 表达式字符串
% variables: 变量名元胞数组，如 {'delta_K','Const_pair_now','f','k','R'}

% 构造函数定义头
argStr = strjoin(variables, ',');
funStr = ['@(', argStr, ') ', expr_str];

% 生成函数句柄
calcFun = str2func(funStr);


end

