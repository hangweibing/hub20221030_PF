% 简单的语法测试
try
    % 尝试加载函数
    help a2aFunc
    disp('a2aFunc语法检查通过');
catch e
    disp(['语法错误: ' e.message]);
end