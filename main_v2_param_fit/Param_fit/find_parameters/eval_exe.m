function [eval_out] = eval_exe(evalstr,xtrain_in,Const_pair_now_in)
%EVAL_EXE 此处显示有关此函数的摘要
   xtrain=xtrain_in;
   Const_pair_now=Const_pair_now_in;
   eval_out=eval([evalstr ';']);
end

