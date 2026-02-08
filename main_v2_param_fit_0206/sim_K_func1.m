function [deltaKSet,K_max] = sim_K_func1(m_name,input,aver_delta_sigma,Smax,testErrSet)
%SIM_K 此处显示有关此函数的摘要
%   此处显示详细说明
% KSet = simlssvm(curModel,input');%simlssvm将122行fprintf('m');注释掉了,lssvm里每一行为一个样本，最终都以这个为准
commd=sprintf('KSet=%s(input);',m_name);
% commd=sprintf('%s(input);',m_name);
eval(commd);
crackIndex=strrep(m_name,'nn_stage','');
if length(crackIndex)==1
    testErr=testErrSet(str2num(crackIndex));
else
    if crackIndex(1)=='3'
        testErr=testErrSet(9);
    elseif crackIndex(1)=='5'
        testErr=testErrSet(10);
    end
end
KSet(KSet<0)=20;
% 
% for i=1:length(KSet)
%     KSet(i)=KSet(i)+normrnd(0,abs(KSet(i)*testErr),1,1);
% end
KSet=KSet+normrnd(0,abs(mean(KSet))*testErr,1,1);
KSet(KSet<0)=20;
KSet= filloutliers(KSet,'linear');
KSet(KSet<0)=20;

% if ~isempty(outliers)
%     for k=1:length(outliers)
%         % 找到两个离该离群点最近的点
%         if outliers(k)<min(normal_points)
%             neighbor_points=normal_points(find(normal_points>1,2));
%             KSet(outliers(k))=KSet(neighbor_points(1));
%         elseif outliers(k)>max(normal_points)
%             tmp_array=find(normal_points<num_K_points);
%             neighbor_point=normal_points(tmp_array(end));
%             KSet(outliers(k))=KSet(neighbor_point);
%         else
%             tmp_array=find(normal_points<outliers(k));
%             left_neighbor=normal_points(tmp_array(end));
%             right_neighbor=normal_points(find(normal_points>outliers(k),1));
%             neighbor_points=[left_neighbor,right_neighbor];
%             KSet(outliers(k))=interp1([KSet(neighbor_points(1)),KSet(neighbor_points(2))],[neighbor_points(1),neighbor_points(2)],outliers(k),'linear','extrap');
%         end
%         
%     end
% end

% for i=1:size(KSet,2)
%     if abs(KSet(i)-averK)>0.5*abs(averK)
%         if i==1
%             KSet(i)=KSet(2);
%         elseif i==size(KSet,2)
%             KSet(i)=KSet(size(KSet,2)-1);
%         else
%             KSet(i)=(KSet(i-1)+KSet(i+1))/2;
%         end
%     end
% end

% fprintf('Kset mean: %f m_name: %s ',mean(KSet),m_name)
deltaKSet=abs(KSet.*aver_delta_sigma/130);
K_max=abs(KSet(end)*Smax/130)/sqrt(1000);  %上表面处的应力强度因子
end