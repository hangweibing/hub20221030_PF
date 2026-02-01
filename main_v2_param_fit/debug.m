% figure
% subplot(2,1,1)
% for i=509
%     [F_ksdensity,XI_ksdensity]=ksdensity(upcrackparticles(:,i));
%     plot(XI_ksdensity,F_ksdensity,'LineWidth',5);
% %     legend(['第' num2str((i-1)*step/1950.70866) '小时'])
%     hold on
% end
% for i=570
%     [F_ksdensity,XI_ksdensity]=ksdensity(upcrackparticles(:,i));
%     plot(XI_ksdensity,F_ksdensity,'LineWidth',5);
% %     legend(['第' num2str((i-1)*step/1950.70866) '小时'])
%     hold on
% end
% 
% for i=609
%     [F_ksdensity,XI_ksdensity]=ksdensity(upcrackparticles(:,i));
%     plot(XI_ksdensity,F_ksdensity,'LineWidth',5);
% %     legend(['第' num2str((i-1)*step/1950.70866) '小时'])
%     hold on
% end
% % for i=509:570
% %     [F_ksdensity,XI_ksdensity]=ksdensity(K_maxSet(:,i));
% %     plot(XI_ksdensity,F_ksdensity);
% % %     legend(['第' num2str((i-1)*step/1950.70866) '小时'])
% %     hold on
% % end
% 
% legend('第520个小时','第583个小时','第623个小时','FontSize',40);
% xlabel('Crack length/mm','FontSize',30);
% ylabel('pdf','FontSize',30);
% % set(gca,'fontname','Times New Roman');
% set(gca,'FontSize',40);
% 
% 
% subplot(2,1,2)
% for i=509
%     [F_ksdensity,XI_ksdensity]=ksdensity(K_maxSet(:,i));
%     plot(XI_ksdensity,F_ksdensity,'LineWidth',5);
% %     legend(['第' num2str((i-1)*step/1950.70866) '小时'])
%     hold on
% end
% for i=570
%     [F_ksdensity,XI_ksdensity]=ksdensity(K_maxSet(:,i));
%     plot(XI_ksdensity,F_ksdensity,'LineWidth',5);
% %     legend(['第' num2str((i-1)*step/1950.70866) '小时'])
%     hold on
% end
% 
% for i=609
%     [F_ksdensity,XI_ksdensity]=ksdensity(K_maxSet(:,i));
%     plot(XI_ksdensity,F_ksdensity,'LineWidth',5);
% %     legend(['第' num2str((i-1)*step/1950.70866) '小时'])
%     hold off
% end
% % for i=509:570
% %     [F_ksdensity,XI_ksdensity]=ksdensity(K_maxSet(:,i));
% %     plot(XI_ksdensity,F_ksdensity);
% % %     legend(['第' num2str((i-1)*step/1950.70866) '小时'])
% %     hold on
% % end
% 
% legend('第520个小时','第583个小时','第623个小时','FontSize',40);
% xlabel('Kmax','FontSize',30);
% ylabel('pdf','FontSize',30);
% % set(gca,'fontname','Times New Roman');
% set(gca,'FontSize',40);

%%
hold off
clc
n_nodes=21;
testErrSet=[0.019803420755871 0.058377623733943 0.013343080193008 0.009221345729625 0.037283991994545 0.011408674471790 0.082711915490345 0.032375406062069 0.041104885753232 0.057544490601307];
%
centers=[13,30];
nRegPoint=n_nodes; % 等间距点数目
thetas=linspace(3/2*pi,2*pi,n_nodes);
load('pod_models.mat'); %加载代理模型数据

K_maxSet_temp=[];
Uinput_splitted_1=Uinput_splitted{1};
averInput_splitted_1=averInput_splitted{1};

upcrackparticles_temp=30.001:0.01:35;
aver_delta_sigma_set=130*ones(1,length(upcrackparticles_temp));
for j=1:length(upcrackparticles_temp)
    aver_delta_sigma=aver_delta_sigma_set(1);
    Smax=130;
    curUinput={};
    curAverInput={};
    SPLITTED=false;
    a_up=upcrackparticles_temp(1,j)-30;
%         a_up=1.4;
    a_down=0;
    stage=[1,2,3,4,5,6,7,0,0,0,0,8]; % 每个阶段在元胞中的位置
    m_index = getModelIndexFunc(a_up,a_down); % 判断裂纹阶段
    if SPLITTED && (m_index==3 || m_index==5)  % 如果前缘已分离，从分离裂纹代理模型集中提取代理模型
        if m_index==3
            curUinput=Uinput_splitted_1;
            curAverInput=averInput_splitted_1;
        elseif m_index==5
            curUinput=Uinput_splitted_2;
            curAverInput=averInput_splitted_2;
        end
        m_name=sprintf('nn_stage%ds',stage(m_index));
    else % 如果前缘未分离，从完整裂纹代理模型集中提取代理模型
        curUinput=Uinput_integrated{stage(m_index)};
        curAverInput=averInput_integrated{stage(m_index)};
        m_name=sprintf('nn_stage%d',stage(m_index));
    end
    yRegSet=centers(1)+a_up*sin(thetas);
    zRegSet=centers(2)+a_up*cos(thetas);
    inputRegSet=[yRegSet,zRegSet]';
    input=curUinput'*(inputRegSet-curAverInput);%POD里每一列表示一个样本
    [~,K_max] = sim_K_func1(m_name,input,aver_delta_sigma,Smax,testErrSet);
    K_maxSet_temp(:,j)=K_max;
%     end
% %     [~,x1]=ksdensity(K_maxSet_temp(:,j));
% %     [y1,x1]=ksdensity(K_maxSet_temp(:,j),linspace(min(x1),max(x1),N));
% %     label=find(x1>0);
% %     xs=x1(label);
% %     ys=y1(label);
% %     yr=normcdf(xs,33.4,3.34);
% %     yrs=yr.*ys;
% %     pofCum=cumtrapz(xs,yrs);
% %     pof=pofCum(end);
% %     truePoFSet_temp=[truePoFSet_temp pof];
%     disp(['计算进展为' num2str(j/k) '; 失效概率为  ' num2str(pof)]);
end
plot(upcrackparticles_temp-30,K_maxSet_temp,'LineWidth',5);
xlabel('Crack length/mm','FontSize',30);
ylabel('Kmax','FontSize',30);
% set(gca,'fontname','Times New Roman');
set(gca,'FontSize',40);