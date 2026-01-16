clc
curUinput={};
curAverInput={};
% SPLITTED=SPLITTE_temp(i);
xparticlei=xparticlek_1(i,:);
a_up=xparticlei(42)-30;
a_down=xparticlei(22)-30;
a_up=3;
a_down=2;
SPLITTED=false;
m_index = getModelIndexFunc(a_up,a_down); % 判断裂纹阶段
stage=[1,2,3,4,5,6,7,0,0,0,0,8];% 每个阶段在元胞中的位�?
if ~SPLITTED % 如果前缘未分离，从完整裂纹代理模型集中提取代理模�?
    curUinput=Uinput_integrated{stage(m_index)};
    curAverInput=averInput_integrated{stage(m_index)};
    m_name=sprintf('nn_stage%d',m_index)
else % 如果前缘已分离，从分离裂纹代理模型集中提取代理模�?
    if m_index==3
        curUinput=Uinput_splitted{1};
        curAverInput=averInput_splitted{1};
    elseif m_index==5
        curUinput=Uinput_splitted{2};
        curAverInput=averInput_splitted{2};
    end
    m_name=sprintf('nn_stage%ds',m_index)
end
yRegSet=xparticlei(1:21);
zRegSet=xparticlei(22:42);
logCstar=xparticlei(43);
gamma=xparticlei(44);
[yRegSet,zRegSet,SPLITTED,logCstar,gamma] = a2aNew(yRegSet,zRegSet,delta_sigma_set,m_name,curUinput,curAverInput,logCstar,gamma,step);
xparticleForCalII(i,:)=[yRegSet,zRegSet,logCstar,gamma];
% SPLITTE_temp(i)=SPLITTED
