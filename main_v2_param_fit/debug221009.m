%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% debug
% date: 2022-10-09
% author: Xuan Zhou
% mail: zhoux@buaa.edu.cn
% description: 裂纹前缘报错
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xparticles=importdata('221009particleSetErr.mat');
last_particle=xparticles(end,1:44);
y_set=last_particle(1:21);
z_set=last_particle(22:42);

plot_geometry_20
hold on
% plot(z_set,y_set)
for i =1:100000
    plot(xparticles(i,22:42),xparticles(i,1:21));
end

plot(xparticles(:,22:42)',xparticles(:,1:21)');