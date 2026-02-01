% 根据pod_model_data文件夹内训练的mat文件生成pod_models.mat，格式与之前相同
stage_intg=[1,2,3,4,5,6,7,8];
stage_split=[3,5];
averInput_integrated=cell(1,8);
averInput_splitted=cell(1,2);
Uinput_integrated=cell(1,8);
Uinput_splitted=cell(1,2);

for i=1:8
    tmp_data=importdata(sprintf('./pod_model_data/model_data_stage%d.mat',stage_intg(i)));
    Uinput_integrated{i}=tmp_data.Uinput;
    averInput_integrated{i}=tmp_data.averInput;
end


for i=1:2
    tmp_data=importdata(sprintf('./pod_model_data/model_data_stage%ds.mat',stage_split(i)));
    Uinput_splitted{i}=tmp_data.Uinput;
    averInput_splitted{i}=tmp_data.averInput;
end

save('pod_models','averInput_integrated','averInput_splitted','Uinput_integrated','Uinput_splitted');
fprintf('finished\n')
