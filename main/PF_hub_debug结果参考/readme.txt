1. pod_models.mat根据pod_model_data文件夹内训练的mat文件生成，格式与之前相同，只包含POD投影矩阵和均值，无模型文件
2. 每个模型的误差在pod_model_data文件夹对应的mat文件里，有训练误差和测试误差，可根据需要提取。例如提取完整前缘模型的数据
for i=1:8
    tmp_data=importdata(sprintf('./pod_model_data/model_data_stage%d.mat',stage_intg(i)));
end
tmp_data内为该模型所有相关的训练数据，包含原始数据、训练测试数据和误差数据等。
3. 替换curModel为m_name，传递到sim_K_func中，使用eval函数执行