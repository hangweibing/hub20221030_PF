model_integrated=cell(1,7);
model_datas=importdata('surrogate_model_data.mat');
model_datas_integrated=model_datas.model_integrated;
for i=1:8
    tmp_train_data=model_datas_integrated{i};
    train_x=tmp_train_data.xtrain;
    train_y=tmp_train_data.ytrain;
    tmp_net=feedforwardnet(20);
    train(tmp_net,train_x',train_y');
    genFunction(tmp_net,sprintf('nn_func/nn_integrated_%d',i));
end