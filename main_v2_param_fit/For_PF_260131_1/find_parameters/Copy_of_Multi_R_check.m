function [test_index,output_struct] = Copy_of_Multi_R_check(p_temp,showfigure,Const_pair_now,theta,delta_kth,kc)
 
    %% 初始化
    evalstr_test=p_temp.evalstr_test_fun;
    diff_delta_K_fun=p_temp.diff_delta_K_fun;

    numGenes = p_temp.numGenes; 

    NR = 30;
    R_temp =linspace(-1,0.95,NR);
    material_i_deltaK_min=0.1;
    %material_i_deltaK_max=max(D_i(:,4));
    material_i_deltaK_max=120;
    x1_temp = linspace (log10(material_i_deltaK_min), log10(material_i_deltaK_max),2000); 
    x1_temp=10.^x1_temp;

    %—— 构造网格并向量化计算 y ——%
    [Rg, x1_temp_g] = meshgrid(R_temp, x1_temp);      % 大小 400×40
    KMg       = x1_temp_g ./ (1 - Rg);


    %% 常量和参数
    f        = theta;
    k        = [delta_kth; kc];
    N_slope  = 5;
    c_log10  = log(10);


    %% 构造 xtest
    xtest     = cell(2,1);
    xtest{1}  = x1_temp_g ./ delta_kth;
    xtest{2}  = KMg ./ kc;

    %% 计算 y_grid
    y_grid = theta(1) * ones(size(Rg));
    % 基因部分累加
    for j = 1:numGenes
        gene_temp = evalstr_test{j}(xtest,Const_pair_now);     
        y_grid    = y_grid + theta(j+1) * lg(gene_temp);
    end
    % 加上最后一项 lg(ΔK)
    y_grid = y_grid + theta(2+numGenes) * lg(x1_temp_g);

    gy_grid=diff_delta_K_fun(x1_temp_g, Const_pair_now, f, k, Rg);
    gy_grid = gy_grid * c_log10 .* x1_temp_g;

    % % 找到所有 > 100 的元素下标
    % [row_high, col_high] = find(gy_grid > 100);
    % 
    % % 找到所有 < -100 的元素下标
    % [row_low,  col_low ] = find(gy_grid < -100);

    %为了数值稳定，设置gy的上下限
    mask_temp = isfinite(gy_grid) & (imag(gy_grid) == 0);
    gy_grid(mask_temp & gy_grid > 100)  =  80;
    gy_grid(mask_temp & gy_grid < -100) = -80;


    %—— 生成有效点掩码 ——%
    mask = isfinite(y_grid) & (imag(y_grid)==0);

    %—— 按列做 check 1 (连续性) + check 2 (梯度有效性) ——%
    test_index = 1;
    nmaxPoints = numel(x1_temp_g);


    %% 预分配三段存储数组
    X1_1 = zeros(nmaxPoints,1); Y_1 = zeros(nmaxPoints,1); GY_1 = zeros(nmaxPoints,1); R_1 = zeros(nmaxPoints,1);
    X1_2 = zeros(nmaxPoints,1); Y_2 = zeros(nmaxPoints,1); GY_2 = zeros(nmaxPoints,1); R_2 = zeros(nmaxPoints,1);
    X1_3 = zeros(nmaxPoints,1); Y_3 = zeros(nmaxPoints,1); GY_3 = zeros(nmaxPoints,1); R_3 = zeros(nmaxPoints,1);
    ptr1 = 0; ptr2 = 0; ptr3 = 0;
 
    %% 主循环：每列一次计算梯度，随后切片存储三段
    for i = 1:NR
        col_mask = mask(:,i);
        idx      = find(col_mask);
        % 连续性检查
        if numel(idx)>=2 && ~all(  col_mask( idx(1):idx(end) )  )
            test_index = 0;
            break;
        end
    
        % 一次性提取整列有效点并计算
        dK_full = x1_temp_g(idx,i);
        y_full  = y_grid(idx,i);
        gy_full = gy_grid(idx,i);
        R_val   = R_temp(i);

        if any(~isreal(gy_full)) || any(isnan(gy_full)) || any(isinf(gy_full))
            test_index = 0;
            break;
        end
    
        % 将结果分段存储
        % 分界点位置
        total = numel(idx);
        p = min(N_slope, total);
        if total > 2*p
            split1 = 1:p;
            split2 = (p+1):(total-p);
            split3 = (total-p+1):total;
        else
            split1 = [];
            split2 = 1:total;
            split3 = [];
        end

        n1 = numel(split1);
        n2 = numel(split2);
        n3 = numel(split3);
    
        % 头部
        if n1 > 0
            X1_1(ptr1+1:ptr1+n1) = dK_full(split1);
            Y_1(ptr1+1:ptr1+n1)  = y_full(split1);
            GY_1(ptr1+1:ptr1+n1) = gy_full(split1);
            R_1(ptr1+1:ptr1+n1)  = R_val;
            ptr1 = ptr1 + n1;
        end
        % 中部
        if n2 > 0
            X1_2(ptr2+1:ptr2+n2) = dK_full(split2);
            Y_2(ptr2+1:ptr2+n2)  = y_full(split2);
            GY_2(ptr2+1:ptr2+n2) = gy_full(split2);
            R_2(ptr2+1:ptr2+n2)  = R_val;
            ptr2 = ptr2 + n2;
        end
        % 尾部
        if n3 > 0
            X1_3(ptr3+1:ptr3+n3) = dK_full(split3);
            Y_3(ptr3+1:ptr3+n3)  = y_full(split3);
            GY_3(ptr3+1:ptr3+n3) = gy_full(split3);
            R_3(ptr3+1:ptr3+n3)  = R_val;
            ptr3 = ptr3 + n3;
        end
    end

    %% 截断多余部分
    output_struct = struct();  % 创建输出结构体
    output_struct.X1_1 = X1_1(1:ptr1);
    output_struct.GY_1 = GY_1(1:ptr1);
    output_struct.R_1 = R_1(1:ptr1);
    
    output_struct.X1_2 = X1_2(1:ptr2);
    output_struct.GY_2 = GY_2(1:ptr2);
    output_struct.R_2 = R_2(1:ptr2);
    
    output_struct.X1_3 = X1_3(1:ptr3);
    output_struct.GY_3 = GY_3(1:ptr3);
    output_struct.R_3 = R_3(1:ptr3);

    output_struct.NI   = ptr1;
    output_struct.NII  = ptr2;
    output_struct.NIII = ptr3;

    output_struct.GY_2_mean = mean(GY_2)-1.0;         %为了数值稳定减1



    %% 作图
    if showfigure==1
        % 从结构体中提取数据
        X1 = [output_struct.X1_1; output_struct.X1_2; output_struct.X1_3];
        Y = [Y_1(1:ptr1); Y_2(1:ptr2); Y_3(1:ptr3)];  % 注意：Y_* 变量未存入结构体
        GY = [output_struct.GY_1; output_struct.GY_2; output_struct.GY_3];
        R_EFF = [output_struct.R_1; output_struct.R_2; output_struct.R_3];

        D_i=p_temp.data_i;
        % 2. 在 (x1,R) 平面做网格
        nx = 200;  % x1 方向格点数，可调
        nR = 200;   % R 方向格点数，可调
        x1_vec = linspace(log10(min(X1)), log10(max(X1)), nx);
        x1_vec=10.^x1_vec;
        R_vec  = linspace(min(R_EFF),  max(R_EFF),  nR);
        [ X1g, Rg ] = meshgrid(x1_vec, R_vec);
        
        % 3. 使用 griddata 插值
        Zg = griddata(X1, R_EFF, Y, X1g, Rg, 'cubic');  % 'cubic' 插值，亦可试 'linear'、'natural' 等
        z_data = 10.^Zg;         % 对实部应用 10.^ 操作

        Zg_y = griddata(X1, R_EFF, GY, X1g, Rg, 'cubic');        
        % 4. 绘制曲面
        figure('Visible', showfigure)  
        scatter3((D_i(:,1)),D_i(:,3),(D_i(:,2)),'o')
        hold on;
        scatter3(X1,R_EFF,10.^Y,'o')
        set(gca, 'XScale', 'log', 'YScale', 'linear', 'ZScale', 'log');
        h1=surf(X1g, Rg, z_data);      % 无网格线的平滑曲面
        % 设置曲面为光滑曲面，不显示网格线
        set(h1, 'EdgeColor', 'none');            
        % 设置曲面颜色为蓝色
        set(h1, 'FaceColor', 'blue');
        % 设置曲面透明度为50%
        set(h1, 'FaceAlpha', 0.3);
        % 设置y轴的方向为反向 
        set(gca, 'YDir', 'reverse')
        %axis([0.5 120 -1 1 10^(-14) 10^(-2)]);

        figure('Visible', showfigure)  
        h2=surf(X1g, Rg, Zg_y);      
        % 设置曲面为光滑曲面，不显示网格线
        set(h2, 'EdgeColor', 'none');            
        % 设置曲面颜色为蓝色
        set(h2, 'FaceColor', 'blue');
        % 设置曲面透明度为50%
        set(h2, 'FaceAlpha', 0.3);
        % 设置y轴的方向为反向 
        set(gca, 'YDir', 'reverse')
    end


end

