clear;
clc;

% 读取文件
file_path = './G13Data/HF_G13';
fp = fopen(file_path, 'r');
A_0 = fread(fp, 'uint8');
fclose(fp);
A_0 = double(A_0);

% 检测零点
zero_line = 126.4;
num_zero = 0;
shift_list = zeros(length(A_0),1);
index_list = zeros(length(A_0),1);
for j = 1:length(A_0)-1
    if (A_0(j)-zero_line)*(A_0(j+1)-zero_line)<=0
        num_zero = num_zero + 1;
        index_list(num_zero)=j;
        if A_0(j)<A_0(j+1)
            shift_list(num_zero)=1;       % 1代表峰
        elseif A_0(j)>A_0(j+1)
            shift_list(num_zero)=-1;      % -1代表谷
        end
    end
end
shift_list = shift_list(1:num_zero);
index_list = index_list(1:num_zero);

% 求线速度均值、线速度随时间变化曲线
frame_len_of_T = 1488;
T_len = 400/3*1e-9;
sample_interval = 1e-8;
upper_bound = 60;
lower_bound = 50;
v_list = zeros(length(index_list),1);
time_of_v_list = zeros(length(index_list),1);
v_index = 0;
last_header_index = 0;
for j = 1:length(index_list)-1
    half_wave_len = index_list(j+1)-index_list(j);
    if half_wave_len > lower_bound && half_wave_len < upper_bound
        if last_header_index == 0
            last_header_index = index_list(j);
        else
            v_index = v_index + 1;
            v_list(v_index) = frame_len_of_T*T_len/(sample_interval*(index_list(j)-last_header_index));
            time_of_v_list(v_index) = index_list(j)*sample_interval;
            last_header_index = index_list(j);
        end
    end
end
v_list = v_list(1:v_index);
time_of_v_list = time_of_v_list(1:v_index);
v_mean = mean(v_list);
v_middle = (max(v_list) + min(v_list)) / 2;
figure;
p1 = plot(time_of_v_list, v_list, 'Color', 'red', 'DisplayName', '线速度', 'LineWidth', 1.5);
hold on;
p2 = plot(time_of_v_list, v_middle*ones(length(time_of_v_list),1), 'Color', 'blue', 'DisplayName', '线速度中值', 'LineWidth', 1.5);
text_x_values = [0.031 0.048 0.059 0.075];
text_y_values = [3.3664 3.360 3.366 3.373];
p3 = plot(text_x_values(1)-0.00036, text_y_values(1), 'o', 'MarkerSize', 8,'MarkerFaceColor', 'cyan', 'MarkerEdgeColor', 'black');
text(text_x_values(1), text_y_values(1)+0.0005, "time1 = " + text_x_values(1) + " s", 'FontSize', 12, 'FontName', 'Times New Roman');
p4 = plot(0.0478156, 3.35986, 'o', 'MarkerSize', 8,'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'black');
text(text_x_values(2)+0.0003, text_y_values(2)-0.0003, "v min = 3.360 m/s", 'FontSize', 12, 'FontName', 'Times New Roman');
p5 = plot(0.0590852, 3.36643, 'o', 'MarkerSize', 8,'MarkerFaceColor', 'cyan', 'MarkerEdgeColor', 'black');
text(0.0606, 3.36671, "time2  = " + text_x_values(3) + " s", 'FontSize', 12, 'FontName', 'Times New Roman');
p6 = plot(0.0749856, 3.373, 'o', 'MarkerSize', 8,'MarkerFaceColor', 'green', 'MarkerEdgeColor', 'black');
text(0.075603, text_y_values(4), "v max = " + text_y_values(4) + " m/s", 'FontSize', 12, 'FontName', 'Times New Roman');
legend([p1 p2],{'线速度',"线速度中值"},'FontSize', 16, 'FontName', 'Songti SC');
xlabel('时间 (s)', 'FontSize', 20, 'FontName', 'Songti SC');
ylabel('线速度 (m/s)', 'FontSize', 20, 'FontName', 'Songti SC');
title('线速度随时间变化（HF信号，无均衡器）', 'FontSize', 25, 'FontName', 'Songti SC');