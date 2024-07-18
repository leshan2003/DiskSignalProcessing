clear;
clc;

file_path = './G13Data/Std.RF';
fp = fopen(file_path, 'r');
A_0 = fread(fp, 'uint8');
fclose(fp);
A_0 = double(A_0);

% 均值滤波
ave_num = 1;
A_0_new_temp = zeros(length(A_0)-ave_num+1,1);
for i = 1:length(A_0)-ave_num+1
    A_0_new_temp(i) = sum(A_0(i:i+ave_num-1))/ave_num;
end
% 滤波器
sample_rate = 1e8;
b = fir2(88,[0 0.06 0.14 0.21 0.28 0.4 1],[0.8 1.4 2.0 3.0 1.8 0.01 0]);
% freqz(b,1,1024,sample_rate);
A_0_new = filter(b,1,A_0_new_temp);

% 线性插值
x = 1:length(A_0_new);
interp_rate = 10;
new_x = linspace(x(1), x(end), interp_rate * length(x) - interp_rate + 1);
new_A = interp1(x, A_0_new, new_x, 'linear'); % 使用 interp1 进行线性插值
new_A = double(new_A);
zero_line = 114.4;

% 检测零点
num_zero = 0;
shift_list = zeros(length(new_A),1);
index_list = zeros(length(new_A),1);
for j = 1:length(new_A)-1
    if (new_A(j)-zero_line)*(new_A(j+1)-zero_line)<=0
        num_zero = num_zero + 1;
        index_list(num_zero)=j;
        if new_A(j)<new_A(j+1)
            shift_list(num_zero)=1;       % 1代表峰
        elseif new_A(j)>new_A(j+1)
            shift_list(num_zero)=-1;      % -1代表谷
        end
    end
end
shift_list = shift_list(1:num_zero);
index_list = index_list(1:num_zero);

% 求线速度均值、线速度随时间变化曲线
frame_len_of_T = 1488;
T_len = 400/3*1e-9;
sample_interval = 1e-8/interp_rate;
upper_bound = 900;
lower_bound = 700;
v_list = zeros(length(index_list),1);
time_of_v_list = zeros(length(index_list),1);
v_index = 0;
last_header_index = 0;
for j = 1:length(index_list)-1
    half_wave_len = index_list(j+1)-index_list(j);
    if half_wave_len > lower_bound && half_wave_len < upper_bound
        if last_header_index == 0
            last_header_index = index_list(j);
        elseif (index_list(j)-last_header_index) > 80979 && (index_list(j)-last_header_index) < 84425
            % disp(index_list(j)-last_header_index);
            v_index = v_index + 1;
            v_list(v_index) = frame_len_of_T*T_len/(sample_interval*(index_list(j)-last_header_index));
            time_of_v_list(v_index) = index_list(j)*sample_interval;
            last_header_index = index_list(j);
        else 
            temp = index_list(j)-last_header_index;
            % disp(temp);
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
p3 = plot(0.00767293, 2.40806, 'o', 'MarkerSize', 8,'MarkerFaceColor', 'cyan', 'MarkerEdgeColor', 'black');
text(0.00767293, 2.409, "time = 0.00767 s", 'FontSize', 12, 'FontName', 'Times New Roman');
text(0.00767293, 2.4085, "v = 2.41 m/s", 'FontSize', 12, 'FontName', 'Times New Roman');
p4 = plot(0.0427343, 2.3927, 'o', 'MarkerSize', 8,'MarkerFaceColor', 'cyan', 'MarkerEdgeColor', 'black');
text(0.0427343, 2.392, "time = 0.04273 s", 'FontSize', 12, 'FontName', 'Times New Roman');
text(0.0427343, 2.391, "v = 2.39 m/s", 'FontSize', 12, 'FontName', 'Times New Roman');
p5 = plot(0.0892138, 2.40747, 'o', 'MarkerSize', 8,'MarkerFaceColor', 'cyan', 'MarkerEdgeColor', 'black');
text(0.0892138, 2.409, "time = 0.08921 s", 'FontSize', 12, 'FontName', 'Times New Roman');
text(0.0892138, 2.4085, "v = 2.41 m/s", 'FontSize', 12, 'FontName', 'Times New Roman');
p6 = plot(0.126024, 2.39293, 'o', 'MarkerSize', 8,'MarkerFaceColor', 'cyan', 'MarkerEdgeColor', 'black');
text(0.126024, 2.392, "time = 0.1260 s", 'FontSize', 12, 'FontName', 'Times New Roman');
text(0.126024, 2.391, "v = 2.39 m/s", 'FontSize', 12, 'FontName', 'Times New Roman');
p7 = plot(0.169186, 2.40765, 'o', 'MarkerSize', 8,'MarkerFaceColor', 'cyan', 'MarkerEdgeColor', 'black');
text(0.160186, 2.409, "time = 0.1692 s", 'FontSize', 12, 'FontName', 'Times New Roman');
text(0.160186, 2.4085, "v = 2.41 m/s", 'FontSize', 12, 'FontName', 'Times New Roman');
p8 = plot(0.205491, 2.39316, 'o', 'MarkerSize', 8,'MarkerFaceColor', 'cyan', 'MarkerEdgeColor', 'black');
text(0.19, 2.392, "time = 0.2055 s", 'FontSize', 12, 'FontName', 'Times New Roman');
text(0.19, 2.391, "v = 2.39 m/s", 'FontSize', 12, 'FontName', 'Times New Roman');
legend([p1 p2],{'线速度',"线速度中值"},'FontSize', 16, 'FontName', 'Songti SC');
xlabel('时间 (s)', 'FontSize', 20, 'FontName', 'Songti SC');
xlim([0 0.21]);
ylabel('线速度 (m/s)', 'FontSize', 20, 'FontName', 'Songti SC');
ylim([2.39 2.41]);
title('线速度随时间变化（自采数据，标准光头，标准伺服系统）', 'FontSize', 25, 'FontName', 'Songti SC');