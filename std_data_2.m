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
% sample_rate = 1e8;
% b = fir2(600,[0 0.01 0.02 0.03 0.04 0.05 0.075 1],[0.6 10^(-0.1/20) 10^(3/20) 10^(3/20) 10^(1/20) 10^(-3/20) 10^(-40/20) 0]);
% % freqz(b,1,1024,sample_rate);
% A_0_new = filter(b,1,A_0_new_temp);
A_0_new = A_0_new_temp;

% 线性插值
x = 1:length(A_0_new);
interp_rate = 1;
new_x = linspace(x(1), x(end), interp_rate * length(x) - interp_rate + 1);
new_A = interp1(x, A_0_new, new_x, 'linear'); % 使用 interp1 进行线性插值
new_A = double(new_A);
zero_line = mean(new_A);

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

% 岸、坑半波统计
peak_num = 0;
valley_num = 0;
peak_length_list = zeros(length(shift_list),1);
valley_length_list = zeros(length(shift_list),1);
for j = 1:length(shift_list)-1
    if shift_list(j) == 1
        peak_num = peak_num + 1;
        peak_length_list(peak_num) = index_list(j+1)-index_list(j);
    elseif shift_list(j) == -1
        valley_num = valley_num + 1;
        valley_length_list(valley_num) = index_list(j+1)-index_list(j);
    end
end
peak_length_list = peak_length_list(1:peak_num);
valley_length_list = valley_length_list(1:valley_num);

% max_length = max([max(peak_length_list), max(valley_length_list),240]);
max_length = 250*interp_rate;
min_length = min([min(peak_length_list), min(valley_length_list),0]);
peak_counts = zeros(max_length-min_length+1,1);
valley_counts = zeros(max_length-min_length+1,1);
for j = 1:max_length-min_length+1
    peak_counts(j) = sum(peak_length_list==j);
    valley_counts(j) = sum(valley_length_list==j);
end
% peak_counts(1) = 0;
% valley_counts(1) = 0;

x_values = min_length:max_length;
figure;
plot(x_values, peak_counts, '-*', 'Color', 'red', 'DisplayName', '岸', 'LineWidth', 1.5, 'MarkerSize', 5);
hold on;
plot(x_values, -valley_counts, '-*', 'Color', 'blue', 'DisplayName', '坑', 'LineWidth', 1.5, 'MarkerSize', 5);
hold off;
legend('FontSize', 16, 'FontName', 'Songti SC');
xlabel('记录符长度（采样点个数）', 'FontSize', 20, 'FontName', 'Songti SC');
xlim([0, 100]);
ylabel('记录符个数', 'FontSize', 20, 'FontName', 'Songti SC');
ylim([-25000, 25000]);
title('记录符分布（自采数据，标准光头，标准伺服系统，无滤波器）', 'FontSize', 25, 'FontName', 'Songti SC');
grid on;