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
max_length = 100*interp_rate;
min_length = min([min(peak_length_list), min(valley_length_list),0]);
peak_counts = zeros(max_length-min_length+1,1);
valley_counts = zeros(max_length-min_length+1,1);
for j = 1:max_length-min_length+1
    peak_counts(j) = sum(peak_length_list==j);
    valley_counts(j) = sum(valley_length_list==j);
end

% 统计nT长度方差
peak_split_index = [97 198 248 297 352 409 464 520 581 648 833];
valley_split_index = [97 199 253 308 360 417 475 531 586 648 833];

peak_3T_length_list = zeros(length(peak_length_list),1);
peak_4T_length_list = zeros(length(peak_length_list),1);
peak_5T_length_list = zeros(length(peak_length_list),1);
peak_6T_length_list = zeros(length(peak_length_list),1);
peak_7T_length_list = zeros(length(peak_length_list),1);
peak_8T_length_list = zeros(length(peak_length_list),1);
peak_9T_length_list = zeros(length(peak_length_list),1);
peak_10T_length_list = zeros(length(peak_length_list),1);
peak_11T_length_list = zeros(length(peak_length_list),1);
peak_14T_length_list = zeros(length(peak_length_list),1);
peak_nT_counts = zeros(10,1);
for i = 1:length(peak_length_list)
    if peak_length_list(i) > peak_split_index(1) && peak_length_list(i) < peak_split_index(2)
        peak_nT_counts(1) = peak_nT_counts(1) + 1;
        peak_3T_length_list(peak_nT_counts(1)) = peak_length_list(i);
    elseif peak_length_list(i) > peak_split_index(2) && peak_length_list(i) < peak_split_index(3)
        peak_nT_counts(2) = peak_nT_counts(2) + 1;
        peak_4T_length_list(peak_nT_counts(2)) = peak_length_list(i);
    elseif peak_length_list(i) > peak_split_index(3) && peak_length_list(i) < peak_split_index(4)
        peak_nT_counts(3) = peak_nT_counts(3) + 1;
        peak_5T_length_list(peak_nT_counts(3)) = peak_length_list(i);
    elseif peak_length_list(i) > peak_split_index(4) && peak_length_list(i) < peak_split_index(5)
        peak_nT_counts(4) = peak_nT_counts(4) + 1;
        peak_6T_length_list(peak_nT_counts(4)) = peak_length_list(i);
    elseif peak_length_list(i) > peak_split_index(5) && peak_length_list(i) < peak_split_index(6)
        peak_nT_counts(5) = peak_nT_counts(5) + 1;
        peak_7T_length_list(peak_nT_counts(5)) = peak_length_list(i);
    elseif peak_length_list(i) > peak_split_index(6) && peak_length_list(i) < peak_split_index(7)
        peak_nT_counts(6) = peak_nT_counts(6) + 1;
        peak_8T_length_list(peak_nT_counts(6)) = peak_length_list(i);
    elseif peak_length_list(i) > peak_split_index(7) && peak_length_list(i) < peak_split_index(8)
        peak_nT_counts(7) = peak_nT_counts(7) + 1;
        peak_9T_length_list(peak_nT_counts(7)) = peak_length_list(i);
    elseif peak_length_list(i) > peak_split_index(8) && peak_length_list(i) < peak_split_index(9)
        peak_nT_counts(8) = peak_nT_counts(8) + 1;
        peak_10T_length_list(peak_nT_counts(8)) = peak_length_list(i);
    elseif peak_length_list(i) > peak_split_index(9) && peak_length_list(i) < peak_split_index(10)
        peak_nT_counts(9) = peak_nT_counts(9) + 1;
        peak_11T_length_list(peak_nT_counts(9)) = peak_length_list(i);
    elseif peak_length_list(i) > peak_split_index(10) && peak_length_list(i) < peak_split_index(11)
        peak_nT_counts(10) = peak_nT_counts(10) + 1;
        peak_14T_length_list(peak_nT_counts(10)) = peak_length_list(i);
    end
end
peak_3T_length_list = peak_3T_length_list(1:peak_nT_counts(1));
peak_4T_length_list = peak_4T_length_list(1:peak_nT_counts(2));
peak_5T_length_list = peak_5T_length_list(1:peak_nT_counts(3));
peak_6T_length_list = peak_6T_length_list(1:peak_nT_counts(4));
peak_7T_length_list = peak_7T_length_list(1:peak_nT_counts(5));
peak_8T_length_list = peak_8T_length_list(1:peak_nT_counts(6));
peak_9T_length_list = peak_9T_length_list(1:peak_nT_counts(7));
peak_10T_length_list = peak_10T_length_list(1:peak_nT_counts(8));
peak_11T_length_list = peak_11T_length_list(1:peak_nT_counts(9));
peak_14T_length_list = peak_14T_length_list(1:peak_nT_counts(10));
peak_nT_length_var = zeros(10,1);
peak_nT_length_std = zeros(10,1);
peak_nT_length_var(1) = var(peak_3T_length_list);
peak_nT_length_var(2) = var(peak_4T_length_list);
peak_nT_length_var(3) = var(peak_5T_length_list);
peak_nT_length_var(4) = var(peak_6T_length_list);
peak_nT_length_var(5) = var(peak_7T_length_list);
peak_nT_length_var(6) = var(peak_8T_length_list);
peak_nT_length_var(7) = var(peak_9T_length_list);
peak_nT_length_var(8) = var(peak_10T_length_list);
peak_nT_length_var(9) = var(peak_11T_length_list);
peak_nT_length_var(10) = var(peak_14T_length_list);
peak_nT_length_std(1) = std(peak_3T_length_list);
peak_nT_length_std(2) = std(peak_4T_length_list);
peak_nT_length_std(3) = std(peak_5T_length_list);
peak_nT_length_std(4) = std(peak_6T_length_list);
peak_nT_length_std(5) = std(peak_7T_length_list);
peak_nT_length_std(6) = std(peak_8T_length_list);
peak_nT_length_std(7) = std(peak_9T_length_list);
peak_nT_length_std(8) = std(peak_10T_length_list);
peak_nT_length_std(9) = std(peak_11T_length_list);
peak_nT_length_std(10) = std(peak_14T_length_list);


valley_3T_length_list = zeros(length(valley_length_list),1);
valley_4T_length_list = zeros(length(valley_length_list),1);
valley_5T_length_list = zeros(length(valley_length_list),1);
valley_6T_length_list = zeros(length(valley_length_list),1);
valley_7T_length_list = zeros(length(valley_length_list),1);
valley_8T_length_list = zeros(length(valley_length_list),1);
valley_9T_length_list = zeros(length(valley_length_list),1);
valley_10T_length_list = zeros(length(valley_length_list),1);
valley_11T_length_list = zeros(length(valley_length_list),1);
valley_14T_length_list = zeros(length(valley_length_list),1);
valley_nT_counts = zeros(10,1);
for i = 1:length(valley_length_list)
    if valley_length_list(i) > valley_split_index(1) && valley_length_list(i) < valley_split_index(2)
        valley_nT_counts(1) = valley_nT_counts(1) + 1;
        valley_3T_length_list(valley_nT_counts(1)) = valley_length_list(i);
    elseif valley_length_list(i) > valley_split_index(2) && valley_length_list(i) < valley_split_index(3)
        valley_nT_counts(2) = valley_nT_counts(2) + 1;
        valley_4T_length_list(valley_nT_counts(2)) = valley_length_list(i);
    elseif valley_length_list(i) > valley_split_index(3) && valley_length_list(i) < valley_split_index(4)
        valley_nT_counts(3) = valley_nT_counts(3) + 1;
        valley_5T_length_list(valley_nT_counts(3)) = valley_length_list(i);
    elseif valley_length_list(i) > valley_split_index(4) && valley_length_list(i) < valley_split_index(5)
        valley_nT_counts(4) = valley_nT_counts(4) + 1;
        valley_6T_length_list(valley_nT_counts(4)) = valley_length_list(i);
    elseif valley_length_list(i) > valley_split_index(5) && valley_length_list(i) < valley_split_index(6)
        valley_nT_counts(5) = valley_nT_counts(5) + 1;
        valley_7T_length_list(valley_nT_counts(5)) = valley_length_list(i);
    elseif valley_length_list(i) > valley_split_index(6) && valley_length_list(i) < valley_split_index(7)
        valley_nT_counts(6) = valley_nT_counts(6) + 1;
        valley_8T_length_list(valley_nT_counts(6)) = valley_length_list(i);
    elseif valley_length_list(i) > valley_split_index(7) && valley_length_list(i) < valley_split_index(8)
        valley_nT_counts(7) = valley_nT_counts(7) + 1;
        valley_9T_length_list(valley_nT_counts(7)) = valley_length_list(i);
    elseif valley_length_list(i) > valley_split_index(8) && valley_length_list(i) < valley_split_index(9)
        valley_nT_counts(8) = valley_nT_counts(8) + 1;
        valley_10T_length_list(valley_nT_counts(8)) = valley_length_list(i);
    elseif valley_length_list(i) > valley_split_index(9) && valley_length_list(i) < valley_split_index(10)
        valley_nT_counts(9) = valley_nT_counts(9) + 1;
        valley_11T_length_list(valley_nT_counts(9)) = valley_length_list(i);
    elseif valley_length_list(i) > valley_split_index(10) && valley_length_list(i) < valley_split_index(11)
        valley_nT_counts(10) = valley_nT_counts(10) + 1;
        valley_14T_length_list(valley_nT_counts(10)) = valley_length_list(i);
    end
end
valley_3T_length_list = valley_3T_length_list(1:valley_nT_counts(1));
valley_4T_length_list = valley_4T_length_list(1:valley_nT_counts(2));
valley_5T_length_list = valley_5T_length_list(1:valley_nT_counts(3));
valley_6T_length_list = valley_6T_length_list(1:valley_nT_counts(4));
valley_7T_length_list = valley_7T_length_list(1:valley_nT_counts(5));
valley_8T_length_list = valley_8T_length_list(1:valley_nT_counts(6));
valley_9T_length_list = valley_9T_length_list(1:valley_nT_counts(7));
valley_10T_length_list = valley_10T_length_list(1:valley_nT_counts(8));
valley_11T_length_list = valley_11T_length_list(1:valley_nT_counts(9));
valley_14T_length_list = valley_14T_length_list(1:valley_nT_counts(10));
valley_nT_length_var = zeros(10,1);
valley_nT_length_std = zeros(10,1);
valley_nT_length_var(1) = var(valley_3T_length_list);
valley_nT_length_var(2) = var(valley_4T_length_list);
valley_nT_length_var(3) = var(valley_5T_length_list);
valley_nT_length_var(4) = var(valley_6T_length_list);
valley_nT_length_var(5) = var(valley_7T_length_list);
valley_nT_length_var(6) = var(valley_8T_length_list);
valley_nT_length_var(7) = var(valley_9T_length_list);
valley_nT_length_var(8) = var(valley_10T_length_list);
valley_nT_length_var(9) = var(valley_11T_length_list);
valley_nT_length_var(10) = var(valley_14T_length_list);
valley_nT_length_std(1) = std(valley_3T_length_list);
valley_nT_length_std(2) = std(valley_4T_length_list);
valley_nT_length_std(3) = std(valley_5T_length_list);
valley_nT_length_std(4) = std(valley_6T_length_list);
valley_nT_length_std(5) = std(valley_7T_length_list);
valley_nT_length_std(6) = std(valley_8T_length_list);
valley_nT_length_std(7) = std(valley_9T_length_list);
valley_nT_length_std(8) = std(valley_10T_length_list);
valley_nT_length_std(9) = std(valley_11T_length_list);
valley_nT_length_std(10) = std(valley_14T_length_list);

x_values = min_length:max_length;
figure;
plot(x_values, peak_counts, '-*', 'Color', 'red', 'DisplayName', '岸', 'LineWidth', 1.5, 'MarkerSize', 5);
hold on;
peak_peak_x_values = [170 218 270 321 376 430 488 542 601 774];
peak_peak_y_values = [3304 3559 2782 1698 1099 637 362 206 43 57];
for i = 1:10
    % 保留两位小数
    peak_nT_length_std(i) = peak_nT_length_std(i)/peak_peak_x_values(10)*14*400/3;
    peak_nT_length_std(i) = round(peak_nT_length_std(i));
    text(peak_peak_x_values(i)-5, peak_peak_y_values(i)+200, "std = " + peak_nT_length_std(i) + " nm", 'FontSize', 12, 'FontName', 'Times New Roman');
end
valley_valley_x_values = [169 230 284 335 391 448 505 555 618 783];
valley_valley_y_values = [3702 3430 2427 1588 945 593 330 197 32 61];
for i = 1:10
    % 保留两位小数
    valley_nT_length_std(i) = valley_nT_length_std(i)/peak_peak_x_values(10)*14*400/3;
    valley_nT_length_std(i) = round(valley_nT_length_std(i));
    text(valley_valley_x_values(i)-5, -valley_valley_y_values(i)-200, "std = " + valley_nT_length_std(i) + " nm", 'FontSize', 12, 'FontName', 'Times New Roman');
end
plot(x_values, -valley_counts, '-*', 'Color', 'blue', 'DisplayName', '坑', 'LineWidth', 1.5, 'MarkerSize', 5);
hold off;
legend('FontSize', 16, 'FontName', 'Songti SC');
xlabel('记录符长度（插值后采样点个数）', 'FontSize', 20, 'FontName', 'Songti SC');
ylabel('记录符个数', 'FontSize', 20, 'FontName', 'Songti SC');
ylim([-5000, 5000]);
title('记录符分布（自采数据，标准光头，标准伺服系统）', 'FontSize', 25, 'FontName', 'Songti SC');
grid on;