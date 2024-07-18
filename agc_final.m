clear;
clc;
close all;

file_path = './G13Data/HF_G13';
fp = fopen(file_path, 'r');
A_0 = fread(fp, 'uint8');
fclose(fp);
% A_0 = A_0(1:5000);

% 滤波
sample_rate = 1e8;
b = fir2(88,[0 0.04 0.08 0.12 0.16 0.2 1],[1 10^(1/20) 10^(3/20) 10^(3/20) 10^(1/20) 10^(-3/20) 0]);
% freqz(b,1,1024,sample_rate);
A_0_new = filter(b,1,A_0);
% A_0_new = A_0; 


% 线性插值
x = 1:length(A_0_new);
interp_rate = 6;
new_x = linspace(x(1), x(end), interp_rate * length(x) - interp_rate + 1);
new_A = interp1(x, A_0_new, new_x, 'linear'); % 使用 interp1 进行线性插值
new_A = double(new_A);

% AGC
zero_line = 126.45;
y = new_A' - 126.45;
agc = comm.AGC('AveragingLength',1,'MaxPowerGain',1,'AdaptationStepSize',1e-1);
gain = 35;

r = agc(y);
r = r*gain;
figure
subplot(2,1,1)
plot(y(100:5099))
ylim([-150,150])
% subplot(3,1,2)
% plot(r(100:5099))
ylim([-150,150])
r(r > 128.6) = 128.6;
r(r < -126.4) = -126.4;
subplot(2,1,2)
plot(r(100:5099))
ylim([-150,150])

r = agc(y);
r = r*gain;
figure
subplot(2,1,1)
plot(y)
ylim([-150,150])
% subplot(3,1,2)
% plot(r)
ylim([-150,150])
r(r > 128.6) = 128.6;
r(r < -126.4) = -126.4;
subplot(2,1,2)
plot(r)
ylim([-150,150])
var_origin = var(y(100:end));
var_AGC = var(r(100:end));
new_A = r;

% 检测零点
num_zero = 0;
shift_list = zeros(length(new_A),1);
index_list = zeros(length(new_A),1);
for j = 1:length(new_A)-1
    if new_A(j)*new_A(j+1)<=0
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
max_length = 60*interp_rate;
min_length = min([min(peak_length_list), min(valley_length_list),0]);
peak_counts = zeros(max_length-min_length+1,1);
valley_counts = zeros(max_length-min_length+1,1);
for j = 1:max_length-min_length+1
    peak_counts(j) = sum(peak_length_list==j);
    valley_counts(j) = sum(valley_length_list==j);
end

peak_counts_min = zeros(200,2);
valley_counts_min = zeros(200,2);
peak_counts_min_num = 0;
valley_counts_min_num = 0;
for i = 3:200
    if peak_counts(i)<=peak_counts(i-1) && peak_counts(i)<=peak_counts(i+1) && peak_counts(i-1)<=peak_counts(i-2) && peak_counts(i+1)<=peak_counts(i+2) && peak_counts(i)<peak_counts(i-2) && peak_counts(i)<peak_counts(i+2)
        peak_counts_min_num = peak_counts_min_num + 1;
        peak_counts_min(peak_counts_min_num,1) = i;
        peak_counts_min(peak_counts_min_num,2) = peak_counts(i);
    end
    if valley_counts(i)<=valley_counts(i-1) && valley_counts(i)<=valley_counts(i+1) && valley_counts(i-1)<=valley_counts(i-2) && valley_counts(i+1)<=valley_counts(i+2) && valley_counts(i)<valley_counts(i-2) && valley_counts(i)<valley_counts(i+2)
        valley_counts_min_num = valley_counts_min_num + 1;
        valley_counts_min(valley_counts_min_num,1) = i;
        valley_counts_min(valley_counts_min_num,2) = valley_counts(i);
    end
end
peak_counts_min = peak_counts_min(1:peak_counts_min_num,:);
valley_counts_min = valley_counts_min(1:valley_counts_min_num,:);
error_num = sum(peak_counts_min(:,2))+sum(valley_counts_min(:,2));

% 统计nT长度方差
peak_split_index = [55 84 106 128 152 177 202 226 249 282 351];
valley_split_index = [55 84 106 128 152 177 202 226 249 282 351];

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
peak_peak_x_values = [71 96 118 141 166 190 213 238 260 322];
peak_peak_y_values = [9647 12282 8799 5550 3243 1963 999 559 113 158];
for i = 1:10
    % 保留两位小数
    peak_nT_length_std(i) = peak_nT_length_std(i)/peak_peak_x_values(10)*14*400/3;
    peak_nT_length_std(i) = round(peak_nT_length_std(i));
    text(peak_peak_x_values(i)-5, peak_peak_y_values(i)+300, "std = " + peak_nT_length_std(i) + " nm", 'FontSize', 12, 'FontName', 'Times New Roman');
end
valley_valley_x_values = [71 94 115 138 164 189 213 236 261 331];
valley_valley_y_values = [11590 12871 8661 4932 3078 1829 932 588 112 135];
for i = 1:10
    % 保留两位小数
    valley_nT_length_std(i) = valley_nT_length_std(i)/peak_peak_x_values(10)*14*400/3;
    valley_nT_length_std(i) = round(valley_nT_length_std(i));
    text(valley_valley_x_values(i)-5, -valley_valley_y_values(i)-300, "std = " + valley_nT_length_std(i) + " nm", 'FontSize', 12, 'FontName', 'Times New Roman');
end
plot(x_values, -valley_counts, '-*', 'Color', 'blue', 'DisplayName', '坑', 'LineWidth', 1.5, 'MarkerSize', 5);
hold off;
legend('FontSize', 16, 'FontName', 'Songti SC');
xlabel('记录符长度（插值后采样点个数）', 'FontSize', 20, 'FontName', 'Songti SC');
xlim([0, max_length]);
ylabel('记录符个数', 'FontSize', 20, 'FontName', 'Songti SC');
ylim([-15000, 15000]);
title('记录符分布（HF信号，有均衡器，有AGC）', 'FontSize', 25, 'FontName', 'Songti SC');
grid on;