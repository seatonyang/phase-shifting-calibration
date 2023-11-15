%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 文件:    phase_shifting_calib.m
% 作者:    Seaton 杨旭彤
% 日期:    2023年11月10日
% 版本:    1.0
% 目的:    此脚本分析一系列相位图像，计算各种参数，如最大相位值、质心等，并展示它们在图像序列中的趋势。
%
% 输入:    确保相位图像存储在指定文件夹中，至少有5张。
%
% 输出:    显示带有最大值和质心标记的相位直方图，以及这些值的趋势在一个独立的图中。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

% 设置包含相位图像的文件夹路径
file_path = 'imgs\';          % 使用真实图像
% file_path = 'simulation\';  % 使用仿真图像

% 读取相位图像文件列表
img_path_list = dir(strcat(file_path, '*.bmp'));
img_num = length(img_path_list);

% 初始化一个元胞数组以存储相位图像
I1 = cell(1, img_num);

% 读取和预处理相位图像
if img_num > 0
    for j = 1:img_num
        img_name = [file_path, int2str(j - 1), '.bmp'];
        I1{j} = double(imread(img_name));
        I1{j} = I1{j}(301:1200, 301:1200); % 仅选择中间部分进行分析
    end
end

% 初始化用于分析的变量
[X, Y] = size(I1{1});
shift = zeros(X, Y);
realshift = zeros(X, Y);
number = length(I1) - 4;
S_phase = zeros(number + 1, 180);
sumshift = zeros(X, Y);
maxavg = zeros(1, number);
center_of_mass = zeros(1, number);

% 为分析创建一个图
figure;

% 对每个图像执行相位分析
for k = 1:number
    % 使用Hari-Harilela算法从强度值计算相位
    D = (I1{k + 4} - I1{k}) ./ (2 * (I1{k + 3} - I1{k + 1}));
    shift = acosd(D);
    realshift = real(shift);
    sumshift = sumshift + realshift;
    
    % 生成相位值的直方图
    for i = 0:179
        count = find(realshift > i & realshift < i + 1);
        S_phase(k, i + 1) = size(count, 1);
    end
    
    % 找到直方图中最大值对应的角度
    maxavg(1, k) = find(S_phase(k, :) == max(S_phase(k, :)));
end

% 显示每个图像的最大相位值
disp(["移相量最大值：" maxavg(1, :)]) 

% 计算每个直方图的质心
for k = 1:number
    histogram_values = S_phase(k, :);
    angles = 0:179;

    % 计算直方图的质心（加权平均值）
    center_of_mass(1, k) = sum(histogram_values .* angles) / sum(histogram_values);
end

% 显示每个图像的质心值
disp(["直方图质心：" center_of_mass]);

% 计算平均值和趋势
aver_angle = sum(maxavg) / number;
aver_ = sqrt(sum((maxavg - aver_angle) .^ 2) / (number - 1));
aver_90 = sqrt(sum((maxavg - 90) .^ 2) / (number - 1));

shiftavg = sumshift / number;

% 生成平均相位值在图像序列中的直方图
for i = 1:180
    count = find(shiftavg >= i - 1 & shiftavg < i);
    S_phase(number + 1, i) = size(count, 1);
end

% 找到平均直方图中最大值对应的角度
avgSphase = find(S_phase(number + 1, :) == max(S_phase(number + 1, :)));
angle = (0:179);
B = (angle - avgSphase) .^ 2 .* S_phase(number + 1, :);
SS = sum((angle - avgSphase) .^ 2 .* S_phase(number + 1, :));

% 计算相位值围绕平均角度的标准差
SSavg=sqrt(SS/(sum(S_phase(number+1,:))));
disp(SSavg)

% 在每个子图中创建显示单个直方图的区域，并在图上标记最大值和质心值
for i = 1:number
    subplot(3, 3, i);
    area(S_phase(i, :));
    set(gca, 'XTick', (0:10:180));
    hold on;
    % 在图上标记最大值和质心值的位置
    plot(maxavg(i), S_phase(i, maxavg(i)), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    plot(center_of_mass(i), S_phase(i, round(center_of_mass(i)) + 1), 'bx', 'MarkerSize', 8, 'LineWidth', 2);
end
hold off;

% 创建一个单独的图，显示最大值和质心值的趋势
figure;
groundtruth = 45*ones(9);
plot(1:number, maxavg, 'r-', 'LineWidth', 2, 'DisplayName', 'Max Value');
hold on;
plot(1:number, center_of_mass, 'b-', 'LineWidth', 2, 'DisplayName', 'Center of Mass');
hold on;
plot(1:number, groundtruth, 'g-', 'LineWidth', 2, 'DisplayName', 'Ground truth');
hold off;
xlabel('图像编号');
ylabel('数值');
title('最大值和质心值趋势');
legend('show');
