% stage 识别阶段

close all; clear;clc;
% 导入
% J变量，K时间，I批次
% 修改pollution_events_10para，pollution_events_site_all，bd_site_all
load('pollution_events_all_fenlei_11para.mat', 'Fevent')
load('pollution_events_all_fenlei_11para.mat', 'pollution_events_club1')
load('pollution_events_all_fenlei_11para.mat', 'pollution_events_club2')
load('pollution_events_all_fenlei_11para.mat', 'pollution_events_club3')
load('pollution_events_all_fenlei_11para.mat', 'pollution_events_club4')
load('pollution_events_all_fenlei_11para.mat', 'club_mean')

%% STAGE 识别阶段
% 选取6个变量parameter6_name = {'no2';'so2';'T@2m';'RH';'WS@10m';'HV'};
% parameter6_name = {'no2';'so2';'T@2m';'RH';'WS@10m';'HV'};
%选取原数据11个变量：{'pm25';'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'}
parameter10_name = {'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'};
pollution_events_10para = cell(size(pollution_events_club4,1),1);
for i0 = 1:size(pollution_events_10para,1)
    pollution_events_10para{i0,1} = pollution_events_club4{i0,2}(:,2:end);
end
% 求得 pollution_events_10para

% 选取合适的 pollution_events_10para 分阶段
[I,~] = size(pollution_events_10para);
[~,J] = size(pollution_events_10para{1,1});
[K,~] = size(pollution_events_10para{1,1});
% I J K 
% 基于变量的展开方式
pollution_events_site_all = zeros(I,K,J);
for i1 = 1:I
    for j1 = 1:J
        for k1 = 1:K
            if isnan(pollution_events_10para{i1,1}(k1,j1))
                pollution_events_10para{i1,1}(k1,j1)=0;
            end
            pollution_events_site_all(i1,k1,j1) = pollution_events_10para{i1,1}(k1,j1);% I批次，K时间，J变量
        end
    end
end
% 转化
bd_site_all = zeros(I*K,J);
for i1 = 1:I
    for k1 = 1:K
        for j1 = 1:J
            bd_site_all((i1-1)*K+k1,j1) = pollution_events_site_all(i1,k1,j1);% I批次*K时间，J变量
        end
    end
end
% 二维数据
clear i1 j1 k1 tt
fprintf('Step 1 is OK...>>>>>>\n')

% 自动识别阶段
num_PC = 1 ;%指定主元个数
JD = zeros(K,3); %整个阶段，时间长度为K
% 第一列为PC贡献，第二列为其增值，第三列为增值排序
% 创建子阶段cell
for tt = 1:K  %从1到K每一个分析主元
    bd_site_all = zeros(I*tt,J);%赋初值
    for i1 = 1:I
        for k1 = 1:tt
            for j1 = 1:J
                bd_site_all((i1-1)*tt+k1,j1) = pollution_events_site_all(i1,k1,j1);% I批次*K时间，J变量
            end
        end
    end
    [U, latent, V_t, frac_var] = pca_yjy(bd_site_all);
    %对bd_site_all进行主元分析
    JD(tt,1) = sum(frac_var(1:num_PC));%第一列为num_PC个主元的贡献率
    if tt == 1
        JD(tt,2) = 0;
    else
        JD(tt,2) = JD(tt,1)-JD(tt-1,1);%第二列为第一列增值bd_site_c1
    end
end

[~,tt_2] = sort(-JD(:,2));
JD(tt_2,3) = (1:length(tt_2))';

figure(1)
plot(JD(:,2));
% hold on
figure(2)
plot(club_mean{4, 1}(:,1))
% hold on
fprintf('Step 2 is OK...>>>>>>\n')
   