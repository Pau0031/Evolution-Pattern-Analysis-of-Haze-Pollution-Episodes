% stage 识别阶段

close all; clear;clc;
%%
% 导入
% J变量，K时间，I批次
% 修改pollution_events_10para，pollution_events_site_all，bj_site_all
cd('H:\works in SWU\Haze Early Warning Research\PM2.5 Early Warning System\MPCA_step\重要程序-于君毅\beijing_pollution_event_site')
load('result_step2_3.mat', 'Fevent')
load('result_step2_3.mat', 'pollution_events_club1')
load('result_step2_3.mat', 'pollution_events_club2')
load('result_step2_3.mat', 'pollution_events_club3')
load('result_step2_3.mat', 'pollution_events_club4')
load('result_step2_3.mat', 'pollution_events_club5')
load('result_step2_3.mat', 'club_mean')
%% STAGE 识别阶段
% 选取6个变量parameter6_name = {'no2';'so2';'T@2m';'RH';'WS@10m';'HV'};
% parameter6_name = {'no2';'so2';'T@2m';'RH';'WS@10m';'HV'};
%选取原数据11个变量：{'pm25';'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tar_club = pollution_events_club5;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
parameter10_name = {'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'};
pollution_events_10para = cell(size(tar_club,1),1);
for i0 = 1:size(pollution_events_10para,1)
    pollution_events_10para{i0,1} = tar_club{i0,2}(:,2:end);
end
% 求得 pollution_events_10para

% 选取合适的 pollution_events_10para 分阶段
[I,~] = size(pollution_events_10para);
[K,J] = size(pollution_events_10para{1,1});
% [K,~] = size(pollution_events_10para{1,1});
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
bj_site_all = zeros(I*K,J);%unfold variablely 
for i1 = 1:I
    for k1 = 1:K
        for j1 = 1:J
            bj_site_all((i1-1)*K+k1,j1) = pollution_events_site_all(i1,k1,j1);% I批次*K时间，J变量
        end
    end
end
% 二维数据
fprintf('Step 1 is OK...>>>>>>\n')
%%
% 自动识别阶段
num_PC = 1 ;%指定主元个数
JD = ones(K,2)*-99; %整个阶段，时间长度为K
% 第一列为PC贡献，第二列为其增值，第三列为增值排序
% 创建子阶段cell
for tt = 1:K  %从1到K每一个分析主元
    bj_site_all = zeros(I*tt,J);%赋初值
    for i1 = 1:I
        for k1 = 1:tt
            for j1 = 1:J
                bj_site_all((i1-1)*tt+k1,j1) = pollution_events_site_all(i1,k1,j1);% I批次*K时间，J变量
            end
        end
    end
    [U, latent, V_t, frac_var] = pca_yjy(bj_site_all);
    %对bj_site_all进行主元分析
    JD(tt,1) = sum(frac_var(1:num_PC));%第一列为num_PC个主元的贡献率
    if tt == 1
        JD(tt,2) = 0;
    else
        JD(tt,2) = JD(tt,1)-JD(tt-1,1);%第二列为第一列增值bj_site_c1
    end
end
fprintf('Step 2 is OK...>>>>>>\n')
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
JD_club5 = JD;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% [~,tt_2] = sort(-JD(:,2));%降序排列，返回位置
% JD(tt_2,3) = (1:length(tt_2))';
% 
% figure(1)
% plot(JD(:,2));
% hold on
% figure(2)
% plot(club_mean{1, 1}(:,1))
% hold on
%%
% tar_club(:,4) = [];
% TF = islocalmax(JD(:,2),'MaxNumExtrema',5);
% TF([3,14,52]) = 0;
% TF(find(JD(:,2)== min(JD(:,2)),1)) = 1;
TF = zeros(length(JD),1);
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
TF(pollution_events_club5{1,4}(:,1)) = 1;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
TF = logical(TF);
key_point_stages = find(TF ==1);
tar_pm25 = zeros(length(JD),size(tar_club,1));
%find the start_end point of each stage
for i = 1:size(tar_club,1)
    tar_pm25(:,i) = tar_club{i,2}(:,1);
    tar_club{i,4}(:,1) = key_point_stages;
    for j = 1:length(key_point_stages)
        tmp = find(tar_club{i,3}(key_point_stages(j),:)==1);
        tar_club{i,4}(j,2) = tmp(end);
    end
end
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
pollution_events_club5 = tar_club;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clear i* j* k* tmp K*
%%
clc;
% TF = zeros(length(JD),1);
% TF(pollution_events_club5{1,4}(:,1)) = 1;
% TF = logical(TF);
subplot(3,2,5)
x = 1:length(JD);
m_std = ones(length(JD),2)*-99;
m_std(:,1) = mean(tar_pm25,2);
m_std(:,2) = std(tar_pm25,0,2);
% figure(404)
title('Evolution pattern of haze polution episodes of the 5^{th} cluster in Beijing')
hold on 
yyaxis left
plot(x,JD(:,2),'k-',x(TF),JD(TF,2),'r*')
xline(x(TF),'k--');
xlabel('Duration time (hour)')
ylabel('Change rate of \itCPV\rm_1')
yyaxis right
%plot(x,tar_pm25(:,1:5));
errorbar(x,m_std(:,1),m_std(:,2),'b-.o'); 
ylabel('PM_{2.5} concentration (\mug/m^3)')
yline(115,'r-')
% xlim([10 75])
hold off
%% online identify
clc;
close all
JD = zeros(size(pollution_events_club1{500,1},1),3);
for i = 1:size(pollution_events_club1{500,1},1)

    tar_mat = pollution_events_club1{500,1}(1:i,:);
    if size(tar_mat,1)<=71
        tar_mat(i+1:71,:) = zeros(71-i,size(tar_mat,2));
        %tar_mat(i+1:71,:) = repmat(tar_mat(end,:),71-i,1);
    end
    [U, latent, V_t, frac_var] = pca_yjy(tar_mat);
    JD(i,1) = sum(frac_var(1:num_PC));
    if i == 1
        JD(i,2) = 0;
    else
        JD(i,2) = JD(i,1)-JD(i-1,1);
    end
    [~,tt_2] = sort(-JD(:,2));%降序排列，返回位置
    JD(tt_2,3) = (1:length(tt_2))';
end
x = 1:size(tar_mat,1);
TF = islocalmax(JD(:,2),'MaxNumExtrema',4);
figure(404)
hold on
yyaxis left
plot(x,JD(:,2),x(TF),JD(TF,2),'r*');
yyaxis right
plot(x,tar_mat(:,1));
xline(x(TF));
hold off