%轨道批次同步化
close all;clear ;clc
load('beijing_pollution_event.mat', 'bj_site_name')
load('beijing_pollution_event.mat', 'pollution_events_1')
load('beijing_pollution_event.mat', 'tt')
%% 预处理
num_sm = 5;%采用滤波处理阶数
pollution_events_1_dist = zeros(size(pollution_events_1,1),size(pollution_events_1,1));%初始化
pollution_events_1_7para = cell(size(pollution_events_1,1),1);
%选取原数据7个变量：{'pm25';'no2';'so2';'T@2m';'RH';'WD@10m';'WS@10m'}
parameter7_name = {'pm25';'no2';'so2';'T@2m';'RH';'WD@10m';'WS@10m'};
for i0 = 1:size(pollution_events_1,1)
    pollution_events_1_7para{i0,1} = [pollution_events_1{i0,1}(:,2),pollution_events_1{i0,1}(:,4),pollution_events_1{i0,1}(:,7:8),pollution_events_1{i0,1}(:,10:12)];%7个变量
end
% 得到pollution_events_1_7para

% 采用滤波处理
for ii = 1:size(pollution_events_1_7para,1)
    for jj = 1:7
        pollution_events_1_7para{ii,1}(:,jj) = smooth(pollution_events_1_7para{ii,1}(:,jj),num_sm);
    end
end
% 得到滤波处理 pollution_events_1_7para
% close all;
% plot(pollution_events_1_7para{1, 1}(:,1))
% hold on;
% plot(pollution_events_1_7para{7, 1}(:,1))
% hold on;
% plot(pollution_events_1_7para{10, 1}(:,1))
% hold on;
% plot(pollution_events_1_7para{11, 1}(:,1))
% hold on;
% plot(pollution_events_1_7para{13, 1}(:,1))
% hold on;

% plot(pollution_events_1_7para{5, 1}(:,1))
% hold on;

% 计算模板帧匹配矩阵距离
% i行,j列
% 每一行为一个事件，共i个事件，每行中的每一项数据为以其做模板时的距离
tic
for i1 = 1:size(pollution_events_1_7para,1)
    for j1 = 1:size(pollution_events_1_7para,1)
        pollution_events_1_dist(i1,j1) = dtw_yjy(pollution_events_1_7para{i1,1}(25:end-24,1),pollution_events_1_7para{j1,1}(25:end-24,1));
    end
    i2 = size(pollution_events_1_7para,1)-i1;
    fprintf('pollution_events_1_dist剩余 %.0f\n',i2);
end
toc
%%
load('result_step1.mat', 'pollution_events_all')
pollution_events_all_dist = zeros(size(pollution_events_all,1));
tic
for i1 = 1:size(pollution_events_all_11para,1)
    for j1 = 1:size(pollution_events_all_11para,1)
        pollution_events_all_dist(i1,j1) = dtw_yjy(pollution_events_all_11para{i1,1}(25:end-24,1),pollution_events_all_11para{j1,1}(25:end-24,1));
    end
    i2 = size(pollution_events_all_11para,1)-i1;
    fprintf('pollution_events_all_dist剩余 %.0f\n',i2);
end
toc
% 计算得到距离矩阵 pollution_events_all_dist
%% 事件_自动分类
pollution_events_1_dist_sum = sum(pollution_events_1_dist,2);%求得每一行的总距离
[~,num_sort] = sort(pollution_events_1_dist_sum(:,1));%排序
pollution_events_1_dist_sum(num_sort,2) = (1:length(num_sort))';
%选取模板事件
num_rank = 1;
[num_Fevent,~] = find(pollution_events_1_dist_sum(:,2)==num_rank);%总距离最小对应的事件
Fevent = pollution_events_1_7para{num_Fevent,1};  
plot(Fevent(:,1));

% num_rank = 5;%选取按相对距离排名的事件
% Fevent = cell(num_rank,1);
% for i=1:num_rank
%     [num_Fevent,~] = find(pollution_events_1_dist_sum(:,2)==i);%总距离最小对应的事件
%     Fevent{i,1} = pollution_events_1_7para{num_Fevent,1};   
% end
% figure(1);
% plot(Fevent{1,1}(:,1));
% figure(2);
% plot(Fevent{2,1}(:,1))
% figure(3);
% plot(Fevent{3,1}(:,1))
% figure(4);
% plot(Fevent{4,1}(:,1))
% figure(5);
% plot(Fevent{5,1}(:,1))


%模板事件 Fevent
% 分类事件组 Fevent_club
num_Fevent_club = 15;% 选取模板事件组个数
Fevent_club = cell(num_Fevent_club,1);
Fevent_club_dist = pollution_events_1_dist(num_Fevent,:)';
[~,num_sort1] = sort(Fevent_club_dist(:,1));%排序
Fevent_club_dist(num_sort1,2) = (1:length(num_sort1))';
for i3 = 1:num_Fevent_club
    [num_sort2,~] = find(Fevent_club_dist(:,2)==i3);
    Fevent_club{i3,1} = pollution_events_1_7para{num_sort2,1};
    plot(pollution_events_1_7para{num_sort2,1}(:,1));%作图
    hold on
end

%% 轨道批次同步化DTW
% pollution_events_1_7para{i1,1}(:,2:end),pollution_events_1_7para{j1,1}(:,2:end)
pollution_events_1_dtw = cell(size(pollution_events_1_7para,1),1);%初始化
%
for i2 = 1:size(pollution_events_1_7para,1)
    t = Fevent;
    r = pollution_events_1_7para{i2,1};
    r_new = dtw_yjy_deal(t(:,1),r(:,1));
    pollution_events_1_dtw{i2,1} = r_new;
    %i3 = size(pollution_events_1_7para,1)-i2;
    %fprintf('pollution_events_dtw 剩余 %.0f\n',i3);
end
parameter6_name = {'no2';'so2';'T@2m';'RH';'WS@10m';'HV'};
%求得等长度的事件 pollution_events_dtw

% 分类事件组 Fevent_club
% num_Fevent_club = 20;% 选取模板事件组个数
% Fevent_club = cell(num_Fevent_club,1);
% Fevent_club_dist = pollution_events_1_dist(num_Fevent,:)';
% [~,num_sort1] = sort(Fevent_club_dist(:,1));%排序
% Fevent_club_dist(num_sort1,2) = (1:length(num_sort1))';
for i3 = 1:num_Fevent_club
    [num_sort2,~] = find(Fevent_club_dist(:,2)==i3);
    Fevent_club{i3,1} = pollution_events_1_dtw{num_sort2,1};
    plot(pollution_events_1_dtw{num_sort2,1}(:,1));%作图
    hold on
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fevent_club_dist_s = sort(Fevent_club_dist(:,1));
% % plot(Fevent_club_dist_s(1:end))
% plot(Fevent_club_dist_s(1:end-1))
%作图
%判断选取的事件Fevent，对应其他事件的距离


% plot(pollution_events_1_dtw{1, 1}(:,1))
% hold on;
% plot(pollution_events_1_dtw{7, 1}(:,1))
% hold on;
% plot(pollution_events_1_dtw{10, 1}(:,1))
% hold on;
% plot(pollution_events_1_dtw{11, 1}(:,1))
% hold on;
% plot(pollution_events_1_dtw{13, 1}(:,1))
% hold on;