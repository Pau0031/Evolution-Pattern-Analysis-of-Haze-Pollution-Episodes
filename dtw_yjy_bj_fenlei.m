close all;clear ;clc
load('beijing_pollution_event.mat', 'bj_site_name')
load('beijing_pollution_event.mat', 'pollution_events_all')
load('beijing_pollution_event.mat', 'tt')
%% 预处理
num_sm = 5;%采用滤波处理阶数
% pollution_events_all_7para = cell(size(pollution_events_all,1),1);
pollution_events_all_11para = cell(size(pollution_events_all,1),1);
%选取原数据7个变量：{'pm25';'no2';'so2';'T@2m';'RH';'WD@10m';'WS@10m'}
parameter11_name = {'pm25';'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'};
% parameter7_name = {'pm25';'no2';'so2';'T@2m';'RH';'WD@10m';'WS@10m'};
for i0 = 1:size(pollution_events_all,1)
    pollution_events_all_11para{i0,1} = pollution_events_all{i0,1}(:,[2:5 7:13]);
    %     pollution_events_all_7para{i0,1} = [pollution_events_all{i0,1}(:,2),pollution_events_all{i0,1}(:,4),pollution_events_all{i0,1}(:,7:8),pollution_events_all{i0,1}(:,10:12)];%7个变量
end
% 得到pollution_events_all_7para

% 采用滤波处理
for ii = 1:size(pollution_events_all_11para,1)
    for jj = 1:11
        pollution_events_all_11para{ii,1}(:,jj) = smooth(pollution_events_all_11para{ii,1}(:,jj),num_sm);
    end
end
% plot(pollution_events_all_7para{1, 1}(:,1))
% hold on
% 计算模板帧匹配矩阵距离
% i行,j列
% 每一行为一个事件，共i个事件，每行中的每一项数据为以其做模板时的距离
%% 事件_自动分类 （规定num_Fevent_club）
class = 0;
num_Fevent_club = 200;% 选取模板事件组个数（每一类事件数）
while size(pollution_events_all_11para,1)>num_Fevent_club
    class = class+1;
    disp(num2str(class));
    pollution_events_all_dist = zeros(size(pollution_events_all_11para,1),size(pollution_events_all_11para,1));%初始化
    for i1 = 1:size(pollution_events_all_11para,1)
        for j1 = 1:size(pollution_events_all_11para,1)
            pollution_events_all_dist(i1,j1) = dtw_yjy(pollution_events_all_11para{i1,1}(:,1),pollution_events_all_11para{j1,1}(:,1));
        end
    %     i2 = size(pollution_events_all_7para,1)-i1;
    %     fprintf('pollution_events_all_dist剩余 %.0f\n',i2);
    end
    pollution_events_all_dist_sum = sum(pollution_events_all_dist,2);%求得每一行的总距离
    [~,num_sort] = sort(pollution_events_all_dist_sum(:,1));%排序
    pollution_events_all_dist_sum(num_sort,2) = (1:length(num_sort))';
    num_rank = 1;%选取距离排最小的值
    [num_Fevent,~] = find(pollution_events_all_dist_sum(:,2)==num_rank);%总距离最小对应的事件
    Fevent = pollution_events_all_11para{num_Fevent,1};  
    %模板事件 Fevent
    % 分类事件组 pollution_events_class_club
    % eval(strcat('pollution_events_',num2str(class),'_club = cell(num_Fevent_club,2);'));
    % 根据距离求 Fevent_club_ini
    Fevent_club_ini = cell(num_Fevent_club,1);
    Fevent_club_dist = pollution_events_all_dist(num_Fevent,:)';
    [~,num_sort1] = sort(Fevent_club_dist(:,1));%排序
    Fevent_club_dist(num_sort1,2) = (1:length(num_sort1))';
    for i3 = 1:num_Fevent_club
        [num_sort2,~] = find(Fevent_club_dist(:,2)==i3);
        Fevent_club_ini{i3,1} = pollution_events_all_11para{num_sort2,1};
        eval(strcat('pollution_events_club',num2str(class),'{i3,1} = Fevent_club_ini{i3,1};'));
    end
    % 轨道批次同步化求 Fevent_club_dtw
    Fevent_club_dtw = cell(size(Fevent_club_ini,1),1);%初始化
    for i2 = 1:size(Fevent_club_dtw,1)
        t = Fevent;
        r = Fevent_club_ini{i2,1};
        r_new = dtw_yjy_deal(t(:,1),r(:,1));
        Fevent_club_dtw{i2,1} = r_new;  
        eval(strcat('pollution_events_club',num2str(class),'{i2,2} = Fevent_club_dtw{i2,1};'));  
    end
    % 清空
    for i5 = 1:num_Fevent_club
        [num_sort3,~] = find(Fevent_club_dist(:,2)==i5);
        pollution_events_all_11para{num_sort3,1} = [];
    end
    f = cellfun(@isempty,pollution_events_all_11para);
    pollution_events_all_11para = pollution_events_all_11para(~f);
    % 剩余 pollution_events_all_11para     
end
pollution_events_clubelse = pollution_events_all_11para;% 剩余其他类 pollution_events_all_7para

%% 事件_自动分类 规定dist_max
close all;clear;clc
load('pollution_events_1_dist.mat', 'pollution_events_all_dist')
load('pollution_events_1_dist.mat', 'pollution_events_all_7para')
class = 0;
dist_max = 2*10^5;
%% num_Fevent_club = 200;% 选取模板事件组个数（每一类事件数）
% while size(pollution_events_all_7para,1)>num_Fevent_club
class = class+1;
% disp(num2str(class));
% pollution_events_all_dist = zeros(size(pollution_events_all_7para,1),size(pollution_events_all_7para,1));%初始化
% for i1 = 1:size(pollution_events_all_7para,1)
%     for j1 = 1:size(pollution_events_all_7para,1)
%         pollution_events_all_dist(i1,j1) = dtw_yjy(pollution_events_all_7para{i1,1}(:,1),pollution_events_all_7para{j1,1}(:,1));
%     end
% %     i2 = size(pollution_events_all_7para,1)-i1;
% %     fprintf('pollution_events_all_dist剩余 %.0f\n',i2);
% end

% 选模板事件
pollution_events_all_dist_sum = sum(pollution_events_all_dist,2);%求得每一行的总距离
[~,num_sort] = sort(pollution_events_all_dist_sum(:,1));%排序
pollution_events_all_dist_sum(num_sort,2) = (1:length(num_sort))';
num_rank = 1;%选取距离排最小的值
[num_Fevent,~] = find(pollution_events_all_dist_sum(:,2)==num_rank);%总距离最小对应的事件
Fevent = pollution_events_all_7para{num_Fevent,1};  

% 分类事件组 pollution_events_class_club
% 根据距离求 Fevent_club_ini
dist_num = [];
for j = 1:size(pollution_events_all_dist,2)
    if pollution_events_all_dist(num_Fevent,j) <=dist_max
        dist_num = [dist_num,j];% 求第一类位置
    end
end
Fevent_club_ini = cell(size(dist_num,2),1);%初始化
% 求 Fevent_club_ini
eval(strcat('pollution_events_',num2str(class),'_club = cell(size(dist_num,2),2);'));
for i = 1:size(dist_num,2)
    Fevent_club_ini{i,1} = pollution_events_all_7para{dist_num(1,i),1};  
    eval(strcat('pollution_events_club',num2str(class),'{i,1} = Fevent_club_ini{i,1};'));
end
% 轨道批次同步化求 Fevent_club_dtw
Fevent_club_dtw = cell(size(Fevent_club_ini,1),1);%初始化
for i2 = 1:size(Fevent_club_dtw,1)
    t = Fevent;
    r = Fevent_club_ini{i2,1};
    r_new = dtw_yjy_deal(t(:,1),r(:,1));
    Fevent_club_dtw{i2,1} = r_new;  
    eval(strcat('pollution_events_club',num2str(class),'{i2,2} = Fevent_club_dtw{i2,1};'));  
end
% 清空
for i5 = 1:size(dist_num,2)
    pollution_events_all_7para{dist_num(1,i5),1} = [];
end
f = cellfun(@isempty,pollution_events_all_7para);
pollution_events_all_7para = pollution_events_all_7para(~f);
% 剩余 pollution_events_all_7para  

pollution_events_all_dist(dist_num,:) = [];
pollution_events_all_dist(:,dist_num) = [];
% 剩余 pollution_events_all_dist

% end
% pollution_events_clubelse = pollution_events_all_7para;% 剩余其他类 pollution_events_all_7para
%% 画图 plot

% plot figure123...
num_figure = 0; % 设初值

while num_figure<2*class    
    % before dtw    
    num_figure = num_figure+1; 
    eval(strcat('figure(',num2str(num_figure),');'));
    eval(strcat('i0=size(pollution_events_club',num2str(num_figure*0.5+0.5),',1);'));
    for i = 1:i0
        eval(strcat('plot(pollution_events_club',num2str(num_figure*0.5+0.5),'{i,1}(:,1));'));
        hold on
    end
    % after dtw
    num_figure = num_figure+1;        
    eval(strcat('figure(',num2str(num_figure),');')); 
    eval(strcat('i0=size(pollution_events_club',num2str(num_figure*0.5),',1);'));
    for i = 1:i0
        eval(strcat('plot(pollution_events_club',num2str(num_figure*0.5),'{i,2}(:,1));'));            
        hold on
    end
end
% figure 1,3,5,7,9...before dtw
% figure 2,4,6,8,10...after dtw






