% 计算距离
close all;clear ;clc
load('tianjin_pollution_event.mat', 'tj_site_name')
load('tianjin_pollution_event.mat', 'pollution_events_all')
load('tianjin_pollution_event.mat', 'tt')
%% 预处理
% num_sm = 5;%采用滤波处理阶数
% pollution_events_all_11para = cell(size(pollution_events_all,1),1);
% %选取原数据7个变量：{'pm25';'no2';'so2';'T@2m';'RH';'WD@10m';'WS@10m'}
% parameter7_name = {'pm25';'no2';'so2';'T@2m';'RH';'WD@10m';'WS@10m'};
% for i0 = 1:size(pollution_events_all,1)
%     pollution_events_all_11para{i0,1} = [pollution_events_all{i0,1}(:,2),pollution_events_all{i0,1}(:,4),pollution_events_all{i0,1}(:,7:8),pollution_events_all{i0,1}(:,10:12)];%7个变量
% end
% % 得到pollution_events_all_11para
% % 采用滤波处理
% for ii = 1:size(pollution_events_all_11para,1)
%     for jj = 1:size(pollution_events_all_11para,2)
%         pollution_events_all_11para{ii,1}(:,jj) = smooth(pollution_events_all_11para{ii,1}(:,jj),num_sm);
%     end
% end
%
num_sm = 5;%采用滤波处理阶数
% pollution_events_all_11para = cell(size(pollution_events_all,1),1);
%选取原数据11个变量：{'pm25';'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'}
parameter11_name = {'pm25';'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'};
for i0 = 1:size(pollution_events_all,1)
    pollution_events_all_11para{i0,1} = [pollution_events_all{i0,1}(:,2:5),pollution_events_all{i0,1}(:,7:13)];%7个变量
end
% 得到pollution_events_all_11para
%
for ii = 1:size(pollution_events_all_11para,1)
    for jj = 1:size(pollution_events_all_11para{ii,1},2)
        for i = 1:size(pollution_events_all_11para{ii,1},1)
            if isnan(pollution_events_all_11para{ii,1}(i,jj))
                pollution_events_all_11para{ii,1}(i,jj)=0;% 去除掉nan数据
            end
        end
        pollution_events_all_11para{ii,1}(:,jj) = smooth(pollution_events_all_11para{ii,1}(:,jj),num_sm);
    end
end
fprintf('Step 预处理 is OK...>>>>>>\n')
%%
tic
% 计算模板帧匹配矩阵距离
% i行,j列
% 每一行为一个事件，共i个事件，每行中的每一项数据为以其做模板时的距离
pollution_events_all_dist = zeros(size(pollution_events_all_11para,1),size(pollution_events_all_11para,1));%初始化
for i1 = 1:size(pollution_events_all_11para,1)
    for j1 = 1:size(pollution_events_all_11para,1)
        pollution_events_all_dist(i1,j1) = dtw_yjy(pollution_events_all_11para{i1,1},pollution_events_all_11para{j1,1});
    end
    i2 = size(pollution_events_all_11para,1)-i1;
    fprintf('pollution_events_all_dist剩余 %.0f\n',i2);
end
fprintf('Step 帧匹配矩阵距离计算 is OK...>>>>>>\n')
% %% 树状图--聚类分析
% % 所有事件聚类分析
% clear ;clc
% % A=zscore(A);
% % y=pdist(A,'Euclid');%计算欧氏距离
% % yc=squareform(y);
% %  load('pollution_events_1_dist.mat', 'pollution_events_all_dist')
%  load('pollution_events_all_dist.mat')
%  load('tianjin_pollution_event.mat', 'pollution_events_all_11para')
% z=linkage(pollution_events_all_dist);
% figure(2)
% [h,t]=dendrogram(z);%可视化聚类树
% % plot
% % 第一类事件 Fevent_class_1
% [class_1,~] = find(t(:,1)==1);%总距离最小对应的事件
% Fevent_class_1 = cell(size(class_1,1),1);% 初始化
% for i = 1:size(class_1,1)
%     Fevent_class_1{i,1} = pollution_events_all_11para{class_1(i,1),1};% 赋值第一类
% end
% 
% % figure(2)
% % for i = 1:size(Fevent_class_1,1)
% %     plot(Fevent_class_1{i,1}(:,1));
% %     hold on
% % end
% %% 第1类事件距离
% Fevent_class_1_dist = zeros(size(Fevent_class_1,1),size(Fevent_class_1,1));%初始化
% for i1 = 1:size(Fevent_class_1,1)
%     for j1 = 1:size(Fevent_class_1,1)
%         Fevent_class_1_dist(i1,j1) = dtw_yjy(Fevent_class_1{i1,1}(:,1),Fevent_class_1{j1,1}(:,1));
%     end
%     i2 = size(Fevent_class_1,1)-i1;
%     fprintf('Fevent_class_1_dist %.0f\n',i2);
% end
% %% 第1类事件聚类分析
% clear ;clc
% % A=zscore(A);
% % y=pdist(A,'Euclid');%计算欧氏距离
% % yc=squareform(y);
% % A=zscore(A);
% % y=pdist(A,'Euclid');%计算欧氏距离
% % yc=squareform(y);
% load('pollution_events_1_dist.mat', 'Fevent_class_1_dist')
% load('pollution_events_1_dist.mat', 'Fevent_class_1')
% z=linkage(Fevent_class_1_dist);
% figure(3)
% [h,t]=dendrogram(z);%可视化聚类树
% 
% %
% %% Fevent_class_1_dist_lg 取对数
% clear ;clc
% load('pollution_events_1_dist.mat', 'Fevent_class_1_dist')
% load('pollution_events_1_dist.mat', 'Fevent_class_1')
% 
% Fevent_class_1_dist_lg = zeros(size(Fevent_class_1_dist,1),size(Fevent_class_1_dist,2));
% for i = 1:size(Fevent_class_1_dist,1)
%     for j = 1:size(Fevent_class_1_dist,2)
%         if Fevent_class_1_dist(i,j) ~=0
%             Fevent_class_1_dist_lg(i,j) = log10(Fevent_class_1_dist(i,j));
%         end
%     end
% end
% % Fevent_class_1_dist_lg
% % z=linkage(Fevent_class_1_dist);
% z=linkage(Fevent_class_1_dist_lg);
% figure(4)
% [h,t]=dendrogram(z);
% 
% %% 第1_1类事件聚类分析
% clear ;clc
% % A=zscore(A);
% % y=pdist(A,'Euclid');%计算欧氏距离
% % yc=squareform(y);
% % A=zscore(A);
% % y=pdist(A,'Euclid');%计算欧氏距离
% % yc=squareform(y);
% load('pollution_events_1_dist.mat', 'Fevent_class_1_dist')
% load('pollution_events_1_dist.mat', 'Fevent_class_1')
% z=linkage(Fevent_class_1_dist);
% figure(3)
% [h,t]=dendrogram(z);%可视化聚类树
% % 第一类事件 Fevent_class_1
% [class_1_1,~] = find(t(:,1)==1);%总距离最小对应的事件
% Fevent_class_1_1 = cell(size(class_1_1,1),1);% 初始化
% for i = 1:size(class_1_1,1)
%     Fevent_class_1_1{i,1} = Fevent_class_1{class_1_1(i,1),1};% 赋值第一类
% end
% %% 第1_1类事件距离
% Fevent_class_1_1_dist = zeros(size(Fevent_class_1_1,1),size(Fevent_class_1_1,1));%初始化
% for i1 = 1:size(Fevent_class_1_1_dist,1)
%     for j1 = 1:size(Fevent_class_1_1_dist,1)
%         Fevent_class_1_1_dist(i1,j1) = dtw_yjy(Fevent_class_1_1{i1,1}(:,1),Fevent_class_1_1{j1,1}(:,1));
%     end
%     i2 = size(Fevent_class_1_1_dist,1)-i1;
%     fprintf('Fevent_class_1_1_dist 剩余 %.0f\n',i2);
% end
% % 计算的到 Fevent_class_1_1_dist
% %% 第1_1类事件聚类分析
% clear ;clc
% % A=zscore(A);
% % y=pdist(A,'Euclid');%计算欧氏距离
% % yc=squareform(y);
% % A=zscore(A);
% % y=pdist(A,'Euclid');%计算欧氏距离
% % yc=squareform(y);
% load('pollution_events_1_dist.mat', 'Fevent_class_1_1_dist')
% load('pollution_events_1_dist.mat', 'Fevent_class_1_1')
% z=linkage(Fevent_class_1_1_dist);
% figure(3)
% [h,t]=dendrogram(z);%可视化聚类树
% 
% %
toc
