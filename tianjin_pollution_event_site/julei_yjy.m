% �������
close all;clear ;clc
load('tianjin_pollution_event.mat', 'tj_site_name')
load('tianjin_pollution_event.mat', 'pollution_events_all')
load('tianjin_pollution_event.mat', 'tt')
%% Ԥ����
% num_sm = 5;%�����˲��������
% pollution_events_all_11para = cell(size(pollution_events_all,1),1);
% %ѡȡԭ����7��������{'pm25';'no2';'so2';'T@2m';'RH';'WD@10m';'WS@10m'}
% parameter7_name = {'pm25';'no2';'so2';'T@2m';'RH';'WD@10m';'WS@10m'};
% for i0 = 1:size(pollution_events_all,1)
%     pollution_events_all_11para{i0,1} = [pollution_events_all{i0,1}(:,2),pollution_events_all{i0,1}(:,4),pollution_events_all{i0,1}(:,7:8),pollution_events_all{i0,1}(:,10:12)];%7������
% end
% % �õ�pollution_events_all_11para
% % �����˲�����
% for ii = 1:size(pollution_events_all_11para,1)
%     for jj = 1:size(pollution_events_all_11para,2)
%         pollution_events_all_11para{ii,1}(:,jj) = smooth(pollution_events_all_11para{ii,1}(:,jj),num_sm);
%     end
% end
%
num_sm = 5;%�����˲��������
% pollution_events_all_11para = cell(size(pollution_events_all,1),1);
%ѡȡԭ����11��������{'pm25';'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'}
parameter11_name = {'pm25';'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'};
for i0 = 1:size(pollution_events_all,1)
    pollution_events_all_11para{i0,1} = [pollution_events_all{i0,1}(:,2:5),pollution_events_all{i0,1}(:,7:13)];%7������
end
% �õ�pollution_events_all_11para
%
for ii = 1:size(pollution_events_all_11para,1)
    for jj = 1:size(pollution_events_all_11para{ii,1},2)
        for i = 1:size(pollution_events_all_11para{ii,1},1)
            if isnan(pollution_events_all_11para{ii,1}(i,jj))
                pollution_events_all_11para{ii,1}(i,jj)=0;% ȥ����nan����
            end
        end
        pollution_events_all_11para{ii,1}(:,jj) = smooth(pollution_events_all_11para{ii,1}(:,jj),num_sm);
    end
end
fprintf('Step Ԥ���� is OK...>>>>>>\n')
%%
tic
% ����ģ��֡ƥ��������
% i��,j��
% ÿһ��Ϊһ���¼�����i���¼���ÿ���е�ÿһ������Ϊ������ģ��ʱ�ľ���
pollution_events_all_dist = zeros(size(pollution_events_all_11para,1),size(pollution_events_all_11para,1));%��ʼ��
for i1 = 1:size(pollution_events_all_11para,1)
    for j1 = 1:size(pollution_events_all_11para,1)
        pollution_events_all_dist(i1,j1) = dtw_yjy(pollution_events_all_11para{i1,1},pollution_events_all_11para{j1,1});
    end
    i2 = size(pollution_events_all_11para,1)-i1;
    fprintf('pollution_events_all_distʣ�� %.0f\n',i2);
end
fprintf('Step ֡ƥ����������� is OK...>>>>>>\n')
% %% ��״ͼ--�������
% % �����¼��������
% clear ;clc
% % A=zscore(A);
% % y=pdist(A,'Euclid');%����ŷ�Ͼ���
% % yc=squareform(y);
% %  load('pollution_events_1_dist.mat', 'pollution_events_all_dist')
%  load('pollution_events_all_dist.mat')
%  load('tianjin_pollution_event.mat', 'pollution_events_all_11para')
% z=linkage(pollution_events_all_dist);
% figure(2)
% [h,t]=dendrogram(z);%���ӻ�������
% % plot
% % ��һ���¼� Fevent_class_1
% [class_1,~] = find(t(:,1)==1);%�ܾ�����С��Ӧ���¼�
% Fevent_class_1 = cell(size(class_1,1),1);% ��ʼ��
% for i = 1:size(class_1,1)
%     Fevent_class_1{i,1} = pollution_events_all_11para{class_1(i,1),1};% ��ֵ��һ��
% end
% 
% % figure(2)
% % for i = 1:size(Fevent_class_1,1)
% %     plot(Fevent_class_1{i,1}(:,1));
% %     hold on
% % end
% %% ��1���¼�����
% Fevent_class_1_dist = zeros(size(Fevent_class_1,1),size(Fevent_class_1,1));%��ʼ��
% for i1 = 1:size(Fevent_class_1,1)
%     for j1 = 1:size(Fevent_class_1,1)
%         Fevent_class_1_dist(i1,j1) = dtw_yjy(Fevent_class_1{i1,1}(:,1),Fevent_class_1{j1,1}(:,1));
%     end
%     i2 = size(Fevent_class_1,1)-i1;
%     fprintf('Fevent_class_1_dist %.0f\n',i2);
% end
% %% ��1���¼��������
% clear ;clc
% % A=zscore(A);
% % y=pdist(A,'Euclid');%����ŷ�Ͼ���
% % yc=squareform(y);
% % A=zscore(A);
% % y=pdist(A,'Euclid');%����ŷ�Ͼ���
% % yc=squareform(y);
% load('pollution_events_1_dist.mat', 'Fevent_class_1_dist')
% load('pollution_events_1_dist.mat', 'Fevent_class_1')
% z=linkage(Fevent_class_1_dist);
% figure(3)
% [h,t]=dendrogram(z);%���ӻ�������
% 
% %
% %% Fevent_class_1_dist_lg ȡ����
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
% %% ��1_1���¼��������
% clear ;clc
% % A=zscore(A);
% % y=pdist(A,'Euclid');%����ŷ�Ͼ���
% % yc=squareform(y);
% % A=zscore(A);
% % y=pdist(A,'Euclid');%����ŷ�Ͼ���
% % yc=squareform(y);
% load('pollution_events_1_dist.mat', 'Fevent_class_1_dist')
% load('pollution_events_1_dist.mat', 'Fevent_class_1')
% z=linkage(Fevent_class_1_dist);
% figure(3)
% [h,t]=dendrogram(z);%���ӻ�������
% % ��һ���¼� Fevent_class_1
% [class_1_1,~] = find(t(:,1)==1);%�ܾ�����С��Ӧ���¼�
% Fevent_class_1_1 = cell(size(class_1_1,1),1);% ��ʼ��
% for i = 1:size(class_1_1,1)
%     Fevent_class_1_1{i,1} = Fevent_class_1{class_1_1(i,1),1};% ��ֵ��һ��
% end
% %% ��1_1���¼�����
% Fevent_class_1_1_dist = zeros(size(Fevent_class_1_1,1),size(Fevent_class_1_1,1));%��ʼ��
% for i1 = 1:size(Fevent_class_1_1_dist,1)
%     for j1 = 1:size(Fevent_class_1_1_dist,1)
%         Fevent_class_1_1_dist(i1,j1) = dtw_yjy(Fevent_class_1_1{i1,1}(:,1),Fevent_class_1_1{j1,1}(:,1));
%     end
%     i2 = size(Fevent_class_1_1_dist,1)-i1;
%     fprintf('Fevent_class_1_1_dist ʣ�� %.0f\n',i2);
% end
% % ����ĵ� Fevent_class_1_1_dist
% %% ��1_1���¼��������
% clear ;clc
% % A=zscore(A);
% % y=pdist(A,'Euclid');%����ŷ�Ͼ���
% % yc=squareform(y);
% % A=zscore(A);
% % y=pdist(A,'Euclid');%����ŷ�Ͼ���
% % yc=squareform(y);
% load('pollution_events_1_dist.mat', 'Fevent_class_1_1_dist')
% load('pollution_events_1_dist.mat', 'Fevent_class_1_1')
% z=linkage(Fevent_class_1_1_dist);
% figure(3)
% [h,t]=dendrogram(z);%���ӻ�������
% 
% %
toc
