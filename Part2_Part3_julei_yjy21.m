%% hierarchical clustering-- GUI to select clusters     
% both average and ward linkage give similar solutions (NMI = 0.75)
% cd('E:\mfile')
% ����·��
close all;clear;clc
load('beijing_pollution_event.mat', 'pollution_events_all_dist');
%%
distance_matrix = pollution_events_all_dist;
% cd('E:\mfile')
% distance_matrix 
hierarchical_linkage = 'average';
hierarchical_fig_num = 404;
plot_styles = {'r-'; 'b-'; 'g-'; 'm-'; 'c-'; 'y-'; 'k-'; 'ro-'; 'bo-'; 'go-'; 'yo-'; 'co-'; 'mo-'; 'ko-'; 'rx-'; 'bx-'; 'gx-'; 'yx-'; 'cx-'; 'mx-'; 'kx-'; 'r+-'; 'b+-'; 'g+-'; 'y+-'; 'c+-'; 'm+-'; 'k+-'};
[Z_link, cop_coef] = calc_dendrogram(distance_matrix, hierarchical_linkage);
dend_fig_data = draw_dendrogram('draw', Z_link, hierarchical_fig_num, plot_styles, 'aggregated_clusters', 0.4 );
title(['Dendrogram for clustering of D^{Avg} using ', hierarchical_linkage, ' linkage. Cophenetic Coef. = ', num2str(cop_coef, 2), '.'])
ylabel('Aggregated Distance')

%cd('E:\North_China_Plain_data12_12\beijing_pollution_event_site1')

%% �¼�����
% close all;clear;clc
% load('pollution_events_all_dist.mat') % ����
%% Ѱ���¼� pollution_events_club*
load('beijing_pollution_event.mat', 'pollution_events_all_11para')
for class = 1:size(aggregated_clusters,1)% ����
    eval(strcat('pollution_events_club',num2str(class),'= cell(size(aggregated_clusters{class,1},2),2);')); 
    eval(strcat('num_events_club = size(pollution_events_club',num2str(class),',1);')); 
    for i = 1:num_events_club
        eval(strcat('pollution_events_club',num2str(class),'{i,1}=pollution_events_all_11para{aggregated_clusters{class,1}(1,i),1};'));          
    end
end
% pollution_events_club*
%% Ѱ�Ҿ��� pollution_events_dist*
for class = 1:size(aggregated_clusters,1)% ����
    eval(strcat('pollution_events_dist',num2str(class),'= zeros(size(aggregated_clusters{class,1},2),size(aggregated_clusters{class,1},2));')); 
    eval(strcat('num_events_club = size(pollution_events_club',num2str(class),',1);')); 
    for i = 1:num_events_club
        for j = 1:num_events_club
            eval(strcat('pollution_events_dist',num2str(class),'(i,j)=pollution_events_all_dist(aggregated_clusters{class,1}(1,i),aggregated_clusters{class,1}(1,j));'));    
        end
    end
end
% pollution_events_dist*

%% Ѱ��ģ���¼� Fevent*
Fevent = cell(size(aggregated_clusters,1),1);
for class = 1:size(aggregated_clusters,1)
    eval(strcat('pollution_events_dist_sum = zeros(size(pollution_events_dist',num2str(class),',1),2);'));
    eval(strcat('pollution_events_dist_sum = sum(pollution_events_dist',num2str(class),',2);'));%���ÿһ�е��ܾ���
    [~,num_sort] = sort(pollution_events_dist_sum(:,1));%����
    pollution_events_dist_sum(num_sort,2) = (1:length(num_sort))';
    num_rank = 1;%ѡȡ��������С��ֵ
    [num_Fevent,~] = find(pollution_events_dist_sum(:,2)==num_rank);%�ܾ�����С��Ӧ���¼�
    eval(strcat('Fevent{class,1} = pollution_events_club',num2str(class),'{num_Fevent,1};'));   
end
%ģ���¼� Fevent*

%% �������ͬ����dtw �� Fevent_club**
% Fevent_club* ��һ��Ϊdtwǰ���ڶ���Ϊdtw��
for class = 1:size(aggregated_clusters,1)
    eval(strcat('num_events_club = size(pollution_events_club',num2str(class),',1);'));
    for i = 1:num_events_club
        t = Fevent{class,1};
        eval(strcat('r = pollution_events_club',num2str(class),'{i,1};'));
        %r_new = dtw_yjy_deal(t,r);
        [D_flag,r_new] = dtw_yjy_deal(t,r);
        eval(strcat('pollution_events_club',num2str(class),'{i,2} = r_new;'));
        eval(strcat('pollution_events_club',num2str(class),'{i,3} = D_flag;'));
    end
end

% ʣ���¼�
for i1 = 1:size(aggregated_clusters,1)
    num_events_club = size(aggregated_clusters{i1,1},2);
    for i2 = 1:num_events_club
        pollution_events_all_11para{aggregated_clusters{i1,1}(1,i2),1} = [];
    end
end
f = cellfun(@isempty,pollution_events_all_11para);
pollution_events_all_11para = pollution_events_all_11para(~f);
pollution_events_clubelse = pollution_events_all_11para;
clear pollution_events_all_11para

% ʣ�� pollution_events_all_11para  

%% Ѱ���¼���ֵ club_mean
club_mean = cell(size(aggregated_clusters,1),1);

for class = 1:size(aggregated_clusters,1)
    eval(strcat('num_events_club = size(pollution_events_club',num2str(class),',1);')); 
    eval(strcat('[hang,lie] = size(pollution_events_club',num2str(class),'{1,2});'));    
    for i2 = 1:lie
        for i3 = 1:hang
            num = 0;
            for i = 1:num_events_club
                eval(strcat('num = num + pollution_events_club',num2str(class),'{i,2}(i3,i2);'));
            end
            club_mean{class,1}(i3,i2) = num/num_events_club;
        end
    end
end
 % ��ֵ club_mean
%% ��ͼ plot

% plot figure123...
num_figure = 0; % ���ֵ

while num_figure<2*class
    % before dtw
    num_figure = num_figure+1;
    subplot(class,2,num_figure)
    %     eval(strcat('figure(',num2str(num_figure),');'));
    eval(strcat('i0=size(pollution_events_club',num2str(num_figure*0.5+0.5),',1);'));
    for i = 1:i0
        eval(strcat('plot(pollution_events_club',num2str(num_figure*0.5+0.5),'{i,1}(:,1));'));
        hold on
    end
    xlabel('Duration time (hour)')
    ylabel('PM_{2.5} value (\mug/m^3)')
    eval(strcat(['title(''Original haze polution episodes of cluster ',num2str(ceil(num_figure/2)),' in Shijiazhuang'')']))
    % after dtw
    num_figure = num_figure+1;
    subplot(class,2,num_figure)
    %     eval(strcat('figure(',num2str(num_figure),');'));
    eval(strcat('i0=size(pollution_events_club',num2str(num_figure*0.5),',1);'));
    for i = 1:i0
        eval(strcat('plot(pollution_events_club',num2str(num_figure*0.5),'{i,2}(:,1));'));
        hold on
    end
    xlabel('Duration time (hour)')
    ylabel('PM_{2.5} value (\mug/m^3)')
    eval(strcat(['title(''Aligned haze polution episodes of cluster ',num2str(ceil(num_figure/2)),' in Shijiazhuang'')']))
end
% figure 1,3,5,7,9...before dtw
% figure 2,4,6,8,10...after dtw