% stage ʶ��׶�

close all; clear;clc;
% ����
% J������Kʱ�䣬I����
% �޸�pollution_events_10para��pollution_events_site_all��bd_site_all
load('pollution_events_all_fenlei_11para.mat', 'Fevent')
load('pollution_events_all_fenlei_11para.mat', 'pollution_events_club1')
load('pollution_events_all_fenlei_11para.mat', 'pollution_events_club2')
load('pollution_events_all_fenlei_11para.mat', 'pollution_events_club3')
load('pollution_events_all_fenlei_11para.mat', 'pollution_events_club4')
load('pollution_events_all_fenlei_11para.mat', 'club_mean')

%% STAGE ʶ��׶�
% ѡȡ6������parameter6_name = {'no2';'so2';'T@2m';'RH';'WS@10m';'HV'};
% parameter6_name = {'no2';'so2';'T@2m';'RH';'WS@10m';'HV'};
%ѡȡԭ����11��������{'pm25';'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'}
parameter10_name = {'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'};
pollution_events_10para = cell(size(pollution_events_club4,1),1);
for i0 = 1:size(pollution_events_10para,1)
    pollution_events_10para{i0,1} = pollution_events_club4{i0,2}(:,2:end);
end
% ��� pollution_events_10para

% ѡȡ���ʵ� pollution_events_10para �ֽ׶�
[I,~] = size(pollution_events_10para);
[~,J] = size(pollution_events_10para{1,1});
[K,~] = size(pollution_events_10para{1,1});
% I J K 
% ���ڱ�����չ����ʽ
pollution_events_site_all = zeros(I,K,J);
for i1 = 1:I
    for j1 = 1:J
        for k1 = 1:K
            if isnan(pollution_events_10para{i1,1}(k1,j1))
                pollution_events_10para{i1,1}(k1,j1)=0;
            end
            pollution_events_site_all(i1,k1,j1) = pollution_events_10para{i1,1}(k1,j1);% I���Σ�Kʱ�䣬J����
        end
    end
end
% ת��
bd_site_all = zeros(I*K,J);
for i1 = 1:I
    for k1 = 1:K
        for j1 = 1:J
            bd_site_all((i1-1)*K+k1,j1) = pollution_events_site_all(i1,k1,j1);% I����*Kʱ�䣬J����
        end
    end
end
% ��ά����
clear i1 j1 k1 tt
fprintf('Step 1 is OK...>>>>>>\n')

% �Զ�ʶ��׶�
num_PC = 1 ;%ָ����Ԫ����
JD = zeros(K,3); %�����׶Σ�ʱ�䳤��ΪK
% ��һ��ΪPC���ף��ڶ���Ϊ����ֵ��������Ϊ��ֵ����
% �����ӽ׶�cell
for tt = 1:K  %��1��Kÿһ��������Ԫ
    bd_site_all = zeros(I*tt,J);%����ֵ
    for i1 = 1:I
        for k1 = 1:tt
            for j1 = 1:J
                bd_site_all((i1-1)*tt+k1,j1) = pollution_events_site_all(i1,k1,j1);% I����*Kʱ�䣬J����
            end
        end
    end
    [U, latent, V_t, frac_var] = pca_yjy(bd_site_all);
    %��bd_site_all������Ԫ����
    JD(tt,1) = sum(frac_var(1:num_PC));%��һ��Ϊnum_PC����Ԫ�Ĺ�����
    if tt == 1
        JD(tt,2) = 0;
    else
        JD(tt,2) = JD(tt,1)-JD(tt-1,1);%�ڶ���Ϊ��һ����ֵbd_site_c1
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
   