% stage ʶ��׶�

close all; clear;clc;
%%
% ����
% J������Kʱ�䣬I����
% �޸�pollution_events_10para��pollution_events_site_all��bj_site_all
cd('H:\works in SWU\Haze Early Warning Research\PM2.5 Early Warning System\MPCA_step\��Ҫ����-�ھ���\beijing_pollution_event_site')
load('result_step2_3.mat', 'Fevent')
load('result_step2_3.mat', 'pollution_events_club1')
load('result_step2_3.mat', 'pollution_events_club2')
load('result_step2_3.mat', 'pollution_events_club3')
load('result_step2_3.mat', 'pollution_events_club4')
load('result_step2_3.mat', 'pollution_events_club5')
load('result_step2_3.mat', 'club_mean')
%% STAGE ʶ��׶�
% ѡȡ6������parameter6_name = {'no2';'so2';'T@2m';'RH';'WS@10m';'HV'};
% parameter6_name = {'no2';'so2';'T@2m';'RH';'WS@10m';'HV'};
%ѡȡԭ����11��������{'pm25';'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'}
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
tar_club = pollution_events_club5;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
parameter10_name = {'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'};
pollution_events_10para = cell(size(tar_club,1),1);
for i0 = 1:size(pollution_events_10para,1)
    pollution_events_10para{i0,1} = tar_club{i0,2}(:,2:end);
end
% ��� pollution_events_10para

% ѡȡ���ʵ� pollution_events_10para �ֽ׶�
[I,~] = size(pollution_events_10para);
[K,J] = size(pollution_events_10para{1,1});
% [K,~] = size(pollution_events_10para{1,1});
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
bj_site_all = zeros(I*K,J);%unfold variablely 
for i1 = 1:I
    for k1 = 1:K
        for j1 = 1:J
            bj_site_all((i1-1)*K+k1,j1) = pollution_events_site_all(i1,k1,j1);% I����*Kʱ�䣬J����
        end
    end
end
% ��ά����
fprintf('Step 1 is OK...>>>>>>\n')
%%
% �Զ�ʶ��׶�
num_PC = 1 ;%ָ����Ԫ����
JD = ones(K,2)*-99; %�����׶Σ�ʱ�䳤��ΪK
% ��һ��ΪPC���ף��ڶ���Ϊ����ֵ��������Ϊ��ֵ����
% �����ӽ׶�cell
for tt = 1:K  %��1��Kÿһ��������Ԫ
    bj_site_all = zeros(I*tt,J);%����ֵ
    for i1 = 1:I
        for k1 = 1:tt
            for j1 = 1:J
                bj_site_all((i1-1)*tt+k1,j1) = pollution_events_site_all(i1,k1,j1);% I����*Kʱ�䣬J����
            end
        end
    end
    [U, latent, V_t, frac_var] = pca_yjy(bj_site_all);
    %��bj_site_all������Ԫ����
    JD(tt,1) = sum(frac_var(1:num_PC));%��һ��Ϊnum_PC����Ԫ�Ĺ�����
    if tt == 1
        JD(tt,2) = 0;
    else
        JD(tt,2) = JD(tt,1)-JD(tt-1,1);%�ڶ���Ϊ��һ����ֵbj_site_c1
    end
end
fprintf('Step 2 is OK...>>>>>>\n')
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
JD_club5 = JD;
%>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% [~,tt_2] = sort(-JD(:,2));%�������У�����λ��
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
    [~,tt_2] = sort(-JD(:,2));%�������У�����λ��
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