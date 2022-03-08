%% mean time
clc;
% load('beijing_pollution_event.mat', 'pollution_events_all')
% CO 3; SO2 7; NO2 4 T 8; WS 12' WD 11; RH 10
tar_club = pollution_events_club6;
mean_v_p_m = [];
sta_du = zeros(4,1);
for i = 1:size(tar_club,1)
    mean_v_p_m = [mean_v_p_m;tar_club{i,1}(25:end-24,1)];
    for j = 1:4
        sta_du(j,1) = sta_du(j,1)+(tar_club{i,4}(j+1,2)-tar_club{i,4}(j,2));
    end
end
sta_du = sta_du./size(tar_club,1);
fprintf('mean value of PM2.5 is %6.2f\n',mean(mean_v_p_m));
fprintf('max value of PM2.5 is %6.2f\n',max(mean_v_p_m));
fprintf('mean value of PM2.5 duration is %6.2fh\n',length(mean_v_p_m)/size(tar_club,1));
%%
clc;
non_ep_sjz_pol = [];
non_ep_sjz_mete = [];
for i = 1:8
    eval(strcat('non_ep_sjz_pol = [non_ep_sjz_pol;non_ep_shijiazhuang_pol_2013_2017_',num2str(i),'];'));
    eval(strcat('non_ep_sjz_mete = [non_ep_sjz_mete;non_ep_shijiazhuang_meteo_2013_2017_',num2str(i),'];'));
end
mean_pol = zeros(size(non_ep_sjz_pol,2),1);
std_pol = zeros(size(non_ep_sjz_pol,2),1);
for i = 1:size(non_ep_sjz_pol,2)
    tmp = non_ep_sjz_pol(:,i);
    mean_pol(i,1) = mean(tmp(tmp>0));
    std_pol(i,1) = std(tmp(tmp>0));
end
mean_mete = zeros(size(non_ep_sjz_mete,2),1);
std_mete = zeros(size(non_ep_sjz_mete,2),1);
for i = 1:size(non_ep_sjz_mete,2)
    tmp = non_ep_sjz_mete(:,i);
    mean_mete(i,1) = mean(tmp(tmp>0));
    std_mete(i,1) = std(tmp(tmp>0));
end
% mean_mete = mean(non_ep_sjz_mete,1);
% std_mete = std(non_ep_sjz_mete,1);
clear tmp
%% stage statistics
clc;
club_num = 6;
stage_stat = zeros(5*4*club_num,size(pollution_events_club1{1,1},2)); % mean, beginning, ending, change range, change rate/ four stages/ = 5*4
stage_wind_s_d = cell(club_num,4);
for i = 1:club_num
    eval(strcat(['tar_club = pollution_events_club',num2str(i),';']));
    eval(strcat(['tmp_m_c',num2str(i),' = zeros(size(tar_club,1),5*4*size(tar_club{1,1},2));']));
    for j = 1:4 %stage
        for k = 1:size(tar_club{1,1},2) %1:11 parameter11_name = {'pm25';'co';'no2';'o3';'so2';'T@2m';'Psta';'RH';'WD@10m';'WS@10m';'HV'};
            for l = 1:size(tar_club,1)
                tmp = tar_club{l,1}(tar_club{l,4}(j,2):tar_club{l,4}(j+1,2),k);%l 1
                eval(strcat(['tmp_m_c',num2str(i),'(l,(j-1)*55 + 5*(k-1)+1) = mean(tmp);']));
                eval(strcat(['tmp_m_c',num2str(i),'(l,(j-1)*55 + 5*(k-1)+2) = tmp(1,1);']));
                eval(strcat(['tmp_m_c',num2str(i),'(l,(j-1)*55 + 5*(k-1)+3) = tmp(end,1);']));
                eval(strcat(['tmp_m_c',num2str(i),'(l,(j-1)*55 + 5*(k-1)+4) = tmp(end)-tmp(1);']));
                eval(strcat(['tmp_m_c',num2str(i),'(l,(j-1)*55 + 5*k) = (tmp(end)-tmp(1))/length(tmp);']));
                stage_wind_s_d{i,j} = [stage_wind_s_d{i,j};tar_club{l,1}(tar_club{l,4}(j,2):tar_club{l,4}(j+1,2),:)];
            end
            eval(strcat(['stage_stat((i-1)*20 + (j-1)*5 + 1:(i-1)*20 + j*5,k) = mean(tmp_m_c',num2str(i),'(:,(j-1)*55 + 5*(k-1)+1:(j-1)*55 + 5*k),1)'';']));
            clear tmp
        end
    end
    fprintf('the %dth club has done\n',i)
end
%% wind field picture
k = 0;
for i = 1:6
    for j = 1:4
        subplot(6,4,k+1);
        polarhistogram(stage_wind_s_d{i,j}(:,9)./360*2*pi,12,'FaceColor','red','FaceAlpha',.3);
        eval(strcat(['title(''C_',num2str(i),'-S_',num2str(j),''');']));
        k = k+1;
    end
end