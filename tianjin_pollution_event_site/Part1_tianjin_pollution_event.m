 % pollution_events 重度污染事件，PM2.5>=115,t>=24h(回落<=5h). 7 step.
% Step 1: Import data
%tianjin北京
close all;clear ;clc
load('NCP_pollu_data_2013_2017.mat', 'tj_site_name')
%%
clc;
for site = 1:size(tj_site_name,1)%站点
    tt = 24;%提前和延后记录时间(h)
%     load('NCP_meteo_data.mat', 'var_name')
%     load('NCP_pollu_data_2013_2017.mat', 'tj_pm25_2013_2017')
%     load('NCP_pollu_data_2013_2017.mat', 'tj_so2_2013_2017')
%     load('NCP_pollu_data_2013_2017.mat', 'tj_pm10_2013_2017')
%     load('NCP_pollu_data_2013_2017.mat', 'tj_o3_2013_2017')
%     load('NCP_pollu_data_2013_2017.mat', 'tj_no2_2013_2017')
%     load('NCP_pollu_data_2013_2017.mat', 'tj_co_2013_2017')
%     load('NCP_meteo_data.mat', 'NCP_tianjin_meteo_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'tj_pm25_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'tj_so2_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'tj_pm10_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'tj_o3_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'tj_no2_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'tj_co_2013_2017')
    load('NCP_meteo_data.mat', 'NCP_tianjin_meteo_2013_2017')
    tj_pm25_2013_2017 = Fill5(tj_pm25_2013_2017,-99);
    tj_so2_2013_2017 = Fill5(tj_so2_2013_2017,-99);
    tj_pm10_2013_2017 = Fill5(tj_pm10_2013_2017,-99);
    tj_o3_2013_2017 = Fill5(tj_o3_2013_2017,-99);
    tj_no2_2013_2017 = Fill5(tj_no2_2013_2017,-99);
    tj_co_2013_2017 = Fill5(tj_co_2013_2017,-99);
    % fprintf('Step 1 is OK...>>>>>>\n')
    %数据
    % Step 2: pollution_data of one site
    tj_pm25 = ones(size(tj_pm25_2013_2017,1),3);
    tj_pm25(:,1) = (1:size(tj_pm25_2013_2017,1))';
    tj_pm25(:,2) = tj_pm25_2013_2017(:,1);
    % for i1 = 1:size(tj_pm25_2013_2017,1)
    %     tj_pm25(i1,3) = mean(tj_pm25_2013_2017(i1,2:end));
    % end
    tj_pm25(:,3) = tj_pm25_2013_2017(:,site+1);
    % fprintf('Step 2 is OK...>>>>>>\n')
    %数据 tj_pm25
    % Step 3:pollution_data PM25>=115
    pollution_data = zeros(size(tj_pm25,1),3);
    for i1 = 1:size(tj_pm25,1)
        if (tj_pm25(i1,3) >= 115)
           pollution_data(i1,:) = tj_pm25(i1,:);
        end
    end
    [n1,~] = find(pollution_data == 0);
     pollution_data(n1,:) = [];
    %  fprintf('Step 3 is OK...>>>>>>\n')
    %污染数据 pollution_data

    % Step 4:create a pollution_data cell 
     %dtime 求时间差
    dtime = zeros(size(pollution_data,1),1);
    for i2 = 1:size(pollution_data,1)-2
        dtime(i2) = pollution_data(i2+1,1)-pollution_data(i2,1);
    end
    %n2,时间差dtime>=5,,,,,n3 时间差条件
    [n2,~] = find(dtime >= 6);

    %n3 = zeros(n2,1);
    %for i4 = 1:size(n2,1)-1
    %    n3(i4+1,1) = n2(i4+1,1) - n2(i4,1);
    %end

    %pollution_events 污染事件
    pollution_events = cell(size(n2,1),1);

    if size(pollution_data(1:n2(1,1),:),1) >= 24
        pollution_events{1} = tj_pm25(pollution_data(1,1)-tt:pollution_data(n2(1,1),1)+tt,2:3);
    else
        pollution_events{1}=[];
    end
    for i3 = 1:size(n2,1)-1
        if size(pollution_data(n2(i3,1):n2(i3+1,1),:),1) >= 24&&pollution_data(n2(i3,1)+1,1)-tt>0&&pollution_data(n2(i3+1,1),1)+tt<=size(tj_pm25,1)
            pollution_events{i3+1} = tj_pm25(pollution_data(n2(i3,1)+1,1)-tt:pollution_data(n2(i3+1,1),1)+tt,2:3);    
        else
            pollution_events{i3+1} = [];
        end
    end
    %pollution_events 清空
    f = cellfun(@isempty,pollution_events);
    pollution_events = pollution_events(~f);
    %  fprintf('Step 4 is OK...>>>>>>\n')
    %污染事件 pollution_events

    % Step 5:pretreatment NCP_city_meteo_2013_2017 
    A = -99*ones(3*size(NCP_tianjin_meteo_2013_2017,1),size(NCP_tianjin_meteo_2013_2017,2));
    A(:,1) = tj_pm25_2013_2017(:,1);
    for i5 = 1:size(NCP_tianjin_meteo_2013_2017,1)
        [A0,~] = find(A(:,1) == NCP_tianjin_meteo_2013_2017(i5,1));
        A(A0,2:10) = NCP_tianjin_meteo_2013_2017(i5,2:10);
    end
    NCP_tianjin_meteo_2013_2017 = A;
    NCP_tianjin_meteo_2013_2017 = Fill5(NCP_tianjin_meteo_2013_2017,-99);
    %  fprintf('Step 5 is OK...>>>>>>\n')
    % 
    % Step 6:pollution_events
    for i4 = 1:size(pollution_events,1)
        A = pollution_events{i4,1};
        for i6 = 1:size(A,1)
            [a0,~] = find(tj_pm25_2013_2017(:,1) == A(i6,1));
            A(i6,3) = tj_co_2013_2017(a0,site+1);
            A(i6,4) = tj_no2_2013_2017(a0,site+1);
            A(i6,5) = tj_o3_2013_2017(a0,site+1);
            A(i6,6) = tj_pm10_2013_2017(a0,site+1);
            A(i6,7) = tj_so2_2013_2017(a0,site+1);
            [b0,~] = find(NCP_tianjin_meteo_2013_2017(:,1) == A(i6,1));
            A(i6,8:9) = NCP_tianjin_meteo_2013_2017(b0,2:3);
            A(i6,10:14) = NCP_tianjin_meteo_2013_2017(b0,6:10);
        end
        pollution_events{i4,1} = A;
    end
    clear A A0 a0 b0 f i1 i2 i3 i4 i5 i6 n1 n2;
    %  fprintf('Step 6 is OK...>>>>>>\n')
    % 输出 pollution_events

    % Step 7:清理含无PM2.5,无SO2,无NO2数据的事件

    for i7 = 1:size(pollution_events,1)
        for i8 = 1:size(pollution_events{i7,1})
            if pollution_events{i7,1}(i8,2)<0||pollution_events{i7,1}(i8,4)<0||pollution_events{i7,1}(i8,7)<0 %清理含无PM2.5,无SO2,无NO2数据的事件
                pollution_events{i7,1}=[];
                break                
            end
        end    
    end
    emp = cellfun(@isempty,pollution_events);
    pollution_events = pollution_events(~emp);

    % for i9 = 1:size(pollution_events,1)
    %     for i10 = 1:tt
    %         if pollution_events{i9,1}(i10,2)>150||pollution_events{i9,1}(end-i10+1,2)>150 %清理开始和结束的tt时间内PM2.5>115的事件
    %             pollution_events{i9,1}=[];
    %             break
    %         end
    %     end
    % end

    emp = cellfun(@isempty,pollution_events);
    pollution_events = pollution_events(~emp);
    eval(strcat('pollution_events_',num2str(site),'= pollution_events(~emp);')); 
    non_ep_time = [];
    for ep_time = 1:size(pollution_events,1)
        non_ep_time = [non_ep_time;pollution_events{ep_time,1}(:,1)];
    end
    non_ep_time = unique(non_ep_time);
    is_time = ismember(NCP_tianjin_meteo_2013_2017(:,1),non_ep_time);
    non_ep_tianjin_meteo_2013_2017 = NCP_tianjin_meteo_2013_2017(is_time == 0,:);
    non_ep_tianjin_pol_2013_2017 = [tj_pm25_2013_2017(is_time == 0,[1,site+1]), tj_co_2013_2017(is_time == 0,site+1), tj_no2_2013_2017(is_time == 0,site+1), tj_o3_2013_2017(is_time == 0,site+1),...
        tj_pm10_2013_2017(is_time == 0,site+1), tj_so2_2013_2017(is_time == 0,site+1)];
    eval(strcat('non_ep_tianjin_meteo_2013_2017_',num2str(site),' = non_ep_tianjin_meteo_2013_2017;'));
    eval(strcat('non_ep_tianjin_pol_2013_2017_',num2str(site),' = non_ep_tianjin_pol_2013_2017;'));
    clear i7 i8 i9 i10 emp A f i1 i2 i3 n1 n2 
    % fprintf('Step 7 is OK...>>>>>>\n')
    fprintf('site %2.0f  is OK...>>>>>>\n',site);
end
%% 
pollution_events_all = cell(0,0);
for site = 1:size(tj_site_name,1)
    eval(strcat('pollution_events_all = [pollution_events_all;pollution_events_',num2str(site),'];'));
end
% 求得 pollution_events_all
    
    
    
    
