% pollution_events �ض���Ⱦ�¼���PM2.5>=115,t>=24h(����<=5h). 7 step.
% Step 1: Import data
%beijing����
%close all;clear ;clc
cd('H:\works in SWU\Haze Early Warning Research\PM2.5 Early Warning System\MPCA_step\��Ҫ����-�ھ���\beijing_pollution_event_site')
load('NCP_pollu_data_2013_2017.mat', 'bj_site_name')
load('NCP_meteo_data.mat', 'var_name')
%%
for site = 1:size(bj_site_name,1)%վ��
    tt = 24;%��ǰ���Ӻ��¼ʱ��(h)

    load('NCP_pollu_data_2013_2017.mat', 'bj_pm25_2013_2017')
    load('NCP_pollu_data_2013_2017.mat', 'bj_so2_2013_2017')
    load('NCP_pollu_data_2013_2017.mat', 'bj_pm10_2013_2017')
    load('NCP_pollu_data_2013_2017.mat', 'bj_o3_2013_2017')
    load('NCP_pollu_data_2013_2017.mat', 'bj_no2_2013_2017')
    load('NCP_pollu_data_2013_2017.mat', 'bj_co_2013_2017')
    load('NCP_meteo_data.mat', 'NCP_beijing_meteo_2013_2017')

    % load('NCP_pollu_data_2013_2017_fin.mat', 'bj_pm25_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'bj_so2_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'bj_pm10_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'bj_o3_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'bj_no2_2013_2017')
    % load('NCP_pollu_data_2013_2017_fin.mat', 'bj_co_2013_2017')
    % load('NCP_meteo_data.mat', 'NCP_beijing_meteo_2013_2017')
    bj_pm25_2013_2017 = Fill5(bj_pm25_2013_2017,-99);% fill the missing data
    bj_so2_2013_2017 = Fill5(bj_so2_2013_2017,-99);
    bj_pm10_2013_2017 = Fill5(bj_pm10_2013_2017,-99);
    bj_o3_2013_2017 = Fill5(bj_o3_2013_2017,-99);
    bj_no2_2013_2017 = Fill5(bj_no2_2013_2017,-99);
    bj_co_2013_2017 = Fill5(bj_co_2013_2017,-99);

    % fprintf('Step 1 is OK...>>>>>>\n')
    %����
    % Step 2: pollution_data of one site
    bj_pm25 = ones(size(bj_pm25_2013_2017,1),3);
    bj_pm25(:,1) = (1:size(bj_pm25_2013_2017,1))';
    bj_pm25(:,2) = bj_pm25_2013_2017(:,1);
    % for i1 = 1:size(bj_pm25_2013_2017,1)
    %     bj_pm25(i1,3) = mean(bj_pm25_2013_2017(i1,2:end));
    % end
    bj_pm25(:,3) = bj_pm25_2013_2017(:,site+1);
    % fprintf('Step 2 is OK...>>>>>>\n')
    %���� bj_pm25
    % Step 3:pollution_data PM25>=115
    pollution_data = zeros(size(bj_pm25,1),3);
    for i1 = 1:size(bj_pm25,1)
        if (bj_pm25(i1,3) >= 115)
           pollution_data(i1,:) = bj_pm25(i1,:);
        end
    end
    [n1,~] = find(pollution_data == 0);
     pollution_data(n1,:) = [];
    %  fprintf('Step 3 is OK...>>>>>>\n')
    %��Ⱦ���� pollution_data

    % Step 4:create a pollution_data cell 
     %dtime ��ʱ���
    dtime = zeros(size(pollution_data,1),1);
    for i2 = 1:size(pollution_data,1)-2
        dtime(i2) = pollution_data(i2+1,1)-pollution_data(i2,1);
    end
    %n2,ʱ���dtime>=5,,,,,n3 ʱ�������
    [n2,~] = find(dtime >= 6);

    %n3 = zeros(n2,1);
    %for i4 = 1:size(n2,1)-1
    %    n3(i4+1,1) = n2(i4+1,1) - n2(i4,1);
    %end

    %pollution_events ��Ⱦ�¼�
    pollution_events = cell(size(n2,1),1);

    if size(pollution_data(1:n2(1,1),:),1) >= 24
        pollution_events{1} = bj_pm25(pollution_data(1,1)-tt:pollution_data(n2(1,1),1)+tt,2:3);
    else
        pollution_events{1}=[];
    end
    for i3 = 1:size(n2,1)-1
        if size(pollution_data(n2(i3,1):n2(i3+1,1),:),1) >= 24

            pollution_events{i3+1} = bj_pm25(pollution_data(n2(i3,1)+1,1)-tt:pollution_data(n2(i3+1,1),1)+tt,2:3);    
        else
            pollution_events{i3+1} = [];
        end
    end
    %pollution_events ���
    f = cellfun(@isempty,pollution_events);
    pollution_events = pollution_events(~f);
    %  fprintf('Step 4 is OK...>>>>>>\n')
    %��Ⱦ�¼� pollution_events

    % Step 5:pretreatment NCP_city_meteo_2013_2017 
    A = -99*ones(3*size(NCP_beijing_meteo_2013_2017,1),size(NCP_beijing_meteo_2013_2017,2));
    A(:,1) = bj_pm25_2013_2017(:,1);
    for i5 = 1:size(NCP_beijing_meteo_2013_2017,1)
        [A0,~] = find(A(:,1) == NCP_beijing_meteo_2013_2017(i5,1));
        A(A0,2:10) = NCP_beijing_meteo_2013_2017(i5,2:10);
    end
    NCP_beijing_meteo_2013_2017 = A;
    NCP_beijing_meteo_2013_2017 = Fill5(NCP_beijing_meteo_2013_2017,-99);
    %  fprintf('Step 5 is OK...>>>>>>\n')
    % 
    % Step 6:pollution_events
    for i4 = 1:size(pollution_events,1)
        A = pollution_events{i4,1};
        for i6 = 1:size(A,1)
            [a0,~] = find(bj_pm25_2013_2017(:,1) == A(i6,1));
            A(i6,3) = bj_co_2013_2017(a0,site+1);
            A(i6,4) = bj_no2_2013_2017(a0,site+1);
            A(i6,5) = bj_o3_2013_2017(a0,site+1);
            A(i6,6) = bj_pm10_2013_2017(a0,site+1);
            A(i6,7) = bj_so2_2013_2017(a0,site+1);
            [b0,~] = find(NCP_beijing_meteo_2013_2017(:,1) == A(i6,1));
            A(i6,8:9) = NCP_beijing_meteo_2013_2017(b0,2:3);
            A(i6,10:14) = NCP_beijing_meteo_2013_2017(b0,6:10);
        end
        pollution_events{i4,1} = A;
    end
    clear A A0 a0 b0 f i1 i2 i3 i4 i5 i6 n1 n2;
    %  fprintf('Step 6 is OK...>>>>>>\n')
    % ��� pollution_events

    % Step 7:������PM2.5,��SO2,��NO2���ݵ��¼�

    for i7 = 1:size(pollution_events,1)
        for i8 = 1:size(pollution_events{i7,1})
            if pollution_events{i7,1}(i8,2)<0||pollution_events{i7,1}(i8,4)<0||pollution_events{i7,1}(i8,7)<0 %������PM2.5,��SO2,��NO2���ݵ��¼�
                pollution_events{i7,1}=[];
                break                
            end
        end    
    end
    emp = cellfun(@isempty,pollution_events);
    pollution_events = pollution_events(~emp);

    % for i9 = 1:size(pollution_events,1)
    %     for i10 = 1:tt
    %         if pollution_events{i9,1}(i10,2)>150||pollution_events{i9,1}(end-i10+1,2)>150 %����ʼ�ͽ�����ttʱ����PM2.5>115���¼�
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
    is_time = ismember(NCP_beijing_meteo_2013_2017(:,1),non_ep_time);
    non_ep_beijing_meteo_2013_2017 = NCP_beijing_meteo_2013_2017(is_time == 0,:);
    non_ep_beijing_pol_2013_2017 = [bj_pm25_2013_2017(is_time == 0,[1,site+1]), bj_co_2013_2017(is_time == 0,site+1), bj_no2_2013_2017(is_time == 0,site+1), bj_o3_2013_2017(is_time == 0,site+1),...
        bj_pm10_2013_2017(is_time == 0,site+1), bj_so2_2013_2017(is_time == 0,site+1)];
    eval(strcat('non_ep_beijing_meteo_2013_2017_',num2str(site),' = non_ep_beijing_meteo_2013_2017;'));
    eval(strcat('non_ep_beijing_pol_2013_2017_',num2str(site),' = non_ep_beijing_pol_2013_2017;'));
    clear i7 i8 i9 i10 emp 
    % fprintf('Step 7 is OK...>>>>>>\n')
    fprintf('site %2.0f  is OK...>>>>>>\n',site);
end
%% 
pollution_events_all = cell(0,0);
for site = 1:size(bj_site_name,1)
    eval(strcat('pollution_events_all = [pollution_events_all;pollution_events_',num2str(site),'];'));
end
% ��� pollution_events_all
    
    
    
    
