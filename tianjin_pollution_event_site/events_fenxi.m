% mean time
load('pollution_events_all_dist.mat', 'pollution_events_all')
sum_time = 0;
num_events = size(pollution_events_all,1);
for i1 = 1:num_events
    sum_time = sum_time + size(pollution_events_all{i1,1},1);
end
mean_time = sum_time/num_events;
%% Mean PM2.5 concentration
clear;clc
load('pollution_events_all_dist.mat', 'pollution_events_all')
tt = 24;
sum_con_simevent = 0;
num_events = size(pollution_events_all,1);
mean_con_event = zeros(num_events,1);

for i1 = 1:num_events
    num_time = size(pollution_events_all{i1,1},1) - 2*tt;
    sum_con_events = 0;
    for i2 = 1:num_time
        sum_con_events = sum_con_events + pollution_events_all{i1,1}(i2+24,2);
    end
    mean_con_event(i1,1) = sum_con_events/num_time;
end
Mean_PM25_concentration = mean(mean_con_event);