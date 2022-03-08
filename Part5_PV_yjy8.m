% 主导变量识别
% {'co';'no2';'so2';'T@2m';'RH';'NSS@10m';'WES@10m';'HV'};
%% 风向风速
close all; clear ; clc
load('result_step2_3_4.mat', 'club_mean')
PV_para = cell(8,5);
% PV_para(:,1) = {'co';'no2';'so2';'T@2m';'Psta';'RH';'NSS@10m';'WES@10m';'HV'};
PV_para(:,1) = {'co';'no2';'so2';'T@2m';'RH';'NSS@10m';'WES@10m';'HV'};
%%
tar_mat = club_mean{1,1};
NSS_WES = zeros(size(tar_mat,1),2);
for i = 1:size(tar_mat,1)
    NSS_WES(i,1) = tar_mat(i,10)*sin(tar_mat(i,9)*pi/180);
    NSS_WES(i,2) = tar_mat(i,10)*cos(tar_mat(i,9)*pi/180);
end
club_mean_new = [tar_mat(:,1:8),NSS_WES,tar_mat(:,11)];
% 北京第一类污染事件分阶段
%%
stage1 = club_mean_new(3:16,[2,3,5,6,8:11]);
stage2 = club_mean_new(17:30,[2,3,5,6,8:11]);
stage3 = club_mean_new(31:51,[2,3,5,6,8:11]);
stage4 = club_mean_new(52:59,[2,3,5,6,8:11]);
clear club_mean_new tar_mat
% 找到四个阶段
%% 找第二阶段主导变量
for i = 1:4
    PC = 5; % 分析PC数
    para = 8;
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    eval(['input_matrix = stage',num2str(i),';']);
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %  X = U,S,V
    %标准化:
    m= size(input_matrix,1);
    mean_x =mean(input_matrix);
    std_x=std(input_matrix);
    center_x= input_matrix - mean_x(ones(m,1),:);
    standard_x=center_x/diag(std_x);
    %PCA建模:
    % R = standard_x'*standard_x/(m-1);
    [U,S,~] = svds(input_matrix,para);
    % where svd is used to decompose X into three matrices U, latent and V.
    % latent is each component's contribute rate...
    % T = standard_x*U;
    % latent1 = diag(latent);
    total_variance = sum(diag(S));
    frac_var=diag(S)/total_variance;

    % 变量系数矩阵 Coe_all,系数矩阵矩阵占比Coe_all_frac
    Coe_all = [];
    B = zeros(para,PC);
    for j = 1:PC
        A = (standard_x)\U(:,j);% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Coe_all = [Coe_all,A(:,1)];
    end

    sum_sum = sum(abs(Coe_all),1);
    for i1 = 1:PC
        for i2 = 1:para
            B(i2,i1) = abs(Coe_all(i2,i1))/sum_sum(i1)*frac_var(i1);
        end
    end
    % 求每个变量的贡献率
    C = sum(B,2);
    for i3 = 1:para
        PV_para{i3,i+1} = (C(i3)/sum(abs(C)))*(sum(frac_var(1:5)));
    end
end
clear A B C i* U S j *_x
disp(PV_para)
%%