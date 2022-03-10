% 实现DTW算法的函数Dtw.m
% 输入 t为模板矩阵，r为待识别矩阵
% 得到r_new为参考模板矩阵后的匹配矩阵
function [D_flag,r_new] = dtw_yjy_deal(t,r)
n=size(t,1);
m=size(r,1);
r_new = zeros(n,size(r,2));
% 帧匹配距离矩阵
d=zeros(n,m);
for i=1:n
    for j=1:m
        d(i,j) = sum((t(i,1)-r(j,1)).^2);
    end
end
% 累积距离矩阵
D=ones(n,m)*realmax;
D_flag=zeros(n,m);
% 动态规划
ii = 1;%设初值
jj = 1;
%%
if 0.5*m<=n&&n<=2*m % 判断
    for i=1:n
        for j=1:m
            %         if n*j/(2*m)-i<=0&&n*j/(2*m)-i+0.5*n>=0&&2*n*j/m-i>=0&&2*n*j/m-i-n<=0
            if floor((3*n-m)*j/(3*m+n)-i)<=0&&ceil((3*n-m)*j/(3*m+n)-i+(n^2+m^2)/(3*m+n))>=0&&ceil((3*n+m)*j/(3*m-n)-i)>=0&&floor((3*n+m)*j/(3*m-n)-i+(n^2+m^2)/(n-3*m))<=0
                if i==1&&j==1
                    D(i,j)=d(1,1);
                    D1=0;
                    D2=0;
                    D3=0;
                end
                if i==1&&j>1
                    D1=D(i,j-1);
                    D2=realmax;
                    D3=realmax;
                end
                if j==1&&i>1
                    D1=D(i-1,j);
                    D2=realmax;
                    D3=realmax;
                end
                if i>1&&j>1
                    D1=D(i-1,j);
                    D2=D(i,j-1);
                    D3=D(i-1,j-1);
                end
                D(i,j)=d(i,j)+min([D1,D2,D3]);
            end
        end
    end
end
%%
if n<0.5*m
    for i=1:n
        for j=1:m
            %         if n*j/(2*m)-i<=0&&n*j/(2*m)-i+0.5*n>=0&&2*n*j/m-i>=0&&2*n*j/m-i-n<=0
            if floor(n^2*j/m^2-i)<=0&&ceil(n^2*j/m^2-i+n-n^2/m)>=0&&j-i>=0&&j-i+n-m<=0
                if i==1&&j==1
                    D(i,j)=d(1,1);
                    D1=0;
                    D2=0;
                    D3=0;
                end
                if i==1&&j>1
                    D1=D(i,j-1);
                    D2=realmax;
                    D3=realmax;
                end
                if j==1&&i>1
                    D1=D(i-1,j);
                    D2=realmax;
                    D3=realmax;
                end
                if i>1&&j>1
                    D1=D(i-1,j);
                    D2=D(i,j-1);
                    D3=D(i-1,j-1);
                end
                D(i,j)=d(i,j)+min([D1,D2,D3]);
            end
        end
    end
end
%%
if n>2*m
    for i=1:n
        for j=1:m
            %         if n*j/(2*m)-i<=0&&n*j/(2*m)-i+0.5*n>=0&&2*n*j/m-i>=0&&2*n*j/m-i-n<=0
            if ceil(n^2*j/m^2-i)>=0&&floor(n^2*j/m^2-i+n-n^2/m)<=0&&j-i<=0&&j-i+n-m>=0
                if i==1&&j==1
                    D(i,j)=d(1,1);
                    D1=0;
                    D2=0;
                    D3=0;
                end
                if i==1&&j>1
                    D1=D(i,j-1);
                    D2=realmax;
                    D3=realmax;
                end
                if j==1&&i>1
                    D1=D(i-1,j);
                    D2=realmax;
                    D3=realmax;
                end
                if i>1&&j>1
                    D1=D(i-1,j);
                    D2=D(i,j-1);
                    D3=D(i-1,j-1);
                end
                D(i,j)=d(i,j)+min([D1,D2,D3]);
            end
        end
    end
end
%%
% dist=D(n,m);
% 计算轨迹
while ii<n %边缘
    if jj<m %边缘
        D_flag(ii,jj)=1;
        D1=D(ii+1,jj);
        D2=D(ii,jj+1);
        D3=D(ii+1,jj+1);
        if D1 == min([D1,D2,D3])
            ii = ii+1;
            D_flag(ii,jj)=1;
        end
        if D2 == min([D1,D2,D3])
            jj = jj+1;
            D_flag(ii,jj)=1;
        end
        if D3 == min([D1,D2,D3])
            ii = ii+1;
            jj = jj+1;
            D_flag(ii,jj)=1;
        end
    else
        ii = ii+1;
        D_flag(ii,jj)=1;
    end

end
%%
% r_new(i,:) = r(ii,:);

% r_new 的赋值
%     if i >= 2
%         if D(i,lie)>D(i-1,lie+1)
%         [~,lie] = find(D==min(D(i,lie:lie+1)));
% %     if lie < i
% %         lie = lie + 1;
% %     end
%     r_new(i,:) = r(lie,:);

% r_new 的赋值
for i1 = 1:n
    for j1 = 1:m
        if D_flag(i1,j1) == 1 %判定条件
            r_new(i1,:) = r(j1,:);% r_new 的赋值
        end
    end
end
end

% 程序中，首先申请两个n×m的距阵D和d，分别为累积距离和帧匹配距离。
% 这里n和m为测试模板与参考模板的帧数。然后通过一个循环计算两个模板的帧匹配距离距阵d。
% 接下来进行动态规划，为每个格点（i，j）都计算其三个可能的前续格点的累积距离D1、D2和D3。
% 考虑到边界问题，有些前续格点可能不存在，因此要加入一些判断条件。
% 最后利用最小值函数min，找到三个前续格点的累积距离的最小值作为累积距离，与当前帧的匹配距离d（i，j）相加，作为当前格点的累积距离。
% 该计算过程一直达到格点（n，m），并将D（n，m）输出，作为模板匹配的结果。
%