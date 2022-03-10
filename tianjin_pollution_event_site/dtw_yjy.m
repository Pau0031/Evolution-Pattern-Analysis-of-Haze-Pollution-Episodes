% DTW的一般算法
% 实现DTW算法的函数Dtw.m
function dist = dtw_yjy(t,r)
%输入 t,r 
n=size(t,1);
m=size(r,1);
%%帧匹配距离矩阵
d=zeros(n,m);
for i=1:n
    for j=1:m
        d(i,j) = sum((t(i,:)-r(j,:)).^2);
    end
end
%%累积距离矩阵
D=ones(n,m)*realmax;
%%动态规划
%全局路径限制
%分三类情况讨论：n=m,n<m,n>m
%横坐标m,纵坐标n 
if 0.5*m<=n<=2*m % 判断
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

dist=D(n,m);
end
% 程序中，首先申请两个n×m的距阵D和d，分别为累积距离和帧匹配距离。
% 这里n和m为测试模板与参考模板的帧数。然后通过一个循环计算两个模板的帧匹配距离距阵d。
% 接下来进行动态规划，为每个格点（i，j）都计算其三个可能的前续格点的累积距离D1、D2和D3。
% 考虑到边界问题，有些前续格点可能不存在，因此要加入一些判断条件。
% 最后利用最小值函数min，找到三个前续格点的累积距离的最小值作为累积距离，与当前帧的匹配距离d（i，j）相加，作为当前格点的累积距离。
% 该计算过程一直达到格点（n，m），并将D（n，m）输出，作为模板匹配的结果。