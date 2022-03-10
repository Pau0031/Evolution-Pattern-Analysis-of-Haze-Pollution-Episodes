% DTW��һ���㷨
% ʵ��DTW�㷨�ĺ���Dtw.m
function dist = dtw_yjy(t,r)
%���� t,r 
n=size(t,1);
m=size(r,1);
%%֡ƥ��������
d=zeros(n,m);
for i=1:n
    for j=1:m
        d(i,j) = sum((t(i,:)-r(j,:)).^2);
    end
end
%%�ۻ��������
D=ones(n,m)*realmax;
%%��̬�滮
%ȫ��·������
%������������ۣ�n=m,n<m,n>m
%������m,������n 
if 0.5*m<=n<=2*m % �ж�
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
% �����У�������������n��m�ľ���D��d���ֱ�Ϊ�ۻ������֡ƥ����롣
% ����n��mΪ����ģ����ο�ģ���֡����Ȼ��ͨ��һ��ѭ����������ģ���֡ƥ��������d��
% ���������ж�̬�滮��Ϊÿ����㣨i��j�����������������ܵ�ǰ�������ۻ�����D1��D2��D3��
% ���ǵ��߽����⣬��Щǰ�������ܲ����ڣ����Ҫ����һЩ�ж�������
% ���������Сֵ����min���ҵ�����ǰ�������ۻ��������Сֵ��Ϊ�ۻ����룬�뵱ǰ֡��ƥ�����d��i��j����ӣ���Ϊ��ǰ�����ۻ����롣
% �ü������һֱ�ﵽ��㣨n��m��������D��n��m���������Ϊģ��ƥ��Ľ����