function [U, latent, V_t, frac_var] = pca_yjy(input_matrix)%, contribution_rate )
% pca_yjy is used to get an X = U,S,Vt
% where U is score matrix, V is loading matrix, if X = [M,N] m is sample
% num and N is the variable num. thus, U have N rows and each column is a component 
% V has M rows and each row in Vt is a principal component.
% I use this program to project data into PC sub-space  
%��׼��:
m= size(input_matrix,1);
mean_x =mean(input_matrix);
std_x=std(input_matrix);
center_x= input_matrix - mean_x(ones(m,1),:);
standard_x=center_x/diag(std_x);

%PCA��ģ:
R = standard_x'*standard_x/(m-1);
[U,latent,V] = svd(R,0);
V_t = V';
% where svd is used to decompose X into three matrices U, latent and V.
% latent is each component's contribute rate..
% T = standard_x*U;
% latent1 = diag(latent);
total_variance = sum(diag(latent));
frac_var=diag(latent)/total_variance;
% sum_f = 0;
% 
% for i = 1:size(frac_var,1)
%     if sum_f <= contribution_rate
%         sum_f = sum_f + frac_var(i);
%     else
%         break
%     end
% end
% pc_num = i -1;
% true_contribution_rate = sum_f;
% clear i sum_f
% 
% latent1 = diag(latent);
% latent2 = latent1(1:k,:);
% l=diag(1./latent2);
% Tk=T(:,1:k);
% Uk=U(:,1:k);
% Unk=U(:,k+1:n);
% tmp=sqrt(diag(1./latent2))*Tk';
% tsquare=sum(tmp.*tmp)';
% a=k*(m-1)/(m-k);
% F1=finv(0.99,k,(m-k));
% F2=finv(0.95,k,(m-k));
% UCL1=a*F1;
% UCL2=a*F2;
% U1=[UCL1 UCL1];
% U2=[UCL2 UCL2];
% tt=[0 m];
% t=0:(m-1);


% I=eye(n);
% c=I-Uk*Uk';
% q=xn*c*xn';
% Q=diag(q);
% j=(k+1):n;
% b=latent1(j,1);
% b1=b.^1;
% b2=b.^2;
% b3=b.^3;
% y1=sum(b1);
% y2=sum(b2);
% y3=sum(b3);
% h=1-2*y1*y3/(3*y2*y2);
% z1=2.325;
% z2=1.645;
% A1=z1*h*sqrt(2*y2)/y1;
% A2=z2*h*sqrt(2*y2)/y1;
% B=y2*h*(h-1)/(y1)^2;
% SPE2=y1*(A2+1+B)^(1/h);
% SPE1=y1*(A1+1+B)^(1/h);
% S1=[SPE1 SPE1];
% S2=[SPE2 SPE2];



%��������
%��׼��:
% [m1,n]=size(x1);
% meanx1=mean(x1);
% stdx1=std(x1);
% centerx1=x1-meanx1(ones(m1,1),:);
% xn1=centerx1/diag(stdx1);
% tt=[0 m1];
% t=[0:1:(m1-1)];
% 
% %PCAģ��:
% T1=xn1*U;
% T1k=T1(:,[1:k]);
% tmp1=sqrt(l)*T1k';
% tsquare1=sum(tmp1.*tmp1)';
% q1=xn1*c*xn1';
% Q1=diag(q1);
% 
% figure(3)
% plot(tt,S1)
% hold on
% plot(t,Q1,'r')
% plot(t,Q1,'r*')
% xlabel('sample number');
% ylabel('SPE');

%�������ݵ�Qֵ���ó��������޵�ԭʼ���ݼ�Qֵ%
% zg=[];%��ʼ����zg��ʾ�������ݵ��ԭʼֵ%
% I=[];%��ʼ����I��ʾ���������޵Ĺ������ݵ����%
% qc=[];%��ʼ����qc��ʾ���������޵Ĺ������ݵ�Qֵ%
% for o=1:m1
%   Q11=Q1(o,:);%Qֵ%
%     x11=x1(o,:);%ԭʼ����%
%     if Q11>SPE1
%         I=[I;o];
%         qc=[qc;Q11];
%         zg=[zg;x11];
%     else
%     end
% end



% figure(4)
% plot(tt,U1)
% hold on
% plot(t,tsquare1,'r')
% plot(t,tsquare1,'r*')
% xlabel('sample number');
% ylabel('T2');
% 
% 
% tsquareg=[];%��ʼ����tsquareg��ʾ������ݵ�T2ֵ%
% J=[];%����ʼ������ʾ���������޵Ĺ������ݵ����%
% tc=[];%��ʼ������ʾ���������޵�T2ֵ%
% xc=[];%��ʼ������ʾ���������޵Ĺ������ݵ�%
% for j=1:m1
%      x12=x1(j,:);
%      tsquare11=tsquare1(j,:);
%     tsquareg=[tsquareg;tsquare11];
%     if tsquare11>UCL1
%       J=[J;j];
%       tc=[tc;tsquare11];
%       xc=[xc;x12];
%     else 
%    end
% end
end

