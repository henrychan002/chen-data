clc;
clear;
k=20;
h=1;

load('D:\matlabwork\Fault-tolerant filtering\1009fdata.mat');
load('D:\matlabwork\Fault-tolerant filtering\1009normaldata.mat');


data2=data1;


for i=1:length(data)
    y(i)=data(i);
end
N=length(data);
m=zeros(1,N);
u=zeros(1,N);
low=zeros(1,N);
me=zeros(1,N);
ue=zeros(1,N);
lowe=zeros(1,N);

left=y(1);
right=y(end);
for p=1:k
    y=[left,y];
end

for l=1:k
    y=[y,right];
end

i=k;
for  ti=1+k:N+k
 for j=i-k:i+k
     tj=ti+(j-i)*h;
     yh(j+1)=y(tj);
 end
     yh=sort(yh);
     m(ti-k)=prctile(yh,50);
     u(ti-k)=prctile(yh,75);
     low(ti-k)=prctile(yh,25);
end


for i=1:N
    sum1=0;
    for j=-k:k
      sum1=sum1+f(y(i+k+j),u(i),m(i),low(i));
    end
      yh2(i)=sum1/(2*k+1);
end

e=data-yh2;
for p=1:k
    e=[0,e];
end

for l=1:k
    e=[e,0];
end

figure(22)
plot(e);

i=k;
for  ti=1+k:N+k
 for j=i-k:i+k
     tj=ti+(j-i)*h;
     eh(j+1)=e(tj);
 end
     eh=sort(eh);
     me(ti-k)=prctile(eh,50);
     ue(ti-k)=prctile(eh,75);
     lowe(ti-k)=prctile(eh,25);
end


for i=1:N
    sum2=0;
    for j=-k:k
      sum2=sum2+f(e(i+k+j),ue(i),me(i),lowe(i));
    end
     eh2(i)=sum2/(2*k+1);
end

figure(33)
plot(eh2);

figure(44)
plot(yh2,'linewidth',1.5);
axis([0,2000,8,11]) 

yftf=yh2+eh2;

figure(20)
plot(eh2);
hold on
plot (data)
axis([0,200,8.7,10]) 

e2=yftf-data;
sum3=0;
sum4=0;
for i=1:N
    sum3=yftf(i)^2+sum3;
    sum4=e2(i)^2+sum4;
end    
s1=10*log10(sum3/sum4);

for i=1:N-1
    d(i)=abs(yftf(i+1)-yftf(i))/h; 
end
d1=sort(d,'descend');
d1s=d1(1:5);
d1m=mean(d1s);
dm=mean(d);


my=mean(yftf);
summad=0;
for i=1:N
    summad=summad+abs(yftf(i)-my);
end    
mad=summad/N;

dss=0;
for i=1:N-1
    ds(i)=(yftf(i+1)-yftf(i))^2; 
    dss=ds(i)+dss;
end
dxs=0;
for i=1:N-1
    dx(i)=(data(i+1)-data(i))^2; 
    dxs=dx(i)+dxs;
end
S3=dss/dxs;

for i=1:N-1
    ds(i)=abs(yftf(i+1)-yftf(i))/h; 
    % ds(i)=(yftf(i+1)-yftf(i))/h; 
end
smax=max(ds);
dsmean=mean(ds);
dsstd=std(ds);
figure(55)
plot(ds);

for i=1:N
    xiang(i)=(yftf(i)-data2(i)).^2;
end
sumx=sum(xiang);
Sx=sqrt(sumx/N);

figure(1);
% plot(data,'linewidth',1.5,'linestyle','--');
% hold on;
plot(yftf,'linewidth',1.5);
% xlabel(['\fontname{Times new roman}Sampling point/60s' ...
%     '(d) Sliding window constrained fault tolerant filter']);
xlabel(['\fontname{Times new roman}Sampling point/60s' ...
    '']);
ylabel('\fontname{Times new roman}Vibration value/\fontname{Times new roman}μm');
% ylabel('Pressure Value /Pa');
% axis([0,200,8.7,10]) 
% axis([0,1000,7,11]) 
axis([0,2000,8,11]) 

% legend('\fontname{宋体}输入数据','\fontname{宋体}滤波数据')
set(gca,'Fontsize',15);
% title('Sliding Window Constrained Fault-Tolerant Filtering');

figure(2);
plot(e2);
axis([0,200,-1,1]);
normplot(e2);
% xlabel(['\fontname{Times new roman}Value' ...
%     '(d) Normal probability plot of filter residuals for sliding window constrained fault tolerant filtering']);
xlabel(['\fontname{Times new roman}Value' ...
    '']);
ylabel('\fontname{Times new roman}Probability value');
set(gca,'Fontsize',15);
% title('Normal Distribution Diagram of Sliding Window Constrained Fault-Tolerant Filtering');

error=yftf-data;
error1=yftf-data2';
RMSE = sqrt(mean((error1).^2));
hold on
% plot([1,200],[0,0])
figure(6);
plot(error,'linewidth',1);
axis([0,200,-1,1]) 



figure(7);
histogram(error,30);

load('1009normaldatagen.mat');
er=abs(yftf-gen);
figure(8);
plot(gen);
axis([0,200,8.7,10]) 
sumgen=mean(er);

erstd=std(er);
yuanm=mean(error);
houm=mean(er);


figure(10);
plot(data,'linewidth',1.5);
% xlabel(['\fontname{Times new roman}Sampling point/60s' ...
%     '(b) Sampling data after adding wild values']);
xlabel(['\fontname{Times new roman}Sampling point/60s' ...
    '']);
ylabel('\fontname{Times new roman}Vibration value/\fontname{Times new roman}μm');
% ylabel('Pressure Value /Pa');
% axis([0,200,8.7,10]) 
% axis([0,1000,7,11]) 
axis([0,2000,8,11]) 
% legend('\fontname{宋体}输入数据','\fontname{宋体}滤波数据')
set(gca,'Fontsize',15);

figure(11);
plot(data2,'linewidth',1.5);
% xlabel(['\fontname{Times new roman}Sampling point/60s' ...
%     '(a) Raw sampling data']);
xlabel(['\fontname{Times new roman}Sampling point/60s' ...
    '']);
ylabel('\fontname{Times new roman}Vibration value/\fontname{Times new roman}μm');
% ylabel('Pressure Value /Pa');
% axis([0,200,8.7,10]) 
% axis([0,1000,7,11]) 
axis([0,2000,8,11]) 
% legend('\fontname{宋体}输入数据','\fontname{宋体}滤波数据')
set(gca,'Fontsize',15);


function y=f(z,u,m,low)
 if z>2*u-m
  y=2*u-m;
 elseif (z+m>=2*low)&&(z+m<=2*u)
  y=z;
 elseif z<2*low-m
  y=2*low-m;
 end
end
