function [dx] = guozijie(t,x)
dx=zeros(23,1);
%Paramaeter
c1=10;c2=50;
gamma(1)=1;gamma(2)=1;
sigma(1)=50;sigma(2)=50;
zeta_min=0.6;zeta_max=0.8;
M=1;%连杆总质量kg
g=9.8;%重力加速度m/s^2
l=1;%机械臂连杆的质心距连杆的转动中心的距离m
D=2;%连杆转动的粘性摩擦系数N*m*s/rad
J=1;%连杆转动惯量kg*m^2
Mr=3;Ml=1;%死区斜率
m0=min([Mr Ml]);
ar=1.5;al=3;%死区输入信号断点
tou2=0.01;%一阶命令滤波器时间常数
ux=zeros(2,1);
x2d_x=ux(1,1);

x1=x(1);
x2=x(2);
x3=x(3);w1(1)=x3;
x4=x(4);w1(2)=x4;
x5=x(5);w1(3)=x5;
x6=x(6);w1(4)=x6;
x7=x(7);w1(5)=x7;
x8=x(8);w2(1)=x8;
x9=x(9);w2(2)=x9;
x10=x(10);w2(3)=x10;
x11=x(11);w2(4)=x11;
x12=x(12);w2(5)=x12;
Gamma1=x(13);
Gamma2=x(14);
niu1=x(15);
niu2=x(16);
wc1=x(17);
wc2=x(18);
wc3=x(19);
wc4=x(20);
wc5=x(21);
Zs1=x(22);
Zs2=x(23);

%Destination
x1d=sin(t);
dot_x1d=cos(t);

%% 前馈控制器设计
z1=x1-niu1;
%Gamma1误差补偿信号，ba_z1补偿跟踪误差
ba_z1=z1-Gamma1;
error_z1=x1-x1d;
xi=z1+1/2*log10(zeta_min/zeta_max);
S=(zeta_max*exp(xi)-zeta_min*exp(-xi))/(exp(xi)+exp(-xi));

%指定性能函数
neu=2.5*exp(-0.5*t)+0.05;
dot_neu=2.5*exp(-0.5*t)*(-0.5);
p=1/(2*neu)*(1/(S+zeta_min)-1/(S-zeta_max));

%模糊隶属度函数(x1d)
s11(1)=exp(-((x1d-4)^2)/4);
s11(2)=exp(-((x1d-2)^2)/4);
s11(3)=exp(-((x1d-0)^2)/4);
s11(4)=exp(-((x1d+2)^2)/4);
s11(5)=exp(-((x1d+4)^2)/4);
NN1=w1(1)*s11(1)+w1(2)*s11(2)+w1(3)*s11(3)+w1(4)*s11(4)+w1(5)*s11(5);

%控制器x2d_x
x2d_a=-(c1*z1)/p-(p*ba_z1)/2-NN1+dot_x1d+(dot_neu*error_z1)/neu;

z2=x2-niu2;
ba_z2=z2-Gamma2;

x2d=x2d_a+x2d_x;
%模糊隶属度函数(ba_x2d=[x1d,x2d]T)
s21(1)=exp(-((x1d-4)^2)/4);
s21(2)=exp(-((x1d-2)^2)/4);
s21(3)=exp(-((x1d-0)^2)/4);
s21(4)=exp(-((x1d+2)^2)/4);
s21(5)=exp(-((x1d+4)^2)/4);

s22(1)=exp(-((x2d-4)^2)/4);
s22(2)=exp(-((x2d-2)^2)/4);
s22(3)=exp(-((x2d-0)^2)/4);
s22(4)=exp(-((x2d+2)^2)/4);
s22(5)=exp(-((x2d+4)^2)/4);

for i=1:5
s2(i)=s21(i)*s22(i);
end
NN2=w2(1)*s2(1)+w2(2)*s2(2)+w2(3)*s2(3)+w2(4)*s2(4)+w2(5)*s2(5);

%前馈控制器va
va=1/m0*(-c2*z2-ba_z1-ba_z2/2-NN2+dx(16));

%% 最优反馈控制器设计
beta1=0.01;beta2=0.01;
R=[0.2,0;0,0.01];
P=[p,0;0,1];
G=[1,0;0,m0];
QZ=ba_z1^2+ba_z2^2;
E=P*G*R^(-1)*G'*P';
 
%hat_ba_hZ
%模糊隶属度函数(x1)
s3(1)=exp(-((x1-4)^2)/4);
s3(2)=exp(-((x1-2)^2)/4);
s3(3)=exp(-((x1-0)^2)/4);
s3(4)=exp(-((x1+2)^2)/4);
s3(5)=exp(-((x1+4)^2)/4);
NN4=w1(1)*s3(1)+w1(2)*s3(2)+w1(3)*s3(3)+w1(4)*s3(4)+w1(5)*s3(5);

%模糊隶属度函数([x1,x2]T)
s41(1)=exp(-((x1-4)^2)/4);
s41(2)=exp(-((x1-2)^2)/4);
s41(3)=exp(-((x1-0)^2)/4);
s41(4)=exp(-((x1+2)^2)/4);
s41(5)=exp(-((x1+4)^2)/4);

s42(1)=exp(-((x2-4)^2)/4);
s42(2)=exp(-((x2-2)^2)/4);
s42(3)=exp(-((x2-0)^2)/4);
s42(4)=exp(-((x2+2)^2)/4);
s42(5)=exp(-((x2+4)^2)/4);

for i=1:5
s4(i)=s41(i)*s42(i);
end
NN5=w2(1)*s4(1)+w2(2)*s4(2)+w2(3)*s4(3)+w2(4)*s4(4)+w2(5)*s4(5);

hat_hZw1=NN4-NN1;
hat_hZw2=NN5-NN2;
hat_ba_hZ=[hat_hZw1;hat_hZw2];

%模糊隶属度函数(dot_Z=[ba_z1,ba_z2]T)
sc1(1)=exp(-((ba_z1-2)^2)/3)*((-ba_z1-2)*2/3);
sc1(2)=exp(-((ba_z1-1)^2)/3)*((-ba_z1-1)*2/3);
sc1(3)=exp(-((ba_z1-0)^2)/3)*((-ba_z1+0)*2/3);
sc1(4)=exp(-((ba_z1+1)^2)/3)*((-ba_z1+1)*2/3);
sc1(5)=exp(-((ba_z1+2)^2)/3)*((-ba_z1+2)*2/3);

sc2(1)=exp(-((ba_z2-2)^2)/3)*((-ba_z2-2)*2/3);
sc2(2)=exp(-((ba_z2-1)^2)/3)*((-ba_z2-1)*2/3);
sc2(3)=exp(-((ba_z2-0)^2)/3)*((-ba_z2-0)*2/3);
sc2(4)=exp(-((ba_z2+1)^2)/3)*((-ba_z2+1)*2/3);
sc2(5)=exp(-((ba_z2+2)^2)/3)*((-ba_z2+2)*2/3);
for i=1:5
sc(:,i)=[sc1(i);sc2(i)];
end
%dotz_Z*hat_wc
NN3=sc(:,1).*wc1+sc(:,2).*wc2+sc(:,3).*wc3+sc(:,4).*wc4+sc(:,5).*wc5;

%最优控制器([x2d_x,vx]T)
ux=-1/2*R^(-1)*G'*P'*NN3;

hat_gamma=sc'*P*(hat_ba_hZ+G*ux);

for i=1:5
Pai(i)=sc(:,i)'*E*sc(:,i);
end

dot_Zs=P*(hat_ba_hZ+G*ux);
Zs=[Zs1;Zs2];
%JZs=(Zs'*Zs)^(5/2)/5;
dot_JZs=Zs.^4.*dot_Zs;
if(dot_JZs'*dot_Zs<0)
 sum_Z_hat_ux=0;
else
    sum_Z_hat_ux=1;
end
%% 前馈+最优 自适应律更改所需参数
B1=[1;0];
B2=[0;1];

ba_faZ1=[s3(1)-s11(1),0;0,s4(1)-s2(1)];
ba_faZ2=[s3(2)-s11(2),0;0,s4(2)-s2(2)];
ba_faZ3=[s3(3)-s11(3),0;0,s4(3)-s2(3)];
ba_faZ4=[s3(4)-s11(4),0;0,s4(4)-s2(4)];
ba_faZ5=[s3(5)-s11(5),0;0,s4(5)-s2(5)];

%累加 P*B*gamma*ba_z*s
LeiA11=P*B1*gamma(1)*ba_z1*s11(1)+P*B2*gamma(2)*ba_z2*s2(1);
LeiA12=P*B1*gamma(1)*ba_z1*s11(2)+P*B2*gamma(2)*ba_z2*s2(2);
LeiA13=P*B1*gamma(1)*ba_z1*s11(3)+P*B2*gamma(2)*ba_z2*s2(3);
LeiA14=P*B1*gamma(1)*ba_z1*s11(4)+P*B2*gamma(2)*ba_z2*s2(4);
LeiA15=P*B1*gamma(1)*ba_z1*s11(5)+P*B2*gamma(2)*ba_z2*s2(5);

LeiA21=B1*sigma(1)*w1(1)*norm(w1(1))+B2*sigma(2)*w2(1)*norm(w2(1));
LeiA22=B1*sigma(1)*w1(1)*norm(w1(2))+B2*sigma(2)*w2(1)*norm(w2(2));
LeiA23=B1*sigma(1)*w1(1)*norm(w1(3))+B2*sigma(2)*w2(1)*norm(w2(3));
LeiA24=B1*sigma(1)*w1(1)*norm(w1(4))+B2*sigma(2)*w2(1)*norm(w2(4));
LeiA25=B1*sigma(1)*w1(1)*norm(w1(5))+B2*sigma(2)*w2(1)*norm(w2(5));

%% 控制输入
vx=ux(2,1);
u=va+vx;
%%
%状态方程
dx(1)=x2;
dx(2)=-(M*g*l/J)*sin(x1)-(D/J)*x2+u;

%前馈虚拟自适应律(综合）
dot_hat_w1=LeiA11-LeiA21-beta2*sum_Z_hat_ux*P*ba_faZ1*dot_JZs;
dot_hat_w2=LeiA12-LeiA22-beta2*sum_Z_hat_ux*P*ba_faZ2*dot_JZs;
dot_hat_w3=LeiA13-LeiA23-beta2*sum_Z_hat_ux*P*ba_faZ3*dot_JZs;
dot_hat_w4=LeiA14-LeiA24-beta2*sum_Z_hat_ux*P*ba_faZ4*dot_JZs;
dot_hat_w5=LeiA15-LeiA25-beta2*sum_Z_hat_ux*P*ba_faZ5*dot_JZs;

dx(3)=dot_hat_w1(1,1);
dx(4)=dot_hat_w2(1,1);
dx(5)=dot_hat_w3(1,1);
dx(6)=dot_hat_w4(1,1);
dx(7)=dot_hat_w5(1,1);

dx(8)=dot_hat_w1(2,1);
dx(9)=dot_hat_w2(2,1);
dx(10)=dot_hat_w3(2,1);
dx(11)=dot_hat_w4(2,1);
dx(12)=dot_hat_w5(2,1);

%只有前馈
% dx(3)=gamma(1)*p*ba_z1*s11(1)-sigma(1)*x3*norm(x3);
% dx(4)=gamma(1)*p*ba_z1*s11(2)-sigma(1)*x4*norm(x4);
% dx(5)=gamma(1)*p*ba_z1*s11(3)-sigma(1)*x5*norm(x5);
% dx(6)=gamma(1)*p*ba_z1*s11(4)-sigma(1)*x6*norm(x6);
% dx(7)=gamma(1)*p*ba_z1*s11(5)-sigma(1)*x7*norm(x7);
% 
% dx(8)=gamma(2)*ba_z2*s2(1)-sigma(2)*x8*norm(x8);
% dx(9)=gamma(2)*ba_z2*s2(2)-sigma(2)*x9*norm(x9);
% dx(10)=gamma(2)*ba_z2*s2(3)-sigma(2)*x10*norm(x10);
% dx(11)=gamma(2)*ba_z2*s2(4)-sigma(2)*x11*norm(x11);
% dx(12)=gamma(2)*ba_z2*s2(5)-sigma(2)*x12*norm(x12);

%dot_Gamma1
dx(13)=-c1*Gamma1+p*(niu2-x2d+Gamma2);
%dot_Gamma2
dx(14)=-c2*Gamma2;
%dot_niu1
dx(15)=dot_x1d;
%dot_niu2
dx(16)=-(niu2-x2d)/tou2;
%最优反馈控制自适应律
dx(17)=-beta1*hat_gamma(1,1)*(QZ+NN3'*P*hat_ba_hZ-1/4*wc1*Pai(1)*wc1)+1/2*beta2*sum_Z_hat_ux*sc(:,1)'*E*dot_JZs;
dx(18)=-beta1*hat_gamma(2,1)*(QZ+NN3'*P*hat_ba_hZ-1/4*wc2*Pai(2)*wc2)+1/2*beta2*sum_Z_hat_ux*sc(:,2)'*E*dot_JZs;
dx(19)=-beta1*hat_gamma(3,1)*(QZ+NN3'*P*hat_ba_hZ-1/4*wc3*Pai(3)*wc3)+1/2*beta2*sum_Z_hat_ux*sc(:,3)'*E*dot_JZs;
dx(20)=-beta1*hat_gamma(4,1)*(QZ+NN3'*P*hat_ba_hZ-1/4*wc4*Pai(4)*wc4)+1/2*beta2*sum_Z_hat_ux*sc(:,4)'*E*dot_JZs;
dx(21)=-beta1*hat_gamma(5,1)*(QZ+NN3'*P*hat_ba_hZ-1/4*wc5*Pai(5)*wc5)+1/2*beta2*sum_Z_hat_ux*sc(:,5)'*E*dot_JZs;

dx(22)=dot_Zs(1,1);
dx(23)=dot_Zs(2,1);
end

