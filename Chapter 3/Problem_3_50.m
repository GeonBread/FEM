
clearvars;
clc; close;

%% Given Parameters
E=30*10^(6);
A1=2;
A2=2;
A3=2;
L1=30*sqrt(2);
L2=30;
L3=30;

%% Element Stiffness Equation
% Element 1
Theta1=45;
alpha1=-45;
alpha2=0;
co1_1=cos((Theta1-alpha1)/180*pi);
si1_1=sin((Theta1-alpha1)/180*pi);
co1_2=cos((Theta1-alpha2)/180*pi);
si1_2=sin((Theta1-alpha2)/180*pi);
matT1=[ co1_1 si1_1    0   0;
       -si1_1 co1_1    0   0;
          0   0  co1_2 si1_2;
          0   0 -si1_2 co1_2];
kLocal1=(E*A1)/(L1)*[ 1 0 -1 0;
                      0 0  0 0;
                     -1 0  1 0;
                      0 0  0 0];
kGlobal1=transpose(matT1)*kLocal1*matT1;
% Element 2
Theta2=180;
alpha2=0;
alpha3=0;
co2_2=cos((Theta2-alpha2)/180*pi);
si2_2=sin((Theta2-alpha2)/180*pi);
co2_3=cos((Theta2-alpha3)/180*pi);
si2_3=sin((Theta2-alpha3)/180*pi);
matT2=[ co2_2 si2_2    0   0;
       -si2_2 co2_2    0   0;
          0   0  co2_3 si2_3;
          0   0 -si2_3 co2_3];
kLocal2=(E*A2)/(L2)*[ 1 0 -1 0;
                      0 0  0 0;
                     -1 0  1 0;
                      0 0  0 0];
kGlobal2=transpose(matT2)*kLocal2*matT2;
% Element 3
Theta3=90;
alpha1=-45;
alpha3=0;
co3_1=cos((Theta3-alpha1)/180*pi);
si3_1=sin((Theta3-alpha1)/180*pi);
co3_3=cos((Theta3-alpha3)/180*pi);
si3_3=sin((Theta3-alpha3)/180*pi);
matT3=[ co3_1 si3_1    0   0;
       -si3_1 co3_1    0   0;
          0   0  co3_3 si3_3;
          0   0 -si3_3 co3_3];
kLocal3=(E*A3)/(L3)*[ 1 0 -1 0;
                      0 0  0 0;
                     -1 0  1 0;
                      0 0  0 0];
kGlobal3=transpose(matT3)*kLocal3*matT3;

%% Assemblage
KK_hat=zeros(6,6);
% Element 1
eDoFs1=[2*(1)-1,2*(1),2*(2)-1,2*(2)];
KK_hat(eDoFs1,eDoFs1)=KK_hat(eDoFs1,eDoFs1)+kGlobal1;
% Element 2
eDoFs2=[2*(2)-1,2*(2),2*(3)-1,2*(3)];
KK_hat(eDoFs2,eDoFs2)=KK_hat(eDoFs2,eDoFs2)+kGlobal2;
% Element 3
eDoFs3=[2*(1)-1,2*(1),2*(3)-1,2*(3)];
KK_hat(eDoFs3,eDoFs3)=KK_hat(eDoFs3,eDoFs3)+kGlobal3;
%% Impose Boundary Conditions
FF_hat=zeros(6,1);
FF_hat(2*(2),1)=-2000;
dd_hat=zeros(6,1);
eDoFs = [2*(1)-1, 2*(2)];
dd_hat(eDoFs,1) = KK_hat(eDoFs,eDoFs) \ FF_hat(eDoFs,1);
%% Reaction Force Calculation
FF_hat = KK_hat*dd_hat;
%% Element Stress Calculation
% Element 1
eStress1=(1/A1)*kLocal1(3,:)*matT1*dd_hat(eDoFs1,1)
% Element 2
eStress2=(1/A2)*kLocal2(3,:)*matT2*dd_hat(eDoFs2,1)
% Element 3
eStress3=(1/A3)*kLocal3(3,:)*matT3*dd_hat(eDoFs3,1)
%% Visualization
scale=5*10^(3);
%node 1 좌표계를 글로벌로 바꿔야 됨.
TGlobal=eye(6,6);
co=cos(alpha1/180*pi);
si=sin(alpha1/180*pi);
TGlobal(1:2,1:2)=[co -si;
                  si  co];
dd=TGlobal*dd_hat;
x = 1;
y = 2;
Node1 = [1,2];
Node2 = [3,4];
Node3 = [5,6];
Node1_point = [0, 0];
Node2_point = [30, 30];
Node3_point = [0, 30];
% Element 1
xInitial_1=[Node1_point(1); Node2_point(1)];
yInitial_1=[Node1_point(2); Node2_point(2)];
xDeformed_1 = xInitial_1 + scale * dd([Node1(1), Node2(1)], 1);
yDeformed_1 = yInitial_1 + scale * dd([Node1(2), Node2(2)], 1);
% Element 2
xInitial_2=[Node2_point(1); Node3_point(1)];
yInitial_2=[Node2_point(2); Node3_point(2)];
xDeformed_2 = xInitial_2 + scale * dd([Node2(1), Node3(1)], 1);
yDeformed_2 = yInitial_2 + scale * dd([Node2(2), Node3(2)], 1);
% Element 3
xInitial_3=[Node1_point(1); Node3_point(1)];
yInitial_3=[Node1_point(2); Node3_point(2)];
xDeformed_3 = xInitial_3 + scale * dd([Node1(1), Node3(1)], 1);
yDeformed_3 = yInitial_3 + scale * dd([Node1(2), Node3(2)], 1);
figure(1);
plot(xInitial_1,yInitial_1,'k','LineWidth',1.1);
hold on;
plot(xInitial_2,yInitial_2,'k','LineWidth',1.1);
hold on;
plot(xInitial_3,yInitial_3,'k','LineWidth',1.1);
hold on;
plot(xDeformed_1,yDeformed_1,'b','LineWidth',1.1);
hold on;
plot(xDeformed_2,yDeformed_2,'b','LineWidth',1.1);
hold on;
plot(xDeformed_3,yDeformed_3,'b','LineWidth',1.1);
hold on;
axis padded;
