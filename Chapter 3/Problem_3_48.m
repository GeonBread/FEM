clearvars;
clc; close;

%% Given Parameters
E = 210 * 10^9;
A = 10 * 10^-4;
Theta1 = atand(4/8);
Theta2 = 0;
Theta3 = -atand(4/8);
Theta4 = atand(4/8);
Theta5 = -90;
L1 = sqrt(20);
L2 = 8;
L3 = sqrt(20);
L4 = sqrt(20);
L5 = 4;

%% Element Stiffness Equation
% Element 1
co1=cos(Theta1/180*pi);
si1=sin(Theta1/180*pi);
matT1=[ co1 si1    0   0;
       -si1 co1    0   0;
          0   0  co1 si1;
          0   0 -si1 co1];
kLocal1=(E*A)/(L1)*[ 1 0 -1 0;
                     0 0  0 0;
                    -1 0  1 0;
                     0 0  0 0];
kGlobal1=transpose(matT1)*kLocal1*matT1;
% Element 2
co2=cos(Theta2/180*pi);
si2=sin(Theta2/180*pi);
matT2=[ co2 si2    0   0;
       -si2 co2    0   0;
          0   0  co2 si2;
          0   0 -si2 co2];
kLocal2=(E*A)/(L2)*[ 1 0 -1 0;
                     0 0  0 0;
                    -1 0  1 0;
                     0 0  0 0];
kGlobal2=transpose(matT2)*kLocal2*matT2;
% Element 3
co3=cos(Theta3/180*pi);
si3=sin(Theta3/180*pi);
matT3=[ co3 si3    0   0;
       -si3 co3    0   0;
          0   0  co3 si3;
          0   0 -si3 co3];
kLocal3=(E*A)/(L3)*[ 1 0 -1 0;
                     0 0  0 0;
                    -1 0  1 0;
                     0 0  0 0];
kGlobal3=transpose(matT3)*kLocal3*matT3;
% Element 4
co4=cos(Theta4/180*pi);
si4=sin(Theta4/180*pi);
matT4=[ co4 si4    0   0;
       -si4 co4    0   0;
          0   0  co4 si4;
          0   0 -si4 co4];
kLocal4=(E*A)/(L4)*[ 1 0 -1 0;
                     0 0  0 0;
                    -1 0  1 0;
                     0 0  0 0];
kGlobal4=transpose(matT4)*kLocal4*matT4;
% Element 5
co5=cos(Theta5/180*pi);
si5=sin(Theta5/180*pi);
matT5=[ co5 si5    0   0;
       -si5 co5    0   0;
          0   0  co5 si5;
          0   0 -si5 co5];
kLocal5=(E*A)/(L5)*[ 1 0 -1 0;
                     0 0  0 0;
                    -1 0  1 0;
                     0 0  0 0];
kGlobal5=transpose(matT5)*kLocal5*matT5;

%% Assemblage
KK=(zeros(8,8));
% Element 1
eDoFs1=[2*(1)-1,2*(1),2*(2)-1,2*(2)];
KK(eDoFs1,eDoFs1) = KK(eDoFs1,eDoFs1) + kGlobal1;
% Element 2
eDoFs2=[2*(1)-1,2*(1),2*(4)-1,2*(4)];
KK(eDoFs2,eDoFs2) = KK(eDoFs2,eDoFs2) + kGlobal2;
% Element 3
eDoFs3=[2*(2)-1,2*(2),2*(4)-1,2*(4)];
KK(eDoFs3,eDoFs3) = KK(eDoFs3,eDoFs3) + kGlobal3;
% Element 4
eDoFs4=[2*(2)-1,2*(2),2*(3)-1,2*(3)];
KK(eDoFs4,eDoFs4) = KK(eDoFs4,eDoFs4) + kGlobal4;
% Element 5
eDoFs5=[2*(3)-1,2*(3),2*(4)-1,2*(4)];
KK(eDoFs5,eDoFs5) = KK(eDoFs5,eDoFs5) + kGlobal5;

%% Impose Boundary Conditions
FF=(zeros(8,1));
FF(2*(2), 1) = -20 * 10^3;
FF(2*(3), 1) = -10 * 10^3;
dd=(zeros(8,1));
% 1의 x, 2의 x,y 3의 y, 4의 y
eDoFs_solve = [2*(1)-1, 2*(2)-1, 2*(2), 2*(3), 2*(4)];
dd(eDoFs_solve, 1) = KK(eDoFs_solve, eDoFs_solve) \ FF(eDoFs_solve, 1)

%% Reaction Force Calculation
FF = KK*dd

%% stress calculation
% element 1
stress1 = (1/A) * kLocal1(3, :) * matT1 * dd(eDoFs1, 1)
% element 2
stress2 = (1/A) * kLocal2(3, :) * matT2 * dd(eDoFs2, 1)
% element 3
stress3 = (1/A) * kLocal3(3, :) * matT3 * dd(eDoFs3, 1)
% element 4
stress4 = (1/A) * kLocal4(3, :) * matT4 * dd(eDoFs4, 1)
% element 5
stress5 = (1/A) * kLocal5(3, :) * matT5 * dd(eDoFs5, 1)

%% Visualization
scale=1*10^(2);
dd = dd * scale;
x = 1; y = 2;
Node1_point=[0; 0];
Node2_point=[4; 2];
Node3_point=[8; 4];
Node4_point=[8; 0];
Node1 = [1;2];
Node2 = [3;4];
Node3 = [5;6];
Node4 = [7;8];
Node5 = [9;10];

xInitial_1=[Node1_point(x); Node2_point(x)];
yInitial_1=[Node1_point(y); Node2_point(y)];
xInitial_2=[Node1_point(x); Node4_point(x)];
yInitial_2=[Node1_point(y); Node4_point(y)];
xInitial_3=[Node2_point(x); Node4_point(x)];
yInitial_3=[Node2_point(y); Node4_point(y)];
xInitial_4=[Node2_point(x); Node3_point(x)];
yInitial_4=[Node2_point(y); Node3_point(y)];
xInitial_5=[Node3_point(x); Node4_point(x)];
yInitial_5=[Node3_point(y); Node4_point(y)];
xDeformed_1 = xInitial_1 + [dd([Node1(x), Node2(x)], 1)];
yDeformed_1 = yInitial_1 + [dd([Node1(y), Node2(y)], 1)];
xDeformed_2 = xInitial_2 + [dd([Node1(x), Node4(x)], 1)];
yDeformed_2 = yInitial_2 + [dd([Node1(y), Node4(y)], 1)];
xDeformed_3 = xInitial_3 + [dd([Node2(x), Node4(x)], 1)];
yDeformed_3 = yInitial_3 + [dd([Node2(y), Node4(y)], 1)];
xDeformed_4 = xInitial_4 + [dd([Node2(x), Node3(x)], 1)];
yDeformed_4 = yInitial_4 + [dd([Node2(y), Node3(y)], 1)];
xDeformed_5 = xInitial_5 + [dd([Node3(x), Node4(x)], 1)];
yDeformed_5 = yInitial_5 + [dd([Node3(y), Node4(y)], 1)];
figure(1);
plot(xInitial_1,yInitial_1,'k','LineWidth',1.1);
hold on;
plot(xInitial_2,yInitial_2,'k','LineWidth',1.1);
hold on;
plot(xInitial_3,yInitial_3,'k','LineWidth',1.1);
hold on;
plot(xInitial_4,yInitial_4,'k','LineWidth',1.1);
hold on;
plot(xInitial_5,yInitial_5,'k','LineWidth',1.1);
hold on;
plot(xDeformed_1,yDeformed_1,'b','LineWidth',1.1);
hold on;
plot(xDeformed_2,yDeformed_2,'b','LineWidth',1.1);
hold on;
plot(xDeformed_3,yDeformed_3,'b','LineWidth',1.1);
hold on;
plot(xDeformed_4,yDeformed_4,'b','LineWidth',1.1);
hold on;
plot(xDeformed_5,yDeformed_5,'b','LineWidth',1.1);
hold on;

axis padded;
axis equal;
