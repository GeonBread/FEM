clearvars; close;
clc;

%% Given Parameters
E = 70*10^(9);
A = 0.003125;
Theta = [0;
         atand(3/2);
         -atand(1/2);
         90;
         0;
         45;
         -atand(3/2);
         -atand(1/2);
         90;
         0;
         atand(1/2);
         -atand(1/2);
         90;
         0;
         -atand(1/2)];
L = [2;
     sqrt(13);
     sqrt(5);
     3;
     2;
     sqrt(8);
     sqrt(13);
     sqrt(5);
     2;
     2;
     sqrt(5);
     sqrt(5);
     1;
     2;
     sqrt(5)];
Node1 = [1,2];
Node2 = [3,4];
Node3 = [5,6];
Node4 = [7,8];
Node5 = [9,10];
Node6 = [11,12];
Node7 = [13,14];
Node8 = [15,16];
Node9 = [17,18];

%element 총 15개이고, element당 총 4개의 정보를 가져야 함.
eDoFs = zeros(15, 4);
eDoFs(1, :) = [Node1, Node3];
eDoFs(2, :) = [Node1, Node4];
eDoFs(3, :) = [Node2, Node4];
eDoFs(4, :) = [Node3, Node4];
eDoFs(5, :) = [Node3, Node5];
eDoFs(6, :) = [Node3, Node6];
eDoFs(7, :) = [Node4, Node5];
eDoFs(8, :) = [Node4, Node6];
eDoFs(9, :) = [Node5, Node6];
eDoFs(10, :) = [Node5, Node7];
eDoFs(11, :) = [Node5, Node8];
eDoFs(12, :) = [Node6, Node8];
eDoFs(13, :) = [Node7, Node8];
eDoFs(14, :) = [Node7, Node9];
eDoFs(15, :) = [Node8, Node9];

%% Element Stiffness Equation
KK=zeros(2*9,2*9);
for k = 1 : 15  
    L1 = L(k);
    Theta1 = Theta(k);
    
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
    
    % assemble
    eDoFs1=eDoFs(k, :);
    KK(eDoFs1,eDoFs1)=KK(eDoFs1,eDoFs1)+kGlobal1;
end

%% Impose Boundary Conditions
x = 1; y = 2;
FF=zeros(2*9,1);
FF(Node4(y),1) = -1000;
FF(Node6(y),1) = -1000;
FF(Node8(y),1) = -1000;
FF(Node9(y),1) = -500;
dd=zeros(2*9,1);
eDoFs_solve = [2*(3)-1 : 2*(9)];
dd(eDoFs_solve,1) = KK(eDoFs_solve, eDoFs_solve) \ FF(eDoFs_solve,1)
%% Reaction Force Calculation
KK*dd

%% Visualization
Node1_point = [0, 0];
Node2_point = [0, 4];
Node3_point = [2, 0];
Node4_point = [2, 3];
Node5_point = [4, 0];
Node6_point = [4, 2];
Node7_point = [6, 0];
Node8_point = [6, 1];
Node9_point = [8, 0];

scale=5*10^(3);
dd = dd * scale;

xInitial_1=[Node1_point(x); Node3_point(x)];
yInitial_1=[Node1_point(y); Node3_point(y)];
xInitial_2=[Node1_point(x); Node4_point(x)];
yInitial_2=[Node1_point(y); Node4_point(y)];
xInitial_3=[Node2_point(x); Node4_point(x)];
yInitial_3=[Node2_point(y); Node4_point(y)];
xInitial_4=[Node3_point(x); Node4_point(x)];
yInitial_4=[Node3_point(y); Node4_point(y)];
xInitial_5=[Node3_point(x); Node5_point(x)];
yInitial_5=[Node3_point(y); Node5_point(y)];
xInitial_6=[Node3_point(x); Node6_point(x)];
yInitial_6=[Node3_point(y); Node6_point(y)];
xInitial_7=[Node4_point(x); Node5_point(x)];
yInitial_7=[Node4_point(y); Node5_point(y)];
xInitial_8=[Node4_point(x); Node6_point(x)];
yInitial_8=[Node4_point(y); Node6_point(y)];
xInitial_9=[Node5_point(x); Node6_point(x)];
yInitial_9=[Node5_point(y); Node6_point(y)];
xInitial_10=[Node5_point(x); Node7_point(x)];
yInitial_10=[Node5_point(y); Node7_point(y)];
xInitial_11=[Node5_point(x); Node8_point(x)];
yInitial_11=[Node5_point(y); Node8_point(y)];
xInitial_12=[Node6_point(x); Node8_point(x)];
yInitial_12=[Node6_point(y); Node8_point(y)];
xInitial_13=[Node7_point(x); Node8_point(x)];
yInitial_13=[Node7_point(y); Node8_point(y)];
xInitial_14=[Node7_point(x); Node9_point(x)];
yInitial_14=[Node7_point(y); Node9_point(y)];
xInitial_15=[Node8_point(x); Node9_point(x)];
yInitial_15=[Node8_point(y); Node9_point(y)];
xDeformed_1 = xInitial_1 + dd([Node1(x), Node3(x)]);
yDeformed_1 = yInitial_1 + dd([Node1(y), Node3(y)]);
xDeformed_2 = xInitial_2 + dd([Node1(x), Node4(x)]);
yDeformed_2 = yInitial_2 + dd([Node1(y), Node4(y)]);
xDeformed_3 = xInitial_3 + dd([Node2(x), Node4(x)]);
yDeformed_3 = yInitial_3 + dd([Node2(y), Node4(y)]);
xDeformed_4 = xInitial_4 + dd([Node3(x), Node4(x)]);
yDeformed_4 = yInitial_4 + dd([Node3(y), Node4(y)]);
xDeformed_5 = xInitial_5 + dd([Node3(x), Node5(x)]);
yDeformed_5 = yInitial_5 + dd([Node3(y), Node5(y)]);
xDeformed_6 = xInitial_6 + dd([Node3(x), Node6(x)]);
yDeformed_6 = yInitial_6 + dd([Node3(y), Node6(y)]);
xDeformed_7 = xInitial_7 + dd([Node4(x), Node5(x)]);
yDeformed_7 = yInitial_7 + dd([Node4(y), Node5(y)]);
xDeformed_8 = xInitial_8 + dd([Node4(x), Node6(x)]);
yDeformed_8 = yInitial_8 + dd([Node4(y), Node6(y)]);
xDeformed_9 = xInitial_9 + dd([Node5(x), Node6(x)]);
yDeformed_9 = yInitial_9 + dd([Node5(y), Node6(y)]);
xDeformed_10 = xInitial_10 + dd([Node5(x), Node7(x)]);
yDeformed_10 = yInitial_10 + dd([Node5(y), Node7(y)]);
xDeformed_11 = xInitial_11 + dd([Node5(x), Node8(x)]);
yDeformed_11 = yInitial_11 + dd([Node5(y), Node8(y)]);
xDeformed_12 = xInitial_12 + dd([Node6(x), Node8(x)]);
yDeformed_12 = yInitial_12 + dd([Node6(y), Node8(y)]);
xDeformed_13 = xInitial_13 + dd([Node7(x), Node8(x)]);
yDeformed_13 = yInitial_13 + dd([Node7(y), Node8(y)]);
xDeformed_14 = xInitial_14 + dd([Node7(x), Node9(x)]);
yDeformed_14 = yInitial_14 + dd([Node7(y), Node9(y)]);
xDeformed_15 = xInitial_15 + dd([Node8(x), Node9(x)]);
yDeformed_15 = yInitial_15 + dd([Node8(y), Node9(y)]);

figure(1);
plot(xInitial_1,yInitial_1,'k','LineWidth',1.1);hold on;
plot(xInitial_2,yInitial_2,'k','LineWidth',1.1);hold on;
plot(xInitial_3,yInitial_3,'k','LineWidth',1.1);hold on;
plot(xInitial_4,yInitial_4,'k','LineWidth',1.1);hold on;
plot(xInitial_5,yInitial_5,'k','LineWidth',1.1);hold on;
plot(xInitial_6,yInitial_6,'k','LineWidth',1.1);hold on;
plot(xInitial_7,yInitial_7,'k','LineWidth',1.1);hold on;
plot(xInitial_8,yInitial_8,'k','LineWidth',1.1);hold on;
plot(xInitial_9,yInitial_9,'k','LineWidth',1.1);hold on;
plot(xInitial_10,yInitial_10,'k','LineWidth',1.1);hold on;
plot(xInitial_11,yInitial_11,'k','LineWidth',1.1);hold on;
plot(xInitial_12,yInitial_12,'k','LineWidth',1.1);hold on;
plot(xInitial_13,yInitial_13,'k','LineWidth',1.1);hold on;
plot(xInitial_14,yInitial_14,'k','LineWidth',1.1);hold on;
plot(xInitial_15,yInitial_15,'k','LineWidth',1.1);hold on;
plot(xDeformed_1,yDeformed_1,'b','LineWidth',1.1);hold on;
plot(xDeformed_2,yDeformed_2,'b','LineWidth',1.1);hold on;
plot(xDeformed_3,yDeformed_3,'b','LineWidth',1.1);hold on;
plot(xDeformed_4,yDeformed_4,'b','LineWidth',1.1);hold on;
plot(xDeformed_5,yDeformed_5,'b','LineWidth',1.1);hold on;
plot(xDeformed_6,yDeformed_6,'b','LineWidth',1.1);hold on;
plot(xDeformed_7,yDeformed_7,'b','LineWidth',1.1);hold on;
plot(xDeformed_8,yDeformed_8,'b','LineWidth',1.1);hold on;
plot(xDeformed_9,yDeformed_9,'b','LineWidth',1.1);hold on;
plot(xDeformed_10,yDeformed_10,'b','LineWidth',1.1);hold on;
plot(xDeformed_11,yDeformed_11,'b','LineWidth',1.1);hold on;
plot(xDeformed_12,yDeformed_12,'b','LineWidth',1.1);hold on;
plot(xDeformed_13,yDeformed_13,'b','LineWidth',1.1);hold on;
plot(xDeformed_14,yDeformed_14,'b','LineWidth',1.1);hold on;
plot(xDeformed_15,yDeformed_15,'b','LineWidth',1.1);hold on;

axis padded;
axis equal;
