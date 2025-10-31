%3.34번 풀이

clear; clc; close;
E = 30 * 10^6;
A = 2;
L1 = 8;
L2 = 4;
L3 = 8;
L4 = 4;
L5 = 4*sqrt(5);
theta1 = 0 * (pi/180);
theta2 = 90 * (pi/180);
theta3 = 180 * (pi/180);
theta4 = 90 * (pi/180);
theta5 = pi - atan(4/8);
K = [ 1 0 -1 0;
      0 0  0 0;
     -1 0  1 0;
      0 0  0 0;];
Node1 = [1,2];
Node2 = [3,4];
Node3 = [5,6];
Node4 = [7,8];
Node1_point = [0, 0];
Node2_point = [8, 0];
Node3_point = [8, 4];
Node4_point = [0, 4];

%% element 1
C1 = cos (theta1);
S1 = sin (theta1);
T1 = [C1 S1  0  0;
    -S1 C1  0  0;
      0  0  C1 S1;
      0  0 -S1 C1];
KLocal1 = (E*A/L1) * T1' * K * T1
eDoFs1 = [Node1, Node2]
%% element 2
C2 = cos (theta2);
S2 = sin (theta2);
T2 = [C2 S2  0  0;
    -S2 C2  0  0;
      0  0  C2 S2;
      0  0 -S2 C2];
KLocal2 = (E*A/L2) * T2' * K * T2
eDoFs2 = [Node2, Node3]
%% element 3
C3 = cos (theta3);
S3 = sin (theta3);
T3 = [C3 S3  0  0;
    -S3 C3  0  0;
      0  0  C3 S3;
      0  0 -S3 C3];
KLocal3 = (E*A/L3) * T3' * K * T3
eDoFs3 = [Node3, Node4]
%% element 4
C4 = cos (theta4);
S4 = sin (theta4);
T4 = [C4 S4  0  0;
    -S4 C4  0  0;
      0  0  C4 S4;
      0  0 -S4 C4];
KLocal4 = (E*A/L4) * T4' * K * T4
eDoFs4 = [Node1, Node4]
%% element 5
C5 = cos (theta5);
S5 = sin (theta5);
T5 = [C5 S5  0  0;
    -S5 C5  0  0;
      0  0  C5 S5;
      0  0 -S5 C5];
KLocal5 = (E*A/L5) * T5' * K * T5
eDoFs5 = [Node2, Node4]
KGlobal = zeros (8, 8);
KGlobal(eDoFs1, eDoFs1) = KGlobal(eDoFs1, eDoFs1) + KLocal1;
KGlobal(eDoFs2, eDoFs2) = KGlobal(eDoFs2, eDoFs2) + KLocal2;
KGlobal(eDoFs3, eDoFs3) = KGlobal(eDoFs3, eDoFs3) + KLocal3;
KGlobal(eDoFs4, eDoFs4) = KGlobal(eDoFs4, eDoFs4) + KLocal4;
KGlobal(eDoFs5, eDoFs5) = KGlobal(eDoFs5, eDoFs5) + KLocal5

%B.C
u = zeros(8, 1);

u(4) = -0.02;
eDoFs = [7, 8];
F = 0 - KGlobal([7, 8], 4) * u(4)

u(eDoFs, 1) = u(eDoFs, 1) + KGlobal(eDoFs, eDoFs) \ F
F = KGlobal * u

%visualization
scale = 10^2;
u = u * scale;
x1_init = [Node1_point(1); Node2_point(1)];
y1_init = [Node1_point(2); Node2_point(2)];
x2_init = [Node2_point(1); Node3_point(1)];
y2_init = [Node2_point(2); Node3_point(2)];
x3_init = [Node3_point(1); Node4_point(1)];
y3_init = [Node3_point(2); Node4_point(2)];
x4_init = [Node1_point(1); Node4_point(1)];
y4_init = [Node1_point(2); Node4_point(2)];
x5_init = [Node2_point(1); Node4_point(1)];
y5_init = [Node2_point(2); Node4_point(2)];
x1_deformed = x1_init + u([Node1(1), Node2(1)], 1);
y1_deformed = y1_init + u([Node1(2), Node2(2)], 1);
x2_deformed = x2_init + u([Node2(1), Node3(1)], 1);
y2_deformed = y2_init + u([Node2(2), Node3(2)], 1);
x3_deformed = x3_init + u([Node3(1), Node4(1)], 1);
y3_deformed = y3_init + u([Node3(2), Node4(2)], 1);
x4_deformed = x4_init + u([Node1(1), Node4(1)], 1);
y4_deformed = y4_init + u([Node1(2), Node4(2)], 1);
x5_deformed = x5_init + u([Node2(1), Node4(1)], 1);
y5_deformed = y5_init + u([Node2(2), Node4(2)], 1);
plot(x1_init, y1_init, color='b', LineWidth=1.1); hold on;
plot(x2_init, y2_init, color='b', LineWidth=1.1); hold on;
plot(x3_init, y3_init, color='b', LineWidth=1.1); hold on;
plot(x4_init, y4_init, color='b', LineWidth=1.1); hold on;
plot(x5_init, y5_init, color='b', LineWidth=1.1); hold on;
plot(x1_deformed, y1_deformed, color='r', LineWidth=1.1); hold on;
plot(x2_deformed, y2_deformed, color='r', LineWidth=1.1); hold on;
plot(x3_deformed, y3_deformed, color='r', LineWidth=1.1); hold on;
plot(x4_deformed, y4_deformed, color='r', LineWidth=1.1); hold on;
plot(x5_deformed, y5_deformed, color='r', LineWidth=1.1); hold on;

axis padded;
axis equal
