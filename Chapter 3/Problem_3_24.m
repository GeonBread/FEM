
clear; clc; close;
E = 10^12;
A = 10^-3;
K = [1 0 -1 0;
     0 0 0 0;
     -1 0 1 0;
     0 0 0 0];
Node1 = [1, 2];
Node2 = [3, 4];
Node3 = [5, 6];
Node4 = [7, 8];
%% element 1
theta1 = 0 * pi/180;
C1 = cos (theta1);
S1 = sin (theta1);
Le1 = 20;
T1 = [C1 S1 0 0;
     -S1 C1 0 0;
     0 0 C1 S1;
     0 0 -S1 C1]
kLocal1 = ((E*A)/Le1) * T1' * K * T1;
eDoFs1 = [Node1, Node2];
%% element 2
theta2 = (atan (15/20)) * pi/180;
C2 = cos (theta2);
S2 = sin (theta2);
Le2 = 25;
T2 = [C2 S2 0 0;
     -S2 C2 0 0;
     0 0 C2 S2;
     0 0 -S2 C2];
kLocal2 = ((E*A)/Le2) * T2' * K * T2;
eDoFs2 = [Node1, Node3];
%% element 3
theta3 = (90 + atan (20/15)) * pi/180;
C3 = cos (theta3);
S3 = sin (theta3);
Le3 = 25;
T3 = [C3 S3 0 0;
     -S3 C3 0 0;
     0 0 C3 S3;
     0 0 -S3 C3];
kLocal3 = ((E*A)/Le3) * T3' * K * T3;
eDoFs3 = [Node2, Node4];
%% element 4
theta4 = 90 * pi/180;
C4 = cos (theta4);
S4 = sin (theta4);
Le4 = 15;
T4 = [C4 S4 0 0;
     -S4 C4 0 0;
     0 0 C4 S4;
     0 0 -S4 C4];
kLocal4 = ((E*A)/Le4) * T4' * K * T4;
eDoFs4 = [Node2, Node3];
%% element 5
theta5 = 180 * pi/180;
C5 = cos (theta5)
S5 = sin (theta5)
Le5 = 20;
T5 = [C5 S5 0 0;
     -S5 C5 0 0;
     0 0 C5 S5;
     0 0 -S5 C5];
kLocal5 = ((E*A)/Le5) * T5' * K * T5;
eDoFs5 = [Node3, Node4];
%% element 5
theta6 = 90 * pi/180;
C6 = cos (theta6);
S6 = sin (theta6);
Le6 = 15;
T6 = [C6 S6 0 0;
     -S6 C6 0 0;
     0 0 C6 S6;
     0 0 -S6 C6];
kLocal6 = ((E*A)/Le6) * T6' * K * T6;
eDoFs6 = [Node1, Node4];
%%Assemble
kGlobal = zeros(8, 8);
kGlobal(eDoFs1, eDoFs1) = kGlobal(eDoFs1, eDoFs1) + kLocal1;
kGlobal(eDoFs2, eDoFs2) = kGlobal(eDoFs2, eDoFs2) + kLocal2;
kGlobal(eDoFs3, eDoFs3) = kGlobal(eDoFs3, eDoFs3) + kLocal3;
kGlobal(eDoFs4, eDoFs4) = kGlobal(eDoFs4, eDoFs4) + kLocal4;
kGlobal(eDoFs5, eDoFs5) = kGlobal(eDoFs5, eDoFs5) + kLocal5;
kGlobal(eDoFs6, eDoFs6) = kGlobal(eDoFs6, eDoFs6) + kLocal6
%% Solve
FF = [0; 1000; 0; 1000];
uu = zeros(8, 1);
eDoFs = [[Node2, Node3]]
uu(eDoFs, 1) = uu(eDoFs, 1) + kGlobal(eDoFs, eDoFs) \ FF
FF = kGlobal * uu
%%visualization
initial_x1 = [0, 20];
initial_y1 = [0, 0];
initial_x2 = [0, 20];
initial_y2 = [0, 15];
initial_x3 = [20, 0];
initial_y3 = [0, 15];
initial_x4 = [20, 20];
initial_y4 = [0, 15];
initial_x5 = [20, 0];
initial_y5 = [15, 15];
initial_x6 = [0, 0];
initial_y6 = [0, 15];
scale = 10^5;
uu = uu * scale;
deform_x1 = initial_x1 + uu([Node1(1), Node2(1)])';
deform_x2 = initial_x2 + uu([Node1(1), Node3(1)])';
deform_x3 = initial_x3 + uu([Node2(1), Node4(1)])';
deform_x4 = initial_x4 + uu([Node2(1), Node3(1)])';
deform_x5 = initial_x5 + uu([Node3(1), Node4(1)])';
deform_x6 = initial_x6 + uu([Node1(1), Node4(1)])';
deform_y1 = initial_y1 + uu([Node1(2), Node2(2)])';
deform_y2 = initial_y2 + uu([Node1(2), Node3(2)])';
deform_y3 = initial_y3 + uu([Node2(2), Node4(2)])';
deform_y4 = initial_y4 + uu([Node2(2), Node3(2)])';
deform_y5 = initial_y5 + uu([Node3(2), Node4(2)])';
deform_y6 = initial_y6 + uu([Node1(2), Node4(2)])';
plot(initial_x1, initial_y1, LineWidth=1.1, Color='b');
hold on;
plot(initial_x2, initial_y2, LineWidth=1.1, Color='b');
hold on;
plot(initial_x3, initial_y3, LineWidth=1.1, Color='b');
hold on;
plot(initial_x4, initial_y4, LineWidth=1.1, Color='b');
hold on;
plot(initial_x5, initial_y5, LineWidth=1.1, Color='b');
hold on;
plot(initial_x6, initial_y6, LineWidth=1.1, Color='b');
hold on;
plot(deform_x1, deform_y1, LineWidth=1.1, Color='r');
hold on;
plot(deform_x2, deform_y2, LineWidth=1.1, Color='r');
hold on;
plot(deform_x3, deform_y3, LineWidth=1.1, Color='r');
plot([deform_x3(2), deform_x3(1)], [deform_y3(2), deform_y3(1)], LineWidth=1.1, Color='r');
hold on;
plot(deform_x4, deform_y4, LineWidth=1.1, Color='r');
hold on;
plot(deform_x5, deform_y5, LineWidth=1.1, Color='r');
hold on;
plot(deform_x6, deform_y6, LineWidth=1.1, Color='r');
hold on;

axis padded;
axis equal;