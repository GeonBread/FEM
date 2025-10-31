clearvars;
clc; close;
%% Given Parameters
k1 = 500;
k2 = 300;
k3 = 300;
k4 = 400;
k5 = 400;
%% Element 1
kLocal1=k1*[ 1 -1;
           -1  1];
kGlobal1=kLocal1;
%% Element 2
kLocal2=k2*[ 1 -1;
           -1  1];
kGlobal2=kLocal2;
%% Element 3
kLocal3=k3*[ 1 -1;
           -1  1];
kGlobal3=kLocal3;
%% Element 4
kLocal4=k4*[ 1 -1;
           -1  1];
kGlobal4=kLocal4;
%% Element 5
kLocal5=k5*[ 1 -1;
           -1  1];
kGlobal5=kLocal5;
%% Assembly
KK=zeros(4,4);
% Element 1
eDoFs1=[1, 2];
KK(eDoFs1,eDoFs1)=KK(eDoFs1,eDoFs1)+kGlobal1;
% Element 2
eDoFs2=[2, 3];
KK(eDoFs2,eDoFs2)=KK(eDoFs2,eDoFs2)+kGlobal2;
% Element 3
eDoFs3=[2, 3];
KK(eDoFs3,eDoFs3)=KK(eDoFs3,eDoFs3)+kGlobal3;
% Element 4
eDoFs4=[3, 4];
KK(eDoFs4,eDoFs4)=KK(eDoFs4,eDoFs4)+kGlobal4;
% Element 5
eDoFs5=[2, 4];
KK(eDoFs5,eDoFs5)=KK(eDoFs5,eDoFs5)+kGlobal5;
%% Boundary Conditions
FF=zeros(4,1);
dd=zeros(4,1);
FF(3, 1) = 1000;
eDoFs = [2, 3];
dd(eDoFs,1) = KK(eDoFs, eDoFs) \ FF(eDoFs, 1)
Reaction_Force = KK * dd
