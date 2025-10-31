clearvars;
clc; close;

%% Given Parameters
syms P E k;

%% Element 1
kLocal1=k*[ 1 -1;
           -1  1];
kGlobal1=kLocal1;
%% Element 2
kLocal2=k*[ 1 -1;
           -1  1];
kGlobal2=kLocal2;
%% Element 3
kLocal3=k*[ 1 -1;
           -1  1];
kGlobal3=kLocal3;
%% Element 4
kLocal4=k*[ 1 -1;
           -1  1];
kGlobal4=kLocal4;

%% Assembly
KK=sym(zeros(5,5));
% Element 1
eDoFs1=[1, 2];
KK(eDoFs1,eDoFs1)=KK(eDoFs1,eDoFs1)+kGlobal1;
% Element 2
eDoFs2=[2, 3];
KK(eDoFs2,eDoFs2)=KK(eDoFs2,eDoFs2)+kGlobal2;
% Element 3
eDoFs3=[3, 4];
KK(eDoFs3,eDoFs3)=KK(eDoFs3,eDoFs3)+kGlobal3;
% Element 4
eDoFs4=[4, 5];
KK(eDoFs4,eDoFs4)=KK(eDoFs4,eDoFs4)+kGlobal4

%% Boundary Conditions
FF=sym(zeros(5,1));
FF(3,1) = P; 
dd=sym(zeros(5,1));
eDoFs = [2:4];
% Nodal Displacement
dd(eDoFs, 1) = KK(eDoFs, eDoFs) \ FF(eDoFs, 1)
% Reaction Force
KK * dd
