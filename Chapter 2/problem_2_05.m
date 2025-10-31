clearvars;
clc; close;
%% Given Parameters
k1 = 1;
k2 = 2;
k3 = 3;
k4 = 4;
k5 = 5;
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
kLocal5=k5*[1  -1;
             -1 1];
kGlobal5=kLocal5;
%% Assembly
KK=sym(zeros(4,4));
% Element 1
eDoFs1=[1, 2];
KK(eDoFs1,eDoFs1)=KK(eDoFs1,eDoFs1)+kGlobal1;
% Element 2
eDoFs2=[2, 4];
KK(eDoFs2,eDoFs2)=KK(eDoFs2,eDoFs2)+kGlobal2;
% Element 3
eDoFs3=[2, 4];
KK(eDoFs3,eDoFs3)=KK(eDoFs3,eDoFs3)+kGlobal3;
% Element 4
eDoFs4=[2, 4];
KK(eDoFs4,eDoFs4)=KK(eDoFs4,eDoFs4)+kGlobal4;
% Element 5
eDoFs5=[3, 4];
KK(eDoFs5,eDoFs5)=KK(eDoFs5,eDoFs5)+kGlobal5
