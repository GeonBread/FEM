clearvars;
clc; close;
%% Given Parameters
k1 = 1000;
k2 = 3000;
%% Element 1
kLocal1=k1*[ 1 -1;
           -1  1];
kGlobal1=kLocal1;
%% Element 2
kLocal2=k2*[ 1 -1;
           -1  1];
kGlobal2=kLocal2;

%% Assembly
KK=zeros(3,3);
% Element 1
eDoFs1=[1, 2];
KK(eDoFs1,eDoFs1)=KK(eDoFs1,eDoFs1)+kGlobal1;
% Element 2
eDoFs2=[2,3];
KK(eDoFs2,eDoFs2)=KK(eDoFs2,eDoFs2)+kGlobal2
%% Boundary Conditions
FF=zeros(3,1);
dd=zeros(3,1);

dd(3, 1) = 20 * 10^-3;
eDoFs = [2];
FF(eDoFs, 1) = - KK(eDoFs, 3) * dd(3, 1)
dd(eDoFs,1) = KK(eDoFs, eDoFs) \ FF(eDoFs, 1)
Reaction_Force = KK * dd
eForce1 = kLocal1 * dd(eDoFs1, 1)
eForce2 = kLocal2 * dd(eDoFs2, 1)
