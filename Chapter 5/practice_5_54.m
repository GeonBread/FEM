%% practice_5_5
clearvars;
clc;close
nNod=251 + 13;
nEle=250;
%% Given Parameters
E=30*10^(6);
G=12*10^(6);
% I=50;
% J=20;
matI = zeros(nEle, 1);
matI(1:120, 1) = 200;
matI(121:200, 1) = 50;
matI(201:250, 1) = 200;
% matI(1:250, 1) = 200;
matJ = zeros(nEle, 1);
matJ(1:120, 1) = 100;
matJ(121:200, 1) = 20;
matJ(201:250, 1) = 100;
% matJ(1:250, 1) = 20;
%% matP & matE
matP=zeros(nNod,2);
matP(  1:1: 61, 1) = 0:36/60:36;
matP( 61:1:121, 1) = 0:36/60:36;
matP(121:1:131, 1) = 6:6/10:12;
matP(131:1:141, 1) = 6:6/10:12;
matP(141:1:151, 1) = [6:6/10:12] + 6;
matP(151:1:161, 1) = [6:6/10:12] + 6;
matP(161:1:171, 1) = [6:6/10:12] + 12;
matP(171:1:181, 1) = [6:6/10:12] + 12;
matP(181:1:191, 1) = [6:6/10:12] + 18;
matP(191:1:201, 1) = [6:6/10:12] + 18;
matP(201:1:211, 1) = 6;
matP(211:1:221, 1) = 12;
matP(221:1:231, 1) = 18;
matP(231:1:241, 1) = 24;
matP(241:1:251, 1) = 30;
matP(  1:1: 61, 2) = 0;
matP( 61:1:121, 2) = 10;
matP(121:1:131, 2) = 10:-1:0;
matP(131:1:141, 2) = 0:1:10;
matP(141:1:151, 2) = 10:-1:0;
matP(151:1:161, 2) = 0:1:10;
matP(161:1:171, 2) = 10:-1:0;
matP(171:1:181, 2) = 0:1:10;
matP(181:1:191, 2) = 10:-1:0;
matP(191:1:201, 2) = 0:1:10;
matP(201:1:211, 2) = 10:-1:0;
matP(211:1:221, 2) = 10:-1:0;
matP(221:1:231, 2) = 10:-1:0;
matP(231:1:241, 2) = 10:-1:0;
matP(241:1:251, 2) = 10:-1:0;
matP(252, :) = [6 0];
matP(253, :) = [12 0];
matP(254, :) = [18 0];
matP(255, :) = [24 0];
matP(256, :) = [30 0];
matP(257, :) = [36 0];
matP(258, :) = [0 10];
matP(259, :) = [6 10];
matP(260, :) = [12 10];
matP(261, :) = [18 10];
matP(262, :) = [24 10];
matP(263, :) = [30 10];
matP(264, :) = [36 10];
matE=zeros(nEle,2);
for kk=1:1:nEle
    matE(kk,:)=[kk kk+1];
end
Right = [252 253 254 255 256 257 ...
         259 260 261 262 263 264 ...
         253 260 ...
         254 261 ...
         255 262 ...
         256 263 ...
         252 253 254 255 256];
Left = [  1 252 253 254 255 256 ...
        258 259 260 261 262 263 ...
        259 252 ...
        260 253 ...
        261 254 ...
        262 255 ...
        259 260 261 262 263];
matE(10:10:250, 2) = Right;                                                                                                                                                                                 
matE(1:10:241, 1) = Left;
usedNode = unique(matE);
%% Global Stiffness Equation
KK=zeros(3*nNod,3*nNod);
for kk=1:1:nEle
    p=matE(kk,1);
    q=matE(kk,2);
    x1=matP(p,1);
    y1=matP(p,2);
    x2=matP(q,1);
    y2=matP(q,2);
    J = matJ(kk);
    I = matI(kk);
    Le=sqrt((x2-x1)^2+(y2-y1)^2);
    co=(x2-x1)/Le;
    si=(y2-y1)/Le;
    kLocal=zeros(6,6);
    kLocal([2,5],[2,5])=(G*J)/Le*[ 1 -1;
                                  -1  1];
    kLocal([1,3,4,6],[1,3,4,6])=(E*I)/(Le^3)*[   12  -6*Le   -12  -6*Le;
                                              -6*Le 4*Le^2  6*Le 2*Le^2;
                                                -12   6*Le    12   6*Le;
                                              -6*Le 2*Le^2  6*Le 4*Le^2];
    matT=[1   0  0 0   0  0;
          0  co si 0   0  0;
          0 -si co 0   0  0;
          0   0  0 1   0  0;
          0   0  0 0  co si;
          0   0  0 0 -si co];
    eDoFs=[3*p-2:3*p,3*q-2:3*q];
    KK(eDoFs,eDoFs)=KK(eDoFs,eDoFs)+transpose(matT)*kLocal*matT;
end
%% Boundary Conditions
FF=zeros(3*nNod,1);
FF(3*([252, 253, 254, 255, 256])-2, 1) = -1;
FF(3*([259, 260, 261, 262, 263])-2, 1) = -1;
dd=zeros(3*nNod,1);
FullDoFs=sort([3*(usedNode)-2, 3*(usedNode)-1, 3*(usedNode)]);
FixedDoFs=[3*(1)-2:3*(1)-1,3*(258)-2:3*(258)-1, ...
           3*(257)-2, 3*(264)-2];

FreeDoFs=setdiff(FullDoFs,FixedDoFs);
dd(FreeDoFs,1)=KK(FreeDoFs,FreeDoFs)\FF(FreeDoFs,1);
%% Reaction Force
% KK(3*1-2:3*1,:)*dd
% KK(3*21-2:3*21,:)*dd
%% Visualization
scale=1000000 * 5;
for kk=1:1:nEle
    p=matE(kk,1);
    q=matE(kk,2);
    x1=matP(p,1);
    y1=matP(p,2);
    x2=matP(q,1);
    y2=matP(q,2);
    if(mod(kk, 11) == 0)
        1;
    end
    Le=sqrt((x2-x1)^2+(y2-y1)^2);
    co=(x2-x1)/Le;
    si=(y2-y1)/Le;
    matT=[1   0  0 0   0  0;
          0  co si 0   0  0;
          0 -si co 0   0  0;
          0   0  0 1   0  0;
          0   0  0 0  co si;
          0   0  0 0 -si co];
    xInitial=[x1 x2];
    yInitial=[y1 y2];
    zInitial=[ 0  0];
    eDoFs=[3*p-2:3*p,3*q-2:3*q];
    de=dd(eDoFs,1);
    dLocal=matT*de;
    Sshear=(G/Le)*(2*J/pi)^(1/4)*[-1 1]*[dLocal(2);
                                         dLocal(5)];
    d2N1=-(6*Le - 12*(Le/2))/Le^3;
    d2N2=(- 4*Le^2 + 6*(Le/2)*Le)/Le^3;
    d2N3=(6*Le - 12*(Le/2))/Le^3;
    d2N4=(- 2*Le^2 + 6*(Le/2)*Le)/Le^3;
    Snormal=(E)*(2*J/pi)^(1/4)*[-d2N1 d2N2 -d2N3 d2N4]*[dLocal(1);
                                                        dLocal(3);
                                                        dLocal(4);
                                                        dLocal(6)];
    Svm=sqrt(Snormal^(2)+3*(Sshear)^2);
    SCol=Svm*[1 1];
    xDeformed=[x1 x2];
    yDeformed=[y1 y2];
    zDeformed=[scale*de(1) scale*de(4)];
    figure(1);
    plot3(xInitial,yInitial,zInitial,'k','LineWidth',1.1);
    hold on;
    % plot3(xDeformed,yDeformed,zDeformed,'b','LineWidth',1.1);
    % hold on;
    patch(xDeformed,yDeformed,zDeformed,SCol,'EdgeColor','flat','LineWidth',1.1);
    hold on;
end
colormap jet;
colorbar;
axis equal;
axis padded;
