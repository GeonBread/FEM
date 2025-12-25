%% 연습 문제 5-57번
clearvars;
clc; close;
%% Given Parameters
E=210*10^(9);
G=84*10^(9);
I=2*10^(-4);
J=1*10^(-4);
A = 1 * 10^-2;
nNod=34;
nEle=30;
%% matP & matE
matP=zeros(nNod,2);
matP(1:11,1)=0:0.3:3;
matP(11:21,1)=6:-3/10:3;
matP(21:31,1)=0:0.3:3;
matP(1:11,2)=0:-0.51962:-5.1962;
matP(11:21,2)=0:-0.51962:-5.1962;
matP(21:31,2)=-6.5:1.3038/10:-5.1962;
matP(32, :) = [6 0];
matP(33, :) = [3 -5.1962];
matP(34, :) = [0 -6.5];
matE=zeros(nEle,2);
for kk=1:1:nEle
    matE(kk,:)=[kk kk+1];
end
matE(10, :) = [10 33];
matE(11, :) = [32 12];
matE(20, :) = [20 33];
matE(21, :) = [34 22];
matE(30, :) = [30 33];
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
FF(3*33-2,1)=-100*10^(3);
FF(3*16-2,1)=-100*10^(3);
dd=zeros(3*nNod,1);
FullDoFs=sort([3*(usedNode)-2, 3*(usedNode)-1, 3*(usedNode)]);
FixedDoFs=[3*1-2:3*1,3*32-2:3*32, 3*34-2:3*34];
FreeDoFs=setdiff(FullDoFs,FixedDoFs);
dd(FreeDoFs,1)=KK(FreeDoFs,FreeDoFs)\FF(FreeDoFs,1);
%% Reaction Force
KK(3*1-2:3*1,:)*dd
KK(3*32-2:3*32,:)*dd
KK(3*34-2:3*34,:)*dd
%% Visualization
scale=100/2;
allStress = zeros(nEle, 1);
for kk=1:1:nEle
    p=matE(kk,1);
    q=matE(kk,2);
    x1=matP(p,1);
    y1=matP(p,2);
    x2=matP(q,1);
    y2=matP(q,2);
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
    allStress(kk, 1) = Svm;
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
allStress
colormap jet;
colorbar;
axis equal;
axis padded;
