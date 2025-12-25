%% practice 5_37
clearvars;
clc; close;
%% Given Parameters
E=200*10^(9);
A=0.04;
I=4*10^(-4);
nEle=44;
nNod=52;
SxxBar=zeros(nEle,1);
SxxBeam=zeros(nEle,1);
%% matP & matE
matP=zeros(nNod,2);
matP(1:25, 1) = 0:1:24;
matP(25:35, 1) = 0:6/10:6;
matP(35:45, 1) = (0:6/10:6) + 18;
matP(1:25, 2) = 0;
matP(25:35, 2) = -8:8/10:0;
matP(35:45, 2) = 0:-8/10:-8;
matP(46, :) = [ 6  0];
matP(47, :) = [10  0];
matP(48, :) = [13  0];
matP(49, :) = [18  0];
matP(50, :) = [24  0];
matP(51, :) = [ 0 -8];
matP(52, :) = [24 -8];
matE=zeros(nEle,2);
for kk=1:1:nEle
    matE(kk,:)=[kk kk+1];
end
Right = [46 47 48 49 50 46 52];
Left = [1 46 47 48 49 51 49];
matE([6 10 13 18 24 34 44], 2) = Right;
matE([1 7 11 14 19 25 35], 1) = Left;
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
    kLocal=zeros(3*2,3*2);
    kLocal([1,4],[1,4])=(E*A)/(Le)*[ 1 -1;
                                    -1  1];
    kLocal([2,3,5,6],[2,3,5,6])=(E*I)/(Le^3)*[  12   6*Le   -12   6*Le;
                                              6*Le 4*Le^2 -6*Le 2*Le^2;
                                               -12  -6*Le    12  -6*Le;
                                              6*Le 2*Le^2 -6*Le 4*Le^2];
    matT=[ co si 0   0  0 0;
          -si co 0   0  0 0;
            0  0 1   0  0 0;
            0  0 0  co si 0;
            0  0 0 -si co 0;
            0  0 0   0  0 1];
    eDofs=[3*(p)-2:3*(p),3*(q)-2:3*(q)];
    KK(eDofs,eDofs)=KK(eDofs,eDofs)+transpose(matT)*kLocal*matT;
end
%% Boundary Conditions
FF=zeros(3*nNod,1);
FF(3*(47)-1, 1) = -18*10^3;
FF(3*(48)-1, 1) = -72*10^3;
FF(3*(49)-1, 1) = -72*10^3;
dd=zeros(3*nNod,1);
FullDoFs=sort([3*(usedNode) - 2, 3*(usedNode) - 1, 3*(usedNode)]);
FixedDoFs=[3*(1)-1, 3*(50)-2:3*(50)-1, 3*(51)-2:3*(51),3*(52)-2:3*(52)];
FreeDoFs=setdiff(FullDoFs,FixedDoFs);
dd(FreeDoFs,1)=KK(FreeDoFs,FreeDoFs)\FF(FreeDoFs,1);
%% Visualization
scale=10^(2);
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
    
    matT=[ co si 0   0  0 0;
          -si co 0   0  0 0;
            0  0 1   0  0 0;
            0  0 0  co si 0;
            0  0 0 -si co 0;
            0  0 0   0  0 1];
    xInitial=[x1 x2];
    yInitial=[y1 y2];
    eDofs=[3*(p)-2:3*(p),3*(q)-2:3*(q)];
    de=dd(eDofs,1);
    xDeformed=[x1+scale*de(1) x2+scale*de(4)];
    yDeformed=[y1+scale*de(2) y2+scale*de(5)];
    dLocal=matT*de;
    d2N1=-(6*Le - 12*(Le/2))/Le^3;
    d2N2=(- 4*Le^2 + 6*(Le/2)*Le)/Le^3;
    d2N3=(6*Le - 12*(Le/2))/Le^3;
    d2N4=(- 2*Le^2 + 6*(Le/2)*Le)/Le^3;
    y=sqrt(3*I/A);
    SxxBar(kk)=(E/Le)*[-1 1]*dLocal([1,4],1);
    SxxBeam(kk)=E*(-y)*[d2N1 d2N2 d2N3 d2N4]*dLocal([2,3,5,6],1);
    if abs(SxxBar(kk)+SxxBeam(kk))>abs(SxxBar(kk)-SxxBeam(kk))
        Smax=SxxBar(kk)+SxxBeam(kk);
    else
        Smax=SxxBar(kk)-SxxBeam(kk);
    end
    SCol=Smax*[1 1];
    figure(1);
    plot(xInitial,yInitial,'k','LineWidth',1.1);
    hold on;
    % plot(xDeformed,yDeformed,'b','LineWidth',1.1);
    % hold on;
    patch(xDeformed,yDeformed,SCol,'EdgeColor','flat','LineWidth',1.1);
    hold on;
end
colormap('jet');
colorbar;
axis equal;
axis padded;
