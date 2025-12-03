%% 연습 문제 5-54
clearvars;
clc; close
%% Given Parameters
E=30 * 10^6;
G=12 * 10^6;
I=50;
J=20;
nNod=264;
nEle=250;
%% matP & matE
matP=zeros(nNod,2);
matP(1:11,:)=[(0:(6/10):6)', zeros(11, 1)];
matP(11:21,:)=[(6:(6/10):12)', zeros(11, 1)];
matP(21:31,:)=[(12:(6/10):18)', zeros(11, 1)];
matP(31:41,:)=[(18:(6/10):24)', zeros(11, 1)];
matP(41:51,:)=[(24:(6/10):30)', zeros(11, 1)];
matP(51:61,:)=[(30:(6/10):36)', zeros(11, 1)];
matP(61:71,:)=[(0:(6/10):6)', -10 * ones(11, 1)];
matP(71:81,:)=[(6:(6/10):12)', -10 * ones(11, 1)];
matP(81:91,:)=[(12:(6/10):18)', -10 * ones(11, 1)];
matP(91:101,:)=[(18:(6/10):24)', -10 * ones(11, 1)];
matP(101:111,:)=[(24:(6/10):30)', -10 * ones(11, 1)];
matP(111:121,:)=[(30:(6/10):36)', -10 * ones(11, 1)];
matP(121:131,:)=[6 * ones(11, 1), (-10:(10/10):0)'];
matP(131:141,:)=[(6:(6/10):12)', (0:(-10/10):-10)'];
matP(141:151,:)=[(6:(6/10):12)', (-10:(10/10):0)'];
matP(151:161,:)=[12 * ones(11, 1), (-10:(10/10):0)'];
matP(161:171,:)=[(12:(6/10):18)', (0:(-10/10):-10)'];
matP(171:181,:)=[(12:(6/10):18)', (-10:(10/10):0)'];
matP(181:191,:)=[18 * ones(11, 1), (-10:(10/10):0)'];
matP(191:201,:)=[(18:(6/10):24)', (0:(-10/10):-10)'];
matP(201:211,:)=[(18:(6/10):24)', (-10:(10/10):0)'];
matP(211:221,:)=[24 * ones(11, 1), (-10:(10/10):0)'];
matP(221:231,:)=[(24:(6/10):30)', (0:(-10/10):-10)'];
matP(231:241,:)=[(24:(6/10):30)', (-10:(10/10):0)'];
matP(241:251,:)=[30 * ones(11, 1), (-10:(10/10):0)'];
matP(252, :) = [6, 0];
matP(253, :) = [12, 0];
matP(254, :) = [18, 0];
matP(255, :) = [24, 0];
matP(256, :) = [30, 0];
matP(257, :) = [36, 0];
matP(258, :) = [0, -10];
matP(259, :) = [6, -10];
matP(260, :) = [12, -10];
matP(261, :) = [18, -10];
matP(262, :) = [24, -10];
matP(263, :) = [30, -10];
matP(264, :) = [36, -10];
matE=zeros(nEle,2);
for kk=1:1:nEle
    matE(kk,:)=[kk kk+1];
end
%10, 20, 30 ...
LastElementLeft = (10:10:nEle)';
LastElementRight = [252;253;254;255;256;257;
                     259;260;261;262;263;264;
                     252;
                     260;253;
                     253;
                     261;254;
                     254;
                     262;255;
                     255;
                     263;256;
                     256];
matE(10:10:nEle, :) = [LastElementLeft, LastElementRight];
%11, 21, 31 ...
StartElementRight = (12:10:nEle)';
StartElementLeft = [252;253;254;255;256;
                     258;259;260;261;262;263;
                     259;
                     252;259;
                     260;
                     253;260;
                     261;
                     254;261;
                     262;
                     255;262;
                     263];
matE(11:10:nEle, :) = [StartElementLeft, StartElementRight];
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
FF(3*(252)-2,1)=-1;
FF(3*(253)-2,1)=-1;
FF(3*(254)-2,1)=-1;
FF(3*(255)-2,1)=-1;
FF(3*(256)-2,1)=-1;
FF(3*(259)-2,1)=-1;
FF(3*(260)-2,1)=-1;
FF(3*(261)-2,1)=-1;
FF(3*(262)-2,1)=-1;
FF(3*(263)-2,1)=-1;
dd=zeros(3*nNod,1);
FullDoFs=sort([3*(usedNode)-2, 3*(usedNode)-1, 3*(usedNode)]);
FixedDoFs=[3*(1)-2:3*(1)-1, 3*(258)-2:3*(258)-1,3*(257)-2, 3*(264)-2 ];
FreeDoFs=setdiff(FullDoFs,FixedDoFs);
dd(FreeDoFs,1)=KK(FreeDoFs,FreeDoFs)\FF(FreeDoFs,1);
%% Reaction Force
KK(3*1-2:3*1,:)*dd
KK(3*21-2:3*21,:)*dd
%% Visualization
scale=1000000 * 5;
for kk=1:1:nEle
    if mod(kk, 10) == 0
        ;
    end
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
