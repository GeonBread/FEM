%% 연습문제 5-23
clearvars;
clc; close;
nEle = 100;
nNod = 106;
%% Given Parameters
E=30*10^(6);
matA=zeros(nEle, 1);
matA([1:10, 31:40, 51:60, 81:90]) = 10;
matA([41:50]) = 12;
matA([91:100]) = 12;
matA([11:30, 61:80]) = 2;
matI=zeros(nEle, 1);
matI([1:10, 31:40, 51:60, 81:90]) = 150;
matI([41:50]) = 300;
matI([91:100]) = 200;
matI([11:30, 61:80]) = 1;
SxxBar=zeros(nEle,1);
SxxBeam=zeros(nEle,1);
%% matP & matE
matP=zeros(nNod,2);
matP(1:11, :) = [zeros(11, 1), (0:15/10:15)'];
matP(11:21, :) = [(0:25/10:25)', (15:-15/10:0)'];
matP(21:31, :) = [(0:25/10:25)', (0:15/10:15)'];
matP(31:41, :) = [25*ones(11, 1), (0:15/10:15)'];
matP(41:51, :) = [(0:25/10:25)', 15*ones(11, 1)];
matP(51:61, :) = [zeros(11, 1), (0:15/10:15)'+15];
matP(61:71, :) = [(0:25/10:25)', (15:-15/10:0)'+15];
matP(71:81, :) = [(0:25/10:25)', (0:15/10:15)'+15];
matP(81:91, :) = [25*ones(11, 1), (0:15/10:15)'+15];
matP(91:101, :) = [(0:25/10:25)', 15*ones(11, 1)+15];
matP(102, :) = [0, 15];
matP(103, :) = [0, 30];
matP(104, :) = [25, 30];
matP(105, :) = [25, 15];
matP(106, :) = [25, 0];
matE=zeros(nEle,2);
Right=[102 106 105 105 105 ...
       103 105 104 104 104];
Left=[  1 102   1 106 102 ...
      102 103 102 105 103];
for kk=1:1:nEle
    matE(kk,:)=[kk kk+1];
end
matE(10:10:100, 2) = Right';
matE(1:10:91, 1) = Left';
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
    A = matA(kk, 1);
    I = matI(kk, 1);
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
FF(3*(102)-2, 1) = 2400;
FF(3*(103)-2, 1) = 1200;
FF(3*(103), 1) = 3000;
dd=zeros(3*nNod,1);
FullDoFs=sort([3*(usedNode)-2, 3*(usedNode)-1, 3*(usedNode)]);
FixedDoFs=[3*(1)-2:3*(1)-1,3*(106)-2:3*(106)-1];
FreeDoFs=setdiff(FullDoFs,FixedDoFs);
dd(FreeDoFs,1)=KK(FreeDoFs,FreeDoFs)\FF(FreeDoFs,1);
%% Visualization
scale=5*10^(3);
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
