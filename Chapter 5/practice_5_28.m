%% practice_5_28
clearvars;
clc; close;
%% Given Parameters
E=200*10^(9);
A=0.5 * 10^-2;
I=0.5*10^(-4);
nEle=40;
nNod=45;
SxxBar=zeros(nEle,1);
SxxBeam=zeros(nEle,1);
%% matP & matE
matP=zeros(nNod,2);
matP(1:11, 1) = 0;
matP(11:21, 1) = 0;
matP(21:31, 1) = 0:10/10:10;
matP(31:41, 1) = 0:10/10:10;
matP(1:11, 2) = 0:4/10:4;
matP(11:21, 2) = [0:4/10:4] + 4;
matP(21:31, 2) = 4;
matP(31:41, 2) = 8;
matP(42, :) = [0 4];
matP(43, :) = [10 4];
matP(44, :) = [0 8];
matP(45, :) = [10 8];
matE=zeros(nEle,2);
for kk=1:1:nEle
    matE(kk,:)=[kk kk+1];
end
Right = [42 44 43 45];
Left = [1 42 42 44];
matE([10 20 30 40], 2) = Right';
matE([1 11 21 31], 1) = Left';
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
for kk=1:1:20
    p=matE(kk,1);
    q=matE(kk,2);
    x1=matP(p,1);
    y1=matP(p,2);
    x2=matP(q,1);
    y2=matP(q,2);
    eDoFs=[3*(kk)-2:3*(kk),3*(kk+1)-2:3*(kk+1)];
    Le=sqrt((x2-x1)^2+(y2-y1)^2);
    a = -50*10^3;
    w = a*(kk*Le - 8);
    %rectangle
    FF(eDoFs,1)=FF(eDoFs,1)+[ (w*Le)/2;
                             0;
                             -(w*Le^2)/12;
                              (w*Le)/2;
                             0;
                              (w*Le^2)/12];
    %triangle
    FF(eDoFs,1)=FF(eDoFs,1)+[-(7*Le^2*a)/20;
                             0;
                              (Le^3*a)/20;
                             -(3*Le^2*a)/20;
                             0;
                             -(Le^3*a)/30];

                            [-(w*Le)/2;
                             0;
                              (w*Le^2)/12;
                             -(w*Le)/2;
                             0;
                             -(w*Le^2)/12]
    %triangle
                             [-(7*Le^2*a)/20;
                             0;
                              (Le^3*a)/20;
                             -(3*Le^2*a)/20;
                             0;
                             -(Le^3*a)/30]
end
dd=zeros(3*nNod,1);
FullDoFs=sort([3*(usedNode)-2,3*(usedNode)-1,3*(usedNode)]);
FixedDoFs=[3*(1)-2:3*(1),3*(43)-2:3*(43), 3*(45)-2:3*(45)];
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
