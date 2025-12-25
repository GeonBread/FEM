%% 연습 문제 5-36 Element 수 임의로 변경하기
clearvars;
clc; close;
%% Given Parameters
E=200*10^(9);
A=1*10^-2;
I=1*10^(-4);
w=-13*10^(3);
nEle=400;
nNod=401 + 7;
SxxBar=zeros(nEle,1);
SxxBeam=zeros(nEle,1);
%% matP & matE
matP=zeros(nNod,2);
matP(1:11, 1)=0;
matP(11:21, 1)=0:2/10:2;
matP(21:31, 1)=0:2/10:2;
matP(31:41, 1)=0:2/10:2;
matP(41:51, 1)=2;
matP(51:86, 1)=2:3/35:5;
matP(86:121, 1)=2:3/35:5;
matP(121:156, 1)=2:3/35:5;
matP(156:196, 1)=2:3/40:5;
matP(196:241, 1)=5;
matP(241:281, 1)=5:1.2/40:6.2;
matP(281:321, 1)=5:2/40:7;
matP(321:361, 1)=5:1.2/40:6.2;
matP(361:401, 1)=6.2:0.8/40:7;
matP(1:11, 2)=0:1/10:1;
matP(11:21, 2)=0;
matP(21:31, 2)=1:-1/10:0;
matP(31:41, 2)=1:1.5/10:2.5;
matP(41:51, 2)=0:2.5/10:2.5;
matP(51:86, 2)=0;
matP(86:121, 2)=2.5:-2.5/35:0;
matP(121:156, 2)=0:2.5/35:2.5;
matP(156:196, 2)=2.5;
matP(196:241, 2)=0:2.5/45:2.5;
matP(241:281, 2)=0:1/40:1;
matP(281:321, 2)=0;
matP(321:361, 2)=2.5:-1.5/40:1;
matP(361:401, 2)=1:-1/40:0;
matP(402, :) = [0, 1];
matP(404, :) = [2, 0];
matP(403, :) = [2, 2.5];
matP(405, :) = [5, 2.5];
matP(406, :) = [5, 0];
matP(407, :) = [6.2, 1];
matP(408, :) = [7, 0];
matE=zeros(nEle,2);
Right = [402; 404; 404; 403; 403; ...
         406; 406; 405; 405; 405; ...
         407; 408; 407; 408;];
Left = [  1;   1; 402; 402; 404; ...
        404; 403; 404; 403; 406; ...
        406; 406; 405; 407];
for kk=1:1:nEle
    matE(kk,:)=[kk kk+1];
end
%오른쪽 불규칙
matE([10 20 30 40 50 ...
      85 120 155 195 240 ...
      280 320 360 400], :) = [[10  20   30  40 50 ...
                               85  120 155 195 240 ...
                               280 320 360 400]' Right];
%왼쪽 불규칙
matE([  2;  12;  22;  32;  42;
        52;  87; 122; 157; 197;
        242; 282; 322; 362]-1', :) = [Left, [  2;  12;  22;  32;  42;
                                              52;  87; 122; 157; 197;
                                             242; 282; 322; 362]];
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
FF(3*(403)-1, 1) = -20 * 10^3;
FF(3*(405)-1, 1) = -20 * 10^3;
FF(3*(407)-1, 1) = -20 * 10^3;
dd=zeros(3*nNod,1);
FullDoFs=[3*(usedNode)-2, 3*(usedNode)-1, 3*(usedNode)];
FixedDoFs=[3*(404)-2:3*(404),3*(406)-1];
FreeDoFs=setdiff(FullDoFs,FixedDoFs);
dd(FreeDoFs,1)=KK(FreeDoFs,FreeDoFs)\FF(FreeDoFs,1);
%% Visualization
scale=10^(4);
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
