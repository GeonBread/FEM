%% 왼쪽에서 P로 밀고 가운데에 연직 하중 2P
clearvars;
clc; close;
%% Given Parameters
E=200*10^(9);
G=77.5*10^(9);
v=E/(2*G)-1;
t=0.12;
P=-4000;
nEle = 232;
nNod = 295;
%% matP & matE
matP=zeros(nNod,2);
for k=1:1:23
    matP(1 + 5*(k-1), :) = [0       (1/22)*(k-1)];
    matP(2 + 5*(k-1), :) = [0.025*1 (1/22)*(k-1)];
    matP(3 + 5*(k-1), :) = [0.025*2 (1/22)*(k-1)];
    matP(4 + 5*(k-1), :) = [0.025*3 (1/22)*(k-1)];
    matP(5 + 5*(k-1), :) = [0.025*4 (1/22)*(k-1)];
end
for k=1:1:14
    matP(116 + 5*(k-1), :) = [0.1+(0.8/14)*(k) 1-(1/22)*4];
    matP(117 + 5*(k-1), :) = [0.1+(0.8/14)*(k) 1-(1/22)*3];
    matP(118 + 5*(k-1), :) = [0.1+(0.8/14)*(k) 1-(1/22)*2];
    matP(119 + 5*(k-1), :) = [0.1+(0.8/14)*(k) 1-(1/22)*1];
    matP(120 + 5*(k-1), :) = [0.1+(0.8/14)*(k) 1-(1/22)*0];
end
for k=15:1:18
    matP(186 + 5*(k-15), :) = [0.9+(0.1/4)*(k-14) 1-(1/22)*4];
    matP(187 + 5*(k-15), :) = [0.9+(0.1/4)*(k-14) 1-(1/22)*3];
    matP(188 + 5*(k-15), :) = [0.9+(0.1/4)*(k-14) 1-(1/22)*2];
    matP(189 + 5*(k-15), :) = [0.9+(0.1/4)*(k-14) 1-(1/22)*1];
    matP(190 + 5*(k-15), :) = [0.9+(0.1/4)*(k-14) 1-(1/22)*0];
end
for k = 1:1:18
    matP(206 + 5*(k-1), :) = [0.9 + (0.1/4)*0 1-(1/22)*(5+(k-1))];
    matP(207 + 5*(k-1), :) = [0.9 + (0.1/4)*1 1-(1/22)*(5+(k-1))];
    matP(208 + 5*(k-1), :) = [0.9 + (0.1/4)*2 1-(1/22)*(5+(k-1))];
    matP(209 + 5*(k-1), :) = [0.9 + (0.1/4)*3 1-(1/22)*(5+(k-1))];
    matP(210 + 5*(k-1), :) = [0.9 + (0.1/4)*4 1-(1/22)*(5+(k-1))];
end
matE=zeros(nEle,4);
for k=1:1:22
    matE(1 + 4*(k-1), :) = [1 2  7 6] + 5*(k-1)*ones(1, 4);
    matE(2 + 4*(k-1), :) = [2 3  8 7] + 5*(k-1)*ones(1, 4);
    matE(3 + 4*(k-1), :) = [3 4  9 8] + 5*(k-1)*ones(1, 4);
    matE(4 + 4*(k-1), :) = [4 5 10 9] + 5*(k-1)*ones(1, 4);
end
matE(89, :) = [ 95 116 117 100];
matE(90, :) = [100 117 118 105];
matE(91, :) = [105 118 119 110];
matE(92, :) = [110 119 120 115];
for k = 1:1:17
    matE(93 + 4*(k-1), :) = [116 121 122 117] + 5*(k-1)*ones(1, 4);
    matE(94 + 4*(k-1), :) = [117 122 123 118] + 5*(k-1)*ones(1, 4);
    matE(95 + 4*(k-1), :) = [118 123 124 119] + 5*(k-1)*ones(1, 4);
    matE(96 + 4*(k-1), :) = [119 124 125 120] + 5*(k-1)*ones(1, 4);
end
% matE(161, :) = [181 206 207 186];
% matE(162, :) = [186 207 208 191];
% matE(163, :) = [191 208 209 196];
% matE(164, :) = [196 209 210 201];
matE(161:164, :) = [206 207 186 181;
                    207 208 191 186;
                    208 209 196 191;
                    209 210 201 196];
for k=1:1:17 
    % matE(165 + 4*(k-1), :) = [206 211 212 207] + 5*(k-1)*ones(1, 4);
    % matE(166 + 4*(k-1), :) = [207 212 213 208] + 5*(k-1)*ones(1, 4);
    % matE(167 + 4*(k-1), :) = [208 213 214 209] + 5*(k-1)*ones(1, 4);
    % matE(168 + 4*(k-1), :) = [209 214 215 210] + 5*(k-1)*ones(1, 4);
    matE(165 + 4*(k-1), :) = [ 211 212 207 206 ] + 5*(k-1)*ones(1, 4);
    matE(166 + 4*(k-1), :) = [ 212 213 208 207 ] + 5*(k-1)*ones(1, 4);
    matE(167 + 4*(k-1), :) = [ 213 214 209 208 ] + 5*(k-1)*ones(1, 4);
    matE(168 + 4*(k-1), :) = [ 214 215 210 209 ] + 5*(k-1)*ones(1, 4);
end
%% Element Stiffness Matrix
%% 1
KK=zeros(2*nNod,2*nNod);
b=(0.1/4)/2;
h=(1/22)/2;
syms xx yy
N1= (1/(4*b*h))*(xx-b)*(yy-h);
N2=-(1/(4*b*h))*(xx+b)*(yy-h);
N3= (1/(4*b*h))*(xx+b)*(yy+h);
N4=-(1/(4*b*h))*(xx-b)*(yy+h);
dN1dx=diff(N1,xx);
dN2dx=diff(N2,xx);
dN3dx=diff(N3,xx);
dN4dx=diff(N4,xx);
dN1dy=diff(N1,yy);
dN2dy=diff(N2,yy);
dN3dy=diff(N3,yy);
dN4dy=diff(N4,yy);
matB(xx,yy)=[dN1dx     0 dN2dx     0 dN3dx     0 dN4dx     0;
                 0 dN1dy     0 dN2dy     0 dN3dy     0 dN4dy;
             dN1dy dN1dx dN2dy dN2dx dN3dy dN3dx dN4dy dN4dx];
matC=E/(1-v^2)*[1 v 0;
                v 1 0;
                0 0 (1-v)/2];
kLocal=double(t*int(int(transpose(matB)*matC*matB,xx,-b,b),yy,-h,h));
%% Assembly
for kk=1:1:88
    % kth Element
    p=matE(kk,1);
    q=matE(kk,2);
    r=matE(kk,3);
    s=matE(kk,4);
    eDoFs=[2*p-1:2*p,2*q-1:2*q,2*r-1:2*r,2*s-1:2*s];
    KK(eDoFs,eDoFs)=KK(eDoFs,eDoFs)+kLocal;
end
%% 2
b=(0.8/14)/2;
h=(1/22)/2;
syms xx yy
N1= (1/(4*b*h))*(xx-b)*(yy-h);
N2=-(1/(4*b*h))*(xx+b)*(yy-h);
N3= (1/(4*b*h))*(xx+b)*(yy+h);
N4=-(1/(4*b*h))*(xx-b)*(yy+h);
dN1dx=diff(N1,xx);
dN2dx=diff(N2,xx);
dN3dx=diff(N3,xx);
dN4dx=diff(N4,xx);
dN1dy=diff(N1,yy);
dN2dy=diff(N2,yy);
dN3dy=diff(N3,yy);
dN4dy=diff(N4,yy);
matB(xx,yy)=[dN1dx     0 dN2dx     0 dN3dx     0 dN4dx     0;
                 0 dN1dy     0 dN2dy     0 dN3dy     0 dN4dy;
             dN1dy dN1dx dN2dy dN2dx dN3dy dN3dx dN4dy dN4dx];
matC=E/(1-v^2)*[1 v 0;
                v 1 0;
                0 0 (1-v)/2];
kLocal=double(t*int(int(transpose(matB)*matC*matB,xx,-b,b),yy,-h,h));
%% Assembly
for kk=89:1:144
    % kth Element
    p=matE(kk,1);
    q=matE(kk,2);
    r=matE(kk,3);
    s=matE(kk,4);
    eDoFs=[2*p-1:2*p,2*q-1:2*q,2*r-1:2*r,2*s-1:2*s];
    KK(eDoFs,eDoFs)=KK(eDoFs,eDoFs)+kLocal;
end
%% 3
b=(0.1/4)/2;
h=(1/22)/2;
syms xx yy
N1= (1/(4*b*h))*(xx-b)*(yy-h);
N2=-(1/(4*b*h))*(xx+b)*(yy-h);
N3= (1/(4*b*h))*(xx+b)*(yy+h);
N4=-(1/(4*b*h))*(xx-b)*(yy+h);
dN1dx=diff(N1,xx);
dN2dx=diff(N2,xx);
dN3dx=diff(N3,xx);
dN4dx=diff(N4,xx);
dN1dy=diff(N1,yy);
dN2dy=diff(N2,yy);
dN3dy=diff(N3,yy);
dN4dy=diff(N4,yy);
matB(xx,yy)=[dN1dx     0 dN2dx     0 dN3dx     0 dN4dx     0;
                 0 dN1dy     0 dN2dy     0 dN3dy     0 dN4dy;
             dN1dy dN1dx dN2dy dN2dx dN3dy dN3dx dN4dy dN4dx];
matC=E/(1-v^2)*[1 v 0;
                v 1 0;
                0 0 (1-v)/2];
kLocal=double(t*int(int(transpose(matB)*matC*matB,xx,-b,b),yy,-h,h));
%% Assembly
for kk=145:1:232
    % kth Element
    p=matE(kk,1);
    q=matE(kk,2);
    r=matE(kk,3);
    s=matE(kk,4);
    eDoFs=[2*p-1:2*p,2*q-1:2*q,2*r-1:2*r,2*s-1:2*s];
    KK(eDoFs,eDoFs)=KK(eDoFs,eDoFs)+kLocal;
end
%%%%%%%%%%%%
%% Boundary Conditions
FF=zeros(2*nNod,1);
FF(2*101-1,1)=-P;
FF(2*150,1)=2*P;
% FF(2*106-1,1)=-P;
dd=zeros(2*nNod,1);
FullDoFs=2*1-1:2*nNod;
FixedDoFs=[2*1-1:2*5, 2*291-1:2*(295)];
FreeDoFs=setdiff(FullDoFs,FixedDoFs);
dd(FreeDoFs,1)=KK(FreeDoFs,FreeDoFs)\FF(FreeDoFs,1);
%% Visualization
figure(1);
scale=5*10^3;
% for kk=1:1:nNod
%     scatter(matP(kk,1),matP(kk,2),'r','LineWidth',1.1);
%     hold on;
% end
for kk=1:1:nEle
    if mod(kk, 12) == 0
        a = 1;
    end
    p=matE(kk,1);
    q=matE(kk,2);
    r=matE(kk,3);
    s=matE(kk,4);
    x1=matP(p,1);
    y1=matP(p,2);
    x2=matP(q,1);
    y2=matP(q,2);
    x3=matP(r,1);
    y3=matP(r,2);
    x4=matP(s,1);
    y4=matP(s,2);
    xInitial=[x1 y1;
              x2 y2;
              x3 y3;
              x4 y4];
    f = [1 2 3 4];
    eDoFs=[2*p-1:2*p,2*q-1:2*q,2*r-1:2*r,2*s-1:2*s];
    de=dd(eDoFs,1);
    Se=double(matC*matB(0,0)*de);
    Sxx=Se(2);
    Svm=sqrt((Se(1))^2+(Se(2))^2-(Se(1))*(Se(2))+3*(Se(3))^2);
    xDeformed=[x1+scale*de(1) y1+scale*de(2);
               x2+scale*de(3) y2+scale*de(4);
               x3+scale*de(5) y3+scale*de(6);
               x4+scale*de(7) y4+scale*de(8)];
    % Initial Shape
    patch('Faces',f,'Vertices',xInitial,'EdgeColor','k','FaceColor','none','LineWidth',1.1);
    hold on;
    % Deformed Shape
    % patch('Faces',f,'Vertices',xDeformed,'EdgeColor','b','FaceColor','none','LineWidth',1.1);
    % hold on;
    % Stress Plot
    patch('Faces',f,'Vertices',xDeformed,'FaceVertexCData',Svm,'FaceColor','flat','LineWidth',1.1);
    hold on;
end
colorbar;
axis equal;
axis padded;
