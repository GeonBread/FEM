%% Q4
clearvars;
clc; close;
%% Given Parameters
E=200*10^(9);
G=77.5*10^(9);
v=E/(2*G)-1;
t=0.12;
P=4000;
%% matP & matE
matP=zeros(729,2);
for kk=1:1:81
    matP(1 + 9*(kk-1), :) = [0 * (0.1/8), (1/80)*(kk-1)];
    matP(2 + 9*(kk-1), :) = [1 * (0.1/8), (1/80)*(kk-1)];
    matP(3 + 9*(kk-1), :) = [2 * (0.1/8), (1/80)*(kk-1)];
    matP(4 + 9*(kk-1), :) = [3 * (0.1/8), (1/80)*(kk-1)];
    matP(5 + 9*(kk-1), :) = [4 * (0.1/8), (1/80)*(kk-1)];
    matP(6 + 9*(kk-1), :) = [5 * (0.1/8), (1/80)*(kk-1)];
    matP(7 + 9*(kk-1), :) = [6 * (0.1/8), (1/80)*(kk-1)];
    matP(8 + 9*(kk-1), :) = [7 * (0.1/8), (1/80)*(kk-1)];
    matP(9 + 9*(kk-1), :) = [8 * (0.1/8), (1/80)*(kk-1)];
end
matE=zeros(640,4);
for kk=1:1:80
    matE(1 + 8*(kk - 1), :) = [ 1  2 11 10] + 9*(kk-1);
    matE(2 + 8*(kk - 1), :) = [ 2  3 12 11] + 9*(kk-1);
    matE(3 + 8*(kk - 1), :) = [ 3  4 13 12] + 9*(kk-1);
    matE(4 + 8*(kk - 1), :) = [ 4  5 14 13] + 9*(kk-1);
    matE(5 + 8*(kk - 1), :) = [ 5  6 15 14] + 9*(kk-1);
    matE(6 + 8*(kk - 1), :) = [ 6  7 16 15] + 9*(kk-1);
    matE(7 + 8*(kk - 1), :) = [ 7  8 17 16] + 9*(kk-1);
    matE(8 + 8*(kk - 1), :) = [ 8  9 18 17] + 9*(kk-1);
end
nNod=size(matP,1);
nEle=size(matE,1);
%% Element Stiffness Matrix
b=(0.1/8)/2;
h=(1/80)/2;
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
KK=zeros(2*nNod,2*nNod);
for kk=1:1:nEle
    % kth Element
    p=matE(kk,1);
    q=matE(kk,2);
    r=matE(kk,3);
    s=matE(kk,4);
    eDoFs=[2*p-1:2*p,2*q-1:2*q,2*r-1:2*r,2*s-1:2*s];
    KK(eDoFs,eDoFs)=KK(eDoFs,eDoFs)+kLocal;
end
%% Boundary Conditions
FF=zeros(2*nNod,1);
FF(2*721-1,1)=P;
dd=zeros(2*nNod,1);
FullDoFs=2*1-1:2*nNod;
FixedDoFs=2*1-1:2*9;
FreeDoFs=setdiff(FullDoFs,FixedDoFs);
dd(FreeDoFs,1)=KK(FreeDoFs,FreeDoFs)\FF(FreeDoFs,1);
%% Visualization
figure(1);
scale=5*10^2;
% for kk=1:1:nNod
%     scatter(matP(kk,1),matP(kk,2),'r','LineWidth',1.1);
%     hold on;
% end
for kk=1:1:nEle
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
    patch('Faces',f,'Vertices',xDeformed,'EdgeColor','b','FaceColor','none','LineWidth',1.1);
    hold on;
    % Stress Plot
    patch('Faces',f,'Vertices',xDeformed,'FaceVertexCData',Sxx,'FaceColor','flat','LineWidth',1.1);
    hold on;
end
colorbar;
axis equal;
axis padded;
