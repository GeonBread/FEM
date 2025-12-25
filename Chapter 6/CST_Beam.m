%% CST (Constant Strain Triangle)
clearvars;
clc; close;
%% Given Parameters
E=200*10^(9);
G=77.5*10^(9);
v=E/(2*G)-1;
t=0.12;
P=-4000;
%% matP & matE
matP=zeros(205,2);
for k=1:1:41
    matP(1 + 5*(k-1), :) = [0 + (1/40)*(k-1), 0 * (0.1/4)];
    matP(2 + 5*(k-1), :) = [0 + (1/40)*(k-1), 1 * (0.1/4)];
    matP(3 + 5*(k-1), :) = [0 + (1/40)*(k-1), 2 * (0.1/4)];
    matP(4 + 5*(k-1), :) = [0 + (1/40)*(k-1), 3 * (0.1/4)];
    matP(5 + 5*(k-1), :) = [0 + (1/40)*(k-1), 4 * (0.1/4)]; 
end
matE=zeros(312,3);
for k=1:1:40
    matE(1 + 8*(k-1), : ) = [1 6  7] + 5*(k-1);
    matE(2 + 8*(k-1), : ) = [1 7  2] + 5*(k-1);
    matE(3 + 8*(k-1), : ) = [2 7  8] + 5*(k-1);
    matE(4 + 8*(k-1), : ) = [2 8  3] + 5*(k-1);
    matE(5 + 8*(k-1), : ) = [3 8  4] + 5*(k-1);
    matE(6 + 8*(k-1), : ) = [4 8  9] + 5*(k-1);
    matE(7 + 8*(k-1), : ) = [4 9  5] + 5*(k-1);
    matE(8 + 8*(k-1), : ) = [5 9 10] + 5*(k-1);
end
nEle=size(matE,1);
nNod=size(matP,1);
matC=E/(1-v^2)*[1 v       0;
                v 1       0;
                0 0 (1-v)/2];
%% Assembly
KK=zeros(2*nNod,2*nNod);
for kk=1:1:nEle
    % kth Element
    p=matE(kk,1);
    q=matE(kk,2);
    r=matE(kk,3);
    x1=matP(p,1);
    y1=matP(p,2);
    x2=matP(q,1);
    y2=matP(q,2);
    x3=matP(r,1);
    y3=matP(r,2);
    matXY=[1 x1 y1;
           1 x2 y2;
           1 x3 y3];
    A=(1/2)*det(matXY);
    matABR=(2*A)*inv(matXY);
    matB=(1/(2*A))*[matABR(2,1)           0 matABR(2,2)           0 matABR(2,3)           0;
                              0 matABR(3,1)           0 matABR(3,2)           0 matABR(3,3);
                    matABR(3,1) matABR(2,1) matABR(3,2) matABR(2,2) matABR(3,3) matABR(2,3)];
    kLocal=t*A*transpose(matB)*matC*matB;
    eDoFs=[2*p-1:2*p,2*q-1:2*q,2*r-1:2*r];
    KK(eDoFs,eDoFs)=KK(eDoFs,eDoFs)+kLocal;
end
%% Boundary Conditions
FF=zeros(2*nNod,1);
FF(2*205,1)=P;
dd=zeros(2*nNod,1);
FullDoFs=2*1-1:2*nNod;
FixedDoFs=2*(1)-1:2*(5);
FreeDoFs=setdiff(FullDoFs,FixedDoFs);
dd(FreeDoFs,1)=KK(FreeDoFs,FreeDoFs)\FF(FreeDoFs,1);
%% Visualization
figure(1);
scale=10^(3);
for kk=1:1:nEle
    % kth element
    p=matE(kk,1);
    q=matE(kk,2);
    r=matE(kk,3);
    x1=matP(p,1);
    y1=matP(p,2);
    x2=matP(q,1);
    y2=matP(q,2);
    x3=matP(r,1);
    y3=matP(r,2);
    matXY=[1 x1 y1;
           1 x2 y2;
           1 x3 y3];
    A=(1/2)*det(matXY);
    matABR=(2*A)*inv(matXY);
    matB=(1/(2*A))*[matABR(2,1)           0 matABR(2,2)           0 matABR(2,3)           0;
                              0 matABR(3,1)           0 matABR(3,2)           0 matABR(3,3);
                    matABR(3,1) matABR(2,1) matABR(3,2) matABR(2,2) matABR(3,3) matABR(2,3)];
    xIntial=[x1 y1; 
             x2 y2; 
             x3 y3];
    eDoFs=[2*p-1:2*p,2*q-1:2*q,2*r-1:2*r];
    de=dd(eDoFs,1);
    Se=matC*matB*de;
    Svm=sqrt((Se(1))^2+(Se(2))^2-(Se(1))*(Se(2))+3*(Se(3))^2);
    xDeformed=[x1+scale*de(1) y1+scale*de(2); 
               x2+scale*de(3) y2+scale*de(4); 
               x3+scale*de(5) y3+scale*de(6)];
    f = [1 2 3];
    % Initial Shape
    patch('Faces',f,'Vertices',xIntial,'EdgeColor','k','FaceColor','none','LineWidth',1.1);
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
