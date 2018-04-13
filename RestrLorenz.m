clc,clear
load('Lorenz.mat')%data
A0 = [-10,10,0,0,0,0,0,0,0;28,-1,0,0,0,-1,0,0,0;0,0,-2,1,0,0,0,0,0];%Lorenz actual parameters
N = 10^7;%The amount of data
h=10^(-4);
tau = 0.03;%Sampling step
t = floor(tau/h);

L = zeros(9,N);%base
L(1,:) = DataMat(:,1);
L(2,:) = DataMat(:,2);
L(3,:) = DataMat(:,3);
L(4,:) = DataMat(:,1).*DataMat(:,2);
L(5,:) = DataMat(:,2).*DataMat(:,3);
L(6,:) = DataMat(:,1).*DataMat(:,3);
L(7,:) = DataMat(:,1).^2;
L(8,:) = DataMat(:,2).^2;
L(9,:) = DataMat(:,3).^2;

St = L(:,1+t:N)*L(:,1:N-t)'/(N-t);
S0 = L*L'/N;

A = logm(St*(S0^-1))/(t*h);%Reconstruction results
D = -(A*S0+S0*A');%Noise

A0plus = reshape(A0(:,1:min(length(S0),10)),1,[]);
Aplus = reshape(A(1:3,1:min(length(S0),10)),1,[]);


FS = 20;
FZ = 15;
LW = 1.5;
figure
for i = 1:length(Aplus)
    plot(A0plus(i),Aplus(i),'ro','linesmoothing','on','LineWidth',LW)
    hold on
end

Endpoint1 = floor(min([min(A0plus),min(Aplus)])*2);
Endpoint2 = ceil(max([max(A0plus),max(Aplus)])*1.2);
plot([Endpoint1,Endpoint2],[Endpoint1,Endpoint2],'linesmoothing','on','LineWidth',LW)
axis([Endpoint1,Endpoint2,Endpoint1,Endpoint2])
set(gca,'FontSize',FZ);
xlabel('$A_{ij}$','FontSize',FS,'Interpreter','LaTeX')    
ylabel('$\it{\hat{A}}_{ij}(n=10)$','FontSize',FS,'Interpreter','LaTeX')
title(['$\tau$=',num2str(t*h)],'FontSize',FS,'Interpreter','LaTeX')

figure
for i = 1:3
    for j = 1:3
        if i == j
            plot(100,D(i,j),'ro','linesmoothing','on','LineWidth',LW)
            hold on
        else
            plot(0,D(i,j),'ro','linesmoothing','on','LineWidth',LW)
            hold on
        end
    end
end
plot([-10,110],[-10,110],'linesmoothing','on','LineWidth',LW)
axis([-10,110,-10,110])
set(gca,'FontSize',FZ);
title(['$\tau$=',num2str(t*h)],'FontSize',FS,'Interpreter','LaTeX')
xlabel('${D}_{ij}$','FontSize',FS,'Interpreter','LaTeX')    
ylabel('$\it{\hat{D}}_{ij}$','FontSize',FS,'Interpreter','LaTeX')

