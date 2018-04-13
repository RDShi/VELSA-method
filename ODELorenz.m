function ODELorenz()
clc,clear
h=10^(-4);%Step length
n=10^(7);%The amount of data
D = 100;%Noise intensity
DataMat = zeros(n,3);

for i = 1:n-1
    DataMat(i+1,:) = DataMat(i,:) + h * Lorenz(DataMat(i,:)')'+sqrt(D*h)*randn(1,3);%Euler method
    if mod(i,n/100)==0
        disp(i/(n/100));
    end
end

save Lorenz.mat DataMat

FS = 20;
FZ = 15;
LW = 1.5;
figure
plot(DataMat(:,1),DataMat(:,2)) %Phase Diagrams
xlabel('x','FontSize',FS,'Interpreter','LaTeX')    
ylabel('y','FontSize',FS,'Interpreter','LaTeX')
title('Phase~Diagrams','FontSize',FS,'Interpreter','LaTeX')

function dx = Lorenz(x) %Lorenz equation

dx=zeros(3,1);
o = 10;
% b = 8/3;
b = 2;
p = 28;                                  


dx(1) = o*(x(2)-x(1));
dx(2) = x(1)*(p - x(3)) - x(2);
dx(3) = x(1)*x(2) - b*x(3);