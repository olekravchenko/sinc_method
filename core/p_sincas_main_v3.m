%% Программа по вычислению и решению краевой задачи на неравномерной адаптивной сетке методом КЧК коллокации
% Sinc-collocation method implementation
% Authors: Oleg Kravchenko, Dmitry Churikov
% Date: 19/03/15
% UPD1: 23/03/15 - nice
% UPD2: 25/03/15 - ole
	% use another function

close all;
clear all;
clc;

% make symbolic variable
syms x y w;

%% create symbolic parameters
% f = @(x)(sqrt(x)*(1-x)^(3/4));

%% Прямая замена (**1)
% y = x*exp(-(x)^2/10); % Манипуляция аргумента
% y = (atan(x/100))*(10/x); % Манипуляция аргумента
y = abs(x)*(atan(x/50))*(10); % Манипуляция аргумента
% y = (sqrt(abs(x)))/x; % Манипуляция аргумента

yinv = 100*tan(y/10); % y -> x

% f = @(x)(exp(-y.*x.^2/10).*cos(y.*x));
% f0 = (exp(-(  x)^2/10)*cos(  x));
% f  = (exp(-(  y)^2/10)*cos(  y));

% use another function
f0 = (exp(-(  x)^2/10)*cos(10*x));
f  = (exp(-(  y)^2/10)*cos(10*y));

% f  = (exp(-(y*x)^2/10)*cos(y*x));

f = eval(f);
% figure(10),ezplot(f);

FFTf = fourier(f,x,w);
% pretty(FFTf);

% figure(1), subplot(231), ezplot(f0,[-50,50, -1, 1]);
a1          =   20;
xdiskr      =   -a1:0.1:a1;
xinvdiskr   =   xdiskr;
fdiskr      =   subs(f,     xdiskr);
f0diskr     =   subs(f0,    xdiskr);
ydiskr      =   subs(y,     xdiskr);
yinvdiskr   =   subs(yinv,  xdiskr);

% figure(1), subplot(232), ezplot(f, [-a1, a1]);
figure(1), subplot(231), plot(xdiskr, f0diskr);
title('Исходный сигнал');
axis([-a1 a1 -0.5 1.1]);
figure(1), subplot(232), plot(xdiskr, fdiskr); hold on;
title('Модифицированный сигнал');
axis([-a1 a1 -0.5 1.1]);
figure(1), subplot(233), plot(xdiskr,ydiskr,'b'); hold on;
figure(1), subplot(233), plot(xdiskr,xdiskr,'k--');
% figure(1), subplot(233), plot(xdiskr,yinvdiskr+2,'r.');
title('Закон изменение аргумента');
axis([-a1 a1 -25 25]);

% figure(2), subplot(111), plot(xdiskr,yinvdiskr,'r.');


%% derivatives (производные)
d1f = simplify(diff(f, x));
d2f = simplify(diff(d1f, x));
% coefficients
mu      = sym(' 2*x');
nu      = sym('-2*x');
sigma   = d2f + mu*d1f + nu*f;
% checking
% disp(d2f + mu*d1f + nu*f - sigma)
% disp(simplify(d2f + mu*d1f + nu*f - sigma))

%% create initialization - граничные условия на отрезке [a;b] - нулевые
    a   = -20;
    b   = +20;
NTerms  = 50;
NPoints = 2*NTerms + 1;
dx      = (b-a) / (NPoints - 1);
xarr    = a:0.01*dx:b;
yarr    = subs(f,xarr);
xarr0   = a:dx:b;
yarr0   = subs(f,xarr0);
ind     = -NTerms:NTerms;

% plot discrete  plot
% figure(1), subplot(233), ezplot(FFTf);
% figure(1), subplot(234), plot(xarr, yarr);

%% computation 
% discrete matrixes
delta0 = eye(NPoints);
delta1 = delta0;
delta2 = delta0;
% формирование матриц для вычислений по схеме для чистого Синк-а
for j = ind
    for k = ind
        if j==k
            delta1(j + NTerms + 1, k + NTerms + 1) = 0;
            delta2(j + NTerms + 1, k + NTerms + 1) = -pi^2/3;
        else
            delta1(j + NTerms + 1, k + NTerms + 1) =    (-1)^(k-j)/(k-j);
            delta2(j + NTerms + 1, k + NTerms + 1) = -2*(-1)^(k-j)/(k-j)^2;
        end
    end
end

% initialize additional parameters
h = sqrt(0.5*pi^2/NTerms);
% Поставил такой же шаг, как и был
% не понял, почему изначально он никак не связан с исходной дискретизацией 
% (кроме кол-ва узлов):
% h = dx; 

    xTbl    =   h*meshgrid(ind,1);
   muTbl    =     subs(mu,    xTbl);
   nuTbl    =     subs(nu,    xTbl);
sigmaTbl    =     subs(sigma, xTbl);

%% Обратная замена (**2) - ????
xinvTbl = 10*atan(xTbl/50).*abs(xTbl);

% assembling of global matrix MS
MS = delta0;
dx = h;
for j = ind
    for k = ind
        MS(k + NTerms + 1, j + NTerms + 1) = ...
            (1/dx^2)*delta2(j + NTerms + 1, k + NTerms + 1) + ...
            (1/dx)*muTbl(k + NTerms +1)*delta1(j + NTerms + 1, k + NTerms + 1) + ...
            nuTbl(k + NTerms + 1)*delta0(j + NTerms + 1, k + NTerms + 1);
    end
end

% obtain coefficients
w = MS \ sigmaTbl';

% y approximate
SM = zeros(NPoints, NPoints);
yapproximate = zeros(1, NPoints);
for j = ind
    SM(j + NTerms + 1,:) = sinc_basis('WKS',j, dx, xTbl);
    yapproximate = yapproximate + w(j + NTerms + 1)*SM(j + NTerms + 1,:);
end

% result plot
figure(1), subplot(234),
    plot(xarr, yarr,  xTbl, yapproximate,'r.-')
    axis([-a1 a1 -0.5 1.1]);
    title('Восстановление сигнала');

figure(1), subplot(235), 
%     gp = plot(xdiskr, f0diskr); hold on;
%     set(gp,'Color','b','LineWidth',2,'LineStyle','-');
%     title('Исходный сигнал'); axis([-a1 a1 -0.5 1.1]);    
%                     set(gca,'FontSize', 20);
%                     set(gca,'LineWidth',1.5);    
    
figure(1), subplot(235),
    gp = plot(xinvTbl, yapproximate); hold on;
    set(gp,'Color','r','LineWidth',1.5,'LineStyle','-');
    axis([-a1 a1 -0.5 1.1]);
    title('Восстановление сигнала');  
    
figure(1), subplot(236),
plot(xTbl, xinvTbl,'b')
title('Обратная замена');
% plot(xTbl, abs((yapproximate-yarr0)),'b')
% title('Ошибка восстановления');

% figure(1), subplot(233),
% plot(xTbl, sinc_basis('WKS',0, dx, xTbl),'b.')