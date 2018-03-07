% Rui Yan
% Math 151B
% 02/23/2018
%
%###################################################
% This program solves the ODE IVP
% dy/dt = -20*y+20*t^2+2*t, 0<=t<=1  
% y(t_0) = 1/3
%
% using (a) Euler's method
% (b) Runge-kutta method of orderfour
% (c) Adams fourth-order predictor-corrector method
% (d) Milne-Simpson predictor-corrector method which 
%     combines the explicit Milne’s method
%
% The function defining the ODE is defined in a separate Matlab function.
%###################################################
format long
a = 0.0;
b = 1.0;
y0 = 1.0/3.0;

f = @(t,y) (-20.0)*y + 20.0*t^2 + 2.0*t;

hlist = [0.2, 0.125, 0.1, 0.02];

t1 = linspace(a,b);


% Compare the exact solutions and the solutions 
for i = 1:length(hlist)
    h = hlist(i);
    y = findExact(t1);
    plot(t1,y,'-*', 'Color', 'm');
    hold on;
    Euler(f,a,b,h,y0);
    RungeKutta4(f,a,b,h,y0);
    Adams4(f,a,b,h,y0);
    MilneSimpson(f,a,b,h,y0);
    ylim([-10 10]);
    xlim([0 1]);
    title(['h = ' num2str(h) ' stepsize']);
    legend('exact','Euler','RK-4','Adams-4','Milne-Simpson','Location','northwest');
    filename = ['Figure_' num2str(i) '.fig'];
    savefig(filename);
    hold off;
end


% Euler's method:
function w1 = Euler(f,a,b,h,y0)
N = (b-a)/h;
t(1) = a;
w1(1) = y0;
error(1) = 0;

for n=1:N
    w1(n+1) = w1(n) + h*f(t(n),w1(n));
    t(n+1) = t(n) + h;
    y(n+1) = findExact(t(n+1));
    error(n+1) = abs(y(n+1) - w1(n+1));
end

plot(t, w1,'b--','Color','c');
fprintf('-------------Euler method when h = %.6s-------------\n', h);
fprintf('estimated value : %.6s | error : %.6s\n',[w1;error]);
end


% Runge-kutta method of order four:
function w2 = RungeKutta4(f,a,b,h,y0)
N = (b-a)/h;
t(1) = a;
w2(1) = y0;
error(1) = 0;

for n=1:N
    k1 = h*f(t(n), w2(n));
    k2 = h*f(t(n)+(h/2.0), w2(n)+(k1/2.0));
    k3 = h*f(t(n)+(h/2.0), w2(n)+(k2/2.0));
    k4 = h*f(t(n)+h, w2(n)+k3);
    
    w2(n+1) = w2(n) + (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
    t(n+1) = t(n) + h;
    y(n+1)=findExact(t(n+1));    
    error(n+1) = abs(y(n+1)-w2(n+1));
end
plot(t,w2,'-s','Color','b');
fprintf('-------------RK-4 method when h = %.6s-------------\n', h);
fprintf('estimated value : %.6s | error : %.6s\n',[w2;error]);
end


% Adams fourth-order predictor-corrector method:
function w3 = Adams4(f,a,b,h,y0)
N = (b-a)/h;
t(1) = a;
w3(1) = y0;
error(1) = 0;

for n=1:3
    k1 = h*f(t(n), w3(n));
    k2 = h*f(t(n)+(h/2.0), w3(n)+(k1/2.0));
    k3 = h*f(t(n)+(h/2.0), w3(n)+(k2/2.0));
    k4 = h*f(t(n)+h, w3(n)+k3);
    
    w3(n+1) = w3(n) + (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
    t(n+1) = t(n) + h;
    y(n+1)=findExact(t(n+1));
    error(n+1) = abs(y(n+1)-w3(n+1));
end
    
for n=4:N
    t(n+1) = t(n) + h;
    % Adams-Bashforth predictor step:
    w3_tmp = w3(n) + (h/24) * (55*f(t(n),w3(n)) - 59*f(t(n-1),w3(n-1)) ...
        + 37*f(t(n-2),w3(n-2)) - 9*f(t(n-3),w3(n-3)));
    
    % Adams-Moulton corrector step:
    w3(n+1) = w3(n) + (h/24) * (9*f(t(n+1),w3_tmp) + 19*f(t(n),w3(n)) ...
        - 5*f(t(n-1),w3(n-1)) + f(t(n-2),w3(n-2)));
    
    y(n+1) = findExact(t(n+1));
    error(n+1) = abs(y(n+1)-w3(n+1));
end

plot(t, w3,'-^','Color','r');
fprintf('-------------Adams-4 method when h = %.6s-------------\n', h);
fprintf('estimated value : %.6s | error : %.6s\n',[w3;error]);
end


% Milne-Simpson predictor-corrector method:
function w4 = MilneSimpson(f,a,b,h,y0)
N = (b - a)/h;
t(1) = a;
w4(1) = y0;
error(1) = 0;

for n=1:3
    k1 = h*f(t(n), w4(n));
    k2 = h*f(t(n)+(h/2.0), w4(n)+(k1/2.0));
    k3 = h*f(t(n)+(h/2.0), w4(n)+(k2/2.0));
    k4 = h*f(t(n)+h, w4(n)+k3);
    
    w4(n+1) = w4(n) + (1.0/6.0)*(k1+2.0*k2+2.0*k3+k4);
    t(n+1) = t(n) + h;
    y(n+1)=findExact(t(n+1));
    error(n+1) = abs(y(n+1)-w4(n+1));
end
    
for n=4:N
    t(n+1) = t(n) + h;
    
    % Milne predictor step:
    w4_tmp = w4(n-3) + ((4/3)*h)*(2*f(t(n),w4(n)) - f(t(n-1),w4(n-1))...
        + 2*f(t(n-2),w4(n-2)));
    
    % Simpson corrector step:
    w4(n+1) = w4(n-1) + (h/3)*(f(t(n+1),w4_tmp) + 4.0*f(t(n),w4(n))...
        + f(t(n-1), w4(n-1)));
    
    y(n+1)=findExact(t(n+1));
    error(n+1) = abs(y(n+1)-w4(n+1));
end

plot(t,w4,'-o','Color','g');
fprintf('-------------Milne-Simpson method when h = %.6s-------------\n', h);
fprintf('estimated value : %.6s | error : %.6s\n',[w4;error]);
end


function y = findExact(t)
y = t.^2 + (1.0/3.0)*exp(-20*t);
end