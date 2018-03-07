r = 0.2; K = 4000; a = 0; b = 50; y0 = 1000;
f = @(t,y) r*(1-y/K)*y;
hlist = [10,1,0.1];

% Plot the function y'' to determine M.
syms t;
y = (y0*K)./(y0+(K-y0)*exp(-1*r*t));
z(t) = diff(y,t,2);
t1 = linspace(a,b);
plot(t1,z(t1));
savefig('M.fig');

% Plot the exact solutions and the solutions 
% using Euler's method. Save these 3 figures.
for i=1:length(hlist)
    h = hlist(i);
    y1 = yE(t1);
    plot(t1,y1);
    hold on;
    wE(f,a,b,h,y0);
    filename = ['Figure ' num2str(i) '.fig'];
    savefig(filename);
    hold off;
end

% Approximate the solution using Euler's method
function w=wE(f,a,b,h,y0)
N=(b-a)/h;
errorMax=0; %errorMax stores the maximal actual error 
error=0; %error stores the error bound predicted by Euler's method
M=15;
L=0.2;
t(1)=a;
w(1)=y0;
for n=1:N
    t(n+1)=t(n)+h;
    w(n+1)=w(n)+h*f(t(n),w(n));
    y(n+1)=yE(t(n+1));
    if abs(y(n+1)-w(n+1)) > errorMax
        errorMax = abs(y(n+1)-w(n+1));
    end
    if (h*M)/(2*L)*(exp(L*(t(n+1)-a))-1) > error
        error = (h*M)/(2*L)*(exp(L*(t(n+1)-a))-1);
    end
end
plot(t,w);
title(['Euler Method using h = ' num2str(h) ' stepsize']);
disp(errorMax)
disp(error)
end

function y = yE(t)
r = 0.2; K = 4000; a = 0; b = 50; y0 = 1000;
y = (y0*K)./(y0+(K-y0)*exp(-1*r*t));
end