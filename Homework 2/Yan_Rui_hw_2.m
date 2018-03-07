a=1; b=2; y0=2;
f = @(t,y) 1+y/t;
hlist = [0.2,0.1,0.05];

t1 = linspace(a,b);
% Compare the exact solutions and the solutions 
% approximated by using Taylor's method of order two.
for i = 1:length(hlist)
    h = hlist(i);
    y1 = yE(t1);
    plot(t1,y1,'*');
    hold on;
    WE(f,a,b,h,y0);
    wE(f,a,b,h,y0);
    legend('exact','Midpoint method','Taylor method of order 2','Location','northwest');
    filename = ['Figure_' num2str(i) '.fig'];
    savefig(filename);
    hold off;
end

% Approximate the solution using Taylor method of order two
function w=wE(f,a,b,h,y0)
N=(b-a)/h;
error=0;
t(1)=a;
w(1)=y0;
tic;
for n=1:N
    t(n+1)=t(n)+h;
    w(n+1)=w(n)+h*f(t(n),w(n))+(h^2/2)*zE(t(n),w(n));
    y(n+1)=yE(t(n+1));
    if abs(y(n+1)-w(n+1)) > error
        error = abs(y(n+1)-w(n+1));
    end
end
toc;
plot(t,w,'b--','LineWidth',1.5);
title(['Taylor Method of Order 2 using h = ' num2str(h) ' stepsize']);
e = y(2.0)-w(2.0);
fprintf('Using Taylor method of order 2, for h= %s, the maximum error is: %s\n', h,error)
fprintf('Using Taylor method of order 2, for h= %s, the error for y(2) is: %s\n', h, e)
end

% Approximate the solution using Midpoint method
function W=WE(f,a,b,h,y0)
N=(b-a)/h;
error=0;
t(1)=a;
W(1)=y0;
tic;
for n=1:N
    t(n+1)=t(n)+h;
    W(n+1)=W(n)+h*f(t(n)+(h/2),W(n)+(h/2)*f(t(n),W(n)));
    y(n+1)=yE(t(n+1));
    if abs(y(n+1)-W(n+1)) > error
        error = abs(y(n+1)-W(n+1));
    end
end
toc;
plot(t,W,'r', 'LineWidth',2);
title(['Midpoint method using h = ' num2str(h) ' stepsize']);
e = y(2.0)-W(2.0);
fprintf('Using Midpoint method, for h= %s, the maximum error is: %s\n', h,error)
fprintf('Using Midpoint method, for h= %s, the error for y(2) is: %s\n', h, e)
end

function y = yE(t)
y = t.*log(t)+2*t;
end

function z = zE(t,y)
z = 1/t;
end
