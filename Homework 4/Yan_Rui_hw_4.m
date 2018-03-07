% Rui Yan
% Math 151B
% 03/04/2018
%
%###################################################
% This program solves the system of nonlinear equations
% 		f1 = 15x(1)+x(2)^2-4x(3)-13
%		f2 = x(1)^2+10x(2)-x(3)-11
%		f3 = x(2)^3-25x(3)+22.
% with x(0) = (0,0,0)
% using 
% (a) Newton's method with the stopping criteria x(k)-x(k-1)<10^-6
% (b) Steepest descent method
% (c) Homotopy method with midpoint method, N = 10; 20; 50
% (d) Homotopy method with runge-kutta method of order 4, N = 10; 20; 50
%
% The function defining the ODE is defined in a separate Matlab function.
%###################################################

format long

x0 = [0;0;0];
tol = 1e-6;
fun = @F;
x_exact = fsolve(fun, x0);

N = 100;
%--Newton method--%
fprintf('----------Newton method----------\n')
tic
[x,k] = newton(x0,tol,N);
toc
fprintf('x: (%d, %d, %d)\n', x);
fprintf('The number of iterations is %d.\n', k);
fprintf('Error is %d.\n', norm(x-x_exact, Inf));
%
%--Steepest descent method--%
fprintf('\n----------Steepest descent method----------\n')
tic
[x,k] = descent(x0,tol,N);
toc
fprintf('x: (%d, %d, %d)\n', x);
fprintf('The number of iterations is %d.\n', k);
fprintf('Error is %d.\n', norm(x-x_exact, Inf));

N = [10;20;50];
%--Homotopy method with midpoint method--%
fprintf('\n----------Homotopy method with midpoint method----------\n')
for i = 1:length(N)
    fprintf('N = %d\n', N(i));
    tic
    x = homotopy_midpoint(x0,N(i));
    toc
    fprintf('x: (%d, %d, %d)\n', x);
    fprintf('Error is %s. \n',norm(x-x_exact, Inf));
end
%
%--Homotopy method with RK-4 method--%
fprintf('\n----------Homotopy method with RK-4 method----------\n')
for i = 1:length(N)
    fprintf('N = %d\n', N(i));
    tic
    x = homotopy_RK4(x0, N(i));
    toc
    fprintf('x: (%d, %d, %d)\n', x);
    fprintf('Error is %s. \n', norm(x-x_exact, Inf));
end

%%----------------Newton Method------------------%
function [x,k] = newton(x, tol, N) 
    k = 1;
    while k <= N
        y = -J(x)\F(x);
        x = x + y;
        if norm(y) < tol
            return;
        end
        k = k + 1;
    end
    if k > N
        disp('Maximum number of iterations exceeded');
    end
end

%%----------------Steepest Descent Method------------------%
function [x,k] = descent(x, tol, N)
    k = 1;
    while k <= N
        g1 = G(x);
        delg = 2*J(x)'*F(x);
        z0 = norm(delg);
        if z0 == 0
            disp('Zero gradient');
            disp('x is ');disp(x);
            disp('g1 is ');disp(g1);
            return;
        end
        
        % Unit vector in the steepest descent
        z = delg/z0;
        alpha3 = 1;
        g3 = G(x-alpha3*z);
        while g3 >= g1
          alpha3 = alpha3/2;
          g3 = G(x-alpha3*z);
          if alpha3 < tol/2
              disp('No likely improvement');
              disp('x is ');disp(x);
              disp('g1 is ');disp(g1);
              return;
          end
        end
        alpha2 = alpha3/2; 
        g2 = G(x-alpha2*z);
        
        h1 = (g2-g1)/alpha2;
        h2 = (g3-g2)/(alpha3-alpha2);
        h3 = (h2-h1)/alpha3;
        
        % Choose alpha0 so that g will be minimum in z direction
        alpha0 = 0.5*(alpha2-h1/h3);
        g0 = G(x - alpha0*z);
        if g0 < g3
            alpha = alpha0;g = g0;
        else
            alpha = a3;g = g3;
        end
        x = x - alpha*z;
        
        if abs(g-g1) < tol
            return;
        end
        k = k + 1;
    end
    if k > N
        disp('Maximum iterations exceeded');
    end
end

%%-----------Homotopy method with midpoint and RK-4 method----------%
function x = homotopy_midpoint(x,N)
    h = 1 / N;
    b = -h * F(x);
    for i = 1:N
        k1 = J(x)\b;
        k2 = J(x+ k1/2)\b;
        x = x + k2;
    end
end

function x = homotopy_RK4(x,N)
    h = 1 / N;
    b = -h * F(x);
    for i = 1:N
        k1 = J(x)\b;
        k2 = J(x + k1/2)\b;
        k3 = J(x + k2/2)\b;
        k4 = J(x + k3)\b;
        x = x + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
end

function f = F(x)
f = [15*x(1)+x(2).^2-4*x(3)-13;
    x(1).^2+10*x(2)-x(3)-11;
    x(2).^3-25*x(3)+22];
end

function j = J(x)
j = [15,2*x(2),-4;
    2*x(1),10,-1;
    0,3*x(2).^2,-25];
end

function g = G(x)
g = F(x)' * F(x);
end