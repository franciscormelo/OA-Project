%% Task 2
clear all; close all; clc;
load data1.mat

K = 150;

% Variables s and r
syms r;
s = sym('s', [1 2]);


% Function to minimize in respect to s and r
%F = (1/150)*sum(log10(1+exp(s*X(:,k)-r))-Y(k)*(s*X(:,k)-r));
F = 0;
for k = 1:K
    F = F + (1/K)*(log10(1+exp(s*X(:,k)-r))-Y(k)*(s*X(:,k)-r));
end

a = 1;
gama = 10^-4;
beta = 0.5;

sp = [-1 -1];
rp = 0;
eps = 10^-6;

x(:,1) = [-1; -1; 0];

k = 0;
criteria = 0;

ds1 = diff(F, s(1,1));
ds2 = diff(F, s(1,2));
dr = diff(F, r);

tic
while criteria == 0
    
    g(:,k+1) = [eval(subs(ds1, [s(1,1), s(1,2), r], [x(1, k+1), x(2, k+1), x(3, k+1)])); ...
                eval(subs(ds2, [s(1,1), s(1,2), r], [x(1, k+1), x(2, k+1), x(3, k+1)])); ...
                eval(subs(dr, [s(1,1), s(1,2), r], [x(1, k+1), x(2, k+1), x(3, k+1)]))];
    
    if norm(g(:,k+1)) < eps
        criteria = 1;
        break;
    end
    
    d(:,k+1) = -g(:,k+1);
    
    c = 0;
%     while c == 0
%         subsTerm = a*d(:,k+1);
%         member1 = eval(subs(F, [s(1,1), s(1,2), r], [x(1, end) + subsTerm(1,1), x(2, end) + subsTerm(2,1), x(3, end) + subsTerm(3,1)]));
%         subsTerm2 = eval(subs(F, [s(1,1), s(1,2), r], [x(1, end), x(2, end), x(3, end)]));
%         member2 = subsTerm2 + gama*g(:,k+1)'*subsTerm;
%         
%         if member1 < member2
%             c = 1;
%             break;
%         end
%         a = beta*a;
%     end
    
while eval(subs(F,[s(1) s(2) r],(x(:,k+1)+(a*d(:,k+1)))')) >= eval(subs(F,[s(1) s(2) r],[x(1,k+1) x(2,k+1) x(3,k+1)])) + (gama * g(:,k+1)'*(a*d(:,k+1)))
    a = beta * a;
end

    x(:,k+2) = x(:,k+1) + a*d(:,k+1);
    k = k + 1;
    [k a]
end
toc  
x(:, end)