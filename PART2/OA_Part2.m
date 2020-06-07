%% OA Part 1 - Logistic Regression & Network Localization
% 
% Francisco Melo
%
% Rodrigo Rego
%
%% Task 2
clear all; close all; clc;
load data1.mat

K = length(Y);
[l, c] = size(X);

a = 1;
gama = 10^-4;
beta = 0.5;
eps = 10^-6;

k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));

A = (1/K)*[X; v_ones];

x(:,1) = [-ones(1,l)'; 0];

tic
while criteria == 0
    g(:,k) = sum(A.*((exp(x(:,k)'*[X;v_ones]))./(1 + ...
                                      (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    d(:,k) = -g(:,k);
    
    c = 0;
    a = 1;
    while c == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        z = x(:,k);
        member2 = (1/K)*sum(log(1+exp(z'*[X;v_ones]))-Y.*(z'*[X;v_ones]))...
                                                   + gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            c = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    k = k + 1;
end
toc  
s = [x(1, end); x(2, end)]
r = x(3, end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Gradient Method)',...
                                                    'interpreter','latex')
grid on;

x1 = X(1,:)';
x2 = X(2,:)';
class = Y';

tx = -2:0.001:6;

figure()
plot(x1(class==0), x2(class==0), 'ro', 'linewidth', 1); hold on;
plot(x1(class==1), x2(class==1), 'bo', 'linewidth', 1); hold on;
plot(tx, ((-s(1)*tx)+r)/s(2), 'g--', 'linewidth', 2);
xlabel('$x_{1}$', 'Interpreter', 'latex');
ylabel('$x_{2}$', 'Interpreter', 'latex');
legend('Class 0', 'Class 1', 'Hyperplane', 'location', 'northeast');

%% Task 3
clear all; close all; clc;
load data2.mat

K = length(Y);
[l, c] = size(X);

a = 1;
gama = 10^-4;
beta = 0.5;
eps = 10^-6;

k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));

A = (1/K)*[X; v_ones];

x(:,1) = [-ones(1,l)'; 0];

tic
while criteria == 0
    g(:,k) = sum(A.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    d(:,k) = -g(:,k);
    
    c = 0;
    a = 1;
    while c == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        z = x(:,k);
        member2 = (1/K)*sum(log(1+exp(z'*[X;v_ones]))-Y.*(z'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            c = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    k = k + 1;
end
toc  
s = [x(1, end); x(2, end)]
r = x(3, end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Gradient Method)','interpreter','latex')
grid on;

x1 = X(1,:)';
x2 = X(2,:)';
class = Y';

tx = -3.5:0.001:5.5;

figure()
plot(x1(class==0), x2(class==0), 'ro', 'linewidth', 1); hold on;
plot(x1(class==1), x2(class==1), 'bo', 'linewidth', 1); hold on;
plot(tx, ((-s(1)*tx)+r)/s(2), 'g--', 'linewidth', 2);
xlabel('$x_{1}$', 'Interpreter', 'latex');
ylabel('$x_{2}$', 'Interpreter', 'latex');
legend('Class 0', 'Class 1', 'Hyperplane', 'location', 'northeast');

%% Task 4

%% Task 4 - data3.mat
clear all; close all; clc;
load data3.mat

K = length(Y);
[l, c] = size(X);


a = 1;
gama = 10^-4;
beta = 0.5;
eps = 10^-6;

k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));

A = (1/K)*[X; v_ones];

x(:,1) = [-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    d(:,k) = -g(:,k);
    
    c = 0;
    a = 1;
    while c == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        z = x(:,k);
        member2 = (1/K)*sum(log(1+exp(z'*[X;v_ones]))-Y.*(z'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            c = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    k = k + 1;
end
toc  
s = [x(1, end); x(2, end)]
r = x(3, end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Gradient Method)','interpreter','latex')
grid on;

%% Task 4 - dataset4.m
clear all; close all; clc;
load data4.mat

K = length(Y);
[l, c] = size(X);


a = 1;
gama = 10^-4;
beta = 0.5;
eps = 10^-6;

k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));

A = (1/K)*[X; v_ones];

x(:,1) = [-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    d(:,k) = -g(:,k);
    
    c = 0;
    a = 1;
    while c == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        z = x(:,k);
        member2 = (1/K)*sum(log(1+exp(z'*[X;v_ones]))-Y.*(z'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            c = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    k = k + 1;
    
end
toc  
s = [x(1, end); x(2, end)]
r = x(3, end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Gradient Method)','interpreter','latex')
grid on;

%% Task 6
clear all; close all; clc;
load data1.mat

K = length(Y);
[l, c] = size(X);


a = 1;
gama = 10^-4;
beta = 0.5;
eps = 10^-6;

k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));

A1 = (1/K)*[X; v_ones];

va = [X; v_ones];

x(:,1) = [-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A1.*((exp(x(:,k)'*[X;v_ones]))./(1 + ...
                                    (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    grad = g(:,k);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    H = va*(1/K)*diag(((exp(x(:,k)'*[X;v_ones]))./((1 + ...
                                      (exp(x(:,k)'*[X;v_ones]))).^2)))*va';
    d(:,k) = -(H\grad);

    c = 0;
    a = 1;
    while c == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        z = x(:,k);
        member2 = (1/K)*sum(log(1+exp(z'*[X;v_ones]))-Y.*(z'*[X;v_ones]))... 
                                                   + gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            c = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    step(k) = a;
    k = k + 1;
end
toc  
s = [x(1, end); x(2, end)]
r = x(3, end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Newton Method)',...
                                                    'interpreter','latex');
grid on;

figure()
stem(1:1:k-1, step,'.','MarkerSize', 25, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$\alpha_k$ (Newton Method)','interpreter','latex');

%% Task 6 - data2.mat
clear all; close all; clc;
load data2.mat

K = length(Y);
[l, c] = size(X);


a = 1;
gama = 10^-4;
beta = 0.5;
eps = 10^-6;

k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));

A1 = (1/K)*[X; v_ones];

va = [X; v_ones];

x(:,1) = [-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A1.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    grad = g(:,k);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    H = va*(1/K)*diag(((exp(x(:,k)'*[X;v_ones]))./((1 + (exp(x(:,k)'*[X;v_ones]))).^2)))*va';
    d(:,k) = -(H\grad);

    c = 0;
    a = 1;
    while c == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        z = x(:,k);
        member2 = (1/K)*sum(log(1+exp(z'*[X;v_ones]))-Y.*(z'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            c = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    step(k) = a;
    k = k + 1;
end
toc  
s = [x(1, end); x(2, end)]
r = x(3, end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Newton Method)','interpreter','latex');
grid on;

figure()
stem(1:1:k-1, step,'.','MarkerSize', 25, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$\alpha_k$ (Newton Method)','interpreter','latex');

%% Task 6 - data3.mat
clear all; close all; clc;
load data3.mat

K = length(Y);
[l, c] = size(X);


a = 1;
gama = 10^-4;
beta = 0.5;
eps = 10^-6;

k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));

A1 = (1/K)*[X; v_ones];

va = [X; v_ones];

x(:,1) = [-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A1.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    grad = g(:,k);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    H = va*(1/K)*diag(((exp(x(:,k)'*[X;v_ones]))./((1 + (exp(x(:,k)'*[X;v_ones]))).^2)))*va';
    d(:,k) = -(H\grad);

    c = 0;
    a = 1;
    while c == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        z = x(:,k);
        member2 = (1/K)*sum(log(1+exp(z'*[X;v_ones]))-Y.*(z'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            c = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    step(k) = a;
    k = k + 1;
end
toc  
s = [x(1, end); x(2, end)]
r = x(3, end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlim([1 k]);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Newton Method)','interpreter','latex');
grid on;

figure()
stem(1:1:k-1, step,'.','MarkerSize', 25, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$\alpha_k$ (Newton Method)','interpreter','latex');

%% Task 6 - data4.mat
clear all; close all; clc;
load data4.mat

K = length(Y);
[l, c] = size(X);


a = 1;
gama = 10^-4;
beta = 0.5;
eps = 10^-6;

k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));

A1 = (1/K)*[X; v_ones];

va = [X; v_ones];

x(:,1) = [-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A1.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    grad = g(:,k);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    H = va*(1/K)*diag(((exp(x(:,k)'*[X;v_ones]))./((1 + (exp(x(:,k)'*[X;v_ones]))).^2)))*va';
    d(:,k) = -(H\grad);

    c = 0;
    a = 1;
    while c == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        z = x(:,k);
        member2 = (1/K)*sum(log(1+exp(z'*[X;v_ones]))-Y.*(z'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            c = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    step(k) = a;
    k = k + 1;
end
toc  
s = [x(1, end); x(2, end)]
r = x(3, end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlim([1 k]);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Newton Method)','interpreter','latex');
grid on;

figure()
stem(1:1:k-1, step,'.','MarkerSize', 25, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$\alpha_k$ (Newton Method)','interpreter','latex');

%% Task 8 - lmdata1.mat
clear all; close all; clc;
load lmdata1.mat

x = xinit;

eps = 10^-6;
lambda(1) = 1;

k = 1;
criteria = 0;
tic
while criteria == 0
    previous_f = 0;
    next_f = 0;
    gf1 = zeros(length(x),1);
    gf2 = zeros(length(x),1);
    
    gp1 = zeros(length(x), length(y));
    f1 = 0;
    f2 = 0;
    fp = [];
    for i = 1:length(y)
        ind1 = 2*iA(i,2)-1;
        ind2 = 2*iA(i,2);
        
        aux1 = [gf1(ind1); gf1(ind2)] + (-2)*(norm(A(:,iA(i,1)) -...
            [x(ind1); x(ind2)])-y(i)) * (A(:,iA(i,1)) - ...
            [x(ind1); x(ind2)])/norm(A(:,iA(i,1)) - [x(ind1); x(ind2)]);
        
        gf1(ind1) = aux1(1);
        gf1(ind2) = aux1(2);
        
        gpaux1 = -(A(:,iA(i,1)) - [x(ind1); x(ind2)])/norm(A(:,iA(i,1))...
                                                    - [x(ind1); x(ind2)]);
        
        gp1(ind1,i) = gpaux1(1);
        gp1(ind2,i) = gpaux1(2);
        
        f1 = f1 + (norm(A(:,iA(i,1)) - [x(ind1); x(ind2)])-y(i))^2;
        fp = [fp; (norm(A(:,iA(i,1)) - [x(ind1); x(ind2)])-y(i))];
    end
    
    gp2 = zeros(length(x), length(z));
    for i = 1:length(z)
        ind1 = 2*iS(i,2)-1;
        ind2 = 2*iS(i,2);
        
        idx1 = 2*iS(i,1)-1;
        idx2 = 2*iS(i,1);
        
        aux2 = [gf2(ind1); gf2(ind2); gf2(idx1); gf2(idx2)] +...
            (-2)*(norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)])-z(i)) * ...
            ([x(idx1); x(idx2); x(ind1); x(ind2)] - ...
            [x(ind1); x(ind2);x(idx1); x(idx2)])/norm([x(idx1); x(idx2)]...
            - [x(ind1); x(ind2)]);
        
        gf2(ind1) = aux2(1);
        gf2(ind2) = aux2(2);
        gf2(idx1) = aux2(3);
        gf2(idx2) = aux2(4);
        
        gpaux2 = -([x(idx1); x(idx2); x(ind1); x(ind2)] - [x(ind1);...
            x(ind2);x(idx1); x(idx2)])/norm([x(idx1); x(idx2)] - ...
                                                       [x(ind1); x(ind2)]);
        
        gp2(ind1,i) = gpaux2(1);
        gp2(ind2,i) = gpaux2(2);
        gp2(idx1,i) = gpaux2(3);
        gp2(idx2,i) = gpaux2(4);
        
        f2 = f2 + (norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)])-z(i))^2;
        fp = [fp; (norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)])-z(i))];
    end
    
    previous_f = f1 + f2;
    
    gnorm(k) = norm(gf1+gf2);
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    gp = [gp1 gp2];
    
    Al = [gp'; sqrt(lambda(k))*eye(16)];
    b = [gp'*x - fp; sqrt(lambda(k))*x];
    
    new_x = Al\b;
    
    nf1 = 0;
    nf2 = 0;
    for i = 1:length(y)
        ind1 = 2*iA(i,2)-1;
        ind2 = 2*iA(i,2);
        
        nf1 = nf1 + (norm(A(:,iA(i,1)) - [new_x(ind1); new_x(ind2)])...
                                                                -y(i))^2;
    end
    for i = 1:length(z)
        ind1 = 2*iS(i,2)-1;
        ind2 = 2*iS(i,2);
        
        idx1 = 2*iS(i,1)-1;
        idx2 = 2*iS(i,1);
        
        nf2 = nf2 + (norm([new_x(idx1); new_x(idx2)] - ...
                                       [new_x(ind1); new_x(ind2)])-z(i))^2;
    end
    
    next_f = nf1 + nf2;
    
    if next_f < previous_f
        lambda(k+1) = 0.7*lambda(k);
        x = new_x;
    else
        lambda(k+1) = 2*lambda(k);
    end
    
    k = k + 1;
end
toc

x1 = x(1:2:end)';
x2 = x(2:2:end)';
x = [x1;x2];

xi1 = xinit(1:2:end)';
xi2 = xinit(2:2:end)';
xi = [xi1;xi2];

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
xlim([0 20]);
title('$|| \nabla f(x_{k}) ||$ (LM method)','interpreter','latex')
grid on;

a = A(:,iA(:,1))';
sa = S(:,iA(:,2))';

s1 = S(:,iS(:,1))';
s2 = S(:,iS(:,2))';

figure()
plot(S(1,:), S(2,:), 'b.', 'MarkerSize', 15); hold on;
plot(A(1,:), A(2,:), 'rs', 'LineWidth', 2.5, 'MarkerSize', 10); hold on;
plot(x(1,:), x(2,:), 'bo', 'MarkerSize', 10, 'LineWidth', 1.5); hold on;
plot(xi(1,:), xi(2,:), 'b*', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot([sa(:,1) a(:,1)]',[sa(:,2) a(:,2)]', 'm--'); hold on;
plot([s1(:,1) s2(:,1)]',[s1(:,2) s2(:,2)]', 'm--');
axis([-15 15 -12.25 10.5]);
grid on;

%% Task 8 - lmdata2.mat
clear all; close all; clc;
load lmdata2.mat


for j = 1:25
    x = randi([-15 15], 16, 1);
    xinit = x;
    
    eps = 10^-6;
    lambda(1) = 1;
    
    k = 1;
    criteria = 0;
    gnorm = [];
    tic
    while criteria == 0
        previous_f = 0;
        next_f = 0;
        gf1 = zeros(length(x),1);
        gf2 = zeros(length(x),1);
        
        gp1 = zeros(length(x), length(y));
        f1 = 0;
        f2 = 0;
        fp = [];
        for i = 1:length(y)
            ind1 = 2*iA(i,2)-1;
            ind2 = 2*iA(i,2);
            
            aux1 = [gf1(ind1); gf1(ind2)] + (-2)*(norm(A(:,iA(i,1)) - ...
                [x(ind1); x(ind2)])-y(i)) * ...
                (A(:,iA(i,1)) - [x(ind1); x(ind2)])/norm(A(:,iA(i,1)) -...
                                                    [x(ind1); x(ind2)]);
            
                                                gf1(ind1) = aux1(1);
            gf1(ind2) = aux1(2);
            
            gpaux1 = -(A(:,iA(i,1)) - [x(ind1); x(ind2)])/...
                                norm(A(:,iA(i,1)) - [x(ind1); x(ind2)]);
            
            gp1(ind1,i) = gpaux1(1);
            gp1(ind2,i) = gpaux1(2);
            
            f1 = f1 + (norm(A(:,iA(i,1)) - [x(ind1); x(ind2)])-y(i))^2;
            fp = [fp; (norm(A(:,iA(i,1)) - [x(ind1); x(ind2)])-y(i))];
        end
        
        gp2 = zeros(length(x), length(z));
        for i = 1:length(z)
            ind1 = 2*iS(i,2)-1;
            ind2 = 2*iS(i,2);
            
            idx1 = 2*iS(i,1)-1;
            idx2 = 2*iS(i,1);
            
            aux2 = [gf2(ind1); gf2(ind2); gf2(idx1); gf2(idx2)] +...
                (-2)*(norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)])-z(i))...
                *([x(idx1); x(idx2); x(ind1); x(ind2)] - ...
                [x(ind1); x(ind2);x(idx1); x(idx2)])/...
                norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)]);
            
            gf2(ind1) = aux2(1);
            gf2(ind2) = aux2(2);
            gf2(idx1) = aux2(3);
            gf2(idx2) = aux2(4);
            
            gpaux2 = -([x(idx1); x(idx2); x(ind1); x(ind2)] - ...
                [x(ind1); x(ind2);x(idx1); x(idx2)])/...
                norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)]);
            
            gp2(ind1,i) = gpaux2(1);
            gp2(ind2,i) = gpaux2(2);
            gp2(idx1,i) = gpaux2(3);
            gp2(idx2,i) = gpaux2(4);
            
            f2 = f2 + (norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)])-z(i))^2;
            fp = [fp; (norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)])-z(i))];
        end
        
        previous_f = f1 + f2;
        
        gnorm(k) = norm(gf1+gf2);
        if gnorm(k) < eps
            criteria = 1;
            break;
        end
        
        gp = [gp1 gp2];
        
        Al = [gp'; sqrt(lambda(k))*eye(16)];
        b = [gp'*x - fp; sqrt(lambda(k))*x];
        
        new_x = Al\b;
        
        nf1 = 0;
        nf2 = 0;
        for i = 1:length(y)
            ind1 = 2*iA(i,2)-1;
            ind2 = 2*iA(i,2);
            
            nf1 = nf1 + (norm(A(:,iA(i,1)) - [new_x(ind1);...
                                                    new_x(ind2)])-y(i))^2;
        end
        for i = 1:length(z)
            ind1 = 2*iS(i,2)-1;
            ind2 = 2*iS(i,2);
            
            idx1 = 2*iS(i,1)-1;
            idx2 = 2*iS(i,1);
            
            nf2 = nf2 + (norm([new_x(idx1); new_x(idx2)] - ...
                                       [new_x(ind1); new_x(ind2)])-z(i))^2;
        end
        
        next_f = nf1 + nf2;
        
        if next_f < previous_f
            lambda(k+1) = 0.7*lambda(k);
            x = new_x;
        else
            lambda(k+1) = 2*lambda(k);
        end
        
        k = k + 1;
    end
    eval_f(j) = previous_f;
    num_iter(j) = k;
    eval_x(:,j) = x;
    eval_xinit(:,j) = xinit;
    eval_gnorm(j).norm = gnorm;
end
toc

[min_cost, idx] = min(eval_f);
[min_iter, idxi] = min(num_iter);

if eval_f(idxi) == min_cost
    idx = idxi;
end

x = eval_x(:,idx);
xinit = eval_xinit(:,idx);

x1 = x(1:2:end)';
x2 = x(2:2:end)';
x = [x1;x2];

xi1 = xinit(1:2:end)';
xi2 = xinit(2:2:end)';
xi = [xi1;xi2];

figure()
semilogy(eval_gnorm(idx).norm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
xlim([0 num_iter(idx)+2]);
title('$|| \nabla f(x_{k}) ||$ (LM method)','interpreter','latex')
grid on;

a = A(:,iA(:,1))';
sa = x(:,iA(:,2))';

s1 = x(:,iS(:,1))';
s2 = x(:,iS(:,2))';

figure()
plot(A(1,:), A(2,:), 'rs', 'LineWidth', 2.5, 'MarkerSize', 10); hold on;
plot(x(1,:), x(2,:), 'bo', 'MarkerSize', 10, 'LineWidth', 1.5); hold on;
plot(xi(1,:), xi(2,:), 'b*', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
plot([sa(:,1) a(:,1)]',[sa(:,2) a(:,2)]', 'm--'); hold on;
plot([s1(:,1) s2(:,1)]',[s1(:,2) s2(:,2)]', 'm--');
axis([-15 15 -15 15]);
grid on;