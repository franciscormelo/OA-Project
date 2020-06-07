%% Optimization and Algorithms - Part 2
% Francisco Melo - 84053
%
% Rodrigo Rego - 89213
%
% Group Number - 28
%% Task 2
clear all; close all; clc;
load data1.mat

K = length(Y);
[l c]=size(X);

a = 1;
gama = 10^-4;
beta = 0.5;

eps = 10^-6;

k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));
A = (1/K)*[X; v_ones];

x(:,1)=[-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    d(:,k) = -g(:,k);
    
    z = 0;
    a = 1;
    while z == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        c = x(:,k);
        member2 = (1/K)*sum(log(1+exp(c'*[X;v_ones]))-Y.*(c'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            z = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    k = k + 1;
end
toc  
x(:,end)


figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Gradient Method)','interpreter','latex')
grid on;


x1=X(1,:)';
x2=X(2,:)';
class=Y';

figure();
plot(x1(class==0),x2(class==0),'ro','linewidth',1);hold on;
plot(x1(class==1),x2(class==1),'bo','linewidth',1);hold on;
xlabel('$x_{1}$','interpreter','latex');
ylabel('$x_{2}$','interpreter','latex');


t=-2:0.001:6;


x2=(x(3,end) - x(1,end)*t)/x(2,end);
plot(t,x2,'g--','linewidth',2);
legend('Class 0','Class 1','Hyperplane');

%% Task 3
clear all; close all; clc;
load data2.mat


K = length(Y);
[l c]=size(X);

a = 1;
gama = 10^-4;
beta = 0.5;


eps = 10^-6;


k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));
A = (1/K)*[X; v_ones];

x(:,1)=[-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    d(:,k) = -g(:,k);
    
    z = 0;
    a = 1;
    while z == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        c = x(:,k);
        member2 = (1/K)*sum(log(1+exp(c'*[X;v_ones]))-Y.*(c'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            z = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    k = k + 1;
end
toc 

x(:,end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Gradient Method)','interpreter','latex')
grid on;


x1=X(1,:)';
x2=X(2,:)';
class=Y';

figure();
plot(x1(class==0),x2(class==0),'ro','linewidth',1);hold on;
plot(x1(class==1),x2(class==1),'bo','linewidth',1);hold on;
xlabel('$x_{1}$','interpreter','latex');
ylabel('$x_{2}$','interpreter','latex');


t=-3.5:0.001:5.5;


x2=(x(3,end) - x(1,end)*t)/x(2,end);
plot(t,x2,'g--','linewidth',2);
legend('Class 0','Class 1','Hyperplane');
%% Task 4 data3.mat
clear all; close all; clc;
load data3.mat

K = length(Y);
[l c]=size(X);

a = 1;
gama = 10^-4;
beta = 0.5;



eps = 10^-6;



k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));
A = (1/K)*[X; v_ones];

x(:,1)=[-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    d(:,k) = -g(:,k);
    
    z = 0;
    a = 1;
    while z == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        c = x(:,k);
        member2 = (1/K)*sum(log(1+exp(c'*[X;v_ones]))-Y.*(c'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            z = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    k = k + 1;
end
toc 

x(:,end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Gradient Method)','interpreter','latex')
grid on;


x1=X(1,:)';
x2=X(2,:)';
class=Y';

figure();
plot(x1(class==0),x2(class==0),'ro','linewidth',1);hold on;
plot(x1(class==1),x2(class==1),'bo','linewidth',1);hold on;
xlabel('$x_{1}$','interpreter','latex');
ylabel('$x_{2}$','interpreter','latex');

%% Task 4 data4.mat
clear all; close all; clc;
load task4.mat


K = length(Y);
[l c]=size(X);

a = 1;
gama = 10^-4;
beta = 0.5;

eps = 10^-6;

k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));
A = (1/K)*[X; v_ones];

x(:,1)=[-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    
    gnorm(k) = norm(g(:,k));         
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    d(:,k) = -g(:,k);
    
    z = 0;
    a = 1;
    while z == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        c = x(:,k);
        member2 = (1/K)*sum(log(1+exp(c'*[X;v_ones]))-Y.*(c'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            z = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    k = k + 1;
end
toc 

x(:,end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Gradient Method)','interpreter','latex')
grid on;


x1=X(1,:)';
x2=X(2,:)';
class=Y';

figure();
plot(x1(class==0),x2(class==0),'ro','linewidth',1);hold on;
plot(x1(class==1),x2(class==1),'bo','linewidth',1);hold on;
xlabel('$x_{1}$','interpreter','latex');
ylabel('$x_{2}$','interpreter','latex');
%% Task 6
clear all; close all; clc;
load data1.mat


K = length(Y);
[l c]=size(X);

a = 1;
gama = 10^-4;
beta = 0.5;


eps = 10^-6;


k = 1;
criteria = 0;

v_ones = -1*ones(1, length(X));

va=[X; -v_ones];

A1= (1/K)*[X; v_ones];


x(:,1)=[-ones(1,l)';0];

tic
while criteria == 0
    g(:,k) = sum(A1.*((exp(x(:,k)'*[X;v_ones]))./(1 + (exp(x(:,k)'*[X;v_ones])))  - Y), 2);
    
    gnorm(k) = norm(g(:,k));  
    
    if gnorm(k) < eps
        criteria = 1;
        break;
    end
    
    H =va *((1/K)*diag(((exp(x(:,k)'*[X;v_ones]))./((1 + (exp(x(:,k)'*[X;v_ones]))).^2))))*va';
    
    
    d(:,k) = -inv(H)*g(:,k);
    
    
    
    z = 0;
    a = 1;
    while z == 0
        
        b = x(:,k) + a*d(:,k);
        member1 = (1/K)*sum(log(1+exp(b'*[X;v_ones]))-Y.*(b'*[X;v_ones]));
        c = x(:,k);
        member2 = (1/K)*sum(log(1+exp(c'*[X;v_ones]))-Y.*(c'*[X;v_ones])) + ...
                   gama*g(:,k)'*a*d(:,k);
        
        if member1 < member2
            z = 1;
            break;
        end
        a = beta*a;
    end
    
    x(:,k+1) = x(:,k) + a*d(:,k);
    k = k + 1;
    
    pause();
  
end
toc 

x(:,end)

figure()
semilogy(gnorm, 'linewidth', 1);
xlabel('$k$', 'interpreter','latex');
title('$|| \nabla f(s_{k},r_{k}) ||$ (Gradient Method)','interpreter','latex')
grid on;


x1=X(1,:)';
x2=X(2,:)';
class=Y';

figure();
plot(x1(class==0),x2(class==0),'ro','linewidth',1);hold on;
plot(x1(class==1),x2(class==1),'bo','linewidth',1);hold on;
xlabel('$x_{1}$','interpreter','latex');
ylabel('$x_{2}$','interpreter','latex');


t=-3.5:0.001:5.5;


x2=(x(3,end) - x(1,end)*t)/x(2,end);
plot(t,x2,'g--','linewidth',2);
legend('Class 0','Class 1','Hyperplane');
%% Task 8

close all; clear all; clc;

load lmdata1.mat;