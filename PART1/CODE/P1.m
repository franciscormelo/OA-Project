%% Optimization and Algorithms - Part 1
% Francisco Melo - 84053
%
% Rodrigo Rego - 89213
%
% Group Number - 28

%% Initial Conditions
% Saving the initial conditions
close all; clear all; clc
T=80;
pinitial=[0 5];
pfinal=[15 -15];
K=6;

W=[10 20 30 30 20 10; 10 10 10 0 0 -10];
tau=[10 25 30 40 50 60];
Umax=100;

A=[1 0 0.1 0; 0 1 0 0.1; 0 0 0.9 0;0 0 0 0.9];
B=[0 0; 0 0;0.1 0;0 0.1];
E=[1 0 0 0;0 1 0 0];

save('Initial_Constants');

%% Task1
close all; clear all; clc;

load Initial_Constants;

lambdas=[0.001 0.01 0.1 1 10 100 1000];

for l=1:length(lambdas)
    
    lambda=lambdas(l);
    
    cvx_begin quiet
    variable x(4,T+1);
    variable u(2,T+1);
    
    f1=0;
    for i=1:K
        f1=f1+square_pos(norm(E*x(:,tau(i)+1)-W(:,i)));
    end
    
    f2=0;
    for j=2:T
        f2=f2 + square_pos(norm(u(:,j)-u(:,j-1)));
    end
    
    f= f1+ lambda*f2;
    
    minimize(f);
    
    subject to
    x(:,1)==[pinitial';0;0];
    x(:,T+1)==[pfinal';0;0];
    
    for a=1:T
        norm(u(:,a)) <= Umax;
        x(:,a+1)==A*x(:,a)+ B*u(:,a);
    end
    
    cvx_end;
    
    figure();
    plot(pinitial(1),pinitial(2),'go',pfinal(1), pfinal(2),'go');
    
    hold on;
    plot(x(1,:), x(2,:),'bo','MarkerSize',3);
    hold on;
    
    for q=1:K
        plot(x(1,tau(q)+1), x(2,tau(q)+1),'mo','MarkerSize',10);
        hold on;
    end
    
    hold on;
    plot(W(1,:),W(2,:),'rs','MarkerSize',10);
    grid on;
    title( sprintf('lambda = %g', lambda),'interpreter','latex');
    xlabel('$x$','interpreter','latex');
    ylabel('$y$','interpreter','latex');
    
    
    time=linspace(0,T-1,T);
    
    figure();
    plot(time,u(1,1:end-1),time,u(2,1:end-1));
    grid on;
    title(sprintf('lambda = %g', lambda),'interpreter','latex');
    xlabel('$t$','interpreter','latex');
    ylabel('$u$','interpreter','latex');
    le=legend('$u_{1}(t)$', '$u_{2}(t)$');
    set(le,'interpreter','latex');
    
    counter=0;
    
    for t=2:T
        if norm(u(:,t)-u(:,t-1)) > 10^-6
            counter=counter+1;
        end
    end
    
    fprintf('N of times that optimal control changes (lambda= %g)- %g\n',...
        lambda,counter);
    
    mean_sum=0;
    for r=1:K
        mean_sum=mean_sum + norm(E*x(:,tau(r)+1)-W(:,r));
    end
    
    mean_deviation=(1/K)*mean_sum;
    fprintf('Mean Deviantion (lambda=%g) - %g\n\n', lambda, mean_deviation);
end
%% Task2
close all; clear all; clc;

load Initial_Constants;

lambdas=[0.001 0.01 0.1 1 10 100 1000];

for l=1:length(lambdas)
    
    lambda=lambdas(l);
    cvx_begin quiet
    variable x(4,T+1);
    variable u(2,T+1);
    
    f1=0;
    for i=1:K
        f1=f1+square_pos(norm(E*x(:,tau(i)+1)-W(:,i)));
    end
    
    f2=0;
    for j=2:T
        f2=f2 + norm(u(:,j)-u(:,j-1));
    end
    
    f= f1+ lambda*f2;
    
    minimize(f);
    
    subject to
    x(:,1)==[pinitial';0;0];
    x(:,T+1)==[pfinal';0;0];
    
    for a=1:T
        norm(u(:,a)) <= Umax;
        
        x(:,a+1)==A*x(:,a)+ B*u(:,a);
    end
   
    cvx_end;
    
    figure();
    plot(pinitial(1),pinitial(2),'go',pfinal(1), pfinal(2),'go');
    hold on;
    plot(x(1,:), x(2,:),'bo','MarkerSize',3);
    hold on;
    
    for q=1:K
        plot(x(1,tau(q)+1), x(2,tau(q)+1),'mo','MarkerSize',10);
        hold on;
    end
    
    hold on;
    plot(W(1,:),W(2,:),'rs','MarkerSize',10);
    grid on;
    title( sprintf('lambda = %g', lambda),'interpreter','latex');
    xlabel('$x$','interpreter','latex');
    ylabel('$y$','interpreter','latex');
    
    time=linspace(0,T-1,T);
    
    figure();
    plot(time,u(1,1:end-1),time,u(2,1:end-1));
    grid on;
    title(sprintf('lambda = %g', lambda),'interpreter','latex');
    xlabel('$t$','interpreter','latex');
    ylabel('$u$','interpreter','latex');
    le=legend('$u_{1}(t)$', '$u_{2}(t)$');
    set(le,'interpreter','latex');
    
    counter=0;
    
    for t=2:T
        if norm(u(:,t)-u(:,t-1)) > 10^-6
            counter=counter+1;
        end
    end
    
    fprintf("N of times that optimal control changes (lambda= %g) - %g\n",...
        lambda,counter);
    
    mean_sum=0;
    for r=1:K
        mean_sum=mean_sum + norm(E*x(:,tau(r)+1)-W(:,r));
    end
    
    mean_deviation=(1/K)*mean_sum;
    fprintf('Mean Deviantion (lambda=%g) - %g\n \n', lambda, mean_deviation);
end;

%% Task3
close all; clear all; clc;

load Initial_Constants;

lambdas=[0.001 0.01 0.1 1 10 100 1000];

for l=1:length(lambdas)
    
    lambda=lambdas(l);
    cvx_begin quiet
    variable x(4,T+1);
    variable u(2,T+1);
    
    f1=0;
    for i=1:K
        f1=f1+square_pos(norm(E*x(:,tau(i)+1)-W(:,i)));
    end
    
    f2=0;
    for j=2:T
        f2=f2 + norm(u(:,j)-u(:,j-1),1);%l1 norm
    end
    
    f= f1+ lambda*f2;
    
    minimize(f);
    
    subject to
    x(:,1)==[pinitial';0;0];
    x(:,T+1)==[pfinal';0;0];
    
    for a=1:T
        norm(u(:,a)) <= Umax;
        x(:,a+1)==A*x(:,a)+ B*u(:,a);
    end
    
    
    cvx_end;
    
    figure();
    plot(pinitial(1),pinitial(2),'go',pfinal(1), pfinal(2),'go');
    hold on;
    plot(x(1,:), x(2,:),'bo','MarkerSize',3);
    hold on;
    
    for q=1:K
        plot(x(1,tau(q)+1), x(2,tau(q)+1),'mo','MarkerSize',10);
        hold on;
    end
    
    hold on;
    plot(W(1,:),W(2,:),'rs','MarkerSize',10);
    grid on;
    title( sprintf('lambda = %g', lambda),'interpreter','latex');
    xlabel('$x$','interpreter','latex');
    ylabel('$y$','interpreter','latex');
    
    
    time=linspace(0,T-1,T);
    figure();
    
    plot(time,u(1,1:end-1),time,u(2,1:end-1));
    grid on;
    title(sprintf('lambda = %g', lambda),'interpreter','latex');
    xlabel('$t$','interpreter','latex');
    ylabel('$u$','interpreter','latex');
    le=legend('$u_{1}(t)$', '$u_{2}(t)$');
    set(le,'interpreter','latex');
    
    counter=0;
    
    for t=2:T
        if norm(u(:,t)-u(:,t-1)) > 10^-6
            counter=counter+1;
        end
    end
    
    fprintf("N of times that optimal control changes (lambda= %g) - %g\n",...
        lambda,counter);
    mean_sum=0;
    for r=1:K
        mean_sum=mean_sum + norm(E*x(:,tau(r)+1)-W(:,r));
    end
    
    mean_deviation=(1/K)*mean_sum;
    fprintf('Mean Deviantion (lambda=%g) - %g\n\n', lambda, mean_deviation);
end;
%% Task 6
close all; clear all; clc;

load Initial_Constants;
C=W; %center of the disks
R=2; % radius of the disks

lambdas=[0.001 0.01 0.1 1 10 100 1000];

for l=1:length(lambdas)
    
    lambda=lambdas(l);
    cvx_begin quiet
    variable x(4,T+1);
    variable u(2,T+1);
    
    f1=0;
    for i=1:K
        f1=f1+square_pos(norm(E*x(:,tau(i)+1)- C(:,i))-R);
    end
    
    f2=0;
    for j=2:T
        f2=f2 + norm(u(:,j)-u(:,j-1));
    end
    
    f= f1+ lambda*f2;
    
    minimize(f);
    
    subject to
    x(:,1)==[pinitial';0;0];
    x(:,T+1)==[pfinal';0;0];
    
    for a=1:T
        norm(u(:,a)) <= Umax;
        x(:,a+1)==A*x(:,a)+ B*u(:,a);
    end
    
    
    
    cvx_end;
    figure();
    
    plot(pinitial(1),pinitial(2),'go',pfinal(1), pfinal(2),'go');
    hold on;
    plot(x(1,:), x(2,:),'bo','MarkerSize',3);
    hold on;
    
    for q=1:K
        plot(x(1,tau(q)+1), x(2,tau(q)+1),'mo','MarkerSize',10);
        hold on;
    end
    
    hold on;
    vR = ones(6,1)*R;
    viscircles(C',vR);
    
    grid on;
    title( sprintf('lambda = %g', lambda),'interpreter','latex');
    xlabel('$x$','interpreter','latex');
    ylabel('$y$','interpreter','latex');
    
    time=linspace(0,T-1,T);
    
    figure();
    
    plot(time,u(1,1:end-1),time,u(2,1:end-1));
    grid on;
    title(sprintf('lambda = %g', lambda),'interpreter','latex');
    xlabel('$t$','interpreter','latex');
    ylabel('$u$','interpreter','latex');
    le=legend('$u_{1}(t)$', '$u_{2}(t)$');
    set(le,'interpreter','latex');
    
    counter=0;
    
    for t=2:T
        if norm(u(:,t)-u(:,t-1)) > 10^-6
            counter=counter+1;
        end
    end
    
    fprintf("N of times that optimal control changes (lambda= %g) - %g\n",...
        lambda,counter);
    
    mean_sum=0;
    for r=1:K
        mean_sum=mean_sum + norm(E*x(:,tau(r)+1)-W(:,r));
    end
    
    mean_deviation=(1/K)*mean_sum;
    fprintf('Mean Deviantion (lambda=%g) - %g\n\n', lambda, mean_deviation);
end;
%% Task 7
close all; clear all; clc;

load Initial_Constants;
Umax = 15;

cvx_begin
variable x(4,T+1);
variable u(2,T+1);

f= 0;
minimize(f);

subject to
x(:,1)==[pinitial';0;0];
x(:,T+1)==[pfinal';0;0];

for a=1:T
    norm(u(:,a)) <= Umax;
    x(:,a+1)==A*x(:,a)+ B*u(:,a);
end
for b=1:K
    E*x(:,tau(b)+1) == W(b);
end

cvx_end;
%% Task 9
close all; clear all; clc;

load Initial_Constants;

cvx_begin quiet
variable x(4,T+1);
variable u(2,T+1);

f1=0;
for i=1:K
    f1=f1+square_pos(norm(E*x(:,tau(i)+1)-W(:,i)));
end

f= f1;
minimize(f);

subject to
x(:,1)==[pinitial';0;0];
x(:,T+1)==[pfinal';0;0];

for a=1:T
    norm(u(:,a)) <= Umax;
    x(:,a+1)==A*x(:,a)+ B*u(:,a);
end

cvx_end;

figure();
plot(pinitial(1),pinitial(2),'go',pfinal(1), pfinal(2),'go');

hold on;
plot(x(1,:), x(2,:),'bo','MarkerSize',3);
hold on;

for q=1:K
    plot(x(1,tau(q)+1), x(2,tau(q)+1),'mo','MarkerSize',10);
    hold on;
end

hold on;
plot(W(1,:),W(2,:),'rs','MarkerSize',10);
grid on;
title( '$l_{2}^{2}$ formulation - Optimal Positions','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');


time=linspace(0,T-1,T);
figure();

plot(time,u(1,1:end-1),time,u(2,1:end-1));
grid on;
title('$l_{2}^{2}$ formulation - Optimal Control Signal','interpreter',...
    'latex');
xlabel('$t$','interpreter','latex');
ylabel('$u$','interpreter','latex');
le=legend('$u_{1}(t)$', '$u_{2}(t)$');
set(le,'interpreter','latex');

counter=0;
for t=1:K
    if norm(E*x(:,tau(t)+1)-W(:,t)) <= 10^-6
        counter=counter+1;
    end
end

fprintf("Number of waypoints captured by the robot : %g\n",counter);
%% Task 10
close all; clear all; clc;

load Initial_Constants;

Umax=15;

cvx_begin quiet
variable x(4,T+1);
variable u(2,T+1);

f1=0;
for i=1:K
    f1=f1+norm(E*x(:,tau(i)+1)-W(:,i));
end

f= f1;
minimize(f);

subject to
x(:,1)==[pinitial';0;0];
x(:,T+1)==[pfinal';0;0];

for a=1:T
    norm(u(:,a)) <= Umax;
    x(:,a+1)==A*x(:,a)+ B*u(:,a);
end

cvx_end;

figure();
plot(pinitial(1),pinitial(2),'go',pfinal(1), pfinal(2),'go');

hold on;
plot(x(1,:), x(2,:),'bo','MarkerSize',3);
hold on;

for q=1:K
    plot(x(1,tau(q)+1), x(2,tau(q)+1),'mo','MarkerSize',10);
    hold on;
end

hold on;
plot(W(1,:),W(2,:),'rs','MarkerSize',10);
grid on;
title( '$l_{2}$ formulation - Optimal Positions','interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');


time=linspace(0,T-1,T);
figure();

plot(time,u(1,1:end-1),time,u(2,1:end-1));
grid on;
title('$l_{2}$ formulation - Optimal Control Signal','interpreter','latex');
xlabel('$t$','interpreter','latex');
ylabel('$u$','interpreter','latex');
le=legend('$u_{1}(t)$', '$u_{2}(t)$');
set(le,'interpreter','latex');

counter=0;
for t=1:K
    if norm(E*x(:,tau(t)+1)-W(:,t)) <= 10^-6
        counter=counter+1;
    end
end

fprintf("Number of waypoints captured by the robot : %g\n",counter);
x0 = x;
u0 = u;
    
save('inicialization_task11', 'x0', 'u0');

%% Task 11
clear all; close all; clc;
load Initial_Constants
load inicialization_task11


Umax = 15;
e = 10^-6;
M = 10;

xm = x0;
um=u0;

for l = 1:M
    
    cvx_begin quiet
    variable x(4,T+1);
    variable u(2,T+1);

    
    f1 = 0;
    for i = 1:K
        weight(l,i) = 1/(norm(E*xm(:,tau(i)+1)-W(:,i))+e);
        f1 = f1 + weight(l,i)*norm(E*x(:,tau(i)+1)-W(:,i));
    end
    
    f = f1;
    minimize(f);
    
    subject to
    x(:,1) == [pinitial'; 0; 0];
    x(:,T+1) == [pfinal'; 0; 0];
    for a = 1:T
        norm(u(:,a)) <= Umax;
        x(:,a+1) == A*x(:,a)+B*u(:,a);
    end
    cvx_end;
    
    figure()
    plot(pinitial(1), pinitial(2), 'go', 'MarkerFaceColor', 'g');
    hold on; 
    plot(pfinal(1), pfinal(2), 'go', 'MarkerFaceColor', 'g'); hold on;
    for z = 1:K
        plot(xm(1, tau(z)+1), xm(2, tau(z)+1), 'mo', 'MarkerSize', 10, ...
            'LineWidth', 2); hold on;
    end
    plot(xm(1,:), xm(2,:), 'bo', 'MarkerSize', 3, 'LineWidth', 1.5); 
    hold on; 
    plot(W(1,:), W(2,:), 'rs', 'MarkerSize', 10, 'LineWidth', 2);
    xlabel('$x$', 'interpreter','latex');
    ylabel('$y$','interpreter','latex');
    title(sprintf('m = %d', l),'interpreter','latex');
    grid on;
    
    time = 0:1:T-1;
    
    figure()
    plot(time, um(1, 1:end-1)); hold on;
    plot(time, um(2, 1:end-1));
    xlabel('$t$','interpreter','latex');
    ylabel('$u$','interpreter','latex');
    le=legend('$u_{1}(t)$', '$u_{2}(t)$');
    set(le,'interpreter','latex');

    title(sprintf('m = %d', l), 'interpreter', 'latex');
    grid on;  
        
    cw = 0;
    for r = 1:K
        test = norm(E*xm(:,tau(r)+1)-W(:,r));
        if test <= 10^-6
            cw = cw + 1;
        end
    end
    
    fprintf('Number of Captured Waypoints: %g\n\n', cw);
    
    xm = x;
    um=u;
end