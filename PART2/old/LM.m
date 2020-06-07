%% Task 8
clear all; close all; clc;
load lmdata1.mat

x = xinit;

eps = 10^-6;
lambda(1) = 1;

k = 1;
criteria = 0;
while criteria == 0
    previous_f = 0;
    next_f = 0;
    gf1 = zeros(length(x),1);
    gf2 = zeros(length(x),1);
    
    gp1 = zeros(length(x), length(y));
    
    f1 = 0;
    f2 = 0;
    fp = [];
    
    aux1=0;
    gpaux1=0;
    aux2=0;
    gpaux2=0;
    
    
    for i = 1:length(y)
        ind1 = 2*iA(i,2)-1;
        ind2 = 2*iA(i,2);
        
        aux1 = [gf1(ind1); gf1(ind2)] + (-2)*(norm(A(:,iA(i,1)) - [x(ind1); x(ind2)])-y(i)) * ...
            (A(:,iA(i,1)) - [x(ind1); x(ind2)])/norm(A(:,iA(i,1)) - [x(ind1); x(ind2)]);
        gf1(ind1) = aux1(1);
        gf1(ind2) = aux1(2);
        
        gpaux1 = -(A(:,iA(i,1)) - [x(ind1); x(ind2)])/norm(A(:,iA(i,1)) - [x(ind1); x(ind2)]);
        
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
        
        aux2 = [gf2(ind1); gf2(ind2); gf2(idx1); gf2(idx2)] + (-2)*(norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)])-z(i)) * ...
            ([x(idx1); x(idx2); x(ind1); x(ind2)] - [x(ind1); x(ind2);x(idx1); x(idx2)])/norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)]);
        gf2(ind1) = aux2(1);
        gf2(ind2) = aux2(2);
        gf2(idx1) = aux2(3);
        gf2(idx2) = aux2(4);
        
        gpaux2 = -([x(idx1); x(idx2); x(ind1); x(ind2)] - [x(ind1); x(ind2);x(idx1); x(idx2)])/norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)]);
        
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
        
        nf1 = nf1 + (norm(A(:,iA(i,1)) - [new_x(ind1); new_x(ind2)])-y(i))^2;
    end
    for i = 1:length(z)
        ind1 = 2*iS(i,2)-1;
        ind2 = 2*iS(i,2);
        
        idx1 = 2*iS(i,1)-1;
        idx2 = 2*iS(i,1);
        
        nf2 = nf2 + (norm([new_x(idx1); new_x(idx2)] - [new_x(ind1); new_x(ind2)])-z(i))^2;
    end
    
    next_f = nf1 + nf2;
    
    if next_f < previous_f
        lambda(k+1) = 0.7*lambda(k);
        x = new_x;
    else
        lambda(k+1) = 2*lambda(k);
    end
    
    k = k + 1
end

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
%% Task9
clear all; close all; clc;
load lmdata2.mat

%xv= randi([-15 15],16,5);

xinit= [-11 -10 14 -6 -8 5 -15 -13 -6 -9 6 6 11 2 2 1]'; %best solution
x=xinit;

% for j=1:length(xv)
%  xinit=xv(:,j);
%  x=xinit;

eps = 10^-6;
lambda(1) = 1;

k = 1;
criteria = 0;
while criteria == 0
    previous_f = 0;
    next_f = 0;
    gf1 = zeros(length(x),1);
    gf2 = zeros(length(x),1);
    
    gp1 = zeros(length(x), length(y));
    
    f1 = 0;
    f2 = 0;
    fp = [];
    
    aux1=0;
    gpaux1=0;
    aux2=0;
    gpaux2=0;
    
    
    for i = 1:length(y)
        ind1 = 2*iA(i,2)-1;
        ind2 = 2*iA(i,2);
        
        aux1 = [gf1(ind1); gf1(ind2)] + (-2)*(norm(A(:,iA(i,1)) - [x(ind1); x(ind2)])-y(i)) * ...
            (A(:,iA(i,1)) - [x(ind1); x(ind2)])/norm(A(:,iA(i,1)) - [x(ind1); x(ind2)]);
        gf1(ind1) = aux1(1);
        gf1(ind2) = aux1(2);
        
        gpaux1 = -(A(:,iA(i,1)) - [x(ind1); x(ind2)])/norm(A(:,iA(i,1)) - [x(ind1); x(ind2)]);
        
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
        
        aux2 = [gf2(ind1); gf2(ind2); gf2(idx1); gf2(idx2)] + (-2)*(norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)])-z(i)) * ...
            ([x(idx1); x(idx2); x(ind1); x(ind2)] - [x(ind1); x(ind2);x(idx1); x(idx2)])/norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)]);
        gf2(ind1) = aux2(1);
        gf2(ind2) = aux2(2);
        gf2(idx1) = aux2(3);
        gf2(idx2) = aux2(4);
        
        gpaux2 = -([x(idx1); x(idx2); x(ind1); x(ind2)] - [x(ind1); x(ind2);x(idx1); x(idx2)])/norm([x(idx1); x(idx2)] - [x(ind1); x(ind2)]);
        
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
        
        nf1 = nf1 + (norm(A(:,iA(i,1)) - [new_x(ind1); new_x(ind2)])-y(i))^2;
    end
    for i = 1:length(z)
        ind1 = 2*iS(i,2)-1;
        ind2 = 2*iS(i,2);
        
        idx1 = 2*iS(i,1)-1;
        idx2 = 2*iS(i,1);
        
        nf2 = nf2 + (norm([new_x(idx1); new_x(idx2)] - [new_x(ind1); new_x(ind2)])-z(i))^2;
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
%sa = S(:,iA(:,2))';

% s1 = S(:,iS(:,1))';
% s2 = S(:,iS(:,2))';

figure()
%plot(S(1,:), S(2,:), 'b.', 'MarkerSize', 15); hold on;
plot(A(1,:), A(2,:), 'rs', 'LineWidth', 2.5, 'MarkerSize', 10); hold on;
plot(x(1,:), x(2,:), 'bo', 'MarkerSize', 10, 'LineWidth', 1.5); hold on;
plot(xi(1,:), xi(2,:), 'b*', 'LineWidth', 1.5, 'MarkerSize', 8); hold on;
%plot([sa(:,1) a(:,1)]',[sa(:,2) a(:,2)]', 'm--'); hold on;
%plot([s1(:,1) s2(:,1)]',[s1(:,2) s2(:,2)]', 'm--');
axis([-15 15 -15 15]);
grid on;


fprintf("Cost Function Value %f\n", previous_f);
% min(j)=previous_f;
% k_val(j)=k;
% end