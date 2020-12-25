X = ([0:5:50,0:5:50,0:5:50])'; 
Y = (X - (30*rand + 10)).*(X - (30*rand + 10)).*(X - (30*rand + 10)) + 10000*(rand(size(X))-.5);
wts = 100*rand(size(X)) + .001;
x = 0:.1:50;
a = (wts(1:11).*Y(1:11) + wts(12:22).*Y(12:22) + wts(23:33).*Y(23:33))./(wts(1:11)+wts(12:22)+wts(23:33));

x1 = -5; x2 = 5; y1 = -1000; y2 = 1000;
m = (y2-y1)/(x2-x1); k = y1-m*x1;

coeff2fun = @(coeff, x) arrayfun(@(i) sum(coeff.*(x(i).^(0:(length(coeff)-1)))), 1:length(x));
b = WLS(X,Y,wts,10);
y = coeff2fun(b', x); 
N = 3;
[b2,x0] = WLStangent(X,Y,wts,N, m,k,0,zeros(1,N+1));
x_2 = x0:.1:50;
y_2 = coeff2fun(b2', x_2);
dydx0 = (y_2(2)-y_2(1))/(x_2(2)-x_2(1));

figure; errorbar(X,Y, 10000./wts, 'o');
hold on; plot(x, y);
plot(X(1:11), a, '*')
plot([x1,x2],[y1,y2],'--k','LineWidth',2);
plot(x_2, y_2);

function beta = WLS(x,y,w,N)
    X = cell2mat(arrayfun(@(xi) xi.^(0:N), x, 'UniformOutput', false));
    beta = ((X'*diag(w)*X)^-1)*X'*diag(w)*y;
end

function derrdb = WLSgrad(b, x,y,N,w, m0,c0, coeff2fun, xstart,options)
i = 2:N;

eqnX0 = @(x0, b,i, m0,c0) c0 + x0.*(m0 + 1/m0) - b(1) + sum((i-1).*b(3:end).*(x0.^i));
x0solved = fsolve(@(x0) eqnX0(x0, b,i,m0,c0), xstart, options);
eqnB1 = @(b1, b,i,x0, m0) -b1 - 1/m0 - sum(i.*b(3:end).*(x0.^(i-1)));
b1solved = fsolve(@(b1) eqnB1(b1, b,i,x0solved,m0), b(2), options);
b(2) = b1solved;

derrdb = zeros(size(b));
for j = 1:length(x)
    yj = coeff2fun(b, x(j));
    
    dx0db0 = 1./(m0 + 1/m0 + sum(i.*(i-1).*b(3:end).*(x0solved.^(i-1))));
    dx0dbk = @(k) -((k-1).*x0solved.^k)./(m0 + 1/m0 + sum(i.*(i-1).*b(3:end).*(x0solved.^(i-1))));
    
    dydb0 = 1 - x(j)*dx0db0.*sum(i.*(i-1).*b(3:end).*(x0solved.^(i-2)));
    dydbk = @(k) x(j).^k - x(j).*(...
        sum(i.*(i-1).*b(3:end).*(x0solved.^(i-2))).*dx0dbk(k) + k.*x0solved.^(k-1) );
    
    dydb = [dydb0, 0, arrayfun(@(k) dydbk(k), 2:N)];
    derrdb = derrdb + w(j)*(y(j)-yj)*dydb;
end
end

function y = getYvalue(b,x, N, m0,c0, coeff2fun, xstart,options)
i = 2:N;

eqnX0 = @(x0, b,i, m0,c0) c0 + x0.*(m0 + 1/m0) - b(1) + sum((i-1).*b(3:end).*(x0.^i));
x0solved = fsolve(@(x0) eqnX0(x0, b,i,m0,c0), xstart, options);
eqnB1 = @(b1, b,i,x0, m0) -b1 - 1/m0 - sum(i.*b(3:end).*(x0.^(i-1)));
b1solved = fsolve(@(b1) eqnB1(b1, b,i,x0solved,m0), b(2), options);
b(2) = b1solved;

y = coeff2fun(b,x);
end

function [beta,X0] = WLStangent(x,y,w,N, m0,c0, xstart,bstart)
    options = optimoptions('fsolve','Display','none');%'iter','PlotFcn',@optimplotfirstorderopt);
    coeff2fun = @(coeff, x) arrayfun(@(i) sum(coeff.*(x(i).^(0:(length(coeff)-1)))), 1:length(x));
    
    beta = AdaptGradDesc(@(b) WLSgrad(b, x,y,N,w, m0,c0, coeff2fun, xstart,options),...
        @(b,x) getYvalue(b,x, N, m0,c0, coeff2fun, xstart,options),...
        x, y, bstart, 0);
    
    X0 = 0;
end