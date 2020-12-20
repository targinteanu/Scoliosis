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
b2 = WLStangent(X,Y,wts,10, 1/m,k,mean(x1,x2));
b2(2)=0;
y_2 = coeff2fun(b2', x);

figure; errorbar(X,Y, 10000./wts, 'o');
hold on; plot(x, y);
plot(X(1:11), a, '*')
plot([x1,x2],[y1,y2],'--k','LineWidth',2);
plot(x, y_2);

function beta = WLS(x,y,w,N)
    X = cell2mat(arrayfun(@(xi) xi.^(0:N), x, 'UniformOutput', false));
    beta = ((X'*diag(w)*X)^-1)*X'*diag(w)*y;
end

function beta = WLStangent(x,y,w,N, m0,c0,xstart)
    X = cell2mat(arrayfun(@(xi) xi.^([0,2:N]), x, 'UniformOutput', false));
    b = ((X'*diag(w)*X)^-1)*X'*diag(w)*y;
    
    i = 2:N;
    eqn = @(x0) b(1) - c0 + x0.*(m0 + 1/m0) + sum((i-1).*b(2:end).*(x0.^i));
    options = optimoptions('fsolve','Display','none','PlotFcn',@optimplotfirstorderopt);
    solx = fsolve(eqn, xstart, options);
    b1 = m0 - sum(i.*b(2:end)'.*(solx.^(i-1)));
    beta = [b(1); b1; b(2:end)];
end