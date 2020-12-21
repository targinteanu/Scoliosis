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
[b2,x0] = WLStangent(X,Y,wts,1, m,k,x1);
x_2 = (min(0, x0)):.1:50;
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

function [beta,X0] = WLStangent(x,y,w,N, m0,c0,xstart)
    X = cell2mat(arrayfun(@(xi) xi.^(0:N), x, 'UniformOutput', false));
    
    function eqn = eqnsolver(vars, X,y,W,N, m0,c0)
        
        %x0 = vars(1); b1 = vars(2); b = vars(3:end);
        x0 = vars(1); b = vars(2:end);
    
        eqnsB = ((X'*W*X)^-1)*X'*W*y - b;
        
        b = b';
    
        %{
        i = 2:N;
        eqnX0 = b(1) - c0 + x0.*(m0 + 1/m0) + sum((i-1).*b(2:end).*(x0.^i));
        eqnB1 = m0 - b1 - sum(i.*b(2:end).*(x0.^(i-1)));
                
        eqn = [eqnX0; eqnB1; eqnsB];
        %}
        
        i = 1:N;
        eqn1 = b(1) - m0*x0 - c0 + sum(b.*(x0.^i));
        eqn2 = (1/m0) + sum(i.*b.*(x0.^(i-1)));
        
        eqn = [eqn1; eqn2; eqnsB];
    
    end
    
    options = optimoptions('fsolve','Display','iter','PlotFcn',@optimplotfirstorderopt);
    solvars = fsolve(@(v) eqnsolver(v, X,y,diag(w),N, m0,c0),...
        [xstart; zeros(N+1,1)], options);
    %    [xstart; zeros(N+1,1)], options);
    eqnsolver(solvars, X,y,diag(w),N, m0,c0)
    
    %beta = [solvars(3); solvars(2); solvars(4:end)];
    beta = solvars(2:end);
    X0 = solvars(1);
end