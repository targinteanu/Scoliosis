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
[b2,x0] = WLStangent(X,Y,wts,3, m,k,0,zeros(1,3));
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

function [beta,X0] = WLStangent(x,y,w,N, m0,c0, xstart,bstart)
    options = optimoptions('fsolve','Display','none');%'iter','PlotFcn',@optimplotfirstorderopt);
    coeff2fun = @(coeff, x) arrayfun(@(i) sum(coeff.*(x(i).^(0:(length(coeff)-1)))), 1:length(x));
    
    i = 2:N;
    % 2 solve for x0
    eqnX0 = @(x0, b,i, m0,c0) c0 + x0.*(m0 + 1/m0) - b(1) + sum((i-1).*b(3:end).*(x0.^i));
    % 3 solve for b1
    eqnB1 = @(b1, b,i,x0, m0) -b1 - 1/m0 - sum(i.*b(3:end).*(x0.^(i-1)));
    
    % init vars for loop
    nsteps = 10000;
    learnrt = 0;
    b = bstart;
    plt = [0,0,learnrt,b];
    errs = zeros(nsteps,1); dderrs = zeros(size(errs));
    KI=1e-17; KD=1; KP=1e-15;
    
    figure('Position', [50 100 1300 700]); 
    subplot(1,3,1); errplt1 = plot(plt(:,1)); ht1 = title('');
    hold on; errplt2 = plot([0, diff(plt(:,1))]); errplt3 = plot(plt(:,2));
    subplot(1,3,2); learnplt = plot(plt(:,3)); ht2 = title('');
    subplot(1,3,3); coeffplt=cell(length(b),1);
    for j = 1:length(b)
        coeffplt{j} = plot(plt(:,3+j)); hold on; 
    end
    pause(.01);
    
    for n = 1:nsteps        
        % 2 solve for x0
        x0solved = fsolve(@(x0) eqnX0(x0, b,i,m0,c0), ...
            xstart, options);
        % 3 solve for b1
        b1solved = fsolve(@(b1) eqnB1(b1, b,i,x0solved,m0), ...
            bstart(2), options);
        b(2) = b1solved;
        
        % 4 calculate gradient 
        derrdb = zeros(size(b));
        err=0; 
        for j = 1:length(x)
            yj = coeff2fun(b, x(j)); err = err+(y(j)-yj)^2; 
            
            dx0db0 = 1./(m0 + 1/m0 + sum(i.*(i-1).*b(3:end).*(x0solved.^(i-1))));
            dx0dbk = @(k) -((k-1).*x0solved.^k)./(m0 + 1/m0 + sum(i.*(i-1).*b(3:end).*(x0solved.^(i-1))));
            
            dydb0 = 1 - x(j)*dx0db0.*sum(i.*(i-1).*b(3:end).*(x0solved.^(i-2)));
            dydbk = @(k) x(j).^k - x(j).*(...
                sum(i.*(i-1).*b(3:end).*(x0solved.^(i-2))).*dx0dbk(k) + k.*x0solved.^(k-1) );
            
            dydb = [dydb0, 0, arrayfun(@(k) dydbk(k), 2:(N-1))];
            derrdb = derrdb + w(j)*(y(j)-yj)*dydb;
        end
        
        % error tracking, learning rate adjustment 
        errs(n)=err;
        if n>1
            derr = (errs(n)-errs(n-1));
            if n>2
                dderrs(n) = (errs(n)+errs(n-2)-2*errs(n-1));
            else
                dderrs(n) = 0;
            end
        else
            derr=0;
            dderrs(n)=0;
        end
        dderr = mean(abs(dderrs(max(n-20,1):n)));
        learnrt = (learnrt + KI*abs(err) - KP*(derr))/(1 + KD*abs(dderr));
        
        plt = [plt; err,norm(derrdb),learnrt,b];
        
        ht1.String = num2str(err);
        errplt1.YData = plt(:,1); errplt3.YData = plt(:,2); errplt2.YData = diff(plt(:,1));
        learnplt.YData = plt(:,3); ht2.String = num2str(learnrt);
        for j = 1:length(coeffplt)
            coeffplt{j}.YData = plt(:,3+j);
        end
        pause(.00001);
        
        % gradient descend
        b = b + learnrt*derrdb;
        
    end
    
    % 2 solve for x0
    x0solved = fsolve(@(x0) eqnX0(x0, b,i,m0,c0), ...
        xstart);
    % 3 solve for b1
    b1solved = fsolve(@(b1) eqnB1(b1, b,i,x0solved,m0), ...
        bstart(2));
    b(2) = b1solved;
    
    beta = b';
    X0 = x0solved;
end