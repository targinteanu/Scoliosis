function [coeffs, y_out] = WLSperp(x_in,y_in,w,N, x_perp,y_perp, KPID)

    I = 2:N;
    m0 = diff(y_perp)/diff(x_perp); c0 = y_perp(1) - m0*x_perp(1);
    
    xstart = mean(x_perp); bstart = zeros(1,N+1);

    coeff2fun = @(coeff, x) arrayfun(@(i) sum(coeff.*(x(i).^(0:(length(coeff)-1)))), 1:length(x));
    
    eqnX0 = @(x0,b,SP, i,m0,c0) c0 + x0.*(m0 + 1/m0) - b(1) + sum((i-1).*b(3:end).*(x0.^i));
    eqnB1 = @(b1,b,SP, i,m0) -b1 - 1/m0 - sum(i.*b(3:end).*(SP{1}(1).^(i-1)));
    
    paramsToSolve = @(b) {@(x0,b,SP) eqnX0(x0,b,SP, I,m0,c0), @(b1,b,SP) eqnB1(b1,b,SP, I,m0); ...
        xstart, b(2)};
    
    coeffs = AdaptGradDesc(@(b,SP) WLSgrad(b,SP, x_in,y_in,N,w, m0, coeff2fun),...
        @(b,x,SP) getYvalue(b,x, SP, coeff2fun),...
        x_in, y_in, bstart, 0, KPID, paramsToSolve);
        
    y_out = coeff2fun(coeffs', x_in);

function derrdb = WLSgrad(b,solvedParams, x,y,N,w, m0, coeff2fun)
i = 2:N;

x0solved = solvedParams{1}(1);

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

function y = getYvalue(b,x,solvedParams, coeff2fun)

b1solved = solvedParams{2}(1);

b(2) = b1solved;

y = coeff2fun(b,x);

end

end