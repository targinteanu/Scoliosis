function [coef, func] = EvenPolyFourierSeries(fun, L, N)

A = PolyCosineMatrix(L, N);

M = N/2;
alphaT = zeros(1, M+1);
for m = 1:1:M
    intgrd = @(x) fun(x).*cos(m*x);
    alphaT(m+1) = (1/L)*integral(intgrd, -L, L);
end
alphaT(1) = (1/(2*L))*integral(fun, -L, L);

coef = alphaT*(A^-1);
func = @(x) 0;
n = 0;
for c = coef
    func = @(x) func(x) + c*x.^n;
    n = n+2;
end
figure; fplot(fun, [-L, L]);
hold on; fplot(func, [-L, L]);
legend('original', 'estimate')

end