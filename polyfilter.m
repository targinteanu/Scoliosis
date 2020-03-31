ply_helper = @(p, x) p(1) + p(2)*x + p(3)*x.^2 + p(4)*x.^3;

ply = @(p, x, b) ply_helper(p, x) .* (x>b(1) & x<b(2));

filtply_helper = @(w, p, x, t) (1/(pi*w^3))*(...
    w^3 * sinint(w*(t-x)) .* (ply_helper(p, x)) -...
    cos(w*(t-x)) .* (w^2 * (p(2) + p(3)*(t+x)) + p(4)*(t*t*w*w + t*w*w*x + w*w.*x.*x - 2)) +...
    w * (p(3) + p(4)*(2*t + x)) .* sin(w*(t-x)));

filtply = @(w, p, x, b) filtply_helper(w, p, x, b(2)) - filtply_helper(w, p, x, b(1));

FScoef_helper = @(p, x, n, T) ((exp(-2*j*pi*n*x/T))/(8 * pi^4 * n^4)) .* (...
    4*j*p(1) * pi^3 * n^3 + 2*p(2) * pi^2 * n^2 *(T + 2*j*pi*n*x) + ...
    4*j*p(3) * pi^3 * n^3 * x.^2 + 4*p(3)*T*x * pi^2 * n^2 -...
    2*j*pi*p(3)*n * T^2 + 4*j*p(4) * pi^3 * n^3 * x.^3 +...
    6*p(4)*T * pi^2 * n^2 * x.^2 - 6*j*pi*p(4)*n*x * T^2 - 3*p(4) * T^3);

FScoef0_helper = @(p, x, T) (p(1)*x + (p(2)*x.^2)/2 + (p(3)*x.^3)/3 + (p(4)*x.^4)/4)/T;

FScoef = @(p, n, a, b) FScoef_helper(p, b, n, b-a) - FScoef_helper(p, a, n, b-a);
FScoef0 = @(p, a, b) FScoef0_helper(p, b, b-a) - FScoef0_helper(p, a, b-a);

FSn = @(p, x, n, b) exp(j*2*pi*n*x/(b(2)-b(1))) * FScoef(p, n, b(1), b(2));
FS0 = @(p, x, b) FScoef0(p, b(1), b(2));

FSfilt = @(w, p, x, b) sum(cell2mat(arrayfun(@(n) FSn(p,x,n,b)', ...
    [-ceil(w*(b(2)-b(1))/(2*pi)):-1, 1:ceil(w*(b(2)-b(1))/(2*pi))], ...
    'UniformOutput', false))') + FS0(p, x, b);

P = [.2, -3, -5.8, .6];
B = [0, 10]; 

x = linspace(B(1)-10, B(2)+10, 1000);
PLY = ply(P, x, B);
figure; plot(x, PLY, '-k'); 
hold on; grid on; 
wvals = [2, 5, 8]; 
PLYfilt = zeros(length(wvals), length(PLY));
PLYfiltFS = zeros(size(PLYfilt));
colr = [0,.447,.741; .85,.325,.098; .929,.694,.125; .494,.184,.556; .446,.674,.188; .301,.745,.933; .653,.078,.184];
for i = 1:length(wvals)
    w = wvals(i);
    PLYfilt(i,:) = filtply(w, P, x, B);
    PLYfiltFS(i,:) = real(FSfilt(w, P, x, B));
end
for i = 1:length(wvals)
    plot(x, PLYfilt(i,:), '--', 'Color', colr(i,:))
    plot(x, PLYfiltFS(i,:), ':', 'Color', colr(i,:))
end
PLYfilt = PLYfilt'; PLYfiltFS = PLYfiltFS';
wvals2 = reshape([wvals; wvals], 1, 2*length(wvals));
lgd = arrayfun(@(w) num2str(w), wvals2, 'UniformOutput', false);
legend(['Inf', lgd])