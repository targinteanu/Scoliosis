[num, txt] = xlsread('Writhe-pre-post_new-metrics.csv');
N = 32;
num = num(1:N, :);
XYZ = num(:, 13:end); 
    x = XYZ(idx, 1:3:51); 
    y = XYZ(idx, 2:3:51); 
    z = XYZ(idx, 3:3:51); 
    p = [x;y;z]';

%figure; plot3(p(:,1), p(:,2), p(:,3), 'o'); grid on;

deriv_front = @(cf, n) factorial(n)*cf(:,(end-n));

nplots = 3;
figure; subplot(1,nplots,1); plot(p(:,1), p(:,3), 'ok'); grid on; hold on;
ppx = spline(p(:,3), p(:,1));
%plot(linspace(p(1,3), p(end,3)), ppval(ppx, linspace(p(1,3), p(end,3))));
cfx = ppx.coefs; 
bx = ppx.breaks; 
for i = 1:(length(bx)-1)
    t = linspace(bx(i), bx(i+1)); t2 = t-bx(i);
    cfi = cfx(i,:);
    %spl = cfi(4)*t.^3 + cfi(3)*t.^2 + cfi(2)*t + cfi(1);
    cfi = fliplr(cfi);
    spl = zeros(size(t));
    for j = 1:length(cfi)
        spl = spl + cfi(j)*t2.^(j-1);
    end
    plot(spl, t, 'b');
end

    q = 4; % 2 vertebrae above, 2 vertebrae below, current vertebra -> 5 points to fit each cubic
    vertebrae = (1+q):(size(p,1)-q); 
    taulew = zeros(size(vertebrae)); % local torsion at each point on the spine 
    d2 = zeros(size(vertebrae)); d1 = d2; d3 = d2;
    for vertebra = vertebrae
        [taulew(vertebra-q), d, dd, ddd] = lewinerTorsion(p, vertebra, q);
        d2x(vertebra-q) = dd(1); d1x(vertebra-q) = d(1); d3x(vertebra-1) = ddd(1);
        d2z(vertebra-q) = dd(3); d1z(vertebra-q) = d(3); d3z(vertebra-1) = ddd(3);
    end
    d2 = d2x./d2z; d1 = d1x./d1z; d3 = d3x./d3z;
    
for i = 1:length(taulew)
    vertebra = vertebrae(i);
    t = linspace((vertebra-q), (vertebra+q)); t2 = t-(vertebra);
    cub = x(vertebra) + d1x(i).*t2 + (d2x(i)/2).*t2.^2 + (d3x(i)/6).*t2.^3;
    t = linspace(z(vertebra-q), z(vertebra+q));
    plot(cub, t, 'r');
end
    
df1 = deriv_front(cfx, 1); df2 = deriv_front(cfx, 2);
db1 = deriv_back(cfx, 1, bx); db2 = deriv_back(cfx, 1, bx);
    
subplot(1, nplots, 2); plot(diff(x)./diff(z), (z(1:(end-1))+z(2:end))/2, 'ok');
grid on; hold on;
plot(d1, p(vertebrae,3), '^r'); 
plot(df1, p(1:(end-1),3), '+b');
plot(db1, p(2:end,3), 'xb');

subplot(1, nplots, 3); plot(diff(diff(x)./diff(z))./diff(z(2:end)), z(2:(end-1)), 'ok');
grid on; hold on;
plot(d2, p(vertebrae,3), '^r'); 
plot(df2, p(1:(end-1),3), '+b');
plot(db2, p(2:end,3), 'xb');

function db = deriv_back(cf, n, b)
    b = diff(b)';
    cf = cf(:,1:(end-n)); cf = fliplr(cf); 
    db = zeros(size(b));
    for j = 1:size(cf,2)
        db = db + cf(:,j).*(b.^(j-1));
    end
end