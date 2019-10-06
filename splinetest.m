[num, txt] = xlsread('Writhe-pre-post_new-metrics.csv');
N = 32;
num = num(1:N, :);
XYZ = num(:, 13:end); 
    x = fliplr(XYZ(idx, 1:3:51)); 
    y = fliplr(XYZ(idx, 2:3:51)); 
    z = fliplr(XYZ(idx, 3:3:51)); 
    p = [x;y;z]';
    var = x; nvar = 1;
    
%figure; plot3(p(:,1), p(:,2), p(:,3), 'o'); grid on;

deriv_front = @(cf, n) factorial(n)*cf(:,(end-n));

%pdisp = p(2:end,:)-p(1,:);
pdisp = diff(p);
L = sqrt(diag(pdisp*pdisp')); 
%L = [0;L];
L = [0; cumsum(L)]; 

nplots = 4;
%figure; subplot(1,nplots,1); plot(var, L, 'ok'); grid on; hold on; ylim([L(1), L(end)]);

%pp = [spline(L, x), spline(L, y), spline(L, z)];
pp = spline(L, var);

%plot(linspace(p(1,3), p(end,3)), ppval(ppx, linspace(p(1,3), p(end,3))));
interpfactor = 1;
Ltot = (L(end)-L(1)); interpext = 2/length(L);
Linterp = linspace(L(1)-(interpext*Ltot), L(end)+(interpext*Ltot), ...
    length(L)*(2*interpext + 1)*interpfactor)';
Linterp2 = flipud(interp(flipud(L), interpfactor));

xbrute = ppval(pp, Linterp);
dxbrute = diff(xbrute)./diff(Linterp); 
ddxbrute = diff(dxbrute)./diff(Linterp(2:end)); 
dddxbrute = diff(ddxbrute)./diff(Linterp(2:(end-1)));
%plot(xbrute, Linterp, '--k', 'LineWidth', 1);

cf = pp.coefs; 
b = pp.breaks; 
for i = 1:(length(b)-1)
    t = linspace(b(i), b(i+1)); t2 = t-b(i);
    cfi = cf(i,:);
    %spl = cfi(4)*t.^3 + cfi(3)*t.^2 + cfi(2)*t + cfi(1);
    cfi = fliplr(cfi);
    spl = zeros(size(t));
    for j = 1:length(cfi)
        spl = spl + cfi(j)*t2.^(j-1);
    end
%    plot(spl, t, 'b');
end

    q = 4; % 2 vertebrae above, 2 vertebrae below, current vertebra -> 5 points to fit each cubic
    vertebrae = (1+q):(size(p,1)-q); 
    taulew = zeros(size(vertebrae)); % local torsion at each point on the spine 
    %d2x = zeros(size(vertebrae)); d1x = d2x; d3x = d2x;
    %d2z = zeros(size(vertebrae)); d1z = d2z; d3z = d2z;
    d2 = zeros(size(vertebrae)); d1 = d2; d3 = d2;
    for vertebra = vertebrae
        [taulew(vertebra-q), d, dd, ddd] = lewinerTorsion(p, vertebra, q);
        %d2x(vertebra-q) = dd(nvar); d1x(vertebra-q) = d(nvar); d3x(vertebra-q) = ddd(nvar);
        %d2z(vertebra-q) = dd(3); d1z(vertebra-q) = d(3); d3z(vertebra-1) = ddd(3);
        d1(vertebra-q) = norm(d); d2(vertebra-q) = norm(dd); d3(vertebra-q) = norm(ddd);
    end
    %d2 = d2x./d2z; d1 = d1x./d1z; d3 = d3x./d3z;
    %d2 = d2x; d1 = d1x; d3 = d3x;

    %{
for i = 1:length(taulew)
    vertebra = vertebrae(i);
    t = linspace(L(vertebra-1), L(vertebra+1)); t2 = t-L(vertebra);
    cub = var(vertebra) + d1(i).*t2 + (d2(i)/2).*t2.^2 + (d3(i)/6).*t2.^3;
    %t = linspace(L(vertebra-q), L(vertebra+q));
    plot(cub, t, 'r');
end
    %}
    
df1 = deriv_front(cf, 1); df2 = deriv_front(cf, 2); df3 = deriv_front(cf, 3);
%db1 = deriv_back(cf, 1, b); db2 = deriv_back(cf, 1, b); db3 = deriv_back(cf, 3, b);

%{
tauspline = arrayfun(@(i) -(cross(df1(i), df2(i)) * df3(i)')/(norm(cross(df1(i), df2(i)))^2), ...
    1:length(df3));
taubrut = arrayfun(@(i) ...
    -(cross(dxbrut(i),ddxbrut(i))*dddxbrut(i)')/(norm(cross(dxbrut(i),ddxbrut(i)))^2), ...
    1:length(dddxbrut));
subplot(1, nplots, 5); plot(taubrut, Linterp(3:(end-1)), '--k', 'LineWidth', 1);
grid on; hold on; 
ylim([L(1), L(end)]);
plot(taulew, L(vertebrae), '^r');
plot(tauspline, L(1:(end-1)), '+b');
%}

%%
ppall = [spline(L, x), spline(L, y), spline(L, z)];
cfs = {ppall.coefs}; 

%{
df = arrayfun(@(n) ...
    arrayfun(@(i) deriv_front(cfs{i}, n), 1:3, 'UniformOutput', false), ...
    1:length(cfs), 'UniformOutput', false);
%}

dr = cell2mat(arrayfun(@(var) deriv_front(cfs{var}, 1), 1:3, 'UniformOutput', false)); 
ddr = cell2mat(arrayfun(@(var) deriv_front(cfs{var}, 2), 1:3, 'UniformOutput', false)); 
dddr = cell2mat(arrayfun(@(var) deriv_front(cfs{var}, 3), 1:3, 'UniformOutput', false)); 

tauspline = zeros(size(dr,1), 1); 
for i = 1:length(tauspline)
    cx = cross(dr(i,:), ddr(i,:));
    tauspline(i) = -(cx * dddr(i,:)')/(norm(cx)^2);
end

df1 = arrayfun(@(i) norm(dr(i,:)), 1:size(dr,1)); 
df2 = arrayfun(@(i) norm(ddr(i,:)), 1:size(ddr,1)); 
df3 = arrayfun(@(i) norm(dddr(i,:)), 1:size(dddr,1)); 


    pinterp = cell2mat(arrayfun(@(j) ppval(ppall(j), Linterp), 1:3, 'UniformOutput', false));

    q = q*interpfactor; 
    vinterp = (1+q):(size(pinterp,1)-q); 
    tauinterp = zeros(size(vinterp)); % local torsion at each point on the spine 
    %d2x = zeros(size(vertebrae)); d1x = d2x; d3x = d2x;
    %d2z = zeros(size(vertebrae)); d1z = d2z; d3z = d2z;
    d2int = zeros(size(vinterp)); d1int = d2int; d3int = d2int;
    for vertebra = vinterp
        [tauinterp(vertebra-q), d, dd, ddd] = lewinerTorsion(pinterp, vertebra, q);
        %d2x(vertebra-q) = dd(nvar); d1x(vertebra-q) = d(nvar); d3x(vertebra-q) = ddd(nvar);
        %d2z(vertebra-q) = dd(3); d1z(vertebra-q) = d(3); d3z(vertebra-1) = ddd(3);
        d1int(vertebra-q) = norm(d); d2int(vertebra-q) = norm(dd); d3int(vertebra-q) = norm(ddd);
    end
    
    % do the same thing using interp instead of spline
    pinterp2 = flipud(cell2mat(arrayfun(@(j) interp(flipud(p(:,j)), interpfactor), 1:3,...
        'UniformOutput', false)));
    vinterp2 = (1+q):(size(pinterp2,1)-q); 
    tauinterp2 = zeros(size(vinterp2));
    d2int2 = zeros(size(vinterp2)); d1int2 = d2int2; d3int2 = d2int2;
    for vertebra = vinterp2
        [tauinterp2(vertebra-q), d, dd, ddd] = lewinerTorsion(pinterp2, vertebra, q);
        d1int2(vertebra-q) = norm(d); d2int2(vertebra-q) = norm(dd); d3int2(vertebra-q) = norm(ddd);
    end

figure; subplot(1, nplots, 1);
plot(tauspline, L(1:(end-1)), '-ob'); hold on; grid on; 
plot(tauinterp, Linterp(vinterp), '.c');
plot(tauinterp2, Linterp2(vinterp2), '.m');
plot(taulew, L(vertebrae), '-^r');
ylim([Linterp(1), Linterp(end)]);
title('torsion');
    
%%
subplot(1, nplots, 2); %plot(diff(var)'./diff(L), (L(1:(end-1))+L(2:end))/2, 'ok'); 
grid on; hold on;
%plot(dxbrute, Linterp(2:end), '--k', 'LineWidth', 1);
ylim([Linterp(1), Linterp(end)]);
plot(d1, L(vertebrae), '-^r'); 
plot(d1int, Linterp(vinterp), '.c');
plot(d1int2, Linterp2(vinterp2), '.m');
plot(df1, L(1:(end-1)), '-+b');
%plot(db1, L(2:end), 'xb');
title('d/ds'); 

subplot(1, nplots, 3); %plot(diff(diff(var)'./diff(L))./diff(L(2:end)), L(2:(end-1)), 'ok');
grid on; hold on;
%plot(ddxbrute, Linterp(2:(end-1)), '--k', 'LineWidth', 1);
ylim([Linterp(1), Linterp(end)]);
plot(d2, L(vertebrae), '-^r'); 
plot(d2int, Linterp(vinterp), '.c');
plot(d2int2, Linterp2(vinterp2), '.m');
plot(df2, L(1:(end-1)), '-+b');
%plot(db2, L(2:end), 'xb');
title('d^2/ds^2');

%{
subplot(1, nplots, 4); 
%plot( diff(diff(diff(var)'./diff(L))./diff(L(2:end)))./diff(L(2:(end-1))), L(3:(end-1)), 'ok');
grid on; hold on;
%plot(dddxbrute, Linterp(3:(end-1)), '--k', 'LineWidth', 1);
ylim([L(1), L(end)]);
plot(d3, L(vertebrae), '-^r'); 
plot(df3, L(1:(end-1)), '-+b');
%plot(db3, L(2:end), 'xb');
title('d^3/ds^3');
%}

subplot(1, nplots, 4); 
grid on; hold on;
plot(sqrt(x.^2 + y.^2), L, '-ok');
plot(sqrt(pinterp(:,1).^2 + pinterp(:,2).^2), Linterp, '.c');
plot(sqrt(pinterp2(:,1).^2 + pinterp2(:,2).^2), Linterp2, '.m');
ylim([Linterp(1), Linterp(end)]);
title('(x^2 + y^2)^{1/2}');

%%    

%{
function db = deriv_back(cf, n, b)
    b = diff(b)';
    cf = cf(:,1:(end-n)); cf = fliplr(cf); 
    db = zeros(size(b));
    for j = 1:size(cf,2)
        db = db + cf(:,j).*(b.^(j-1));
    end
end
%}