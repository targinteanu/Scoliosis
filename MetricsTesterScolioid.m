% test metrics on a scoliosis-like 3D polynomial 

%% build a random scoliotic spine-like curve 

tmax = 58; max_kyphotic_imbalance = .05; max_lordotic_imbalance = .01;
coeffs_cor = [-.0123; .00015; 0] + [.0001;.000001;.0000001].*randn(3,1); 
    c1 = -(tmax.^(2:4))*coeffs_cor/tmax; 
    coeffs_cor = [max_lordotic_imbalance*randn + c1; coeffs_cor];
%coeffs_cor = zeros(4,1);
coeffs_sag = [.074592; -.002004; .000016] + [.001;.00001;.0000001].*randn(3,1);
    c1 = -(tmax.^(2:4))*coeffs_sag/tmax; 
    coeffs_sag = [max_kyphotic_imbalance*randn + c1; coeffs_sag];

r = @(t) [(t.^(1:4))*coeffs_cor, (t.^(1:4))*coeffs_sag, t];
dr = @(t) [((1:4).*t.^(0:3))*coeffs_cor, ((1:4).*t.^(0:3))*coeffs_sag, ones(size(t))];
ddr = @(t) [((0:3).*(1:4).*t.^[0 0:2])*coeffs_cor, ...
    ((0:3).*(1:4).*t.^[0 0:2])*coeffs_sag, zeros(size(t))];
dddr = @(t) [([0 0:2].*(0:3).*(1:4).*t.^[0 0 0:1])*coeffs_cor, ...
    ([0 0:2].*(0:3).*(1:4).*t.^[0 0 0:1])*coeffs_sag, zeros(size(t))];

ds = @(t) arrayfun(@(tt) norm(dr(tt)), t);
smax = integral(ds, 0, tmax);

numvertebrae = 24;
coarse_s = linspace(0, smax, numvertebrae)'; 
coarse = linspace(0, tmax, numvertebrae)';
fine = (0:.1:tmax)';
discrete_approx_coarse = r(coarse);
discrete_approx_fine = r(fine);

figure('Position', [50 100 1400 800]);
ax(1) = subplot(1, 2, 1);
plot3(discrete_approx_fine(:,1), discrete_approx_fine(:,2), discrete_approx_fine(:,3)); hold on;
plot3(discrete_approx_coarse(:,1), discrete_approx_coarse(:,2), discrete_approx_coarse(:,3), '-o');
grid on;
xlabel('cor'); ylabel('sag'); zlabel('ax');
ax(2) = subplot(1, 2, 2);
plot3(discrete_approx_fine(:,2), -discrete_approx_fine(:,1), discrete_approx_fine(:,3)); hold on;
plot3(discrete_approx_coarse(:,2), -discrete_approx_coarse(:,1), discrete_approx_coarse(:,3), '-o');
grid on;
xlabel('sag'); ylabel('cor'); zlabel('ax');
linkprop(ax, {'CameraPosition', 'CameraUpVector', 'PlotBoxAspectRatio'});
clear ax;

%% get Writhe estimates 

Writhe = integral2(@(x,y) ddWr(x, y, r, r, dr, dr), 0,tmax,0,tmax,...
    'Method', 'iterated')/(4*pi)
botWrithe = integral2(@(x,y) ddWr(x, y, r, r, dr, dr), 0,tmax/2,0,tmax/2,...
    'Method', 'iterated')/(4*pi)
topWrithe = integral2(@(x,y) ddWr(x, y, r, r, dr, dr), tmax/2,tmax,tmax/2,tmax,...
    'Method', 'iterated')/(4*pi)

Writhe_diff_est_coarse = getWrithe(discrete_approx_coarse);
Writhe_levitt_est_coarse = levittWrithe(discrete_approx_coarse);
Writhe_diff_est_fine = getWrithe(discrete_approx_fine);
Writhe_levitt_est_fine = levittWrithe(discrete_approx_fine);

figure; bar([Writhe, Writhe_diff_est_coarse, Writhe_levitt_est_coarse, ...
    Writhe_diff_est_fine, Writhe_levitt_est_fine]);
xticklabels({'integral2', 'coarse diff', 'coarse levitt', 'fine diff', 'fine levitt'});
grid on;
xlabel('Estimation Method'); ylabel('Writhe Estimate'); 
title('Comparing methods of estimating Writhe');

%% get Torsion estimates 

total_Torsion = integral(@(t) torsion(t, dr, ddr, dddr), 0,tmax)
Torsion_coarse = torsion(coarse, dr, ddr, dddr);
Torsion_fine = torsion(fine, dr, ddr, dddr);
Torsion_diff_est_coarse = getTorsion(discrete_approx_coarse);
Torsion_diff_est_fine = getTorsion(discrete_approx_fine);

q = 3; 
vertebrae = (1+q):(numvertebrae-q);
concavities = zeros(size(vertebrae));
concavity_true = ddr(fine);
concavity_true = diag(concavity_true*concavity_true').^.5;
    Torsion_lewiner_coarse = zeros(size(vertebrae));
    for vertebra = vertebrae
        [Torsion_lewiner_coarse(vertebra-q),~,concavity] = ...
            lewinerTorsion(discrete_approx_coarse, vertebra, q);
        concavities(vertebra-q) = norm(concavity);
    end
    
figure('Position', [50 100 700 800]); 
ax(1) = subplot(2,1,1);
plot(vertebrae, Torsion_lewiner_coarse);
hold on; plot(fine*numvertebrae/tmax, Torsion_fine); 
plot(fine*numvertebrae/tmax, Torsion_diff_est_fine);
plot(coarse*numvertebrae/tmax, Torsion_diff_est_coarse); grid on;
legend('lewiner', 'actual', 'fine diff est', 'coarse diff est', 'Location', 'southwest');
xlabel('vertebra'); ylabel('torsion');
title('Comparing methods of estimating Torsion');
ax(2) = subplot(2,1,2); 
plot(vertebrae, concavities); grid on; hold on; 
plot(fine*numvertebrae/tmax, concavity_true);
xlabel('vertebra'); ylabel('|| d^2r / dt^2 ||'); 
title('Second derivative');
legend('Lewiner', 'actual');
linkaxes(ax, 'x'); clear ax;

%[~, neutral] = min(concavities); neutral = vertebrae(neutral);
%[~, neutral_true] = min(concavity_true); neutral_true = fine(neutral_true);
neutral = (concavities(2:(end-1)) < concavities(1:(end-2))) & ...
    (concavities(2:(end-1)) < concavities(3:end)); 
neutral = find([0 neutral]); neutral = vertebrae(neutral);
ct = concavity_true'; 
neutral_true = (ct(2:(end-1)) < ct(1:(end-2))) & (ct(2:(end-1)) < ct(3:end));
neutral_true = find([0 neutral_true]);
while (length(neutral_true) > length(neutral))
    [~, elim] = max(ct(neutral_true)); elim = neutral_true(elim);
    ct = ct([1:(elim-1) (elim+1):end]);
    neutral_true = (ct(2:17) < ct(1:16)) & (ct(2:17) < ct(3:18));
    neutral_true = find([0 neutral_true]);
end
neutral_true = fine(neutral_true)';

tdc = Torsion_diff_est_coarse; 
neutral_diff = (tdc(2:(end-1)) < tdc(1:(end-2))) & (tdc(2:(end-1)) < tdc(3:end));
neutral_diff = find([0 neutral_diff]);
while (length(neutral_diff) > length(neutral))
    [~, elim] = max(tdc(neutral_coarse)); elim = neutral_coarse(elim);
    tdc = tdc([1:(elim-1) (elim+1):end]);
    neutral_diff = (tdc(2:(end-1)) < tdc(1:(end-2))) & (tdc(2:(end-1)) < tdc(3:end));
    neutral_diff = find([0 neutral_diff]);
end

[~, apex] = max(concavities); apex = vertebrae(apex);
q = abs(apex - neutral); q = min(q, numvertebrae-neutral);

torsion_kadoury = arrayfun(@(i) lewinerTorsion(discrete_approx_coarse, ...
    neutral(i), q(i)), 1:length(q));
torsion_actual = arrayfun(@(v) torsion(v, dr, ddr, dddr), neutral_true);
torsion_diff = Torsion_diff_est_coarse(neutral_diff);

figure; bar([torsion_actual; torsion_diff; torsion_kadoury]);
xticklabels({'actual', 'coarse diff', 'kadoury'});
grid on;
xlabel('Estimation Method'); ylabel('Torsion Estimates'); 
title('Comparing methods of estimating Torsion');

%% get accuracies 


%% supporting functions 
function tau = torsion(t, dr, ddr, dddr)
    tau = zeros(size(t));
    for j = 1:length(t)
        tt = t(j);
        cx = cross(dr(tt), ddr(tt));
        if norm(cx)
            tau(j) = -cx*dddr(tt)' / (norm(cx)^2);
        end
    end
end

    function ddWr_ = ddWr(t1, t2, r1, r2, dr1, dr2)
        t1 = t1'; t2 = t2';
        r12 = r1(t1) - r2(t2);
        if length(t1) <= length(t2)
            sz = size(t1); imax = length(t1);
        else
            sz = size(t2); imax = length(t2);
        end
        ddWr_ = zeros(sz);
        for i = 1:imax
            tt1 = t1(i); tt2 = t2(i); rr12 = r12(i,:);
            if norm(rr12)
                ddWr_(i) = (cross(dr1(tt1), dr2(tt2)) * rr12')/( norm(rr12)^3 );
            else
                ddWr_(i) = 0;
            end
        end
        ddWr_ = ddWr_';
        %figure; plot3(t1, t2, ddWr_, '.'); grid on; xlabel('t1'); ylabel('t2');
    end