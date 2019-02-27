r = @(s) [cos(s/sqrt(2)), sin(s/sqrt(2)), s/sqrt(2)];
dr = @(s) [-sin(s/sqrt(2)), cos(s/sqrt(2)), 1]/sqrt(2);
ddr = @(s) [-cos(s/sqrt(2)), -sin(s/sqrt(2)), 0]/2;
dddr = @(s) [sin(s/sqrt(2)), -cos(s/sqrt(2)), 0]/(2*sqrt(2));
smax = 2*pi*sqrt(2);

coarse = linspace(0, smax, 17)'; fine = (0:.01:smax)';
discrete_approx_coarse = r(coarse);
discrete_approx_fine = r(fine);

figure; plot3(discrete_approx_fine(:,1), discrete_approx_fine(:,2), discrete_approx_fine(:,3));
hold on; plot3(discrete_approx_coarse(:,1), discrete_approx_coarse(:,2), discrete_approx_coarse(:,3));
grid on;
xlabel('x'); ylabel('y'); zlabel('z');
legend('fine', 'coarse'); 
title('Discrete curves with "coarse" and "fine" spacing');

Writhe = integral2(@(x,y) ddWr(x, y, r, r, dr, dr), 0,smax,0,smax,...
    'Method', 'iterated')/(4*pi)
botWrithe = integral2(@(x,y) ddWr(x, y, r, r, dr, dr), 0,smax/2,0,smax/2,...
    'Method', 'iterated')/(4*pi)
topWrithe = integral2(@(x,y) ddWr(x, y, r, r, dr, dr), smax/2,smax,smax/2,smax,...
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

total_Torsion = integral(@(t) torsion(t, dr, ddr, dddr), 0,smax)
Torsion_coarse = torsion(coarse, dr, ddr, dddr);
Torsion_fine = torsion(fine, dr, ddr, dddr);
Torsion_diff_est_coarse = getTorsion(discrete_approx_coarse);
Torsion_diff_est_fine = getTorsion(discrete_approx_fine);

q = 3; 
vertebrae = (1+q):(17-q);
    Torsion_lewiner_coarse = zeros(size(vertebrae));
    for vertebra = vertebrae
        Torsion_lewiner_coarse(vertebra-q) = ...
            lewinerTorsion(discrete_approx_coarse, vertebra, q);
    end
    
figure; 
plot(coarse(vertebrae), Torsion_lewiner_coarse);
hold on; plot(fine, Torsion_fine); plot(fine, Torsion_diff_est_fine);
plot(coarse, Torsion_diff_est_coarse); grid on;
legend('lewiner', 'actual', 'fine diff est', 'coarse diff est');
xlabel('arc length'); ylabel('torsion');
title('Comparing methods of estimating Torsion');

function tau = torsion(t, dr, ddr, dddr)
    tau = zeros(size(t));
    for j = 1:length(t)
        tt = t(j);
        cx = cross(dr(tt), ddr(tt));
        if norm(cx)
            tau(j) = cx*dddr(tt)' / (norm(cx)^2);
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