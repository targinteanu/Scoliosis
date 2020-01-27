%% input variables: ---------------------------------------------------
r = @(s) [0*s, 0*s, s];
dr = @(s) [0, 0, 1];
v = @(s) [cos(pi/2 + s*pi/(2*pi*sqrt(2))), sin(pi/2 + s*pi/(2*pi*sqrt(2))), 0];
npts = 50;
Smin = 0; Smax = 2*4*pi*sqrt(2);
% --------------------------------------------------------------------

%% operations

S = linspace(Smin, Smax, npts)';
vr = @(s) v(s) - proj(v(s), r(s)); ur = @(s) vr(s)/norm(vr(s));
T = @(s) dr(s)/norm(dr(s)); % (unit) tangent vector 

R = r(S); VR = cell2mat( arrayfun(vr, S, 'UniformOutput', false) );

Twist_integral = integral( @(t) dTw(t, T, ur, dr), S(1), S(end)) / (2*pi)
Twist_estimate = deturckTwist2(R, VR)
Writhe_integral = integral2( @(t1,t2) ddWr(t1,t2, r, r, dr, dr), ...
    S(1), S(end), S(1), S(end), 'Method', 'iterated') / (4*pi)
Writhe_estimate = levittWrithe(R)

figure; 
plot3dSpine(R, VR); view([-30,10]);
zmax = r(Smax); zmax = zmax(3);
xlim([-2, 2]); ylim([-2, 2]); zlim([0, zmax]);

%% twist and writhe functions 

function proj_ = proj(u, x)
    denom = (x*x');
    if denom
        proj_ = ((u*x')/denom)*x;
    else
        proj_ = 0*x;
    end
end

function ds_ = ds(t, dr)
% arc length from ti to tf
t = t';
ds_ = zeros(size(t)); imax = length(t);
    for i = 1:imax
        tt = t(i);
        ds_(i) = norm(dr(tt));
    end
    ds_ = ds_';
end

    function dTw_ = dTw(t, T, u, DR)
    % u, du must be entered s.t. u perp dr
    % dr must be a unit vector 
        t = t';
        dTw_ = zeros(size(t)); AL = zeros(size(t));
        U = zeros(length(t), length(u(0)));
        imax = length(t);
        for i = 1:imax
            tt = t(i);
            U(i,:) = u(tt);
            AL(i) = integral( @(ttt) ds(ttt, DR), 0, tt );
        end
        du = diff(U)./diff(AL);
        du = [zeros(size(u(0))); du];
        for i = 1:imax
            tt = t(i);
            dTw_(i) = cross(T(tt), u(tt)) * du(i,:)';
        end
        dTw_ = dTw_';
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