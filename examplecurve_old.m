%% input variables: ---------------------------------------------------
r = @(s) [sin(s/sqrt(2)), cos(s/sqrt(2)), s/sqrt(2)];
dr = @(s) (1/sqrt(2))*[cos(s/sqrt(2)), -sin(s/sqrt(2)), 1];
v = @(s) [cos(pi/2 + s*pi/(2*pi*sqrt(2))), sin(pi/2 + s*pi/(2*pi*sqrt(2))), 0];
npts = 50;
Smin = 0; Smax = .5*4*pi*sqrt(2);
% --------------------------------------------------------------------

%% operations

S = linspace(Smin, Smax, npts)';
vr = @(s) v(s) - proj(v(s), dr(s)); ur = @(s) vr(s)/norm(vr(s));
T = @(s) dr(s)/norm(dr(s)); % (unit) tangent vector 

R = r(S); VR = cell2mat( arrayfun(vr, S, 'UniformOutput', false) );

%Twist_integral = integral( @(t) dTw(t, T, ur, r), S(1), S(end)) / (2*pi)
dt = .001; t = Smin:dt:Smax; dTW = zeros(size(t));
for i = 2:length(t)
    du = ur(t(i)) - ur(t(i-1));
    ds = norm( r(t(i)) - r(t(i-1)) );
    dTW(i) = cross(T(t(i)), ur(t(i))) * (du)';
end
Twist_integral = sum(dTW) / (2*pi)

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

    function dTw_ = dTw(t, T, u, R)
    % u, du must be entered s.t. u perp dr
    % dr must be a unit vector 
        t = t';
        dTw_ = zeros(size(t)); AL = zeros(size(t));
        U = zeros(length(t), length(u(0)));
        imax = length(t);
        for i = 1:imax
            tt = t(i);
            U(i,:) = u(tt);
        end
        for i = 2:imax
            AL(i) = norm(R(t(i))-R(t(i-1)));
        end
        du = diff(U)./AL(2:end);
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