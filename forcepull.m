function F = forcepull(a, b, alpha, beta, gamma, mech, app)
% inputs: 
%  -- a, b, alpha, beta, gamma: arrays of coefficients for polynomials for r
%  and u, functions of parameter t. Row vectors. 
%  -- mech: struct containing mechanical properties. 
%  -- app: struct containing information about applied loads. 
%  -- a, b, alpha, beta, and gamma must be found; mech and app must be
%  input
% outputs: 
%  -- F(1-N): coordinates of force-bending equation 
%  -- F(N+1): internal torque equation 
%  -- F(N+2): CWF formula
%  -- all elements of F must = 0. 

% extract values 
sF = app.sF; % scalar 
Fapp = app.Fapp; % R3 vector 
r_initial = mech.r_initial; % R3 vector function of t
c = mech.c; % scalar constant
I = mech.I; J = mech.J; E = mech.E; G = mech.G; % column vectors 
dt = mech.dt; t = mech.t; t_final = t(end); t_initial = t(1);
% t = 0:dt:58; t = t'; % t and s are column vectors

% r and u
function r_ = r(t)
    N = min([length(a), length(b)]);
    t_ = t.^(0:N); t_ = t_';
    r_ = zeros(1,3); 
    r_(1) = a*t_; r_(2) = b*t_; r_(3) = c*t;
end
    function q_ = q(t)
        q_ = r(t) - r_initial(t);
    end
function u_ = u(t)
    % u is a unit vector. 
    N = min([length(alpha), length(beta), length(gamma)]);
    t_ = t.^(0:N); t_ = t_';
    u_ = zeros(1,3); 
    u_(1) = sum(alpha.*t_); u_(2) = sum(beta.*t_); u_(3) = sum(gamma.*t_);
    if norm(u_)
        u_ = u_/norm(u_);
    end
end

    function ds_ = ds(t)
        dr_ = r(t) - r(t-dt); 
        dx = dr_(1); dy = dr_(2); dz = dr_(3);
        ds_ = dt*sqrt( (dx/dt)^2 + (dy/dt)^2 + (dz/dt)^2 );
    end

s = zeros(size(t)); 
for idx = 2:length(s)
    s(idx) = s(idx-1) + ds(t(idx));
end
s_final = s(end); s_initial = s(1);

% force bending 
M = zeros(length(t), 3);
for idx = 1:length(s)
    si = s(idx); ti = t(idx);
    if si <= sF
        M(idx,:) = cross(q(ti), Fapp);
        qF = q(ti);
    else
        M(idx,:) = cross(qF, Fapp);
    end
end
forcebend = zeros(size(M));
dq0 = q(0-dt) - q(0-2*dt);
for idx = 1:length(s)
    ti = t(idx);
    dq = q(ti) - q(ti-dt);
    ddq = dq-dq0;
    forcebend(idx,:) = cross(dq, ddq)/( norm(dq)^3 );
    dq0 = dq;
end
forcebend = forcebend + M./(I.*E);

% numerical approx of internal torque and CWF: 
%%{ 
% internal torque and Twist
dphi_ds = zeros(size(s));
dTw_ds = zeros(size(s));
dr00 = r(0-2*dt) - r(0-3*dt);
dr01 = r(0-dt) - r(0-2*dt);
ddr0 = dr01 - dr00;
for idx = 1:length(s)
    ti = t(idx);
    dr = r(ti) - r(ti-dt);
    ddr = dr-dr01;
    dddr = ddr-ddr0;
    
    cx = cross(dr, ddr);
    dphi_ds(idx) = ( cx * dddr' )/( norm(cx)^2 );
    
    du = u(ti) - u(ti-dt);
    dTw_ds(idx) = ((dt/ds(ti))^2) * ( cross(du/dt, u) * (dr/dt)' );
        
    %dr00 = dr01; 
    dr01 = dr;
    ddr0 = ddr;
end
dphi_ds = dTw_ds - dphi_ds;
dT_ds = J.*G.*dphi_ds;
 
T = 0; Tw = 0; 
for idx = 1:length(s)
    ti = t(idx); dsi = ds(ti);
    T = T + dT_ds(idx)*dsi;
    Tw = Tw + dTw_ds(idx)*dsi;
end
Tw = Tw/(2*pi);

% Writhe 
Wr_spine = 0;
for idx1 = 1:length(s)
    t1 = t(idx1);
    r1 = r(t1);
    dr1 = r1 - r(t1-dt);
    for idx2 = 1:length(s)
        t2 = t(idx2);
        r2 = r(t2);
        dr2 = r2 - r(t2-dt);
        r12 = r1 - r2;
        if norm(r12)
            ddW = ( cross(dr1, dr2)*r12' )/( norm(r12)^3 );
            Wr_spine = Wr_spine + ddW;
        end
    end
end
Wr_spine = Wr_spine/(4*pi);
%}

end