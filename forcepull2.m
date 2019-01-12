function F = forcepull2(a, b, c, alpha, mech, app)
% inputs: 
%  -- a, b, alpha, beta: arrays of coefficients for polynomials for r
%  and u, functions of parameter t. Row vectors. 
%  -- mech: struct containing mechanical properties. 
%  -- app: struct containing information about applied loads. 
%  -- a, b, alpha, and beta must be found; mech and app must be input
% outputs: 
%  -- F(1-N): coordinates of force-bending equation 
%  -- F(N+1): internal torque equation 
%  -- F(N+2): CWF formula
%  -- all elements of F must = 0. 

% extract values 
sF = app.sF; % scalar 
Fapp = app.Fapp; % R3 vector 
r_initial = mech.r_initial; % R3 vector function of t
%c = mech.c; % scalar constant
I = mech.I; J = mech.J; E = mech.E; G = mech.G; % functions of s 
dt = mech.dt; %t = mech.t; t_final = t(end); t_initial = t(1);
% t = 0:dt:58; t = t'; % t and s are column vectors
t_initial = 0; %t_initial_topbound = mech.t_initial_topbound;
s_final = mech.s_final;

% r and u
function r_ = r(t)
    %N = min([length(a), length(b)]);
    %t_ = t.^(0:N); t_ = t_';
    r_ = zeros(1,3); 
    r_(1) = a*(t.^(0:length(a)))'; 
    r_(2) = b*(t.^(0:length(b)))'; 
    r_(3) = c*t;
end
    function q_ = q(t)
        q_ = r(t) - r_initial(t);
    end
function u_ = u(t)
    % u is a unit vector. 
    % u is in the plane perp to dr. 
    
    %N = min([length(alpha), length(beta), length(gamma)]);
    %t_ = t.^(0:N); t_ = t_';
    %{
    u_ = zeros(1,3); 
    u_(1) = alpha*(t.^(0:length(alpha)))'; 
    u_(2) = beta*(t.^(0:length(beta)))'; 
    u_(3) = gamma*(t.^(0:length(gamma)))';
    if norm(u_)
        u_ = u_/norm(u_);
    end
    %}
    theta = alpha*(t.^(0:length(alpha)))'; 
    % theta is projected onto the xy (axial) plane from the *y-axis*
    % (theta=0 is no torsion)
    %{
    phi = beta*(t.^(0:length(beta)))';
    u_ = [sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)];
    %}
    
    % u is perpendicular to dr. 
    dr_ = r(t) - r(t-dt); 
    dx = dr_(1); dy = dr_(2); dz = dr_(3);
    uax = dz/sqrt( ( dx*sin(theta) + dy*cos(theta) )^2 + dz^2 );
    u_ = [uax*sin(theta), uax*cos(theta), sqrt(1 - uax^2)];
    if ~norm(u_) | isnan(norm(u_)) | isinf(norm(u_))
        u_ = zeros(size(3));
    end
end

    function ds_ = ds(t)
        dr_ = r(t) - r(t-dt); 
        dx = dr_(1); dy = dr_(2); dz = dr_(3);
        ds_ = dt*sqrt( (dx/dt)^2 + (dy/dt)^2 + (dz/dt)^2 );
    end
    s = @(t) integral(ds, 0, t);
    %s_final = s(t_initial_topbound);
    t_final = fsolve(@(t) s(t)-s_final, s_final);

% internal torque and twist: 
dr = @(t) r(t) - r(t-dt);
ddr = @(t) dr(t) - dr(t-dt);
dddr = @(t) ddr(t) - ddr(t-dt);
cx = @(t) cross(dr(t), ddr(t));
dTau = @(t) (( cx(t) * dddr(t)' )/( norm(cx(t))^2 )); % internal torsion
%Tau = integral(dTau, t_initial, t_final); % internal torsion 

du = @(t) u(t) - u(t-dt);
dTw = @(t) ( cross(du(t), u(t)) * (dr(t))' )/ds(t);  
Tw = integral(dTw, t_initial, t_final); Tw = Tw/(2*pi); % Twist

dT = @(t) J(s(t))*G(s(t))*(dTw(t) - dTau(t));
T = integral(dT, t_initial, t_final); % internal torque 

% writhe: 

r_inferior = @(t) [a(2)*t + a(1), b(2)*t + b(1), c*t];
dr_inferior = [a(2), b(2), c];
r_superior = @(t) r(t_final) + (t-t_final)*dr(t_final)/dt;
dr_superior = dr(t_final);

    function ddWr_ = ddWr(t1, t2, r1, r2, dr1, dr2)
        r12 = r1(t1) - r2(t2);
        ddWr_ = (cross(dr1(t1), dr2(t2)) * r12')/( norm(r12)^3 );
    end

Wr_spine = integral2( @(t1,t2) ddWr(t1, t2, r, r, dr, dr), ...
    t_initial, t_final, t_initial, t_final);
Wr_spine = Wr_spine/(4*pi);

Wr_inferior = integral2( @(t1,t2) ddWr(t1,t2, r, r_inferior, dr, dr_inferior),...
    t_initial, t_final, -Inf, t_initial);
Wr_inferior = Wr_inferior/(2*pi);

Wr_superior = integral2( @(t1,t2) ddWr(t1,t2, r, r_superior, dr, dr_superior),...
    t_initial, t_final, t_final, Inf);
Wr_superior = Wr_superior/(2*pi);

Wr_dec = integral2( @(t1,t2) ddWr(t1,t2, ...
    r_inferior, r_superior, dr_inferior, dr_superior), ...
    -Inf, t_initial, t_final, Inf);
Wr_dec = Wr_dec/(2*pi);

F(1) = Tw + Wr_spine + Wr_inferior + Wr_superior + Wr_dec; % CWF
F(2) = T; % internal torque balance 

% bending moment 

dq = @(t) q(t) - q(t-dt);
ddq = @(t) dq(t) - dq(t-dt);
tF = fsolve(@(t) s(t)-sF, sF); qF = q(tF);
    function M_ = M(t)
        if t < tF
            M = cross(q(t), Fapp);
        else
            M = cross(qF, Fapp);
        end
    end
forcebend = @(t) cross(dq(t), ddq(t))/( norm(dq(t))^3 ) + M(t)/(I(s(t))*E(s(t)));

num_to_select = length(a) + length(b) - 3; 
num_rows_to_select = ceil(num_to_select/3);
t_select = linspace(t_initial, t_final, num_rows_to_select);
forcebendeqn = zeros(num_rows_to_select, 3);
idx = 1; 
for ti = t_select
    forcebendeqn(idx,:) = forcebend(ti);
    idx = idx+1;
end
forcebendeqn = forcebendeqn(:); forcebendeqn = forcebendeqn(1:num_to_select);

F = [F, forcebendeqn'];

end