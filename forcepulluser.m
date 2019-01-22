r_undeformed = @(t) [0, .00001643*t*(t-18)*(t-46)*(t-58), t];
b_undeformed = [-.789, .0746, -.002, .0000164]; % starting at n=1; b0=0. 
% 0 = a0 = b0 fixes the origin as the inferior-most point of the spine
% for the spine to be smooth, dr/dt = dr_undeformed/dt at t=0 
% -> [0, -.789, 1] = [a1, b1, c]
% and 0 = theta = dtheta/dt at t=0 -> 0 = alpha0 = alpha1

global a0 b0 a1 b1 c alpha0 alpha1 dt discE discG VBE VBG %a b alpha %mech app
global count
count = 0;

a0=0; b0=0; 
a1=0; b1=-.78903432; c=1;
alpha0=0; alpha1=0;

s_final = 59.795979; % cm
dt = .1; % cm

discE = 1700; % Newton per sq cm
discG = 10.5; % Newton per sq cm
VBE = 37400; % Newton per sq cm
VBG = 500000; % Newton per sq cm

Rvert = @(s) 2.5 - 1.5*s/s_final; % vertebral body radius, cm
I = @(s) (pi/4)*(Rvert(s).^4); % assuming circular cross-section, cm^4
J = @(s) (pi/2)*(Rvert(s).^4); % assuming circular cross-section, cm^4

mech = struct; app = struct;

mech.r_initial = r_undeformed; mech.dt = dt; mech.s_final = s_final; 
mech.I = I; mech.J = J; mech.E = @(s) E(s); mech.G = @(s) G(s);

app.sF = 11.5; % center of L2, cm
app.Fapp = [50, 0, 0]; % N

%{
soln = fsolve( @([ ...
    a2,a3,a4,a5,a6,a7,a8,a9,a10, ...
    b2,b3,b4,b5,b6,b7,b8,b9,b10, ...
    alpha2,alpha3,alpha4]) ...
    forcepull2(...
    [a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10], ...
    [b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10], c, ...
    [alpha0,alpha1,alpha2,alpha3,alpha4], ...
    mech, app), ...
    [0,0,0,0,0,0,0,0,0, .0746, -.002, .0000164,0,0,0,0,0,0, ...
    0,0,0] );
%}

soln = fsolve(@(x) callfun(x, mech, app), ...
    [0,0,0,0,0,0,0,0,0, .0746, -.002, .0000164,0,0,0,0,0,0, 0,0,0]);
%    [.0746, 0, 0,0,0,0,0,0,0, .0746, -.002, .0000164,0,0,0,0,0,0, 0,0,0]);

a = [a0, a1, soln(1:9)];
b = [b0, b1, soln(10:18)];
alpha = [alpha0, alpha1, soln(19:21)];

[F, Torque, Tw, Wr_spine, Wr_inferior, Wr_superior, Wr_dec] = ...
    forcepull2(a, b, c, alpha, mech, app);
F

t = 0:dt:58; t=t'; xyz = zeros(length(t), 3); dir = xyz; 
T = zeros(size(t));
for idx = 1:length(t)
    xyz(idx,:) = r(t(idx));
    dir(idx,:) = u(t(idx));
    T(idx) = Torque(t(idx));
end
figure; plot3(xyz(:,1), xyz(:,2), xyz(:,3));
hold on; plot3(dir(:,1), dir(:,2), dir(:,3));
hold on; plot3(T, zeros(size(t)), t); 
legend('centerline', 'direction', 'ax torque');
title([num2str(Tw), ' | ', num2str(Wr_spine), ' | ', ...
    num2str(Wr_inferior), ' | ', num2str(Wr_superior), ' | ', ...
    num2str(Wr_dec)]);

function r_ = r(t)
    global a b c

    %N = min([length(a), length(b)]);
    %t_ = t.^(0:N); t_ = t_';
    r_ = zeros(1,3); 
    r_(1) = a*(t.^(0:length(a)))'; 
    r_(2) = b*(t.^(0:length(b)))'; 
    r_(3) = c*t;
end
function u_ = u(t)
    global alpha dt

    % u is a unit vector. 
    % u is in the plane perp to dr. 

    theta = alpha*(t.^(0:length(alpha)))'; 
    % theta is projected onto the xy (axial) plane from the *y-axis*
    % (theta=0 is no torsion)
    
    dr_ = r(t) - r(t-dt); 
    dx = dr_(1); dy = dr_(2); dz = dr_(3);
    uax = dz/sqrt( ( dx*sin(theta) + dy*cos(theta) )^2 + dz^2 );
    u_ = [uax*sin(theta), uax*cos(theta), sqrt(1 - uax^2)];
    if ~norm(u_) | isnan(norm(u_)) | isinf(norm(u_))
        u_ = zeros(size(3));
    end
end

function tf = inbone(s)
    if s <= 17
        % lumbar 
        tf = mod(s, 3.4) < 2.7;
    elseif s <= 48.2
        % thoracic 
        tf = mod(s, 2.6) < 1.9;
    elseif s <= 57.2
        % cervical below C2
        tf = mod(s, 1.7) < 1.2;
    else
        % C2 
        tf = true;
    end
end
function E_ = E(s)
    global VBE discE
    if inbone(s)
        E_ = VBE;
    else
        E_ = discE;
    end
end
function G_ = G(s)
    global VBG discG
    if inbone(s)
        G_ = VBG;
    else
        G_ = discG;
    end
end

function F = callfun(x, mech, app)
    global a0 a1 b0 b1 alpha0 alpha1 c
    A = [a0, a1, x(1:9)] 
    B = [b0, b1, x(10:18)]
    alpha = [alpha0, alpha1, x(19:21)];
    F = forcepull2(A, B, c, alpha, mech, app)
    global count 
    count = count+1
end