r_undeformed = @(t) [0, .00001643*t*(t-18)*(t-46)*(t-58), t];
b_undeformed = [-.789, .0746, -.002, .0000164]; % starting at n=1; b0=0. 
% 0 = a0 = b0 fixes the origin as the inferior-most point of the spine
% for the spine to be smooth, dr/dt = dr_undeformed/dt at t=0 
% -> [0, -.789, 1] = [a1, b1, c]
% and 0 = theta = dtheta/dt at t=0 -> 0 = alpha0 = alpha1
a0=0; b0=0; 
a1=0; b1=-.78903432; c=1;
alpha0=0; alpha1=0;

s_final = 59.795979; 

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

soln = fsolve(callfun, ...
    [0,0,0,0,0,0,0,0,0, .0746, -.002, .0000164,0,0,0,0,0,0, 0,0,0]);

a = [a0, a1, soln(1:9)];
b = [b0, b1, soln(10:18)];
alpha = [alpha0, alpha1, soln(19:21)];

function r_ = r(t)
    %N = min([length(a), length(b)]);
    %t_ = t.^(0:N); t_ = t_';
    r_ = zeros(1,3); 
    r_(1) = a*(t.^(0:length(a)))'; 
    r_(2) = b*(t.^(0:length(b)))'; 
    r_(3) = c*t;
end
function u_ = u(t)
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

function F = callfun(x)
    a = [a0, a1, x(1:9)]; 
    b = [b0, b1, x(10:18)];
    alpha = [alpha0, alpha1, x(19:21)];
    F = forcepull2(a, b, c, alpha, mech, app);
end