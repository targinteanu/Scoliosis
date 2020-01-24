%% input variables: ---------------------------------------------------
r = @(s) [cos(s), sin(s), 0*s];
dr = @(s) (1/sqrt(2))*[-sin(s), cos(s), 0*s];
u = @(s) [cos(s).*cos(s), cos(s).*sin(s), sin(s)];
du = @(s) [-cos(s).*sin(s)-sin(s).*cos(s), cos(s).*cos(s)-sin(s).*sin(s), cos(s)];
npts = 100;
Smin = 0; Smax = 2*pi;
% --------------------------------------------------------------------

%% operations

S = linspace(Smin, Smax, npts)';
vr = @(s) v(s) - proj(v(s), r(s)); ur = @(s) vr(s)/norm(vr(s));
dr_norm = @(s) dr(s)/norm(dr(s)); 

R = r(S); 
%VR = cell2mat( arrayfun(vr, S, 'UniformOutput', false) );
VR = u(S);

Twist_integral = integral( @(t) dTw(t, dr_norm, u, du), S(1), S(end)) / (2*pi)
Twist_estimate = getTwist(R, VR)
Writhe_integral = integral2( @(t1,t2) ddWr(t1,t2, r, r, dr, dr), ...
    S(1), S(end), S(1), S(end), 'Method', 'iterated') / (4*pi)
Writhe_estimate = levittWrithe(R)

figure; 
plot3dSpine(R, VR); view([-30,10]);
zmax = r(Smax); zmax = zmax(3);
%xlim([-2, 2]); ylim([-2, 2]); zlim([0, zmax]);

%% twist and writhe functions 

function proj_ = proj(u, x)
    denom = (x*x');
    if denom
        proj_ = ((u*x')/denom)*x;
    else
        proj_ = 0*x;
    end
end

    function dTw_ = dTw(t, dr, u, du)
    % u, du must be entered s.t. u perp dr
        t = t';
        dTw_ = zeros(size(t)); imax = length(t);
        for i = 1:imax
            tt = t(i);
            if (i == 1) 
                if (i+1 <= imax)
                    dtt = t(i+1)-t(i);
                else
                    dtt = .0001;
                end
            else
                dtt = t(i)-t(i-1);
            end
            %du = (u(tt + dtt) - u(tt))/dtt;
            dTw_(i) = cross(du(tt), u(tt)) * dr(tt)';
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