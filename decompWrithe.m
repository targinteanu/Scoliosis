function Writhe = decompWrithe(cm, range)
% gets Writhe of a curve whose points are given by cm, using approximation
% derived by Levitt (https://doi.org/10.1016%2Fs0022-2836%2883%2980129-6)
% Coordinates are rows of cm and sp.
% range: an optional input vector that specifies which rows of cm to
% use. If unspecified, Writhe will be evaluated from the start to end of cm.

if nargin == 1
    range = 1:length(cm(:,1));
end

%Writhe = 0;
%for i = range(2:(end-2))
%    for j = (i+2):range(end)
I = range(2); 
J = range(end);
        
        p1 = cm(I-1,:); p2 = cm(I,:); 
        p3 = cm(J-1,:); p4 = cm(J,:);
        u = p1-p2; u=u/norm(u); v = p4-p3; v = v/norm(v);
        
%        u = p1-p5; v = p3-p6; r13 = p3-p1;
%        ang = @(a,b) acos(a*b'/(norm(a)*norm(b)));
        
r_inferior = @(t) u*t + p1;
dr_inferior = @(t) u;
r_superior = @(t) v*t + p4;
dr_superior = @(t) v;

    function ddWr_ = ddWr(t1, t2, r1, r2, dr1, dr2)
        t1 = t1'; t2 = t2';
%        r12 = r1(t1) - r2(t2);
        if length(t1) <= length(t2)
            sz = size(t1); imax = length(t1);
        else
            sz = size(t2); imax = length(t2);
        end
        ddWr_ = zeros(sz);
        for i = 1:imax
            tt1 = t1(i); tt2 = t2(i); %rr12 = r12(i,:);
            rr12 = r1(tt1) - r2(tt2);
            if norm(rr12)
                ddWr_(i) = (cross(dr1(tt1), dr2(tt2)) * rr12')/( norm(rr12)^3 );
            else
                ddWr_(i) = 0;
            end
        end
        ddWr_ = ddWr_';
        %figure; plot3(t1, t2, ddWr_, '.'); grid on; xlabel('t1'); ylabel('t2');
    end

Wr_dec = integral2( @(t1,t2) ddWr(t1,t2, ...
    r_inferior, r_superior, dr_inferior, dr_superior), ...
    -Inf, 0, 0, Inf, 'Method', 'iterated');
Writhe = Wr_dec/(2*pi);

end