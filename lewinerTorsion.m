function [Torsion, ddr] = lewinerTorsion(p, t, q, w)
% estimates Torsion of a curve at a point (t) using method described by
% Lewiner et al (doi:10.1016). 
% p: matrix of x,y,z coordinates of all vertebrae 
% t: index of point/vertebra at which to evaluate torsion 
% q: width of window over which to estimate/curve fit
% w: weights corresponding to each vertebra; defaults to 1 for all 
%    vertebrae (unweighted) if unspecified

if nargin < 4
    w = ones(length(p), 1);
end

dl = p(2:end,:) - p(1:(end-1),:);
dL = zeros(length(dl), 1);
for i = 1:length(dL)
    dL(i) = norm(dl(i,:));
end
dL = [dL(1:(t-1)); 0; dL(t:end)];
L = zeros(size(dL));
for i = (t+1):length(L)
    L(i) = L(i-1) + dL(i);
end
for i = fliplr(1:(t-1))
    L(i) = L(i+1) + dL(i);
end

range = (t-q):(t+q); 
w = w(range); L = L(range);
x = p(range,1); y = p(range,2); z = p(range,3);

a1 = sum(w.*(L.^2)); a2 = .5*sum(w.*(L.^3)); a3 = .25*sum(w.*(L.^4));
a4 = (1/6)*sum(w.*(L.^4)); a5 = (1/12)*sum(w.*(L.^5)); a6 = (1/36)*sum(w.*(L.^6));
bx1 = sum(w.*L.*x); bx2 = .5*sum(w.*(L.^2).*x); bx3 = (1/6)*sum(w.*(L.^3).*x);
by1 = sum(w.*L.*y); by2 = .5*sum(w.*(L.^2).*y); by3 = (1/6)*sum(w.*(L.^3).*y);
bz1 = sum(w.*L.*z); bz2 = .5*sum(w.*(L.^2).*z); bz3 = (1/6)*sum(w.*(L.^3).*z);

a = [a1 a2 a4; a2 a3 a5; a4 a5 a6];
Dx = (a^-1)*[bx1; bx2; bx3];
Dy = (a^-1)*[by1; by2; by3];
Dz = (a^-1)*[bz1; bz2; bz3];

dr = [Dx(1), Dy(1), Dz(1)];
ddr = [Dx(2), Dy(2), Dz(2)];
dddr = [Dx(3), Dy(3), Dz(3)];

cx = cross(dr, ddr);
Torsion = -(cx * dddr')/(norm(cx)^2);

end