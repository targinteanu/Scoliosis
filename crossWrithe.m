function [Writhe, Wrx] = crossWrithe(cm, range)
% gets Writhe of a curve whose points are given by cm, using approximation
% derived by Levitt (https://doi.org/10.1016%2Fs0022-2836%2883%2980129-6)
% Coordinates are rows of cm and sp.
% range: an optional input vector that specifies which rows of cm to
% use. If unspecified, Writhe will be evaluated from the start to end of cm.

if nargin == 1
    range = 1:length(cm(:,1));
end

ang = @(a,b) acos(a*b'/(norm(a)*norm(b)));

U = zeros(2,3); p1 = zeros(2,3); Wrx = zeros(2,1);
i = range(2); 
p1(1,:) = cm(i,:); p5 = cm(i-1,:); U(1,:) = p5-p1(1,:);
i = range(end); 
p1(2,:) = cm(i,:); p6 = cm(i-1,:); U(2,:) = p1(2,:)-p6;
        

for i = 1:2
    for j = (3 + (i==1)):(range(end) - (i==2))
        
        p3 = cm(j-1,:); p4 = cm(j,:);
        v = p4-p3; r13 = p3-p1(i,:);
        u = U(i,:);
        
        angles = zeros(1,4); 
        angles(1) = ang( cross(r13,v), cross(u,v)-cross(r13,u) );
        angles(2) = ang( cross(u,v)-cross(r13,u), cross(u,v) );
        angles(3) = ang( cross(u,v), cross(r13,u) );
        angles(4) = ang( cross(r13,u), cross(r13,v) );
        
        Omega = sum(angles);
        Omega = Omega*sign( (cross(v,u)) *r13');
        
        Wrx(i) = Wrx(i) + Omega;
    end
end
Wrx = Wrx/(4*pi);
Writhe = sum(Wrx);
if abs(imag(Writhe)) > 9e-8
    disp('Error: non-real Writhe')
else
    Writhe = real(Writhe);
end

end