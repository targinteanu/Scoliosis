function Writhe = decompWritheLevitt(cm, range)
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
i = range(2); 
j = range(end);
        
        p1 = cm(i-1,:); p5 = cm(i,:); 
        p6 = cm(j-1,:); p3 = cm(j,:);
        
        u = p1-p5; v = p3-p6; r13 = p3-p1;
        ang = @(a,b) acos(a*b'/(norm(a)*norm(b)));
        
        angles = zeros(1,3);
        angles(1) = ang(cross(u,v), cross(r13,v)); 
        angles(2) = ang(cross(u,v), cross(r13,u)); 
        angles(3) = ang(cross(r13,v), cross(r13,u));

%        Omega = sum(angles) - 2*pi;
%        Omega = Omega*sign( (cross(r34,r12)) *r13');
        
%        Writhe = Writhe + Omega;
%    end
%end
Writhe = sum(angles)*sign( cross(v,u)*r13' )/(2*pi);

%{
if abs(imag(Writhe)) > 9e-8
    disp('Error: non-real Writhe')
else
    Writhe = real(Writhe);
end
%}

end