function Writhe = levittWritheOld(cm, range)
% gets Writhe of a curve whose points are given by cm, using approximation
% derived by Levitt (https://doi.org/10.1016%2Fs0022-2836%2883%2980129-6)
% Coordinates are rows of cm and sp.
% range: an optional input vector that specifies which rows of cm to
% use. If unspecified, Writhe will be evaluated from the start to end of cm.

if nargin == 1
    range = 2:length(cm(:,1));
end

range = range(2:end);

Writhe = 0;
for i = range
%    for j = range
    for j = range(1):(i-1)
        
        p1 = cm(i-1,:); p2 = cm(i,:); 
        p3 = cm(j-1,:); p4 = cm(j,:);
        
        r34 = p4-p3; r12 = p2-p1;
        r13 = p3-p1; r14 = p4-p1; r23 = p3-p2; r24 = p4-p2;
        
        n = zeros(4, 3);
        n(1,:) = cross(r13, r14);
        n(2,:) = cross(r14, r24);
        n(3,:) = cross(r24, r23);
        n(4,:) = cross(r23, r13);
        
        rows_to_remove = zeros(4, 1);
        for k = 1:4
            if norm(n(k,:))
                % normalize 
                n(k,:) = n(k,:)/norm(n(k,:));
            else
                % mark for removal 
                rows_to_remove(k) = k;
            end
        end
        n = n([1:(min(rows_to_remove)-1), (max(rows_to_remove)+1):end], :);
        
        N = n([2:end 1],:);
        as = diag(n*N'); angles = asin(as) + pi/2;

        Omega = sum(angles) - 2*pi;
        Omega = Omega*sign( (cross(r34,r12)) *r13');
        
        Writhe = Writhe + Omega;
    end
end
%Writhe = Writhe/(4*pi);
Writhe = Writhe/(2*pi);
if abs(imag(Writhe)) > 9e-8
    disp('Error: non-real Writhe')
else
    Writhe = real(Writhe);
end

end