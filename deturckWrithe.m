function Writhe = deturckWrithe(cm, ep, range)
% gets Writhe of a chain of segments with centers given by cm and endpoints 
% given by ep, using approximation derived by DeTurck. 
% Coordinates are rows of cm and ep.
% range: an optional input vector that specifies which rows of cm to
% use. If unspecified, Writhe will be evaluated from the start to end.
% cm is N-by-3. 
% ep is 2N-by-3, where each pair of rows are the two endpoints surrounding
% one row of cm. 

if nargin < 3
    range = 1:length(cm(:,1));
end

Writhe = 0; 
for i = range
    dr1 = ep(2*i, :) - ep(2*i - 1, :);
    for j = (i+1):range(end)
        dr2 = ep(2*j, :) - ep(2*j - 1, :);
        dr12 = cm(i,:) - cm(j,:);
        
        if norm(dr12)
            dW = ( cross(dr1, dr2)*dr12' )/( norm(dr12)^3 );
            Writhe = Writhe + dW;
        end
    end
end
Writhe = Writhe/(2*pi);

end