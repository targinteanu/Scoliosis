function Twist = deturckTwist(cm, ep, sp, range)
% gets Twist of chain of segments whose centers are given by cm and 
% endpoints by ep, with a vector field whose coordinates are given by sp, 
% using approximation derived by DeTurck. 
% It is assumed the vector field is constant along the segments, and
% nonzero Twist occurs only between segments. 
% Coordinates are rows of cm, ep, and sp.
% Each vector in sp corresponds to the point in the same row of cm. 
% range: an optional input vector that specifies which rows of cm and sp to
% use. If unspecified, Twist will be evaluated from the start to end.
% cm and sp are N-by-3. 
% ep is 2N-by-3, where each pair of rows are the two endpoints surrounding
% one row of cm. 

if nargin < 4
    range = 1:length(cm(:,1));
end
range = range(1:(end-1));

U = sp-cm; 
dX = ep(2:end,:) - ep(1:end-1,:);
dX1 = dX(1:2:end, :); dX2 = dX(2:2:end);

% make U perp dX
Uparl = diag((U*dX1')./(dX1*dX1')) .* dX1;
U = U - Uparl;
% make U unit
normU = sqrt(diag(U*U'));
U = U./normU;

% average U and change in U over the gap between endpoints
Ugap = U(2:end,:) - U(1:end-1,:); 
normUgap = sqrt(diag(Ugap*Ugap'));
Ugap = Ugap./normUgap;
Ugap = Ugap(1:2:end, :);
dU = U(2:end,:) - U(1:end-1,:);

Twist = 0; 
for s = range
    ds = norm(dX2(s,:));
    dT = cross((dU(s,:)/ds), Ugap(s,:)) * (dX2(s,:)/ds)' * ds;
    Twist = Twist + dT; 
end
Twist = Twist/(2*pi);

end