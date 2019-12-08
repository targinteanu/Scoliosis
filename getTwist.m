function Twist = getTwist(cm, sp, range)
% gets Twist of a curve whose points are given by cm with a vector field
% whose coordinates are given by sp, using discrete approximations of
% integrals and derivatives. 
% Coordinates are rows of cm and sp.
% Each vector in sp corresponds to the point in the same row of cm. 
% range: an optional input vector that specifies which rows of cm and sp to
% use. If unspecified, Twist will be evaluated from the start to end of cm.

if nargin == 2
    range = 2:length(cm(:,1));
end

U = sp-cm; U = U(2:end,:);
dX = cm(2:end,:) - cm(1:end-1,:);
    % make U perp dX
    Uparl = diag((U*dX')./(dX*dX')) .* dX;
    U = U - Uparl;
    % make U unit
    normU = sqrt(diag(U*U'));
    U = U./normU; 
U = [0 0 0; U];
clear dX;

Twist = 0; 

% debugging ----------------------
tw = zeros(size(range)); i=1;
% --------------------------------

for s = range
    ds = norm(cm(s,:)-cm((s-1),:)); 
    dX = cm(s,:) - cm((s-1),:);
    
    %{
    U = sp([(s-1),s],:) - cm([(s-1),s],:); % 2 vectors
    % make U perp dX
    Uparl = ((U*dX')/(dX*dX')) * dX;
    U = U - Uparl;
    % make U unit 
    U1 = U(1,:)/norm(U(1,:));
    U2 = U(2,:)/norm(U(2,:));
    dU = U2 - U1;
    %}
    
    dU = U(s,:) - U((s-1),:);
    U2 = U(s,:);
    
    % twist 
    dT = cross((dU/ds), U2) * (dX/ds)' * ds;
    Twist = Twist + dT;   
    
    % debugging ----------------------
    tw(i) = dT; i = i + 1;
    % --------------------------------
end
Twist = Twist/(2*pi);

% debugging ----------------------
tw = tw/(2*pi);
figure; plot(tw); hold on; plot(cumsum(tw));
% --------------------------------

end