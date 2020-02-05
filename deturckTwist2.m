function Twist = deturckTwist2(cm, dir, range)

if nargin == 2
    range = 2:size(cm,1);
end

if sum(size(dir)==1) % sp is xy-plane angle in deg
    dir = dir * pi / 180; % deg to rad 
    dir = [cos(dir), sin(dir), zeros(size(dir))];
end

Twist = 0;

DX = (cm(2:end,:) - cm(1:(end-1),:));
U = .5 * (dir(2:end,:) + dir(1:(end-1),:)); 
U = [0, 0, 0; U]; DX = [0, 0, 0; DX];
for i = range
    u1 = U(i,:); 
    u1 = u1 - proj(u1, DX(i,:)); 
    u1 = u1/norm(u1); 
    U(i,:) = u1;
end

for i = range
    
    dX = DX(i,:);
    u2 = U(i,:); u1 = U(i-1,:);
    
    theta = acos( (u1*u2') );
    
    du = u2-u1;
    sgn = cross(dX, u2) * du';
    
    Twist = Twist + theta * sign(sgn);
    
end

Twist = Twist/(2*pi);

function proj_ = proj(u, x)
    % project u onto x
    denom = (x*x');
    if denom
        proj_ = ((u*x')/denom)*x;
    else
        proj_ = 0*x;
    end
end

end