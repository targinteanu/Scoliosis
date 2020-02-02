function Twist = deturckTwist2(cm, dir, range)

if nargin == 2
    range = 2:size(cm,1);
end

if sum(size(dir)==1) % sp is xy-plane angle in deg
    dir = dir * pi / 180; % deg to rad 
    dir = [cos(dir), sin(dir), zeros(size(dir))];
end

Twist = 0;

for i = range
    
    dX = cm(i,:) - cm((i-1),:);
    u1 = dir((i-1),:); u2 = dir(i,:);
    
    u1 = u1 - proj(u1,dX); u2 = u2 - proj(u2,dX);
    u1 = u1/norm(u1); u2 = u2/norm(u2);
    theta = acos( (u1*u2') );
    
    du = u2-u1;
    sgn = cross(dX, u2) * du';
    
    Twist = Twist + theta * sign(sgn);
    theta * sign(sgn)
    
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