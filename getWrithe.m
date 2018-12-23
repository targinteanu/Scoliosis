function Writhe = getWrithe(cm, range)

if nargin == 1
    range = 2:length(cm(:,1));
end

Writhe = 0;
for i = range
    for j = range
        dr1 = cm(i,:) - cm((i-1),:);
        dr2 = cm(j,:) - cm((j-1),:);
        dr12 = cm(i,:) - cm(j,:);
        if norm(dr12)
            dW = ( cross(dr1, dr2)*dr12' )/( norm(dr12)^3 );
            Writhe = Writhe + dW;
        end
    end
end
Writhe = Writhe/(4*pi);

end