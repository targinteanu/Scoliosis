function R = CenterOfMass3(img)
% Finds center of mass of a 3D volume matrix img where mass at each point
% is given by the corresponding element of img. Uses the centroid formula
% in 3 dimensions. 

    rbound = size(img); dbound = rbound(1); hbound = rbound(3); 
    rbound = rbound(2);
    lbound = 1; ubound = 1;

totalmass=sum(sum(sum(img)));

R = [0, 0, 0];
for l = 1:hbound
    for n=ubound:dbound
       for m=lbound:rbound
           R(1)=R(1)+img(n,m,l)*n;
           R(2)=R(2)+img(n,m,l)*m;
           R(3)=R(3)+img(n,m,l)*l;
       end
    end
end

R=R/totalmass;

end