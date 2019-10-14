function R = CenterOfMass(img, lbound, rbound, ubound, dbound)

if nargin == 1
    rbound = size(img); dbound = rbound(1); rbound = rbound(2);
    lbound = 1; ubound = 1;
end

totalmass=sum(sum(img(:,:)));

R = [0, 0];

    for n=ubound:dbound;
       for m=lbound:rbound;
           R(1)=R(1)+img(n,m)*n;
           R(2)=R(2)+img(n,m)*m;
       end
    end
    R=R/totalmass;

end