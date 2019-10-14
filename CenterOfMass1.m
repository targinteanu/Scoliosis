function R = CenterOfMass1(img, lbound, rbound)

if sum(size(img)==1)

    if nargin == 1
        rbound = length(img); 
        lbound = 1; 
    end

    totalmass=sum(img(:));

    R = 0;

       for m=lbound:rbound
           R=R+img(m)*m;
       end
    
       if totalmass
            R=R/totalmass;
       end
       
else
    R = arrayfun(@(r) CenterOfMass1(img(r,:)), 1:size(img,1));
end

end