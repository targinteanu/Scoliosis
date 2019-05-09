range = 2:length(cm);
%range = 1:length(cm);
M = zeros(length(cm));

for i = range
    for j = range
%for i = range(2:(end-2))
%    for j = (i+2):range(end)
        
        if abs(i-j) > 1
        
        p1 = cm(i-1,:); p2 = cm(i,:); 
        p3 = cm(j-1,:); p4 = cm(j,:);
        
        r34 = p4-p3; r12 = p2-p1;
        r13 = p3-p1; r14 = p4-p1; r23 = p3-p2; r24 = p4-p2;
        
        n = zeros(4, 3);
        n(1,:) = cross(r13, r14);
        n(2,:) = cross(r14, r24);
        n(3,:) = cross(r24, r23);
        n(4,:) = cross(r23, r13);
        
        for k = 1:4
            if norm(n(k,:))
                % normalize 
                n(k,:) = n(k,:)/norm(n(k,:));
            end
        end
        
        N = n([2:end 1],:);
        as = diag(n*N'); angles = asin(as) + pi/2;

        Omega = sum(angles) - 2*pi;
        Omega = Omega*sign( (cross(r34,r12)) *r13');
        
        M(i,j) = Omega;
        end
    end
end

%zerolim = 1e-5; M = M.*(abs(M) > zerolim);