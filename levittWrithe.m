function Writhe = levittWrithe(cm, range)

if nargin == 1
    range = 2:length(cm(:,1));
end

Writhe = 0;
for i = range
    for j = range
        r13 = cm(j-1,:) - cm(i-1,:);
        r14 = cm(j,:)   - cm(i-1,:);
        r23 = cm(j-1,:) - cm(i,:);
        r24 = cm(j,:)   - cm(i,:);
        
        r34 = cm(j,:)   - cm(j-1,:);
        r12 = cm(i,:)   - cm(i-1,:);
        
        n = zeros(4, 3);
        n(1,:) = cross(r13, r14);
        n(2,:) = cross(r14, r24);
        n(3,:) = cross(r24, r23);
        n(4,:) = cross(r23, r13);
        for k = 1:4
            if norm(n(k,:))
                n(k,:) = n(k,:)/norm(n(k,:));
            end
        end
        
        %{
        Omega = 0;
        for k = [1 2 3 4]
            for l = [2 3 4 1]
                Omega = Omega + asin(n(k,:) * n(l,:)');
            end
        end
        %}
        
        %%{
        N = n([2 3 4 1],:);
        as = diag(n*N');
        Omega = sum(asin(as));
        Omega = Omega*sign( (cross(r34,r12)) *r13');
        %}
        
        Writhe = Writhe + Omega;
    end
end
Writhe = Writhe/(4*pi);
if abs(imag(Writhe)) > 4e-9
    disp('Error: non-real Writhe')
%else
%    Writhe = real(Writhe);
end

end