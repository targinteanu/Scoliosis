function A = PolyCosineMatrix(L, N)

M = N/2;
A = zeros(M+1);

for n = 0:2:N
    for m = 1:1:M
        row = n/2 + 1; col = m+1;
        
        S = 0;
        for nu = 0:2:n
            Sp = 1;
            for k = 0:(nu-1)
                Sp = Sp*(n-k);
            end
            S = S +( (-1)^(nu/2) * L^(n-nu) * Sp / m^(nu+1) );
        end
        
        C = 0;
        for nu = 1:2:n
            Cp = 1; 
            for k = 0:(nu-1)
                Cp = Cp*(n-k);
            end
            C = C +( (-1)^((nu+1)/2) * L^(n-nu) * Cp / m^(nu+1) );
        end
        
        A(row,col) = (2/L) *( S*sin(m*L) + C*cos(m*L) );
    end
end

n0 = 0:2:N; A0 = (L.^n0) ./ (n0 + 1);
A(:,1) = A0;

end