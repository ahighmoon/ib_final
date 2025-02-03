function U = particle_interp(u, P)

global h N;
np = size(P, 1);
U = zeros(np, 2);
s = P/h;          
i = floor(s);     
r = s - i;        

for k = 1:np
    w_temp = vec_phi1(r(k,1)) .* vec_phi2(r(k,2));
    w = squeeze(permute(w_temp, [1, 3, 2]));
    i1 = mod((i(k,1)-1):(i(k,1)+2), N) + 1;
    i2 = mod((i(k,2)-1):(i(k,2)+2), N) + 1;
    U(k,1) = sum(sum( w .* u(i1, i2, 1) ));
    U(k,2) = sum(sum( w .* u(i1, i2, 2) ));
end
end
