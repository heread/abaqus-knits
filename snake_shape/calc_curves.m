function ret = calc_curves(t_set, x_set,r,L,num_rows)
p_len = length(x_set); t_len = length(t_set);
if p_len > 1 && t_len > 1
    error("one of these must be a single value")
end
T = 2*pi*r; %period for one cycloid
charWholeLen = num_rows*T; %natural L for (num_rows - 1) turn
% L <==> charWholeLen if nothing weird is going on
added_L = (L - charWholeLen)/num_rows;
T_adj = T + added_L;

a = @(p) 3*pi/2 - p; %p \in [-pi/2, pi/2]
g = @(t,p) t*(a(p) - 2*pi) + 2*pi;
%p = -pi/2: [0,1] --> [2pi, 2pi]
%p = pi/2: [0,1] --> [2pi, pi]

x = @(t,p) r*(p + g(t,p) - sin(g(t,p)));
%todo: start staring here
%r\cdot (- (g (t )-p-\ \sin (g (t ) )-\pi ),1-\cos (g (t ) )+2 )
y = @(t,p) r*(1 - cos(g(t,p))); 
n_p = @(t,p) r*(t.*cos(t*(p + pi/2)) - t + 1); np_t = @(t) n_p(t,pi/2);
m_p = @(t,p) r*t.*sin(t*(p + pi/2)); mp_t = @(t) m_p(t,pi/2);

nn_p = @(t,p) r + r.*(t.*cos(pi.*t) - t + 1) - r.*(t.*cos(t.*(p + pi/2)) - t + 1);
nnp_t = @(t) nn_p(t,pi/2);
mn_p = @(t,p) r.*t.*sin(t.*(p + pi/2)) + r*t.*sin(pi*t);
mnp_t = @(t) mn_p(t,pi/2);

%repmat([1,2,3],1,2)
circ_angles = [-pi/2, pi/2];
for i = 1:(num_rows-1)
    idx1 = 2*(i-1) + 1; idx2 = idx1 + 1;
    bd_mid(idx1:idx2) = i*T_adj + r*circ_angles;
end
%bd_mid = T + r*circ_angles + added_L; %[left bd, right bd]

fxns = struct; fxns.x = x; fxns.y = y; fxns.np_t = np_t; fxns.mp_t = mp_t;
fxns.nnp_t = nnp_t; fxns.mnp_t = mnp_t;
info = struct; info.addL = added_L; info.midBounds = bd_mid; info.r = r;
info.T_adj = T_adj; info.L = L; info.num_rows = num_rows;

%BLUE LINES- all t handled together
if p_len == 1
    ret = calc_things(t_set,x_set,info,fxns);
else %PINK CASE- each x handled by itself
    ret = zeros(2,p_len);
    for i=1:length(x_set)
        ret(:,i) = calc_things(t_set,x_set(i),info,fxns);
    end
end

end

