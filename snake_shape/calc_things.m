function ret = calc_things(t_set,x_val,info,fxns)
bd_mid = info.midBounds; added_L = info.addL; r = info.r; T = 2*pi*r;
bd_mid(end+1) = info.L; %append L to end so min works
T_adj = info.T_adj;
ret = zeros(2,length(t_set));
if info.num_rows > 2
    d_max = bd_mid(4) - bd_mid(2);
end

idx_min = min(nonzeros((1:length(bd_mid)).*(x_val<=bd_mid)));
idx_min_mod4 = mod(idx_min,4);
U_num = mod(idx_min,2);
d_fin = bd_mid(end) - bd_mid(end-1);

y_bd = @(t) fxns.y(t,pi/2);

%idx_min --> even: cyc (2 mod 4 means spin up)
%idx_min --> odd: d case (3 mod 4 means spin up)
if idx_min == 1
    %trivial case
    ret(1,:) = x_val*ones(1,length(t_set));
    ret(2,:) = zeros(1,length(t_set));
elseif idx_min_mod4 == 0
    %cyc spin down case
    x_bd = @(t) fxns.x(t,pi/2) + added_L;
    d = x_val - bd_mid(idx_min-2);
    x_adj = x_val - 2*T_adj; p_val = x_adj / r;
    x_t = @(t) fxns.x(t,p_val);
    y_t = @(t) fxns.y(t,p_val);
    ret(1,:) = -x_t(t_set) + r*(2*pi + p_val) + d/r * fxns.np_t(t_set) + x_bd(t_set);
    ret(2,:) = y_t(t_set) + d/r * fxns.mp_t(t_set) + y_bd(t_set);
elseif idx_min_mod4 == 1
    %d spin down case
    d = x_val - bd_mid(idx_min-1);
    ret(1,:) = d/r * fxns.nnp_t(t_set) + r*5*pi/2 + T_adj/r * fxns.np_t(t_set) + added_L;
    ret(2,:) = d/r * fxns.mnp_t(t_set) + T_adj/r * fxns.mp_t(t_set) + 2*y_bd(t_set);
    %todo: add support (check?) for arbitrary L
elseif idx_min_mod4 == 2
    %cyc spin up case- todo: mult rows??
    x_adj = x_val - T_adj ; p_val = x_adj / r;
    x_t = @(t) fxns.x(t,p_val) + added_L;
    y_t = @(t) fxns.y(t,p_val);
    ret(1,:) = x_t(t_set);
    ret(2,:) = y_t(t_set);
elseif idx_min_mod4 == 3
    %d spin up case- todo: mult rows???
    x_bd = @(t) fxns.x(t,pi/2) + added_L; y_bd = @(t) fxns.y(t,pi/2);
    d = x_val - bd_mid(idx_min-1);
    ret(1,:) = d/r * fxns.np_t(t_set) + x_bd(t_set);
    ret(2,:) = d/r * fxns.mp_t(t_set) + y_bd(t_set);
end

end