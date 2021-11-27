function Diff_A_x = Diff_Para_Azi_Ele(N_dim_sub, x_value)
    n_dim = (0:(N_dim_sub-1)).';
	Lp = length(x_value);
	syms x_diff;
    a_x = exp(1i*n_dim*x_diff);
    diff_a_x = diff(a_x);
    Diff_A_x = zeros(N_dim_sub, Lp);
    for lp = 1:Lp
        x_lp = x_value(lp);
        Diff_A_x(:,lp) = double(subs(diff_a_x,x_diff,x_lp));
    end
end