function Diff_A_x = Diff_Para_tau(x_value, k_index)
	
	Lp = length(x_value);
	syms x_diff;
    a_x = exp(1i*(k_index-1)*x_diff);
    diff_a_x = diff(a_x);
    Diff_A_x = zeros(length(k_index), Lp);
    for lp = 1:Lp
        x_lp = x_value(lp);
        Diff_A_x(:,lp) = double(subs(diff_a_x,x_diff,x_lp));
    end
end