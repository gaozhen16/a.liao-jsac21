function Diff_A_x = Diff_Para_Doppler(x_value, nn)
	
	Lp = length(x_value);
	syms x_diff;
    a_x = exp(1i*(nn-1)*x_diff);
    diff_a_x = diff(a_x);
    Diff_A_x = zeros(length(nn), Lp);
    for lp = 1:Lp
        x_lp = x_value(lp);
        Diff_A_x(:,lp) = double(subs(diff_a_x,x_diff,x_lp));
    end
end