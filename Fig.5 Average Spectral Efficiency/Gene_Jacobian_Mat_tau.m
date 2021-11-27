function Jacobian_Mat_tau = Gene_Jacobian_Mat_tau(miu_tau, K, fs, D_subc, L)

    Dim = length(miu_tau)/L;
    
    Jacobian_Mat_tau = zeros(Dim,Dim,L);
    for ll_Jac = 1:L
        miu_tau_ll = miu_tau(ll_Jac);
        syms miu_tau_diff;
        g_miu_tau = - K*miu_tau_diff/(2*pi*D_subc*fs);	% fs = 1 --> tau*fs = Normalized delay
        diff_g_miu_tau = diff(g_miu_tau,miu_tau_diff,1);
        for dd = 1:Dim
            Jacobian_Mat_tau(dd,dd,ll_Jac) = double(subs(diff_g_miu_tau,miu_tau_diff,miu_tau_ll(dd)));
        end
    end
    
end