function Jacobian_Mat_tau = Gene_Jacobian_Mat_Doppler(niu_Doppler, T_sym, L)

    Dim = length(niu_Doppler)/L;
    
    Jacobian_Mat_tau = zeros(Dim,Dim,L);
    for ll_Jac = 1:L
        niu_Doppler_ll = niu_Doppler(ll_Jac);
        syms niu_Doppler_diff;
        g_niu_Doppler = niu_Doppler_diff/(2*pi*T_sym);
        diff_g_niu_Doppler = diff(g_niu_Doppler,niu_Doppler_diff,1);
        for dd = 1:Dim
            Jacobian_Mat_tau(dd,dd,ll_Jac) = double(subs(diff_g_niu_Doppler,niu_Doppler_diff,niu_Doppler_ll(dd)));
        end
    end
    
end