function Jacobian_Mat = Gene_Jacobian_Mat_Azi_Ele(miu, niu, Delta, lambda, d_ant, L)

    Dim = length(miu)/L;
    
    Jacobian_Mat = zeros(2*Dim,2*Dim,L);
    for ll_Jac = 1:L
        miu_ll = miu(ll_Jac);
        niu_ll = niu(ll_Jac);
        syms miu_diff niu_diff
        g_1_niu = asin(lambda*niu_diff/(2*pi*d_ant*Delta)); % d_ant = lambda/2;
        diff_g_1_niu = diff(g_1_niu,niu_diff,1);
        phi_MS = asin(lambda*niu_ll/(2*pi*d_ant*Delta));
        g_2_miu = asin(lambda*miu_diff./(2*pi*d_ant*Delta*cos(phi_MS)));
        diff_g_2_miu = diff(g_2_miu,miu_diff,1);
        for ii = 1:2
            if ii == 1
                for nn = 1:Dim
                    Jacobian_Mat(nn,nn,ll_Jac) = double(subs(diff_g_1_niu,niu_diff,niu_ll(nn)));
                end
            else
                for mm = 1:Dim
                    Jacobian_Mat(Dim+mm,Dim+mm,ll_Jac) = double(subs(diff_g_2_miu(mm),miu_diff,miu_ll(mm)));
                end
            end
        end
    end
    
end