function Index_W_RF_Delay = A0_Gene_Index_W_RF_2Sub_Array(N_AC_h,N_AC_v,M_AC,M_AC_h,M_AC_v,I_AC,I_AC_h)
% Generate the Index matrix W_RF for delay estimation
    
    Index_W_RF_Delay = zeros(M_AC,I_AC);
    for ii_AC = 1:I_AC
        i_AC_v_rf = ceil(ii_AC/I_AC_h);
        i_AC_h_rf = ii_AC-(i_AC_v_rf-1)*I_AC_h;
        Array_full = zeros(N_AC_h,N_AC_v);
        Array_full((i_AC_h_rf-1)*M_AC_h+1:i_AC_h_rf*M_AC_h,(i_AC_v_rf-1)*M_AC_v+1:i_AC_v_rf*M_AC_v) = ones(M_AC_h,M_AC_v);
        Index_W_RF_Delay(:,ii_AC) = find(Array_full == 1);
    end
end