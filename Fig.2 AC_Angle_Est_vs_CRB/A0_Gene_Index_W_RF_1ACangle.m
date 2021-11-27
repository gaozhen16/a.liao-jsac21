function Index_W_RF_AC_Angle = A0_Gene_Index_W_RF_1ACangle(M_AC_bar,I_AC_bar,I_AC,I_AC_h,I_AC_h_bar,M_AC_h_bar,M_AC_v_bar,M_AC_h,M_AC_v,Delta_AC,N_AC_h,N_AC_v)
% Generate the Index matrix W_RF for angle estimation
    
    Index_W_RF_AC_Angle = zeros(M_AC_bar,I_AC_bar,I_AC);
    Subarray_used = ones(M_AC_h_bar,M_AC_v_bar);
    for ii_AC = 1:I_AC
        i_AC_v_rf = ceil(ii_AC/I_AC_h);
        i_AC_h_rf = ii_AC-(i_AC_v_rf-1)*I_AC_h;
        for ii_AC_bar = 1:I_AC_bar
            i_AC_v_bar = ceil(ii_AC_bar/I_AC_h_bar);
            i_AC_h_bar = ii_AC_bar-(i_AC_v_bar-1)*I_AC_h_bar;
            Subarray_full = zeros(M_AC_h,M_AC_v);
            Subarray_full((i_AC_h_bar-1)*Delta_AC+1:(i_AC_h_bar-1)*Delta_AC+M_AC_h_bar,(i_AC_v_bar-1)*Delta_AC+1:(i_AC_v_bar-1)*Delta_AC+M_AC_v_bar)...
                = Subarray_used;
            Array_full = zeros(N_AC_h,N_AC_v);
            Array_full((i_AC_h_rf-1)*M_AC_h+1:i_AC_h_rf*M_AC_h,(i_AC_v_rf-1)*M_AC_v+1:i_AC_v_rf*M_AC_v) = Subarray_full;
            Index_W_RF_AC_Angle(:,ii_AC_bar,ii_AC) = find(Array_full == 1);
        end
    end
end