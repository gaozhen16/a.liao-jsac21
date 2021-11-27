function Index_W_RF_BS_Angle = A0_Gene_Index_W_RF_1BSangle(N_BS_bar,I_BS_bar,N_BS_v_bar,N_BS_h_bar,I_BS_h_bar,N_BS_v,N_BS_h,Delta_BS)
% Generate the Index matrix W_RF for angle estimation
    
    Index_W_RF_BS_Angle = zeros(N_BS_bar,I_BS_bar);
    Subarray_used = ones(N_BS_h_bar,N_BS_v_bar);
    for ii_BS_bar = 1:I_BS_bar
        i_BS_v_bar = ceil(ii_BS_bar/I_BS_h_bar);
        i_BS_h_bar = ii_BS_bar-(i_BS_v_bar-1)*I_BS_h_bar;
        Array_full = zeros(N_BS_v,N_BS_h);
        Array_full((i_BS_h_bar-1)*Delta_BS+1:(i_BS_h_bar-1)*Delta_BS+N_BS_h_bar,(i_BS_v_bar-1)*Delta_BS+1:(i_BS_v_bar-1)*Delta_BS+N_BS_v_bar)...
            = Subarray_used;
        Index_W_RF_BS_Angle(:,ii_BS_bar) = find(Array_full == 1);
    end
    
end