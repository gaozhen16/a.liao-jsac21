function radial_vt_BS_set = Calculate_Radial_Velocity(vt_vec, Coo_Air_C_init, Coo_BS1_A, Coo_BS2_B, L)
% Tanslation coordinates system to generate Azimuth/Elevation angles at the BS1, BS2, and UPA (fixed at Aircraft)
% according to real coordinate of G-point

    vt = norm(vt_vec); radial_vt_BS_set = zeros(L,1);
    %% Calculate azimuth & elevation angles at BS1 under coordinate system E 
    Coo_BS1_A_to_C = Coo_Air_C_init - Coo_BS1_A;	% A --> C
    angle_BS1 = acos(dot(Coo_BS1_A_to_C,vt_vec)/(norm(Coo_BS1_A_to_C)*norm(vt_vec)));
    if angle_BS1 >= pi/2
        symbol_BS1 = 1;	% 夹角≥pi/2，说明是靠近的，速度为正 --> 多普勒为正
    else
        symbol_BS1 = -1;
    end
    radial_vt_BS_set(1) = symbol_BS1*abs(vt*cos(angle_BS1));	% radial velocity vt_BS1
    
    %% Calculate azimuth & elevation angles at BS2 under coordinate system F
    Coo_BS2_B_to_C = Coo_Air_C_init - Coo_BS2_B;	% B --> C
    angle_BS2 = acos(dot(Coo_BS2_B_to_C,vt_vec)/(norm(Coo_BS2_B_to_C)*norm(vt_vec)));
    if angle_BS2 >= pi/2
        symbol_BS2 = 1;	% 夹角≥pi/2，说明是靠近的，速度为正 --> 多普勒为正
    else
        symbol_BS2 = -1;
    end
    radial_vt_BS_set(2) = symbol_BS2*abs(vt*cos(angle_BS2));	% radial velocity vt_BS2
    
end