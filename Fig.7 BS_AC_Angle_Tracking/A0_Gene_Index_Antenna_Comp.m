function Index_Ant_Comp = A0_Gene_Index_Antenna_Comp(N_Comp,M_Comp,N_Comp_h,N_Comp_v,M_Comp_h,N_h,N_v,I_AC,I_AC_h,M_AC_h,M_AC_v)
% Generate the Index matrix for antenna compensation
    
    if nargin < 8
        Index_Ant_Comp = zeros(N_Comp,M_Comp);
        Subarray_comp_BS = ones(N_Comp_h,N_Comp_v);
        for m_BS_comp = 1:M_Comp
            m_BS_Comp_v = ceil(m_BS_comp/M_Comp_h);
            m_BS_Comp_h = m_BS_comp - (m_BS_Comp_v-1)*M_Comp_h;
            Array_full_BS = zeros(N_h,N_v);
            Array_full_BS((m_BS_Comp_h-1)*N_Comp_h+1:m_BS_Comp_h*N_Comp_h,(m_BS_Comp_v-1)*N_Comp_v+1:m_BS_Comp_v*N_Comp_v)...
                = Subarray_comp_BS;
            Index_Ant_Comp(:,m_BS_comp) = find(Array_full_BS == 1);
        end
    else
        Index_Ant_Comp = zeros(N_Comp,M_Comp,I_AC);
        Subarray_comp_AC = ones(N_Comp_h,N_Comp_v);
        for ii_AC_comp = 1:I_AC
            i_AC_v_comp = ceil(ii_AC_comp/I_AC_h);
            i_AC_h_comp = ii_AC_comp-(i_AC_v_comp-1)*I_AC_h;
            for m_AC_comp = 1:M_Comp
                m_AC_Comp_v = ceil(m_AC_comp/M_Comp_h);
                m_AC_Comp_h = m_AC_comp - (m_AC_Comp_v-1)*M_Comp_h;
                Subarray_AC = zeros(M_AC_h,M_AC_v);
                Subarray_AC((m_AC_Comp_h-1)*N_Comp_h+1:m_AC_Comp_h*N_Comp_h,(m_AC_Comp_v-1)*N_Comp_v+1:m_AC_Comp_v*N_Comp_v)...
                    = Subarray_comp_AC;
                Array_full_AC = zeros(N_h,N_v);
                Array_full_AC((i_AC_h_comp-1)*M_AC_h+1:i_AC_h_comp*M_AC_h,(i_AC_v_comp-1)*M_AC_v+1:i_AC_v_comp*M_AC_v) = Subarray_AC;
                Index_Ant_Comp(:,m_AC_comp,ii_AC_comp) = find(Array_full_AC == 1);
            end
        end
    end
    
end