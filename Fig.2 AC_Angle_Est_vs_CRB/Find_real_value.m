function Real_Spatial_Fre = Find_real_value(Spatial_Fre_est, Spatial_Fre_real, Multi_d)
if Multi_d > 1
    if size(Spatial_Fre_est,1) == 1
        Spatial_Fre_est = Spatial_Fre_est.';
    end
    if size(Spatial_Fre_real,1) == 1
        Spatial_Fre_real = Spatial_Fre_real.';
    end
    N_Spatial_Fre = length(Spatial_Fre_real);
    
    Poss = (-Multi_d:Multi_d).';
    
    Real_Spatial_Fre = zeros(N_Spatial_Fre,1);
    for nn = 1:N_Spatial_Fre
        Spatial_Fre_nn = Spatial_Fre_est(nn) + Poss*pi/Multi_d;
        [~, Index] = min(abs(Spatial_Fre_nn - Spatial_Fre_real(nn)*ones(length(Poss),1)));
        Real_Spatial_Fre(nn) = Spatial_Fre_nn(Index);
    end
else
    Real_Spatial_Fre = Spatial_Fre_est;
end

end