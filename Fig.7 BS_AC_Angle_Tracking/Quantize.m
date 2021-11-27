function Quantized_Matrix = Quantize(Matix, N_bits)
% Quantizing the analog precoding/conmbining matrix

    if N_bits == 0
        Quantized_Matrix = Matix;
    else
        N_Bits = 2^N_bits;
        N_T = size(Matix,1);
        Matix_quan_phase = floor((angle(Matix)+pi)*N_Bits/(2*pi)) *2*pi/N_Bits;% - pi;
        Quantized_Matrix = exp(1i*Matix_quan_phase);%/sqrt(N_T);
    end

end