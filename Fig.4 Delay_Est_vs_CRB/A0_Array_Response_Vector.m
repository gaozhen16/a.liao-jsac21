function a_vec_fianl = A0_Array_Response_Vector(N_h, N_v, Azi, Ele, k, K, fs, fc)
% Generate array response vector of UPA

    n_h = (0:1:N_h-1)'; n_v = (0:1:N_v-1)';
    if nargin == 4
        a_vec_h = exp(1i*pi*n_h*sin(Azi)*cos(Ele)); % /sqrt(N_h)
        a_vec_v = exp(1i*pi*n_v*sin(Ele)); % /sqrt(N_v)
        a_vec_fianl = kron(a_vec_v,a_vec_h);
        
    elseif nargin == 8
        a_vec_h_k = exp(1i*pi*(1+((k-1)/K-0.5)*fs/fc)*n_h*sin(Azi)*cos(Ele)); % /sqrt(N_h)
        a_vec_v_k = exp(1i*pi*(1+((k-1)/K-0.5)*fs/fc)*n_v*sin(Ele)); % /sqrt(N_v)
        a_vec_fianl = kron(a_vec_v_k,a_vec_h_k);
        
    else
        error('Error: invalid nargin!');
    end
    
end