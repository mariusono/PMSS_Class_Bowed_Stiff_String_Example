function [I_grid] = generate_interpolation_grid_1D(N,conn_point,order)

% % Inputs should be in ratio of the length of the element ! 

if strcmp(order,'linear')
    
    frac = floor(conn_point * N);
    alpha = conn_point * N - frac;

    I_grid = zeros(N,1);
    I_grid(frac) = 1-alpha;
    I_grid(frac+1) = alpha;
 
    
elseif strcmp(order,'nearest')
    
    frac = floor(conn_point * N);
    I_grid = zeros(N,1);
    I_grid(frac) = 1;
    
elseif strcmp(order,'cubic')

    conn_point_idx = floor(conn_point * N);
    I_grid = zeros(N,1);
    alph = conn_point * N - conn_point_idx;
    interPolVec = [(alph * (alph - 1) * (alph - 2)) / -6.0;
        ((alph - 1) * (alph + 1) * (alph - 2)) / 2.0;
        (alph * (alph + 1) * (alph - 2)) / -2.0;
        (alph * (alph + 1) * (alph - 1)) / 6.0];

    I_grid(conn_point_idx-1:conn_point_idx+2) = interPolVec;

end

end