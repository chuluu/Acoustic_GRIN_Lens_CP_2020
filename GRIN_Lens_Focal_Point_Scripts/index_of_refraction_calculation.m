function n_y = index_of_refraction_calculation(y,a,n_o)
    n_y = n_o.*sech(a*y);
end