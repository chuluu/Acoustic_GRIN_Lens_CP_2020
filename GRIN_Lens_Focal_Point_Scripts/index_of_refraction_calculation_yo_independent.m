function n_y = index_of_refraction_calculation_yo_independent(y,a,n_o,B1,B2)
    z = a*y;
    g = z/(1 + (B1.*(z.^2)) + (B2.*(z.^4)));
    n_y = n_o.*sech(g);
end