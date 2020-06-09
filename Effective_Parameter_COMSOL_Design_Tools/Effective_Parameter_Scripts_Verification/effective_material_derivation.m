function [Z,n,B_eff,p_eff] = effective_material_derivation(T,R,d,freq)
    cair =343;
    m = 0;
    k0 =  (2.*pi.*freq)./cair;
    cos_expression = (1-(R.^2)+(T.^2))./(2.*T);
    n = acos(cos_expression)./(k0.*d) + (2*pi*m)./(k0.*d);
    Z = sqrt(((1+R).^2 - (T.^2))./((1-R).^2 - (T.^2)));
    B_eff = Z./n;
    p_eff = n.*Z;
end