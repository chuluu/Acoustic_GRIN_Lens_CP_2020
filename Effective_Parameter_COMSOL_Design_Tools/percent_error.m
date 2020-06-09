function P_E = percent_error(Approx,Experimental)
    P_E = (abs(Approx - Experimental)./Experimental).*100;
    P_E = P_E( ~any( isnan( P_E ) | isinf( P_E ), 2 ),: );
end