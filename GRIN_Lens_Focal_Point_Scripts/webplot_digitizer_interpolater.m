function value_interp = webplot_digitizer_interpolater(data, new_x)
    old_x = data(:,1);
    phase = data(:,2);
    value_interp = interp1(old_x,phase,new_x);
end