function plot_function(x,y,head,x_lab,y_lab,FontSize)
figure;
C = linspecer(1);
axes('NextPlot','replacechildren', 'ColorOrder',C);
plot(x,y,'Linewidth',2.5);
title(head,'Fontsize',FontSize);
xlabel(x_lab,'Fontsize',FontSize);
ylabel(y_lab,'Fontsize',FontSize);
grid on;

end