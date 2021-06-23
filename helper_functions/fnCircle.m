function h = fnCircle(x,y,r)
    hold on
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    h = plot(xunit, yunit,'--','Color',[0.5, 0.5, 0.5], 'linewidth', 0.3);
    hold off
end