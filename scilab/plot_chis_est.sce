
sd_chis = stdev(Chi, 'c');
mean_chis = mean(Chi, 'c');

titles = ["$m_{%d}$","$mx_{c%d}$","$my_{c%d}$","$mz_{c%d}$","$I_{%d,xx}$","$I_{%d,yy}$","$I_{%d,zz}$","$I_{%d,xy}$","$I_{%d,xz}$","$I_{%d,yz}$","$I_{a, %d}$","$f_{v,%d}$","$f_{c,%d}$","$f_{off,%d}$"];

s = scf();
for i = 1:70 do
    subplot(10,7,i);
    
    mn = mean_chis(i);
    psd = mn + sd_chis(i);
    nsd = mn - sd_chis(i)
    
    plot2d(1:10, Chi(i,:), 2);
    plot2d(1:10, ones(1,10) * psd, 5);
    plot2d(1:10, ones(1,10) * nsd, 5);
    a = gca();
    a.x_ticks = tlist(["ticks", "locations", "labels"],1:1:10, string(["","","","","","","","","",""])); 

    a.font_size = 0;
    format('e', 8);
    a.y_ticks = tlist(["ticks", "locations", "labels"],[nsd, mn, psd], string([nsd, mn, psd])); 
    a.margins(1) = 0.3;
    a.margins(2) = 0.0;
    
    format('v', 5);
    j = modulo(i,14)
    if j == 0 then
        j = 14;
    end;
    a.title.text = string(sprintf(titles(j), ceil(i/14)));
    a.title.position = [5, psd];
    a.title.font_size = 2;
end
