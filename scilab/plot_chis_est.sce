kolor = 5;
Chi = Chis
sd_chis = stdev(Chi, 'c');
mean_chis = mean(Chi, 'c');

LINEAR_COLS = [2,3,4,6,7,8,9,0,1,5,17,18,20,31,10,14,32,15,16,19,21,22,23,61]+1;
titles = ["$m_{%d}","$mx_{c%d}","$my_{c%d}","$mz_{c%d}","$I_{%d,xx}","$I_{%d,yy}","$I_{%d,zz}","$I_{%d,xy}","$I_{%d,xz}","$I_{%d,yz}","$I_{a, %d}","$f_{v,%d}","$f_{c,%d}","$f_{off,%d}"];

function ans=in(vec, el)
    ans = %f;
    n = length(vec);
    for i = 1:n do
        if el == vec(i) then
            ans = %t;
            return;
        end
    end
endfunction


function init = converter(cols, no)
    //cols in scilab index numeration, i.e. starting with 1
    //init in python index numeration, i.e. starting with 0
    ind = 0;
    for i = 1:70 do
        in = %t;
        for j = 1:length(no)
            if no(j) == i then
                in = %f;
            end
        end
        if in then 
            ind  = ind + 1;
        end
        if ind == cols
            init = i-1
            return 
        end
    end
endfunction
//s = scf();
p = 1;
for k = 1:70 do
    subplot(10,7,k);
    format('v', 10);
    if in(LINEAR_COLS, k) == %f then
        i = p;
        mn = mean_chis(i);
        psd = mn + sd_chis(i);
        nsd = mn - sd_chis(i);
        
        plot2d(1:10, Chi(i,:), kolor);
        plot2d(1:10, ones(1,10) * psd, kolor);
        plot2d(1:10, ones(1,10) * nsd, kolor);
       
        a = gca();
        a.x_ticks = tlist(["ticks", "locations", "labels"],1:1:10, string(["","","","","","","","","",""])); 
    
        a.font_size = 0;
        format('e', 8);
        a.y_ticks = tlist(["ticks", "locations", "labels"],[nsd, mn, psd], string([nsd, mn, psd])); 
        a.margins(1) = 0.3;
        a.margins(2) = 0.0;
        
        format('v');
        j = modulo(k,14)
        if j == 0 then
            j = 14;
        end;
//        a.title.text = string(sprintf(titles(j), ceil(k/14))) + " sd:" + string(sd_chis(i) / abs(mean_chis(i)) * 100) + "\%$";
        //a.title.text = "$sd_{filt}^{blue}:" + string(sd_chis(i) / abs(mean_chis(i)) * 100) + "\%";
        a.title.text = a.title.text + "; sd_{raw}^{red}:" + string(sd_chis(i) / abs(mean_chis(i)) * 100) + "\%$";
        a.title.position = [1, psd];
        a.title.font_size = 2;
        p = p + 1;
    else
        plot(0,0);
        a = gca();
        a.background = 35;
        a.x_ticks = tlist(["ticks", "locations", "labels"],1:1:10, string(["","","","","","","","","",""])); 
    
        a.font_size = 0;
        format('e', 8);
        a.margins(1) = 0.3;
        a.margins(2) = 0.0;
    end
end
