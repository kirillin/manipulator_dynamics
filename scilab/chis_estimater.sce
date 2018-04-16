////////////////////////////////////////////////////////////////////////
//////// SCRIPT FOR ESTIMATION AND PLOTING PARAMETERS ...///////////////
////////////////////////////////////////////////////////////////////////
clc();
//clear;
stacksize('max');

////////////////////////////////////////////////////////////////////////
////////// INITIALIZATION //////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

path = get_absolute_file_path("chis_estimater.sce");
path = path + "../data_for_identification_Youbot/";
QTY_FILES = 10;
QTY_JOINTS = 5;
QTY_COLS = 70;
TOTAL_COLS = 70;

LINEAR_COLS = [2,3,4,6,7,8,9,0,1,5,17,18,20,31,10,14,32,15,16,19,21,22,23,61]+1;
TITLES = ["m_{%d}","mx_{c%d}","my_{c%d}","mz_{c%d}","I_{%d,xx}","I_{%d,yy}","I_{%d,zz}","I_{%d,xy}","I_{%d,xz}","I_{%d,yz}","I_{a, %d}","f_{v,%d}","f_{c,%d}","f_{off,%d}"];

function [xi,tau]=read_measuarments(_path, _data_type, number)
    // !!! need global xi and tau
    // _type: filt or raw
    // number: number of file
    // Example for path to file "/home/data/bigs/filt/big_xi1.txt"
    // [xi,tau]=read_measuarments("/home/data/", "filt", 1);
    fn_xi = _path + sprintf("bigs_full/%s/big_xi%d.txt", _data_type, number);
    fn_tau = _path + sprintf("bigs_full/%s/big_tau%d.txt", _data_type, number);
    xi = read(fn_xi, -1, QTY_COLS);
    tau = read(fn_tau, -1, 1);
    printf("Was read: xi[%d, %d], tau[%d, %d]\n", length(xi(:,1)),length(xi(1,:)),length(tau(:,1)), length(tau(1,:)))
endfunction

function [xi,tau]=normalization(xi, tau)
    //NORMALIZING DATA
    sz = size(xi);
    m = sz(1) / QTY_JOINTS;
    for i = 1:QTY_JOINTS
        taus(:,i) = tau(i:QTY_JOINTS:sz(1),:);
    end
    for i = 1:QTY_JOINTS
        norma(i) = norm(taus(:, i))
    end
    winH=waitbar('Нормировка данных...');
    for i = 1:sz(1)
        waitbar(i/sz(1),winH);
        tau(i) = tau(i) / norma(pmodulo(i-1, QTY_JOINTS)+1);
        xi(i, :) = xi(i, :) / norma(pmodulo(i-1, QTY_JOINTS)+1);
    end
    close(winH);
    clear taus;
endfunction

function Chis=estimate_chis(_data_type, method_id)
    Chi = zeros(QTY_COLS, QTY_FILES);
    est_Chi_k_1 = zeros(QTY_COLS, 1);
    for i = 1:QTY_FILES
        [xi,tau]=read_measuarments(path, _data_type, i);
//        [xi,tau]=normalization(xi, tau);

        select method_id
            case 1 then
                Chi(:, i) = xi \ tau;
            case 2 then 
                Chi(:, i) = lsq(xi, tau);
            case 3 then
                // Sigma = 0.0001 * eye(QTY_COLS, QTY_COLS);
                // Chi(:,l) = inv(xi' * xi + Sigma) * xi' * tau;               
                Chi(:, i) = inv(xi' * xi) * xi' * tau;
            case 4 then
                [U,S,V] = svd([xi]);
                Splus = inv(S' * S) * S';
                Chi(:, i) = V * Splus * U' * Tau;
            case 5 then
                [U,S,V] = svd([xi, tau]);
                Chi(:, i) = - V(1:$-1, $) / V($, $);
            case 6 then
                p = QTY_JOINTS;
                sz = size(xi);
                m = sz(1) / QTY_JOINTS;
                r = 1;
                xi_k = xi((r*p-p+1):(r*p), :)';
                P_k_1 = eye(QTY_COLS, QTY_COLS);
                g_k = P_k_1 * xi_k;
                eps_k_1 = zeros(p, 1);
                for k = 2:m do
                    xi_k = xi(r*p-p+1:r*p, :)';
                    y_k = tau(r*p-p+1:r*p,:);
                    r = r + 1;
                    eps_k = y_k - xi_k' * est_Chi_k_1;
                    g_k = P_k_1 * xi_k * inv(1*eye(p,p) + xi_k' * P_k_1 * xi_k);
                    P_k = P_k_1 - g_k * xi_k' * P_k_1;
                    est_Chi_k = est_Chi_k_1 + g_k * eps_k;
                    ///// Plot eps
                    // plot(k-1:k, [eps_k_1, eps_k]);
                    ///// Plot parameters
                    //s = est_Chi_k_1;
                    //e = est_Chi_k;
                    //plot(k-1:k, [s,e]);
                    //disp(k);
                    eps_k_1 = eps_k;
                    P_k_1 = P_k;
                    est_Chi_k_1 = est_Chi_k;
                end;
                Chi(:,i) = est_Chi_k;
        end
        // CLEAR MEMORY
        clear tau;
        clear xi;
    end
    Chis = Chi;
endfunction

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

function [mn,_sd,psd,nsd]=sd(Chis)
    mn = mean(Chis, 'c');
    _sd = stdev(Chis, 'c');
    psd = mn + _sd;
    nsd = mn - _sd;
endfunction

function plot_chis(Chis, kolor, _data_type)
    n = QTY_FILES;
    params = 14;
    lables = []
    for l = 1:n do
        lables(l) = "";
    end
    k = 1;
    for i = 1:TOTAL_COLS do
        subplot(n, TOTAL_COLS/n, i);
        if in(LINEAR_COLS, i) == %f then
            legenda = [];
            sz = size(_data_type);
            for j = 1:sz(2) do
                Chi = Chis(:, (j-1)*n+1:j*n);
                disp(size(Chi))
                [mn,_sd,psd,nsd] = sd(Chi)
                plot2d(1:n, Chi(k, :), kolor(j));
                plot2d(1:n, ones(1, n) * psd(k), kolor(j));
                plot2d(1:n, ones(1, n) * nsd(k), kolor(j));
                a = gca();
                a.x_ticks = tlist(["ticks", "locations", "labels"],1:1:n, string(lables')); 
                a.font_size = 0;
                format('e', 8);
                a.y_ticks = tlist(["ticks", "locations", "labels"],[nsd(j), mn(j), psd(j)], string([nsd(j), mn(j), psd(j)])); 
                format('v');
                a.margins(1) = 0.3;
                a.margins(2) = 0.0;
                l = modulo(i, params)
                if l == 0 then
                    l = params;
                end;
                a.title.text = "$"+string(sprintf(TITLES(l), ceil(i / params)))+"$";
                a.title.position = [3, psd(k)];
                a.title.font_size = 2;
                legenda(j) = string(sprintf("$sd_{%s}: {%.2f}", _data_type(j), _sd(k) / abs(mn(k)) * 100)) + "\%$"
            end
            a = gca()
            legend(legenda(1), legenda(2));
            aa = a.children(1);
            aa.fill_mode = "off";
            k = k + 1;
        else
            plot(0,0);
            a = gca();
            a.font_size = 0;
            a.x_ticks = tlist(["ticks", "locations", "labels"],1:1:n, string(lables')); 
            a.background = 35;
            a.margins(1) = 0.3;
            a.margins(2) = 0.0; 
        end
    end
endfunction

function plot_taus(Chis, _data_type)
    p = QTY_JOINTS;
    for i = 1:QTY_FILES do
        [xi,tau]=read_measuarments(path, _data_type, i);
//        [xi,tau]=normalization(xi, tau);
        sz = size(xi);
        m = sz(1) / QTY_JOINTS;
        tau_calc = xi * Chis(:,i);
        scf();
        for j = 1:QTY_JOINTS do
            subplot(2, 3, j);
            //subplot(QTY_FILES, QTY_JOINTS, (i-1)*QTY_JOINTS+j);
            plot2d(1:m, tau(j:p:sz(1)), 3);
            plot2d(1:m, tau_calc(j:p:sz(1)), 1);
            legend("raw data", 'calc data')
            a = gca()
            a.title.text = 'Link ' + string(j);
        end
        xs2pdf(gcf(), path + string(i) +'.pdf');
        xs2png(gcf(), path + string(i) +'.png');
    end
endfunction

//Chis_filt = estimate_chis("filt", 1)
//mean_Chi_filt = mean(Chis_filt,'c')
//plot_taus(mean_Chi_filt, "filt");

//Chis_raw = estimate_chis("raw", 1)
//plot_taus(Chis_raw, "raw");



plot_chis([Chis_filt, Chis_raw], [2, 5], ["filt", "raw"])




                                        



