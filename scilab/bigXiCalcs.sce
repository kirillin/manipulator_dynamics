a = [0.5, 0.4];
d = [0, 0];
delta_1 = 60 * %pi / 180;
delta_2 = 40 * %pi / 180;

delta = [delta_1, delta_2]

function [Xi]=getXi(q, dq, ddq)
    theta = delta - q;
    xi111 = a(1)*(a(1)*ddq(1) - 9.82*cos(theta(1)));
    xi112 = 2.0*a(1)*ddq(1) - 9.82*cos(theta(1));
    xi113 = 9.82*sin(theta(1));
    xi114 = 0;
    xi115 = 0;
    xi116 = 0;
    xi117 = ddq(1);
    xi118 = 0;
    xi119 = 0;
    xi1110 = 0;
    xi1111 = ddq(1);
    xi1112 = dq(1);
    xi1113 = sign(dq(1));
    xi1114 = 1;
    
    xi121 = 1.0*a(1)**2*ddq(1) + 2.0*a(1)*a(2)*sin(theta(2))*dq(1)*dq(2) + 1.0*a(1)*a(2)*sin(theta(2))*dq(2)**2 + 2.0*a(1)*a(2)*cos(theta(2))*ddq(1) + 1.0*a(1)*a(2)*cos(theta(2))*ddq(2) - 9.82*a(1)*cos(theta(1)) + 1.0*a(2)**2*ddq(1) + 1.0*a(2)**2*ddq(2) + 9.82*a(2)*sin(theta(1))*sin(theta(2)) - 9.82*a(2)*cos(theta(1))*cos(theta(2));
    xi122 = 2.0*a(1)*sin(theta(2))*dq(1)*dq(2) + 1.0*a(1)*sin(theta(2))*dq(2)**2 + 2.0*a(1)*cos(theta(2))*ddq(1) + 1.0*a(1)*cos(theta(2))*ddq(2) + 2.0*a(2)*ddq(1) + 2.0*a(2)*ddq(2) - 9.82*cos(theta(1) + theta(2));
    xi123 = 0;
    xi124 = 2.0*a(1)*sin(theta(2))*ddq(1) + 1.0*a(1)*sin(theta(2))*ddq(2) - 2.0*a(1)*cos(theta(2))*dq(1)*dq(2) - 1.0*a(1)*cos(theta(2))*dq(2)**2 - 9.82*sin(theta(1) + theta(2));
    xi125 = 0;
    xi126 = ddq(1) + ddq(2);
    xi127 = 0;
    xi128 = 0;
    xi129 = 0;
    xi1210 = 0;
    xi1211 = 0;
    xi1212 = 0;
    xi1213 = 0;
    xi1214 = 0;

    xi221 = a(2)*(-1.0*a(1)*sin(theta(2))*dq(1)**2 + 1.0*a(1)*cos(theta(2))*ddq(1) + 1.0*a(2)*ddq(1) + 1.0*a(2)*ddq(2) - 9.82*cos(theta(1) + theta(2)));
    xi222 = -1.0*a(1)*sin(theta(2))*dq(1)**2 + 1.0*a(1)*cos(theta(2))*ddq(1) + 2.0*a(2)*ddq(1) + 2.0*a(2)*ddq(2) - 9.82*cos(theta(1) + theta(2));
    xi223 = 0;
    xi224 = 1.0*a(1)*sin(theta(2))*ddq(1) + 1.0*a(1)*cos(theta(2))*dq(1)**2 - 9.82*sin(theta(1) + theta(2));
    xi225 = 0;
    xi226 = ddq(1) + ddq(2);
    xi227 = 0;
    xi228 = 0;
    xi229 = 0;
    xi2210 = 0;
    xi2211 = ddq(2);
    xi2212 = dq(2);
    xi2213 = sign(dq(2));
    xi2214 = 1;

    Xi = [xi111, xi112, xi113, xi114, xi115, xi116, xi117, xi118,xi119, xi1110,xi1111,xi1112,xi1113,xi1114, xi121, xi122, xi123, xi124, xi125, xi126, xi127, xi128,xi129, xi1210,xi1211,xi1212,xi1213,xi1214;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0, xi221, xi222, xi223, xi224, xi225, xi226, xi227, xi228,xi229, xi2210,xi2211,xi2212,xi2213,xi2214]    
endfunction

// read ident. data
path = "/media/data/evo/manipulator_dynamics/data_for_identification_2dof/";
data_path = path + "data/data_1_filt.txt";
data = read(data_path, -1, 9);
Q = data(:,1:2);
dQ = data(:,3:4);
ddQ = data(:,5:6);
Tau = data(:,7:8);
T = data(:,8);


// make big Xi and big Tau
bigXiFull = [];
bigTau = [];

k = 1;
for i = 1:length(Q)/2
    Xi = getXi([Q(i,1), Q(i,2)], [dQ(i,1), dQ(i,2)], [ddQ(i,1), ddQ(i,2)]);
    bigXiFull(k,:) = Xi(1,:);
    bigXiFull(k+1,:) = Xi(2,:);
    bigTau(k) = Tau(i,1);
    bigTau(k+1) = Tau(i,2);
    k = k + 2;
end

// compress big Xi with well columns
bigXi = []
wellCols = [1, 2, 10, 11, 12, 13, 15, 17, 19, 24, 25, 26, 27] + 1;
k = 1;
for i = 1:length(wellCols)
    bigXi(:, k) = bigXiFull(:, wellCols(i));
    k = k + 1;
end


sz = size(bigXi);
m = sz(1) / 2;

// calculate parameters Chi
Chi = inv(bigXi' * bigXi) * bigXi' * bigTau;
ChiReal = [4.5, 0, -2.5833333333333335, 0, 0, 0, 1.0, 0, -0.5333333333333334, 0.0, 0, 0, 0]'


// plot graphics
scf()
tau_calc = bigXi * Chi;
tau_calc_real = bigXi * ChiReal;
for i = 1:2
    subplot(2, 2, i);
    plot2d(1:m, bigTau(i:2:sz(1)), 3);
    plot2d(1:m, tau_calc(i:2:sz(1)), 1);
    plot2d(1:m, tau_calc_real(i:2:sz(1)), 5);
    
    a = gca()
    a.title.text = 'Usual calculations. Link ' + string(i);
end
legend("raw data", "ident data", 'calc data')


//###################################//
// normalize data
for i = 1:2
    norma(i) = norm(Tau(:, i))
end

winH=waitbar('Нормировка данных...');
for i = 1:sz(1)
    waitbar(i/sz(1),winH);
    bigTau(i) = bigTau(i) / norma(pmodulo(i-1, 2)+1);
    bigXi(i, :) = bigXi(i, :) / norma(pmodulo(i-1, 2)+1);
end
close(winH);


// calculate parameters Chi
Chi = inv(bigXi' * bigXi) * bigXi' * bigTau;
ChiReal = [4.5, 0, -2.5833333333333335, 0, 0, 0, 1.0, 0, -0.5333333333333334, 0.0, 0, 0, 0]'


// plot graphics
tau_calc = bigXi * Chi;
tau_calc_real = bigXi * ChiReal;
for i = 1:2
    subplot(2, 2, i+2);
    plot2d(1:m, bigTau(i:2:sz(1)), 3);
    plot2d(1:m, tau_calc(i:2:sz(1)), 1);
    plot2d(1:m, tau_calc_real(i:2:sz(1)), 5);
    
    a = gca()
    a.title.text = 'Usual calculations. Link ' + string(i);
end
legend("raw data", "ident data", 'calc data')
