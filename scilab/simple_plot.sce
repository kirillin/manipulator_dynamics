clear;

// INPUT DATA
links_num = 5;
there_are_accels = %t;

// calculation number of columns
if there_are_accels then
    cols_num = 4 * links_num + 1;
else
    cols_num = 3 * links_num + 1;
end

// file reading
data = read("~/Desktop/ident_data/data_1_filt.txt", -1, cols_num)
time = data(:, cols_num);

// graphs plotting
for i = 1:links_num
    kolor = i;
    // angles
    subplot(2, 2, 1);
    plot2d(time, data(:, i), kolor);
    // speeds
    subplot(2, 2, 2);
    plot2d(time, data(:, links_num + i), kolor);
    
    if there_are_accels then
        // accelerations
        subplot(2, 2, 3);
        plot2d(time, data(:, 2*links_num + i), kolor);
        // torques
        subplot(2, 2, 4);
        plot2d(time, data(:,3*links_num+i), kolor);
    else
        // torques
        subplot(2, 1, 2);
        plot2d(time, data(:,2*links_num+i), kolor);
    end
end
