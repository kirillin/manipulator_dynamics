clc();
clear;
N = 50;
P = 70;
for i = 1:N
    name = "/media/data/evo/robotics_report/ros_packages/youbot_arm_control/calculations/data_for_identification/ee/EE" + string(i-1) + ".txt";
    r = read(name, -1, 5);
    E((i-1)*5+1:i*5,1:P) = r'
end


function v=isZeroVec(vec)
    n = max(size(vec))
    flag = %t;
    for i = 1:n
        clear_vec = clean(vec, 10^-20);
        if sum(clear_vec == zeros(vec)) <> n then 
            flag = %f;
            break;
        end    
    end
    v = flag;
endfunction

printf("For !python! use these indexes:\n");
for j = 1:P
    if isZeroVec(E(:, j)) then
        printf('%d, ',j-1)
    end
end


