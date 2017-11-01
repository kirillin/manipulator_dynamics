clc();
clear;
global no
no =  [36,37]+1;
function init = converter(cols)
    global no
    //cols in scilab index numeration, i.e. starting with 1
    //init in python index numeration, i.e. starting with 0
    ind = 0;
    for i = 1:70
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

function answer=in(num, vec)
    for i = 1:length(vec)
        if num == vec(i) then
            answer = %t;
            return;
        end
    end
    answer = %f;
endfunction

N = 50;
for i = 1:N
    name = "/media/data/evo/robotics_report/ros_packages/youbot_arm_control/calculations/data_for_identification/ee/EE" + string(i-1) + ".txt";
    r = read(name, -1, 5);
    EE1(:,:,i) = r'
end
s = size(EE1(:,:,N));
printf("Succuessfully read; size of matrices: %dx%d\n", s(1), s(2))


sn = size(no)
arr = zeros(1,s(2) - sn(2));
ind = 1
for i = 1:s(2)
    if ~in(i, no) then
        arr(ind) = i
        ind = ind + 1    
    end
end

ROW_NUM = 5;
for i = 1:N
    for j = 1:s(2) - sn(2)
        E(i,j) = EE1(ROW_NUM:5:s(1),arr(j), i);
    end
end

sz = size(E);
nc = sz(2);

for i = 1:nc
    printf("checking of %d column:\n", i)
    B = E(:,i);
    if sum(B == zeros(sz(1),1)) == sz(1) then
        continue;
    end        
    A = E(:, [1:i-1, i+1:nc]);
    if rank(A) >= rank([A,B]) then
        k = linsolve(A,-B);
        k = clean(k);
        printf("Column No %2d (%d) is dependent on other cols with koefficients", i, converter(i));
        disp(k');
        printf("\n")
    end 
end
