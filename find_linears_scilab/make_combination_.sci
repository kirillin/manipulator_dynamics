global count flag

count = 0;
flag = %f;
function make_combination(n,m, num)
    global count flag
    pos = pos + 1
    for i = num:n

        if m == 1 then
            count = count + 1;        
            flag = %t
            printf("finish! %d\n", count)
        else
            make_combination(n, m-1, i+1);
            if flag == %t then
                return;
            end
        end
    end
    pos = pos - 1
endfunction
make_combination(54,54,1)
