function params = permute_params(table, args, permutation_type)

if permutation_type == 1
    params = combinatorial_make_array(table, args);
else
    params = iterative_make_array(table, args);
end
end

function result = collapse_map(map)
result = [];
keyset = keys(map);

for i = 1:length(keyset)
    key    = keyset{i};
    result = [result; map(key)];
end
end

function result = combinatorial_make_array(table, args)
result = [table{1,:}];
names  = table.Properties.VariableNames;
keyset = keys(args);

for i = 1:length(keyset)
    key           = keyset{i};
    k             = args(key); 
    size_result   = size(result);
    height_result = size_result(1);
    
    for ii = 1:height_result
        row    = result(ii,:);
        new_table = array2table(row, 'VariableNames', names);
        temp   = make_array(new_table, key, k{1}, k{2}, k{3}, 0);
        result = [result; temp];
    end
end

end

function result = iterative_make_array(table, args)
result = containers.Map();
keyset = keys(args);
result('original') = table{1,:};
for i = 1:length(keyset)
    key  = keyset{i};
    k    = args(key); 
    temp = make_array(table, key, k{1}, k{2}, k{3}, 0);
    result(key) = temp;
end

end

function result = make_array(table, vary, func, interval, range, inclusive)

names = table.Properties.VariableNames;
num   = 1:length(names);
name_to_num = containers.Map(names, num);

p_0 = table{1,:};

variable = name_to_num(vary);
if func == 1 
    arr = log_permute(p_0(variable), interval, range, inclusive);
else
    arr= linear_permute(p_0(variable), interval, inclusive);
end

result = zeros(length(arr), length(p_0));

for i = 1:length(arr)
    result(i,:)        = p_0;
    result(i,variable) = arr(i);
end

end

function result = log_permute(x, interval, range, inclusive)
if inclusive == 1
    result = zeros(range*4+1, 1);
else
    result = zeros(range*4, 1);
end

index = 1;
for i = -range:range
    value = x*10^i;
    for ii = 1:length(interval)
        if inclusive ~= 1 && 10^i*interval(ii) == 1
            continue
        else
            result(index)   = value*interval(ii);
            index = index + 1;
        end
    end
    
end

end

function result = linear_permute(x, interval, inclusive)

if inclusive == 1
    result = x*ones(length(interval), 1) + interval;
else
    interval2 = interval(interval ~= 0);
    result = x*ones(length(interval2), 1) + interval2;
end

end