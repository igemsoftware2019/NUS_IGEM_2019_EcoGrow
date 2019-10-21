function result = format_axis(axis_ticks)
data = zeros(length(axis_ticks),1);

for i=1:length(axis_ticks)
    str = axis_ticks{i};
    num = str2num(str);
    data(i) = num;
end

result = cell(length(data),1);
for i=1:length(data)
    num = data(i);
    str = num2str(num, '%2.2e');
    result{i} = str;
end
end