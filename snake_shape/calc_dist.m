function ret = calc_dist(x_set,y_set)
n = length(x_set);
if length(y_set) ~= n
    error("x set and y set are not the same length")
end
dist_cur = 0;
for i = 1:n-1
    dist_added = sqrt((x_set(i) - x_set(i+1))^2 + (y_set(i) - y_set(i+1))^2);
    dist_cur = dist_cur + dist_added;
end

ret = dist_cur;
end