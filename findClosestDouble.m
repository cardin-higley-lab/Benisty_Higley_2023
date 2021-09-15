function ind = findClosestDouble(vec, val)

[~, ind] = min(abs(vec-val));