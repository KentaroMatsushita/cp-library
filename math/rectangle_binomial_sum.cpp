mint rectangle_binomial_sum(int lx, int ly, int rx, int ry) { return C(rx + ry, rx) - C(rx + ly, rx) - C(lx + ry, ry) + C(lx + ly, lx); }
