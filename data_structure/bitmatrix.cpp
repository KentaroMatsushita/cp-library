using B = bitset<500>;
struct Matrix {
    array<B, 500> a;
    int n, m;
    Matrix(int n) : n(n), m(n) {}
    Matrix() : n(500), m(500) {}
    B &operator[](int k) { return a[k]; }
    Matrix operator*(Matrix &rhs) {
        Matrix mat(n);
        rep(i, n) rep(j, n) mat[j][i] = rhs[i][j];
        Matrix res(n);
        rep(i, n) rep(j, n) { res[i][j] = (a[i] & mat[j]).any(); }
        return res;
    }
    Matrix operator^(int k) const {
        Matrix res(n);
        rep(i, n) res[i][i] = 1;
        auto now = (*this);
        while(k) {
            if(k & 1) res = res * now;
            k >>= 1;
            now = now * now;
        }
        return res;
    }
    B operator*(const B &rhs) const {
        B res;
        rep(i, n) res[i] = (a[i] & rhs).any();
        return res;
    }
    bool any() {
        rep(i, n) if(a[i].any()) return true;
        return false;
    }
};