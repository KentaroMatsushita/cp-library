template <const int H = 2000, const int W = 2000> struct Matrix {
    int h, w;
    using B = bitset<W>;
    array<B, H> v;
    inline const B &operator[](int k) const { return v[k]; }
    inline B &operator[](int k) { return v[k]; }
    Matrix() {}
    Matrix(int h) : h(h), w(h) {}
    Matrix(int h, int w) : h(h), w(w) {}
    static Matrix I(int n) {
        Matrix mat(n);
        for(int i = 0; i < n; ++i) mat[i][i] = 1;
        return mat;
    }
    Matrix &operator+=(const Matrix &B) { rep(i, h) v[i] ^= B[i]; return *this; }
    Matrix &operator^=(const Matrix &B) { rep(i, h) v[i] ^= B[i]; return *this; }
    Matrix operator+(const Matrix &B) const { return (Matrix(*this) += B); }
    Matrix operator^(const Matrix &B) const { return (Matrix(*this) += B); }

    Matrix<H, W> inv() {
        if(h != w) throw std::invalid_argument("inverse requires square matrix");
        Matrix res = I(h);
        auto w = v;
        for(int i = 0; i < h; ++i) {
            if(!w[i][i]) {
                for(int j = i + 1; j < h; ++j) {
                    if(w[j][i]) {
                        swap(w[j], w[i]), swap(res[j], res[i]);
                        break;
                    }
                }
            }
            if(!w[i][i]) throw std::invalid_argument("singular GF(2) matrix");
            for(int j = 0; j < i; ++j)
                if(w[j][i]) w[j] ^= w[i], res[j] ^= res[i];
            for(int j = i + 1; j < h; ++j)
                if(w[j][i]) w[j] ^= w[i], res[j] ^= res[i];
        }
        return res;
    }
};
