template <class T, const int N> struct Matrix {
    int n = N;
    array<array<T, N>, N> A{};

    Matrix() {}

    Matrix(size_t n) : n(n) {}

    inline const array<T, N> &operator[](int k) const { return A[k]; }

    inline array<T, N> &operator[](int k) { return A[k]; }

    static Matrix I(size_t n) {
        Matrix mat(n);
        for(int i = 0; i < n; i++) mat[i][i] = 1;
        return (mat);
    }

    Matrix &operator+=(const Matrix &B) {
        assert(n == B.n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++) (*this)[i][j] += B[i][j];
        return (*this);
    }

    Matrix &operator-=(const Matrix &B) {
        assert(n == B.n);
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++) (*this)[i][j] -= B[i][j];
        return (*this);
    }

    Matrix &operator*=(const Matrix &B) {
        assert(n == B.n);
        array<array<T, N>, N> res{};
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                for(int k = 0; k < n; k++) res[i][j] = (res[i][j] + (*this)[i][k] * B[k][j]);
        A.swap(res);
        return (*this);
    }

    Matrix &operator^=(long long k) {
        Matrix B = Matrix::I(n);
        while(k > 0) {
            if(k & 1) B *= *this;
            *this *= *this;
            k >>= 1LL;
        }
        A.swap(B.A);
        return (*this);
    }

    Matrix operator+(const Matrix &B) const { return (Matrix(*this) += B); }

    Matrix operator-(const Matrix &B) const { return (Matrix(*this) -= B); }

    Matrix operator*(const Matrix &B) const { return (Matrix(*this) *= B); }

    Matrix operator^(const long long k) const { return (Matrix(*this) ^= k); }

    friend ostream &operator<<(ostream &os, Matrix &p) {
        size_t n = p.n, m = p.n;
        os << '\n';
        for(int i = 0; i < n; i++) {
            os << "[";
            for(int j = 0; j < m; j++) { os << p[i][j] << (j + 1 == m ? "]\n" : ","); }
        }
        return (os);
    }

    T determinant() {
        Matrix B(*this);
        T ret = 1;
        for(int i = 0; i < n; i++) {
            int idx = -1;
            for(int j = i; j < n; j++) {
                if(B[j][i] != 0) idx = j;
            }
            if(idx == -1) return (0);
            if(i != idx) {
                ret *= -1;
                swap(B[i], B[idx]);
            }
            ret *= B[i][i];
            T vv = B[i][i];
            for(int j = 0; j < n; j++) { B[i][j] /= vv; }
            for(int j = i + 1; j < n; j++) {
                T a = B[j][i];
                for(int k = 0; k < n; k++) { B[j][k] -= B[i][k] * a; }
            }
        }
        return (ret);
    }
    friend vector<T> operator*(const Matrix &A, const vector<T> &b) {
        assert(A.n == b.size());
        vector<T> res(A.n);
        rep(i, A.n) rep(j, A.n) res[i] += A[i][j] * b[j];
        return res;
    }
    friend vector<T> operator*(const vector<T> &b, const Matrix &A) {
        assert(A.n == b.size());
        vector<T> res(A.n);
        rep(i, A.n) rep(j, A.n) res[j] += b[i] * A[i][j];
        return res;
    }
};
