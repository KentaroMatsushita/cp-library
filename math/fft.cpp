// Iterative radix-2 FFT. Forward uses exp(-2*pi*i*k/n); inverse divides by n.
struct FFT {
    void fft(vector<complex<double>>& a,bool inv=false) const {
        size_t n=a.size();
        if(n==0) return;
        if(n&(n-1)) throw std::invalid_argument("FFT length must be a power of two");
        for(size_t i=1,j=0;i<n;++i) {
            size_t bit=n>>1;
            for(;j&bit;bit>>=1) j^=bit;
            j^=bit;
            if(i<j) swap(a[i],a[j]);
        }
        double pi=acos(-1.0);
        for(size_t len=2;len<=n;len<<=1) {
            complex<double> step=polar(1.0,(inv?2:-2)*pi/len);
            for(size_t i=0;i<n;i+=len) {
                complex<double> w=1;
                for(size_t j=0;j<len/2;++j,w*=step) {
                    auto u=a[i+j],v=w*a[i+j+len/2];
                    a[i+j]=u+v;a[i+j+len/2]=u-v;
                }
            }
        }
        if(inv) for(auto& v:a) v/=double(n);
    }
};
