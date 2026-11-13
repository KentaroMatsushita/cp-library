template <typename T = long long, const int sz = 5 * 10000000> struct IO {
    char reader[sz + 1]{}, writer[sz];
    char *now, *now2 = writer;
    IO() {
        reader[fread(reader, sizeof(char), sz, stdin)] = '\0';
        now = reader;
    }
    inline T read() {
        while(*now && *now <= 32) now++;
        if(*now == '-') {
            now++;
            T res = 0;
            while('0' <= *now and *now <= '9') { res = res * 10 + *now++ - '0'; }
            return -res;
        } else {
            T res = 0;
            while('0' <= *now and *now <= '9') { res = res * 10 + *now++ - '0'; }
            return res;
        }
    }
    inline void read(T &res) {
        while(*now && *now <= 32) now++;
        if(*now == '-') {
            now++;
            res = 0;
            while('0' <= *now and *now <= '9') { res = res * 10 + *now++ - '0'; }
            res = -res;
        } else {
            res = 0;
            while('0' <= *now and *now <= '9') { res = res * 10 + *now++ - '0'; }
        }
    }
    inline string read_str() {
        while(*now && *now <= 32) now++;
        string res;
        while(*now && *now > 32) res += *now++;
        if(*now) now++;
        return res;
    }
    template <class S> inline void write(S x, char margin = ' ') {
        if(x == 0) {
            putchar('0');
            putchar(margin);
            return;
        }
        using U = std::make_unsigned_t<S>;
        U value = static_cast<U>(x);
        if(x < 0) {
            putchar('-');
            value = U(0) - value;
        }
        while(value) {
            *now2 = '0' + value % 10;
            now2++;
            value /= 10;
        }
        do {
            now2--;
            putchar(*now2);
        } while(now2 != writer);
        putchar(margin);
    }
    inline void write_str(string s, char margin = ' ') {
        for(auto c : s) putchar(c);
        putchar(margin);
    }
};
IO<int> io;
inline int read() { return io.read(); }
inline void write(ll x) { io.write(x); }
