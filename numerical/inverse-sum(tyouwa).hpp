// calc (1 / l) + (1 / (l + 1)) + ... + (1 / (r - 1))
long double inv_sum(long double l, long double r) { return boost::math::digamma(r) - boost::math::digamma(l); }

// boost ないとき

long double digamma(long double x) {
    // 定数の定義
    static const long double EulerMascheroni = 0.57721566490153286060651209;
    static const long double pi = 3.14159265358979323846;
    static const long double epsilon = std::numeric_limits<long double>::epsilon();

    // ポール（負の整数やゼロ）の処理
    if(x <= 0.0) {
        if(std::floor(x) == x) {
            // 負の整数ではポールが存在するため無限大を返す
            return std::numeric_limits<long double>::infinity();
        } else {
            // 反射公式を適用
            return digamma(1 - x) - pi / std::tan(pi * x);
        }
    }

    // 小さい x に対する補正
    long double result = 0.0;
    while(x < 10.0) {
        result -= 1.0 / x;
        x += 1.0;
    }

    // ラクンシュ式の級数展開
    long double xx = 1.0 / x;
    long double xx2 = xx * xx;

    result += std::log(x) - 0.5 * xx - xx2 * (1.0 / 12.0 - xx2 * (1.0 / 120.0 - xx2 * (1.0 / 252.0 - xx2 * (1.0 / 240.0 - xx2 * 1.0 / 132.0))));

    return result;
}
long double inv_sum2(long double l, long double r) { return digamma(r) - digamma(l); }