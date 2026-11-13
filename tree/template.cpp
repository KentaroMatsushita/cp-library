template< typename T >
struct edge{
    int from, to;
    T cost;
    int id;
    edge(int to,T cost) : from(-1), to(to), cost(cost){}
    edge(int from,int to,T cost) : from(from), to(to), cost(cost){}
    // edge(int from,int to,T cost,int id):from(from),to(to),cost(cost),id(id){}
    edge &operator=(const int &x) {
        to = x;
        return *this;
    }
    operator int() const{ return to;}
};

template< typename T >
using Edges = vector< edge< T > >;
template< typename T >
using WeightedTree = vector< Edges<T>>;
using tree = vector< vector<int> >;

tree make(int n,int offset = 1){
    tree res(n);
    for(int i = 0;i < n-1; i++){
        int a,b; cin >> a >> b;
        a -= offset,b -= offset;
        res[a].emplace_back(b);
        res[b].emplace_back(a);
    }
    return res;
}
template< typename T >
WeightedTree<T> make2(int n, int offset = 1){
    WeightedTree<T> res(n);
    for(int i = 0;i < n-1 ; i++){
        int a,b ; cin >> a >> b;
        a -= offset, b -= offset;
        T c; cin >> c;
        res[a].emplace_back(b,c);
        res[b].emplace_back(a,c);
    }
    return res;
}