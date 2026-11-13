
template< typename T >
struct edge {
  int src, to;
  T cost;

  edge(int to, T cost) : src(-1), to(to), cost(cost) {}

  edge(int src, int to, T cost) : src(src), to(to), cost(cost) {}

  edge &operator=(const int &x) {
    to = x;
    return *this;
  }


  operator int() const { return to; }
};

template< typename T >
using Edges = vector< edge< T > >;
template< typename T >
using WeightedGraph = vector< Edges< T > >;
using UnWeightedGraph = vector< vector< int > >;
template< typename T >
using Matrix = vector< vector< T > >;

UnWeightedGraph make(int n,int m,int offset=1){
    UnWeightedGraph res(n);
    for(int i = 0;i<m;++i){
        int a,b; cin>>a>>b;
        a-=offset,b-=offset;
        res[a].emplace_back(b);
        res[b].emplace_back(a);
    }
    return res;
}

template< typename T >
WeightedGraph<T> make2(int n,int m,int offset=1){
    WeightedGraph<T> res(n);
    for(int i = 0;i<m;++i){
        int a,b;cin>>a>>b;
        T c;cin>>c;
        a-=offset,b-=offset;
        res[a].emplace_back(b,c);
        res[b].emplace_back(a,c);
    }
    return res;
}
