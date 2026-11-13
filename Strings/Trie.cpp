namespace TrieTree {
constexpr int sz = 26;
constexpr char start = 'a';
struct Node {
    int par;
    int c;
    int nxt[sz];
    // Node に他に持たせる情報

    //
    Node() : c(-1) { memset(nxt, -1, sizeof(nxt)); }
    Node(int par, int c) : c(c), par(par) { memset(nxt, -1, sizeof(nxt)); }
};
struct Trie {
    int num;
    vector<Node> Nodes;
    Trie() : num(0), Nodes(1, Node()) {}
    void add(const string &s, int i, int x) {
        if(i == s.size()) { return; }
        const int c = s[i] - start;
        if(Nodes[x].nxt[c] == -1) Nodes[x].nxt[c] = Nodes.size(), Nodes.emplace_back(Node(x, c));
        add(s, i + 1, Nodes[x].nxt[c]);
    }
    void add(const string &s) { add(s, 0, 0); }
};
} // namespace TrieTree
using namespace TrieTree;

// a ~ z の文字列で書いているので適宜変更！