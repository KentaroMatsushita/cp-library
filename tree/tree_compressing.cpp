#pragma once
#include "HLD.cpp"
#include <algorithm>
#include <utility>
#include <vector>

// The input is an undirected forest; make() includes each selected component's root.
struct LCA_Tree {
    using Graph = std::vector<std::vector<int>>;
    HLDecomposition<Graph> hld;
    Graph &g;

    LCA_Tree(Graph &graph, std::vector<int> roots = {})
        : hld(graph, roots), g(graph) { hld.build(); }

    int lca(int x, int y) { return hld.lca(x, y); }
    int d(int x) { return hld.d[x]; }

    // Mutates nodes into the Euler-ordered vertex set of the returned forest.
    std::pair<std::vector<int>, Graph> make(std::vector<int> &nodes) {
        auto comp = [&](int i, int j) { return hld.in[i] < hld.in[j]; };
        std::sort(nodes.begin(), nodes.end(), comp);
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        const int selected_count = nodes.size();
        for(int i = 1; i < selected_count; ++i) {
            if(hld.root[nodes[i - 1]] == hld.root[nodes[i]])
                nodes.emplace_back(lca(nodes[i - 1], nodes[i]));
        }
        for(int i = 0; i < selected_count; ++i) nodes.emplace_back(hld.root[nodes[i]]);
        std::sort(nodes.begin(), nodes.end(), comp);
        nodes.erase(std::unique(nodes.begin(), nodes.end()), nodes.end());
        Graph compressed(nodes.size());
        std::vector<int> stack;
        for(int i = 0; i < int(nodes.size()); ++i) {
            if(hld.root[nodes[i]] == nodes[i]) {
                stack.clear();
            } else {
                while(hld.out[nodes[stack.back()]] <= hld.in[nodes[i]]) stack.pop_back();
                compressed[stack.back()].emplace_back(i);
                compressed[i].emplace_back(stack.back());
            }
            stack.emplace_back(i);
        }
        return {nodes, compressed};
    }
};
