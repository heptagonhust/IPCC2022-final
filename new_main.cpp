#include "timer.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <queue>
#include <stack>
#include <string>
#include <sys/time.h>
#include <vector>

using namespace std;
struct Edge {
  int a, b;
  double weight;
  bool operator<(const Edge &rhs) const { return this->weight > rhs.weight; }
};

int largest_volume_node(const vector<double> &volume) {
  ScopeTimer __t("largest_volume_node");
  return max_element(volume.begin(), volume.end()) - volume.begin();
}

vector<int> get_unweighted_distance_bfs(const vector<Edge> &edges,
                                        const vector<vector<int>> &G,
                                        int start) {
  ScopeTimer __t("get_unweighted_distance_bfs");
  vector<int> res(G.size());
  vector<bool> vis(G.size());
  queue<int> q;
  q.push(start);
  int cnt = 0;
  while (!q.empty()) {
    cnt++;
    int top = q.front();
    q.pop();
    for (int i = 0; i < G[top].size(); ++i) {
      Edge e = edges[G[top][i]];
      if (e.b == top) {
        swap(e.a, e.b);
      }
      if (!vis[e.b]) {
        res[e.b] = res[e.a] + 1;
        q.push(e.b);
        vis[e.b] = 1;
      }
    }
  }
  return res;
}

vector<Edge> get_new_edges(const vector<Edge> &edges, const vector<int> &deg,
                           const vector<int> &unweighted_distance) {
  ScopeTimer __t("get_new_edges");
  vector<Edge> res(edges.size());
  for (int i = 0; i < edges.size(); ++i) {
    res[i].a = edges[i].a;
    res[i].b = edges[i].b;
    res[i].weight =
        log(1.0 * max(deg[edges[i].a], deg[edges[i].b])) /
        (unweighted_distance[edges[i].a] + unweighted_distance[edges[i].b]);
  }
  return res;
}

struct UnionFindSet {
  vector<int> fa;
  UnionFindSet(size_t sz) : fa(sz) {
    for (int i = 0; i < sz; ++i) {
      fa[i] = i;
    }
  }
  int find_fa(int x) {
    // 寻找x的祖先
    if (fa[x] == x) // 如果 x 是祖先则返回
      return x;
    else
      return fa[x] = find_fa(fa[x]); // 如果不是则 x 的爸爸问 x 的爷爷
  }

  void merge(int a, int b) { fa[a] = fa[b]; }
};

void kruskal(int node_cnt, vector<Edge> &edges, vector<Edge> &tree_edges,
             vector<Edge> &off_tree_edges) {
  ScopeTimer __t("kruskal");
  sort(edges.begin(), edges.end());
  tree_edges.reserve(node_cnt - 1);
  off_tree_edges.reserve(edges.size() - (node_cnt - 1));
  UnionFindSet ufs(node_cnt + 1);
  int edge_cnt = 0;
  for (; tree_edges.size() != node_cnt - 1; ++edge_cnt) {
    auto &edge = edges[edge_cnt];
    if (ufs.find_fa(edge.a) != ufs.find_fa(edge.b)) {
      ufs.merge(edge.a, edge.b);
      tree_edges.push_back(edge);
    } else {
      off_tree_edges.push_back(edge);
    }
  }
  int off_tree_edges_cur = off_tree_edges.size();
  off_tree_edges.resize(off_tree_edges.size() + edges.size() - edge_cnt);
  for (int i = edge_cnt, j = off_tree_edges_cur; i < edges.size(); ++i, ++j) {
    off_tree_edges[j] = edges[j];
  }
}

vector<vector<int>> rebuild_tree(int node_cnt, const vector<Edge> &tree) {
  ScopeTimer __t("rebuild_tree");
  vector<vector<int>> res(node_cnt + 1);
  for (int i = 0; i < tree.size(); ++i) {
    auto &e = tree[i];
    res[e.a].push_back(i);
    res[e.b].push_back(i);
  }
  return res;
}

void tarjan_lca_impl(const vector<vector<int>> &tree,
                     const vector<Edge> &tree_edges,
                     const vector<vector<int>> &query_indices,
                     const vector<Edge> &query_info, int cur, vector<int> &lca,
                     UnionFindSet &ufs, vector<bool> &vis,
                     vector<double> &weighted_depth,
                     vector<int> &unweighted_depth) {
  ufs.fa[cur] = cur;
  vis[cur] = 1;
  for (int i = 0; i < tree[cur].size(); ++i) {
    Edge e = tree_edges[tree[cur][i]];
    if (e.b == cur) {
      swap(e.a, e.b);
    }
    if (!vis[e.b]) {
      weighted_depth[e.b] = e.weight + weighted_depth[e.a];
      unweighted_depth[e.b] = 1 + unweighted_depth[e.a];
      tarjan_lca_impl(tree, tree_edges, query_indices, query_info, e.b, lca,
                      ufs, vis, weighted_depth, unweighted_depth);
      ufs.fa[e.b] = cur;
    }
  }
  for (int i = 0; i < query_indices[cur].size(); ++i) {
    Edge e = query_info[query_indices[cur][i]];
    if (e.b == cur) {
      swap(e.a, e.b);
    }
    if (vis[e.b]) {
      lca[query_indices[cur][i]] = ufs.find_fa(e.b);
    }
  }
}

vector<int> tarjan_lca(const vector<vector<int>> &tree,
                       const vector<Edge> &tree_edges,
                       const vector<Edge> &query_info, int root, int node_cnt,
                       vector<double> &weigthed_depth,
                       vector<int> &unweighted_depth) {
  ScopeTimer __t("tarjan_lca");
  vector<int> res(query_info.size());
  UnionFindSet ufs(node_cnt + 1);
  vector<bool> vis(node_cnt + 1);
  vector<vector<int>> query_indices(node_cnt + 1);
  weigthed_depth.resize(node_cnt + 1);
  for (int i = 0; i < query_info.size(); ++i) {
    auto &e = query_info[i];
    query_indices[e.a].push_back(i);
    query_indices[e.b].push_back(i);
  }
  tarjan_lca_impl(tree, tree_edges, query_indices, query_info, root, res, ufs,
                  vis, weigthed_depth, unweighted_depth);
  return res;
}

void sort_off_tree_edges(vector<Edge> &edges, const vector<int> &lca,
                         const vector<double> &depth) {
  ScopeTimer __t("sort_off_tree_edges");
  for (int i = 0; i < edges.size(); ++i) {
    auto &e = edges[i];
    e.weight *= depth[e.a] + depth[e.b] - 2 * depth[lca[i]];
  }
  sort(edges.begin(), edges.end());
}

int main(int argc, const char *argv[]) {
  // read input file
  const char *file = "byn1.mtx";
  if (argc > 2) {
    printf("Usage : ./main <filename>");
    exit(-1);
  } else if (argc == 2) {
    file = argv[1];
  }
  // the matrix you read must be a adjacency matrix
  ifstream fin(file);
  int M, N, L;
  // Ignore headers and comments
  while (fin.peek() == '%')
    fin.ignore(2048, '\n');
  // declare matrix vector, volume and degree
  fin >> M >> N >> L;
  // M == N, point cnt
  // L, edge cnt
  vector<vector<int>> G(M + 1);
  vector<Edge> origin_edges;

  vector<int> degree(M + 1);
  vector<double> volume(M + 1);
  for (int i = 0; i < L; ++i) {
    int f, t;
    double w;
    fin >> f >> t >> w;
    bool exists = false;
    for (int j = 0; j < G[f].size(); ++j) {
      auto &e = origin_edges[G[f][j]];
      if ((e.a == f && e.b == t) || (e.a == t && e.b == f)) {
        exists = true;
        break;
      }
    }
    if (exists) {
      continue;
    }
    G[f].push_back(origin_edges.size());
    G[t].push_back(origin_edges.size());
    volume[f] += w;
    volume[t] += w;
    degree[f]++;
    degree[t]++;
    origin_edges.push_back(Edge{f, t, w});
  }
  fin.close();
  printf("edge_cnt: %ld\n", origin_edges.size());
  /**************************************************/
  /***************** Start timing *******************/
  /**************************************************/
  struct timeval start, end;
  gettimeofday(&start, NULL);
  int r_node = largest_volume_node(volume);
  auto unweighted_distance =
      get_unweighted_distance_bfs(origin_edges, G, r_node);
  auto new_edges = get_new_edges(origin_edges, degree, unweighted_distance);
  vector<Edge> tree_edges, off_tree_edges;
  kruskal(M, new_edges, tree_edges, off_tree_edges);
  auto new_tree = rebuild_tree(M, tree_edges);
  vector<double> tree_weighted_depth;
  vector<int> tree_unweighted_depth;
  vector<int> lca = tarjan_lca(new_tree, tree_edges, off_tree_edges, r_node, M,
                               tree_weighted_depth, tree_unweighted_depth);
  sort_off_tree_edges(off_tree_edges, lca, tree_weighted_depth);
  gettimeofday(&end, NULL);
  printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000 +
                                     (end.tv_usec - start.tv_usec) / 1000.0);
  /**************************************************/
  /******************* End timing *******************/
  /**************************************************/
}
