#include "parallel_hashmap/phmap.h"
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

// #define DEBUG
struct Edge {
  int a, b;
  double weight, origin_weight;
  int lca;
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
  vis[start] = 1;
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
    res[i].origin_weight = edges[i].origin_weight;
    res[i].weight =
        edges[i].weight * log(1.0 * max(deg[edges[i].a], deg[edges[i].b])) /
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
    if (x != fa[x]) // x 不是自身的父亲，即 x 不是该集合的代表
      fa[x] = find_fa(fa[x]); // 查找 x 的祖先直到找到代表，于是顺手路径压缩
    return fa[x];
  }

  void merge(int a, int b) { fa[fa[a]] = fa[b]; }
};

void kruskal(int node_cnt, vector<Edge> &edges, vector<Edge> &tree_edges,
             vector<Edge> &off_tree_edges) {
  ScopeTimer __t("kruskal");
  stable_sort(edges.begin(), edges.end());
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
    off_tree_edges[j] = edges[i];
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
                     vector<Edge> &query_info, int cur, vector<int> &lca,
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
      weighted_depth[e.b] = 1.0 / e.origin_weight + weighted_depth[e.a];
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
      query_info[query_indices[cur][i]].lca = ufs.find_fa(e.b);
    }
  }
}

vector<int> tarjan_lca(const vector<vector<int>> &tree,
                       const vector<Edge> &tree_edges,
                       vector<Edge> &query_info, int root, int node_cnt,
                       vector<double> &weigthed_depth,
                       vector<int> &unweighted_depth) {
  ScopeTimer __t("tarjan_lca");
  vector<int> res(query_info.size());
  UnionFindSet ufs(node_cnt + 1);
  vector<bool> vis(node_cnt + 1);
  vector<vector<int>> query_indices(node_cnt + 1);
  weigthed_depth.resize(node_cnt + 1, 0);
  unweighted_depth.resize(node_cnt + 1, 0);
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
    e.weight = e.origin_weight * (depth[e.a] + depth[e.b] - 2 * depth[e.lca]);
  }
#ifdef DEBUG
  puts("unsorted_off_tree_edges");
  for (auto &x : edges) {
    printf("%d %d %lf %lf\n", x.a, x.b, x.weight, x.origin_weight);
  }
#endif
  stable_sort(edges.begin(), edges.end());
}

vector<Edge> add_off_tree_edges(const int node_cnt,
                                const vector<vector<int>> &tree,
                                const vector<Edge> &tree_edges,
                                const vector<Edge> &off_tree_edges,
                                const vector<int> &lca,
                                const vector<int> &depth) {
  ScopeTimer __t("add_off_tree_edges");
  vector<vector<int>> rebuilt_off_tree_graph(node_cnt + 1);
  for (int i = 0; i < off_tree_edges.size(); ++i) {
    auto &e = off_tree_edges[i];
    rebuilt_off_tree_graph[e.a].push_back(i);
    rebuilt_off_tree_graph[e.b].push_back(i);
  }
  vector<Edge> edges_to_be_add;
  vector<bool> ban(off_tree_edges.size());
  struct UnweightedEdge {
    int u, v;
    bool operator==(const UnweightedEdge &rhs) const {
      return (u == rhs.u && v == rhs.v) || (u == rhs.v && v == rhs.u);
    }
  };
  struct UnweightedEdgeHash {
    size_t operator()(const UnweightedEdge &other) const {
      return hash<int>()(other.u) ^ hash<int>()(other.v);
    }
  };
  for (int i = 0; i < off_tree_edges.size(); ++i) {
    if (edges_to_be_add.size() == max(int(off_tree_edges.size() / 25), 2)) {
      break;
    }
    auto &e = off_tree_edges[i];
    if (ban[i] == 0) {
      edges_to_be_add.push_back(e);
      int beta = min(depth[e.a], depth[e.b]) - depth[e.lca];
#ifdef DEBUG
      printf("beta: %d, (%d, %d)\n", beta, e.a, e.b);
#endif
      struct QueueEntry {
        int node, layer;
      };
      auto beta_layer_bfs_1 = [&tree, &tree_edges,
                               beta](int start, queue<QueueEntry> q,
                                     vector<int> &black_list) {
        vector<bool> vis(tree.size(), false);
        vis[start] = true;
        black_list.push_back(start);
        while (!q.empty()) {
          int cur_node = q.front().node;
          int cur_layer = q.front().layer;
          q.pop();
          for (int j = 0; j < tree[cur_node].size(); ++j) {
            Edge e = tree_edges[tree[cur_node][j]];
            if (e.b == cur_node) {
              swap(e.a, e.b);
            }
            if (!vis[e.b] && cur_layer + 1 <= beta) {
              vis[e.b] = true;
              black_list.push_back(e.b);
              q.push({e.b, cur_layer + 1});
            }
          }
        }
      };
      auto beta_layer_bfs_2 = [&tree, &tree_edges,
                               beta](int start, queue<QueueEntry> q,
                                     phmap::flat_hash_set<int> &black_list) {
        vector<bool> vis(tree.size(), false);
        vis[start] = true;
        black_list.insert(start);
        while (!q.empty()) {
          int cur_node = q.front().node;
          int cur_layer = q.front().layer;
          q.pop();
          for (int j = 0; j < tree[cur_node].size(); ++j) {
            Edge e = tree_edges[tree[cur_node][j]];
            if (e.b == cur_node) {
              swap(e.a, e.b);
            }
            if (!vis[e.b] && cur_layer + 1 <= beta) {
              vis[e.b] = true;
              black_list.insert(e.b);
              q.push({e.b, cur_layer + 1});
            }
          }
        }
      };
      queue<QueueEntry> q1, q2;
      q1.push({e.a, 0}), q2.push({e.b, 0});
      vector<int> black_list1;
      phmap::flat_hash_set<int> black_list2;
      beta_layer_bfs_1(e.a, q1, black_list1);
      beta_layer_bfs_2(e.b, q2, black_list2);
      for (auto u : black_list1) {
        for (int j = 0; j < rebuilt_off_tree_graph[u].size(); ++j) {
          Edge e = off_tree_edges[rebuilt_off_tree_graph[u][j]];
          if (e.b == u) {
            swap(e.a, e.b);
          }
          if (black_list2.count(e.b) > 0) {
            ban[rebuilt_off_tree_graph[u][j]] = 1;
          }
        }
      }
    }
  }
  return edges_to_be_add;
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
    if (f > t) {
      continue;
    }
    G[f].push_back(origin_edges.size());
    G[t].push_back(origin_edges.size());
    volume[f] += w;
    volume[t] += w;
    degree[f]++;
    degree[t]++;
    origin_edges.push_back(Edge{f, t, w, w});
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

#ifdef DEBUG
  puts("kruscal results: ");
  for (auto &x : tree_edges) {
    printf("%d %d %lf %lf\n", x.a, x.b, x.weight, x.origin_weight);
  }
#endif

  auto new_tree = rebuild_tree(M, tree_edges);
  vector<double> tree_weighted_depth;
  vector<int> tree_unweighted_depth;
  vector<int> lca = tarjan_lca(new_tree, tree_edges, off_tree_edges, r_node, M,
                               tree_weighted_depth, tree_unweighted_depth);
  sort_off_tree_edges(off_tree_edges, lca, tree_weighted_depth);

#ifdef DEBUG
  puts("sorted off_tree_edges: ");
  for (auto &x : off_tree_edges) {
    printf("%d %d %.16lf %.16lf\n", x.a, x.b, x.weight, x.origin_weight);
  }
#endif

  vector<Edge> res = add_off_tree_edges(M, new_tree, tree_edges, off_tree_edges,
                                        lca, tree_unweighted_depth);
  gettimeofday(&end, NULL);
  printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000 +
                                     (end.tv_usec - start.tv_usec) / 1000.0);
  /**************************************************/
  /******************* End timing *******************/
  /**************************************************/
  FILE *out = fopen("result.ng.txt", "w");
  for (int i = 0; i < tree_edges.size(); i++) {
    fprintf(out, "%d %d\n", int(tree_edges[i].a), int(tree_edges[i].b));
  }
  for (int i = 0; i < res.size(); i++) {
    fprintf(out, "%d %d\n", int(res[i].a), int(res[i].b));
  }
  fclose(out);
}
