#include "timer.hpp"
#include <algorithm>
#include <boost/sort/sort.hpp>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <memory>
#include <queue>
#include <stack>
#include <string>
#include <sys/time.h>
#include <unordered_set>
#include <vector>
using namespace std;

// #define DEBUG
struct Edge {
  int a, b;
  double weight, origin_weight;
  int lca;
  bool operator<(const Edge &rhs) const { return this->weight > rhs.weight; }
};

extern "C" __attribute__((noinline)) void magic_trace_stop_indicator() {
  asm volatile("" ::: "memory");
}
template <typename T> struct CSRMatrix {
  std::unique_ptr<T[]> row_indices;
  std::unique_ptr<T[]> neighbors;

  CSRMatrix(const int &matrix_dim, const int &nonzero_count)
      : row_indices(new T[matrix_dim + 1]), neighbors(new T[nonzero_count]) {}
};

template <bool use_edge_list = false>
CSRMatrix<int> build_csr_matrix(const int &node_cnt, const int &edge_cnt,
                                Edge *edges) {
  vector<int> deg(node_cnt + 1);
  for (int i = 0; i < edge_cnt; ++i) {
    deg[edges[i].a]++;
    deg[edges[i].b]++;
  }
  auto mat = CSRMatrix<int>(node_cnt + 1, edge_cnt * 2);
  int start_idx = 0;
  for (int i = 1; i <= node_cnt; i++) {
    mat.row_indices[i] = start_idx;
    start_idx += deg[i];
  }
  mat.row_indices[0] = 0;
  mat.row_indices[node_cnt + 1] = start_idx;
  for (int i = 0; i < edge_cnt; ++i) {
    const auto &e = edges[i];
    if constexpr (use_edge_list) {
      mat.neighbors[mat.row_indices[e.a]] = i;
      mat.neighbors[mat.row_indices[e.b]] = i;
    } else {
      mat.neighbors[mat.row_indices[e.a]] = e.b;
      mat.neighbors[mat.row_indices[e.b]] = e.a;
    }
    mat.row_indices[e.a]++;
    mat.row_indices[e.b]++;
  }
  for (int i = node_cnt; i >= 1; i--) {
    mat.row_indices[i] = mat.row_indices[i - 1];
  }
  return mat;
}

int largest_volume_node(const int &node_cnt, double *volume) {
  ScopeTimer t_("largest_volume_node");
  double max_volume = 0.0;
  int max_index = 0;
  for (int i = 1; i <= node_cnt; i++) {
    if (volume[i] > max_volume) {
      max_index = i;
      max_volume = volume[i];
    }
  }

  return max_index;
}

int *get_unweighted_distance_bfs(Edge *edges, const CSRMatrix<int> &G,
                                 int start, int node_cnt) {
  ScopeTimer t_("get_unweighted_distance_bfs");
  int *res(new int[node_cnt + 1]());
  vector<bool> vis(node_cnt + 1);
  vis[start] = 1;
  queue<int> q;
  q.push(start);
  while (!q.empty()) {
    int top = q.front();
    q.pop();
    for (int i = G.row_indices[top]; i < G.row_indices[top + 1]; ++i) {
      const Edge &e = edges[G.neighbors[i]];
      int v = top ^ e.a ^ e.b;
      if (!vis[v]) {
        res[v] = res[top] + 1;
        q.push(v);
        vis[v] = 1;
      }
    }
  }
  return res;
}

void get_new_edges(const int &edge_cnt, Edge *edges, int *deg,
                   int *unweighted_distance) {
  ScopeTimer t_("get_new_edges");
#pragma omp parallel for
  for (int i = 0; i < edge_cnt; ++i) {
    edges[i].weight =
        edges[i].weight * log(1.0 * max(deg[edges[i].a], deg[edges[i].b])) /
        (unweighted_distance[edges[i].a] + unweighted_distance[edges[i].b]);
  }
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

void kruskal(int node_cnt, int edge_cnt, Edge *edges, Edge *tree_edges,
             Edge *off_tree_edges) {
  ScopeTimer t_("kruskal");
  boost::sort::parallel_stable_sort(edges, edges + edge_cnt);
  t_.tick("sort edges");
  // tree_edges.reserve(node_cnt - 1);
  // off_tree_edges.reserve(edge_cnt - (node_cnt - 1));
  UnionFindSet ufs(node_cnt + 1);
  int edge_index = 0, tree_edges_size = 0, off_tree_edges_size = 0;
  for (; tree_edges_size != node_cnt - 1; ++edge_index) {
    auto edge = edges[edge_index];
    if (ufs.find_fa(edge.a) != ufs.find_fa(edge.b)) {
      ufs.merge(edge.a, edge.b);
      tree_edges[tree_edges_size++] = edge;
    } else {
      off_tree_edges[off_tree_edges_size++] = edge;
    }
  }
  t_.tick("collect mst edges");
  for (int i = edge_index, j = off_tree_edges_size; i < edge_cnt; ++i, ++j) {
    off_tree_edges[j] = edges[i];
  }
  t_.tick("copy off tree edges");
}

CSRMatrix<int> rebuild_tree(int node_cnt, int edge_cnt, Edge *tree) {
  ScopeTimer t_("rebuild_tree");
  return build_csr_matrix</* use_edge_list */ true>(node_cnt, edge_cnt, tree);
}

void euler_tour(const CSRMatrix<int> &tree, Edge *tree_edges, int cur, int fa,
                int &dfn, double *weighted_depth, int *unweighted_depth,
                int *euler_series, int *pos) {
  euler_series[dfn++] = cur;
  pos[cur] = dfn - 1;
  for (int i = tree.row_indices[cur]; i < tree.row_indices[cur + 1]; ++i) {
    const Edge &e = tree_edges[tree.neighbors[i]];
    int v = cur ^ e.a ^ e.b;
    if (v != fa) {
      weighted_depth[v] = 1.0 / e.origin_weight + weighted_depth[cur];
      unweighted_depth[v] = 1 + unweighted_depth[cur];
      euler_tour(tree, tree_edges, v, cur, dfn, weighted_depth,
                 unweighted_depth, euler_series, pos);
      euler_series[dfn++] = cur;
    }
  }
}

constexpr int ilog2(int x) {
  return (
      (unsigned)(8 * sizeof(unsigned long long) - __builtin_clzll((x)) - 1));
}

void rmq_lca(const CSRMatrix<int> &tree, Edge *tree_edges, int query_size,
             Edge *query_info, int root, int node_cnt, double *weighted_depth,
             int *unweighted_depth) {
  ScopeTimer t_("rmq_lca");
  int *euler_series(new int[(node_cnt + 1) * 2 - 1]);
  int *pos(new int[node_cnt + 1]);
  int dfn = 0;
  euler_tour(tree, tree_edges, root, root, dfn, weighted_depth,
             unweighted_depth, euler_series, pos);
  t_.tick("euler tour");
  const int block_size = sqrt(dfn);
  const int block_count = (dfn + block_size - 1) / block_size;
  std::unique_ptr<int[]> prefix_min_per_block(new int[dfn]);
  std::unique_ptr<int[]> postfix_min_per_block(new int[dfn]);
  std::unique_ptr<int[]> min_per_block(new int[block_count]);
  std::unique_ptr<int[]> contiguous_block_min(
      new int[block_count * (block_count + 1) / 2]);
  const auto idx_map = [&block_count](int i, int j) {
    const int col_idx = j - i;
    const int res = (block_count + block_count - i + 1) * i / 2 + col_idx;
    // printf("(%d, %d/%d): %d\n", i, j, block_count, res);
    return res;
  };
  t_.tick("vec init");
#pragma omp parallel for
  for (int i = 0; i < block_count; ++i) {
    const int block_start = block_size * i;
    const int block_end = min(block_start + block_size, dfn);
    prefix_min_per_block[block_start] = pos[euler_series[block_start]];
    for (int j = block_start + 1; j < block_end; ++j) {
      prefix_min_per_block[j] =
          min(pos[euler_series[j]], prefix_min_per_block[j - 1]);
    }
    postfix_min_per_block[block_end - 1] = pos[euler_series[block_end - 1]];
    for (int j = block_end - 1; j > block_start; --j) {
      postfix_min_per_block[j - 1] =
          min(pos[euler_series[j - 1]], postfix_min_per_block[j]);
    }
    min_per_block[i] = postfix_min_per_block[block_start];
  }
#pragma omp parallel for
  for (int i = 0; i < block_count; ++i) {
    contiguous_block_min[idx_map(i, i)] = min_per_block[i];
    for (int j = i + 1; j < block_count; ++j) {
      contiguous_block_min[idx_map(i, j)] =
          min(min_per_block[j], contiguous_block_min[idx_map(i, j - 1)]);
    }
  }
  t_.tick("lca preprocess");
#pragma omp parallel for
  for (int i = 0; i < query_size; ++i) {
    Edge &e = query_info[i];
    const int l = min(pos[e.a], pos[e.b]);
    const int r = max(pos[e.a], pos[e.b]);
    const int l_block = l / block_size;
    const int r_block = r / block_size;
    int lca_pos = dfn;
    if (l_block == r_block) {
      for (int j = l; j <= r; ++j) {
        if (pos[euler_series[j]] < lca_pos) {
          lca_pos = pos[euler_series[j]];
        }
      }
    } else if (l_block + 1 == r_block) {
      lca_pos = min(postfix_min_per_block[l], prefix_min_per_block[r]);
    } else {
      lca_pos = min(contiguous_block_min[idx_map(l_block + 1, r_block - 1)],
                    min(postfix_min_per_block[l], prefix_min_per_block[r]));
    }
    const int lca = euler_series[lca_pos];
    e.lca = lca;
  }
  t_.tick("lca query");
  delete[] euler_series;
  delete[] pos;
}

void sort_off_tree_edges(int edges_size, Edge *edges, double *depth) {
  ScopeTimer t_("sort_off_tree_edges");
#pragma omp parallel for
  for (int i = 0; i < edges_size; ++i) {
    auto &e = edges[i];
    e.weight = e.origin_weight * (depth[e.a] + depth[e.b] - 2 * depth[e.lca]);
  }
  t_.tick("map edges weight");
#ifdef DEBUG
  puts("unsorted_off_tree_edges");
  for (auto &x : edges) {
    printf("%d %d %lf %lf\n", x.a, x.b, x.weight, x.origin_weight);
  }
#endif
  boost::sort::parallel_stable_sort(edges, edges + edges_size);
  t_.tick("sort edges");
}

void mark_ban_edges(vector<bool> &ban, const vector<int> &ban_edges) {
  for (auto &x : ban_edges) {
    ban[x] = true;
  }
}
struct QueueEntry {
  int node, layer, predecessor;
};

int beta_layer_bfs_1(int start, vector<QueueEntry> &q,
                     const CSRMatrix<int> &tree, vector<bool> &black_list1,
                     int beta) {
  int rear = 0;
  q[rear++] = {start, 0, -1};
  for (int idx = 0; idx < rear; idx++) {
    int cur_node = q[idx].node;
    int cur_layer = q[idx].layer;
    int cur_pre = q[idx].predecessor;
    black_list1[cur_node] = true;
    if (cur_layer == beta) {
      continue;
    }
    for (int j = tree.row_indices[cur_node]; j < tree.row_indices[cur_node + 1];
         ++j) {
      int v = tree.neighbors[j];
      if (cur_pre != v) {
        q[rear++] = {v, cur_layer + 1, cur_node};
      }
    }
  }
  return rear;
}

const auto pair_hash = [](const std::pair<int, int> &p) {
  return p.first * 31 + p.second;
};

int beta_layer_bfs_2(
    int start, vector<QueueEntry> &q, const CSRMatrix<int> &tree,
    const CSRMatrix<int> &off_tree_graph, vector<bool> &black_list1, int beta,
    unordered_set<pair<int, int>, decltype(pair_hash)> &banned_edges,
    vector<unsigned short> &banned_nodes) {
  static unsigned short id = 1;
  if (id == 0) {
    id++;
  }
  int rear = 0;
  q[rear++] = {start, 0, -1};
  for (int idx = 0; idx < rear; idx++) {
    int cur_node = q[idx].node;
    int cur_layer = q[idx].layer;
    int cur_pre = q[idx].predecessor;
    for (int j = off_tree_graph.row_indices[cur_node];
         j < off_tree_graph.row_indices[cur_node + 1]; ++j) {
      int v = off_tree_graph.neighbors[j];
      if (black_list1[v]) {
        if (banned_nodes[cur_node] == 0 &&
            banned_nodes[cur_node] == banned_nodes[v]) {
          banned_nodes[cur_node] = banned_nodes[v] = id++;
        } else {
          banned_edges.insert({min(cur_node, v), max(cur_node, v)});
        }
      }
    }
    if (cur_layer == beta) {
      continue;
    }
    for (int j = tree.row_indices[cur_node]; j < tree.row_indices[cur_node + 1];
         ++j) {
      int v = tree.neighbors[j];
      if (v != cur_pre) {
        q[rear++] = {v, cur_layer + 1, cur_node};
      }
    }
  }
  return rear;
}

vector<int> add_off_tree_edges(const int node_cnt, const int tree_edges_size,
                               Edge *tree_edges, const int off_tree_edges_size,
                               Edge *off_tree_edges, int *depth) {

  ScopeTimer t_("add_off_tree_edges");
  // vector<vector<int>> rebuilt_off_tree_graph(node_cnt + 1);
  // vector<vector<int>> rebuilt_tree_graph(node_cnt + 1);
  // for (int i = 0; i < off_tree_edges.size(); ++i) {
  //   auto &e = off_tree_edges[i];
  //   rebuilt_off_tree_graph[e.a].push_back(e.b);
  //   rebuilt_off_tree_graph[e.b].push_back(e.a);
  // }
  // for (int i = 0; i < tree_edges.size(); ++i) {
  //   auto &e = tree_edges[i];
  //   rebuilt_tree_graph[e.a].push_back(e.b);
  //   rebuilt_tree_graph[e.b].push_back(e.a);
  // }
  auto off_tree_graph =
      build_csr_matrix(node_cnt, off_tree_edges_size, off_tree_edges);
  auto tree_graph = build_csr_matrix(node_cnt, tree_edges_size, tree_edges);
  vector<int> edges_to_be_add;
  vector<bool> black_list1(node_cnt + 1, false);
  vector<unsigned short> banned_nodes(node_cnt + 1, 0);
  unordered_set<pair<int, int>, decltype(pair_hash)> banned_edges(1024,
                                                                  pair_hash);
  vector<QueueEntry> q1(node_cnt), q2(node_cnt);

  int alpha = max(int(off_tree_edges_size / 25), 2);
  for (int i = 0; i < off_tree_edges_size; ++i) {
    if (edges_to_be_add.size() == alpha) {
      break;
    }
    auto &e = off_tree_edges[i];
    if (banned_nodes[e.a] != 0 && banned_nodes[e.a] == banned_nodes[e.b]) {
      continue;
    }
    if ((banned_nodes[e.a] != 0 || banned_nodes[e.b] != 0) &&
        banned_edges.find({min(e.a, e.b), max(e.a, e.b)}) !=
            banned_edges.end()) {
      continue;
    }
    edges_to_be_add.push_back(i);
    int beta = min(depth[e.a], depth[e.b]) - depth[e.lca];
#ifdef DEBUG
    printf("beta: %d, (%d, %d)\n", beta, e.a, e.b);
#endif
    int size1 = beta_layer_bfs_1(e.a, q1, tree_graph, black_list1, beta);
    beta_layer_bfs_2(e.b, q2, tree_graph, off_tree_graph, black_list1, beta,
                     banned_edges, banned_nodes);
    for (int j = 0; j < size1; j++) {
      black_list1[q1[j].node] = false;
    }
  }
  // magic_trace_stop_indicator();
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
  assert(!(L & 1));

  Edge *origin_edges(new Edge[L >> 1]);
  int *degree(new int[M + 1]());
  double *volume(new double[M + 1]());
  int edges_cnt = 0;
  for (int i = 0; i < L; ++i) {
    int f, t;
    double w;
    fin >> f >> t >> w;
    if (f > t) {
      continue;
    }
    volume[f] += w;
    volume[t] += w;
    degree[f]++;
    degree[t]++;
    origin_edges[edges_cnt++] = Edge{f, t, w, w};
  }
  auto G =
      build_csr_matrix</* use_edge_list */ true>(M, edges_cnt, origin_edges);
  fin.close();
  printf("edge_cnt: %d\n", edges_cnt);
  /**************************************************/
  /***************** Start timing *******************/
  /**************************************************/
  struct timeval start, end;
  gettimeofday(&start, NULL);
  int r_node = largest_volume_node(M, volume);
  auto unweighted_distance =
      get_unweighted_distance_bfs(origin_edges, G, r_node, M);
  get_new_edges(edges_cnt, origin_edges, degree, unweighted_distance);

  auto new_edges = origin_edges;

  const int tree_edges_size = M - 1;
  const int off_tree_edges_size = edges_cnt - tree_edges_size;
  Edge *tree_edges(new Edge[tree_edges_size]),
      *off_tree_edges(new Edge[off_tree_edges_size]);
  kruskal(M, edges_cnt, new_edges, tree_edges, off_tree_edges);

#ifdef DEBUG
  puts("kruscal results: ");
  for (auto &x : tree_edges) {
    printf("%d %d %lf %lf\n", x.a, x.b, x.weight, x.origin_weight);
  }
#endif

  auto new_tree = rebuild_tree(M, tree_edges_size, tree_edges);
  double *tree_weighted_depth(new double[M + 1]);
  int *tree_unweighted_depth(new int[M + 1]);
  rmq_lca(new_tree, tree_edges, off_tree_edges_size, off_tree_edges, r_node, M,
          tree_weighted_depth, tree_unweighted_depth);
  // tarjan_lca(new_tree, tree_edges, off_tree_edges, r_node, M,
  //            tree_weighted_depth, tree_unweighted_depth);

  sort_off_tree_edges(off_tree_edges_size, off_tree_edges, tree_weighted_depth);

#ifdef DEBUG
  puts("sorted off_tree_edges: ");
  for (auto &x : off_tree_edges) {
    printf("%d %d %.16lf %.16lf\n", x.a, x.b, x.weight, x.origin_weight);
  }
#endif

  vector<int> res =
      add_off_tree_edges(M, tree_edges_size, tree_edges, off_tree_edges_size,
                         off_tree_edges, tree_unweighted_depth);
  gettimeofday(&end, NULL);
  printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000 +
                                     (end.tv_usec - start.tv_usec) / 1000.0);
  /**************************************************/
  /******************* End timing *******************/
  /**************************************************/
  FILE *out = fopen("result.ng.txt", "w");
  for (int i = 0; i < tree_edges_size; i++) {
    fprintf(out, "%d %d\n", int(tree_edges[i].a), int(tree_edges[i].b));
  }
  for (int i = 0; i < res.size(); i++) {
    fprintf(out, "%d %d\n", int(off_tree_edges[res[i]].a),
            int(off_tree_edges[res[i]].b));
  }
  fclose(out);
  delete[] origin_edges;
  delete[] volume;
  delete[] degree;
  delete[] tree_edges;
  delete[] off_tree_edges;
  delete[] unweighted_distance;
  delete[] tree_weighted_depth;
  delete[] tree_unweighted_depth;
}
