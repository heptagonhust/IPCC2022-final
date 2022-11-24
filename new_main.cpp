#include "SPSCQueue.h"
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
#include <oneapi/tbb.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/rw_mutex.h>
#include <oneapi/tbb/spin_mutex.h>
#include <queue>
#include <stack>
#include <string>
#include <sys/time.h>
#include <unordered_set>
#include <vector>
using namespace std;
using namespace rigtorp;

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
  std::unique_ptr<int[]> row_indices;
  std::unique_ptr<T[]> neighbors;

  CSRMatrix(const int &matrix_dim, const int &nonzero_count)
      : row_indices(new int[matrix_dim + 1]), neighbors(new T[nonzero_count]) {}
};

using vec2 = pair<int, int>;
using vec3 = struct Triple { int first, second, third; };

enum class CSRValueType {
  index,
  neighbor,
  neighbor_and_index,
};

template <CSRValueType value_type = CSRValueType::index, class T = int>
CSRMatrix<T> build_csr_matrix(const int &node_cnt, const int &edge_cnt,
                              const Edge *edges) {
  vector<int> deg(node_cnt + 1);
  for (int i = 0; i < edge_cnt; ++i) {
    deg[edges[i].a]++;
    deg[edges[i].b]++;
  }
  auto mat = CSRMatrix<T>(node_cnt + 1, edge_cnt * 2);
  int start_idx = 0;
  for (int i = 1; i <= node_cnt; i++) {
    mat.row_indices[i] = start_idx;
    start_idx += deg[i];
  }
  mat.row_indices[0] = 0;
  mat.row_indices[node_cnt + 1] = start_idx;
  for (int i = 0; i < edge_cnt; ++i) {
    const auto &e = edges[i];
    if constexpr (value_type == CSRValueType::index) {
      mat.neighbors[mat.row_indices[e.a]] = T{i};
      mat.neighbors[mat.row_indices[e.b]] = T{i};
    } else if constexpr (value_type == CSRValueType::neighbor) {
      mat.neighbors[mat.row_indices[e.a]] = T{e.b};
      mat.neighbors[mat.row_indices[e.b]] = T{e.a};
    } else {
      mat.neighbors[mat.row_indices[e.a]] = T{e.b, i};
      mat.neighbors[mat.row_indices[e.b]] = T{e.a, i};
    }
    mat.row_indices[e.a]++;
    mat.row_indices[e.b]++;
  }
  for (int i = node_cnt; i >= 1; i--) {
    mat.row_indices[i] = mat.row_indices[i - 1];
  }
  return mat;
}

int largest_volume_node(const int &node_cnt, const double *volume) {
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

unique_ptr<int[]> get_unweighted_distance_bfs(const Edge *edges,
                                              const CSRMatrix<int> &G,
                                              int start, int node_cnt) {
  ScopeTimer t_("get_unweighted_distance_bfs");
  unique_ptr<int[]> res(new int[node_cnt + 1]{});
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

void get_new_edges(const int &edge_cnt, Edge *edges, const int *deg,
                   const int *unweighted_distance) {
  ScopeTimer t_("get_new_edges");
  tbb::parallel_for(0, edge_cnt, [edges, deg, unweighted_distance](auto i) {
    edges[i].weight =
        edges[i].weight * log(1.0 * max(deg[edges[i].a], deg[edges[i].b])) /
        (unweighted_distance[edges[i].a] + unweighted_distance[edges[i].b]);
  });
}

struct UnionFindSet {
  unique_ptr<int[]> fa;
  UnionFindSet(size_t sz) : fa(new int[sz]) {
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

CSRMatrix<vec2> rebuild_tree(int node_cnt, int edge_cnt, const Edge *tree) {
  ScopeTimer t_("rebuild_tree");
  return build_csr_matrix</* use_edge_list */ CSRValueType::neighbor_and_index,
                          vec2>(node_cnt, edge_cnt, tree);
}

constexpr int parallel_depth = 7;
void get_subtree_size(const CSRMatrix<vec2> &tree, int cur, int fa, int depth,
                      int *unweighted_depth, int *subtree_size) {
  for (int i = tree.row_indices[cur]; i < tree.row_indices[cur + 1]; ++i) {
    int v = tree.neighbors[i].first;
    if (v != fa) {
      unweighted_depth[v] = 1 + unweighted_depth[cur];
      get_subtree_size(tree, v, cur, depth + 1, unweighted_depth, subtree_size);
      subtree_size[cur] += subtree_size[v];
    }
  }
}

void parallel_get_subtree_size(const CSRMatrix<vec2> &tree, int cur, int fa,
                               int depth, int *unweighted_depth,
                               int *subtree_size) {
  unweighted_depth[cur] = depth;
  if (depth < parallel_depth) {
    oneapi::tbb::task_group g;
    oneapi::tbb::spin_mutex lock;
    for (int i = tree.row_indices[cur]; i < tree.row_indices[cur + 1]; ++i) {
      int v = tree.neighbors[i].first;
      if (v != fa) {
        g.run([&tree, v, cur, unweighted_depth, subtree_size, depth, &lock] {
          parallel_get_subtree_size(tree, v, cur, depth + 1, unweighted_depth,
                                    subtree_size);
          lock.lock();
          subtree_size[cur] += subtree_size[v];
          lock.unlock();
        });
      }
    }
    g.wait();
  } else {
    get_subtree_size(tree, cur, fa, depth, unweighted_depth, subtree_size);
  }
}

void euler_tour(const CSRMatrix<vec2> &tree, const Edge *tree_edges, int cur,
                int fa, int &dfn, int *euler_series, int *pos,
                double *weighted_depth) {
  euler_series[dfn] = cur;
  pos[cur] = dfn;
  dfn++;
  for (int i = tree.row_indices[cur]; i < tree.row_indices[cur + 1]; ++i) {
    const Edge &e = tree_edges[tree.neighbors[i].second];
    int v = tree.neighbors[i].first;
    if (v != fa) {
      weighted_depth[v] = 1.0 / e.origin_weight + weighted_depth[cur];
      euler_tour(tree, tree_edges, v, cur, dfn, euler_series, pos,
                 weighted_depth);
      euler_series[dfn++] = cur;
    }
  }
}

void parallel_euler_tour(const CSRMatrix<vec2> &tree, const Edge *tree_edges,
                         int cur, int fa, const int *subtree_size,
                         const int *unweighted_depth, int dfn,
                         int *euler_series, int *pos, double *weighted_depth) {
  if (subtree_size[cur] >= 512 && subtree_size[cur] > subtree_size[fa] * 0.5) {
    euler_series[dfn] = cur;
    pos[cur] = dfn;
    dfn++;
    oneapi::tbb::task_group g;
    for (int i = tree.row_indices[cur]; i < tree.row_indices[cur + 1]; ++i) {
      const Edge &e = tree_edges[tree.neighbors[i].second];
      int v = tree.neighbors[i].first;
      if (v != fa) {
        g.run([&tree, tree_edges, v, cur, subtree_size, unweighted_depth, dfn,
               euler_series, pos, weighted_depth, e]() {
          weighted_depth[v] = 1.0 / e.origin_weight + weighted_depth[cur];
          parallel_euler_tour(tree, tree_edges, v, cur, subtree_size,
                              unweighted_depth, dfn, euler_series, pos,
                              weighted_depth);
        });
        dfn += subtree_size[v] * 2 - 1;
        euler_series[dfn++] = cur;
      }
    }
    g.wait();
  } else {
    euler_tour(tree, tree_edges, cur, fa, dfn, euler_series, pos,
               weighted_depth);
  }
}

constexpr int ilog2(int x) {
  return (
      (unsigned)(8 * sizeof(unsigned long long) - __builtin_clzll((x)) - 1));
}

void rmq_lca(const CSRMatrix<vec2> &tree, Edge *tree_edges, int query_size,
             Edge *query_info, int root, int node_cnt, double *weighted_depth,
             int *unweighted_depth) {
  ScopeTimer t_("rmq_lca");
  const int euler_series_len = 2 * node_cnt - 1;
  std::unique_ptr<int[]> euler_series(new int[euler_series_len]);
  std::fill_n(euler_series.get(), euler_series_len, -1);
  std::unique_ptr<int[]> pos(new int[node_cnt + 1]);
  std::unique_ptr<int[]> subtree_size(new int[node_cnt + 1]);
  std::fill_n(subtree_size.get(), node_cnt + 1, 1);
  parallel_get_subtree_size(tree, root, root, 0, unweighted_depth,
                            subtree_size.get());
  t_.tick("get subtree size");
  parallel_euler_tour(tree, tree_edges, root, root, subtree_size.get(),
                      unweighted_depth, 0, euler_series.get(), pos.get(),
                      weighted_depth);
  t_.tick("parallel euler tour");
  const int block_size = sqrt(euler_series_len);
  const int block_count = (euler_series_len + block_size - 1) / block_size;
  std::unique_ptr<int[]> prefix_min_per_block(new int[euler_series_len]);
  std::unique_ptr<int[]> postfix_min_per_block(new int[euler_series_len]);
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
  tbb::parallel_for(
      0, block_count,
      [block_size, &prefix_min_per_block, euler_series_len, &pos, &euler_series,
       &postfix_min_per_block, &min_per_block](auto i) {
        const int block_start = block_size * i;
        const int block_end = min(block_start + block_size, euler_series_len);
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
      });
  tbb::parallel_for(
      0, block_count,
      [&contiguous_block_min, &min_per_block, idx_map, block_count](auto i) {
        contiguous_block_min[idx_map(i, i)] = min_per_block[i];
        for (int j = i + 1; j < block_count; ++j) {
          contiguous_block_min[idx_map(i, j)] =
              min(min_per_block[j], contiguous_block_min[idx_map(i, j - 1)]);
        }
      });
  t_.tick("lca preprocess");
  tbb::parallel_for(
      0, query_size,
      [query_info, &pos, block_size, euler_series_len, &euler_series,
       &postfix_min_per_block, &prefix_min_per_block, &contiguous_block_min,
       idx_map](auto i) {
        Edge &e = query_info[i];
        const int l = min(pos[e.a], pos[e.b]);
        const int r = max(pos[e.a], pos[e.b]);
        const int l_block = l / block_size;
        const int r_block = r / block_size;
        int lca_pos = euler_series_len;
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
      });
  t_.tick("lca query");
}

void sort_off_tree_edges(int edges_cnt, Edge *edges, const double *depth) {
  ScopeTimer t_("sort_off_tree_edges");
  tbb::parallel_for(0, edges_cnt, [edges, depth](auto i) {
    auto &e = edges[i];
    e.weight = e.origin_weight * (depth[e.a] + depth[e.b] - 2 * depth[e.lca]);
  });
  t_.tick("map edges weight");
#ifdef DEBUG
  puts("unsorted_off_tree_edges");
  for (auto &x : edges) {
    printf("%d %d %lf %lf\n", x.a, x.b, x.weight, x.origin_weight);
  }
#endif
  boost::sort::parallel_stable_sort(edges, edges + edges_cnt);
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

int beta_layer_bfs_1(int start, unique_ptr<QueueEntry[]> &q,
                     const CSRMatrix<int> &tree,
                     unique_ptr<bool[]> &black_list1, int beta) {
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

vector<vec3> beta_layer_bfs_2(int start, unique_ptr<QueueEntry[]> &q,
                              const CSRMatrix<int> &tree,
                              const CSRMatrix<vec2> &off_tree_graph,
                              unique_ptr<bool[]> &black_list1, int beta) {
  vector<vec3> res{};
  int rear = 0;
  q[rear++] = {start, 0, -1};
  for (int idx = 0; idx < rear; idx++) {
    int cur_node = q[idx].node;
    int cur_layer = q[idx].layer;
    int cur_pre = q[idx].predecessor;
    for (int j = off_tree_graph.row_indices[cur_node];
         j < off_tree_graph.row_indices[cur_node + 1]; ++j) {
      int v = off_tree_graph.neighbors[j].first;
      int edge_idx = off_tree_graph.neighbors[j].second;
      if (black_list1[v]) {
        // if (banned_nodes[cur_node] == 0 &&
        //     banned_nodes[cur_node] == banned_nodes[v]) {
        //   banned_nodes[cur_node] = banned_nodes[v] = id++;
        // } else {
        //   banned_edges.insert({min(cur_node, v), max(cur_node, v)});
        // }
        res.push_back({cur_node, v, edge_idx});
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
  return res;
}

const int num_producer = 16;

void produce_ban_off_tree_edges(
    int thread_id, SPSCQueue<vector<vec3>> &spsc, int node_cnt,
    CSRMatrix<int> &tree_graph, CSRMatrix<vec2> &off_tree_graph,
    const Edge *tree_edges, const int off_tree_edges_size,
    const Edge *off_tree_edges, const int *depth, atomic<bool> &done,
    const atomic<int> *banned_nodes, const atomic<bool> *edge_is_banned) {
  unique_ptr<QueueEntry[]> q1(new QueueEntry[node_cnt]);
  unique_ptr<QueueEntry[]> q2(new QueueEntry[node_cnt]);
  unique_ptr<bool[]> black_list1(new bool[node_cnt + 1]());
  for (int i = thread_id; i < off_tree_edges_size && !done; i += num_producer) {
    auto &e = off_tree_edges[i];
    if (banned_nodes[e.a] != 0 && banned_nodes[e.a] == banned_nodes[e.b]) {
      spsc.emplace(vector<vec3>{});
      continue;
    }

    if (edge_is_banned[i]) {
      spsc.emplace(vector<vec3>{});
      continue;
    }
    int beta = min(depth[e.a], depth[e.b]) - depth[e.lca];
    int size1 = beta_layer_bfs_1(e.a, q1, tree_graph, black_list1, beta);
    auto ban_list = beta_layer_bfs_2(e.b, q2, tree_graph, off_tree_graph,
                                     black_list1, beta);
    spsc.emplace(ban_list);

    for (int j = 0; j < size1; j++) {
      black_list1[q1[j].node] = false;
    }
  }
}

vector<int> add_off_tree_edges(const int node_cnt, const int tree_edges_size,
                               const Edge *tree_edges,
                               const int off_tree_edges_size,
                               const Edge *off_tree_edges, const int *depth) {

  ScopeTimer t_("add_off_tree_edges");
  auto off_tree_graph =
      build_csr_matrix<CSRValueType::neighbor_and_index, vec2>(
          node_cnt, off_tree_edges_size, off_tree_edges);
  auto tree_graph = build_csr_matrix<CSRValueType::neighbor>(
      node_cnt, tree_edges_size, tree_edges);
  t_.tick("build graph");
  int alpha = max(int(off_tree_edges_size / 25), 2);

  vector<int> edges_to_be_add(alpha);
  unique_ptr<atomic<int>[]> banned_nodes(new atomic<int>[node_cnt + 1]());
  unique_ptr<atomic<bool>[]> edge_is_banned(
      new atomic<bool>[off_tree_edges_size] {});

  std::thread threads[num_producer];
  SPSCQueue<vector<vec3>> spscs[num_producer];
  int mark_color = 1;
  atomic<bool> threads_done{false};
  for (int i = 0; i < num_producer; i++) {
    threads[i] = std::thread([&, i] {
      produce_ban_off_tree_edges(
          i, spscs[i], node_cnt, tree_graph, off_tree_graph, tree_edges,
          off_tree_edges_size, off_tree_edges, depth, threads_done,
          banned_nodes.get(), edge_is_banned.get());
    });
  }

  int edges_added_size = 0;
  for (int i = 0; i < off_tree_edges_size; ++i) {
    if (edges_added_size == alpha) {
      break;
    }
    auto &e = off_tree_edges[i];
    // if (banned_nodes[e.a] != 0 && banned_nodes[e.a] == banned_nodes[e.b]) {
    //   continue;
    // }
    // if (banned_edges.find({min(e.a, e.b), max(e.a, e.b)}) !=
    //     banned_edges.end()) {
    //   continue;
    // }
    while (!spscs[i % num_producer].front())
      ;

    bool is_banned =
        (banned_nodes[e.a] != 0 && banned_nodes[e.a] == banned_nodes[e.b]) ||
        edge_is_banned[i];

    if (is_banned) {
      spscs[i % num_producer].pop();
      continue;
    }

    auto ban_list = spscs[i % num_producer].front();

    edges_to_be_add[edges_added_size++] = i;
    for (auto &[u, v, edge_idx] : *ban_list) {
      if (banned_nodes[u] == 0 &&
          banned_nodes[v] == 0) { // both remain unpainted
        banned_nodes[u] = banned_nodes[v] = mark_color;
        mark_color++;
      } else { // each is painted
        edge_is_banned[edge_idx] = true;
      }
    }
    spscs[i % num_producer].pop();
  }
  threads_done = true;
  for (auto &q : spscs) {
    if (q.front() != nullptr)
      q.pop();
  }
  for (int i = 0; i < num_producer; ++i) {
    threads[i].join();
  }

  // magic_trace_stop_indicator();
  edges_to_be_add.resize(edges_added_size);
  return edges_to_be_add;
}

int main(int argc, const char *argv[]) {
  oneapi::tbb::global_control global_limit(
      tbb::global_control::thread_stack_size, 16 * 1024 * 1024);
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

  unique_ptr<Edge[]> origin_edges(new Edge[L >> 1]);
  unique_ptr<int[]> degree(new int[M + 1]());
  unique_ptr<double[]> volume(new double[M + 1]());
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
      build_csr_matrix<CSRValueType::index>(M, edges_cnt, origin_edges.get());
  fin.close();
  printf("edge_cnt: %d\n", edges_cnt);
  /**************************************************/
  /***************** Start timing *******************/
  /**************************************************/
  struct timeval start, end;
  gettimeofday(&start, NULL);
  int r_node = largest_volume_node(M, volume.get());
  auto unweighted_distance =
      get_unweighted_distance_bfs(origin_edges.get(), G, r_node, M);
  get_new_edges(edges_cnt, origin_edges.get(), degree.get(),
                unweighted_distance.get());

  auto new_edges = std::move(origin_edges);

  const int tree_edges_size = M - 1;
  const int off_tree_edges_size = edges_cnt - tree_edges_size;
  unique_ptr<Edge[]> tree_edges(new Edge[tree_edges_size]),
      off_tree_edges(new Edge[off_tree_edges_size]);
  kruskal(M, edges_cnt, new_edges.get(), tree_edges.get(),
          off_tree_edges.get());

#ifdef DEBUG
  puts("kruscal results: ");
  for (auto &x : tree_edges) {
    printf("%d %d %lf %lf\n", x.a, x.b, x.weight, x.origin_weight);
  }
#endif

  auto new_tree = rebuild_tree(M, tree_edges_size, tree_edges.get());
  unique_ptr<double[]> tree_weighted_depth(new double[M + 1]);
  unique_ptr<int[]> tree_unweighted_depth(new int[M + 1]);
  rmq_lca(new_tree, tree_edges.get(), off_tree_edges_size, off_tree_edges.get(),
          r_node, M, tree_weighted_depth.get(), tree_unweighted_depth.get());
  // tarjan_lca(new_tree, tree_edges, off_tree_edges, r_node, M,
  //            tree_weighted_depth, tree_unweighted_depth);

  sort_off_tree_edges(off_tree_edges_size, off_tree_edges.get(),
                      tree_weighted_depth.get());

#ifdef DEBUG
  puts("sorted off_tree_edges: ");
  for (auto &x : off_tree_edges) {
    printf("%d %d %.16lf %.16lf\n", x.a, x.b, x.weight, x.origin_weight);
  }
#endif

  vector<int> res = add_off_tree_edges(
      M, tree_edges_size, tree_edges.get(), off_tree_edges_size,
      off_tree_edges.get(), tree_unweighted_depth.get());
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
}
