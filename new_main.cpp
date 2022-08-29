#include <algorithm>
#include <cstddef>
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
  bool operator<(const Edge &rhs) const { return weight >= rhs.weight; }
};

int largest_volume_node(const vector<int> &volume) {
  return max_element(volume.begin(), volume.end()) - volume.begin();
}

vector<int> get_unweighted_distance_bfs(const vector<Edge> &edges,
                                        const vector<vector<int>> &G,
                                        int start) {
  vector<int> res(G.size());
  vector<bool> vis(G.size());
  queue<int> q;
  q.push(start);
  while (!q.empty()) {
    int top = q.front();
    q.pop();
    for (int i = 0; i < G[top].size(); ++i) {
      Edge e = edges[G[top][i]];
      if (e.b == i) {
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
  vector<Edge> res(edges.size());
  for (int i = 0; i < edges.size(); ++i) {
    res[i].a = edges[i].a;
    res[i].b = edges[i].b;
    res[i].weight =
        log(max(deg[edges[i].a], deg[edges[i].b])) /
        (unweighted_distance[edges[i].a], unweighted_distance[edges[i].b]);
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
      return find_fa(fa[x]); // 如果不是则 x 的爸爸问 x 的爷爷
  }

  void merge(int a, int b) { fa[a] = fa[b]; }
};

void kruskal(int node_cnt, vector<Edge> &edges, vector<Edge> &tree_edges,
             vector<Edge> &off_tree_edges) {
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

  vector<int> degree(M + 1), volume(M + 1);
  for (int i = 0; i < L; ++i) {
    int f, t;
    double w;
    fin >> f >> t >> w;
    G[f].push_back(origin_edges.size());
    G[t].push_back(origin_edges.size());
    volume[f] += w;
    volume[t] += w;
    degree[f]++;
    degree[t]++;
    origin_edges.push_back(Edge{f, t, w});
  }
  fin.close();

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

  gettimeofday(&end, NULL);
  printf("Using time : %f ms\n", (end.tv_sec - start.tv_sec) * 1000 +
                                     (end.tv_usec - start.tv_usec) / 1000.0);
  /**************************************************/
  /******************* End timing *******************/
  /**************************************************/
}
