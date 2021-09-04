#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <ctime>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "read_data.cpp"
#include "motif_id.cpp"

using namespace std;

class hmotif{
public:
	int e_a, e_b, e_c;
	int size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc;
	hmotif(int e_a, int e_b, int e_c,
	       int size_a, int size_b, int size_c, int C_ab, int C_bc, int C_ca, int g_abc){
		this->e_a = e_a;
		this->e_b = e_b;
		this->e_c = e_c;
		this->size_a = size_a;
		this->size_b = size_b;
		this->size_c = size_c;
		this->C_ab = C_ab;
		this->C_bc = C_bc;
		this->C_ca = C_ca;
		this->g_abc = g_abc;
	}
};

string hyperedge2hash(vector<int>& nodes)
{
	stringstream ss;
	sort(nodes.begin(), nodes.end());
	for (int i = 0; i < nodes.size(); i++){
		ss << nodes[i];
		if (i < nodes.size() - 1)
			ss << ",";
	}
	return ss.str();
}

bool count_temporal_hmotif(hmotif& instance, vector<long long>& motif2count, vector< vector<double> >& E2times, double delta)
{
	int comb2[6][3] = {{0,0,1}, {0,1,0}, {1,0,0}, {1,1,0}, {1,0,1}, {0,1,1}};
	int comb2_C[6][3] = {{0,1,1}, {1,1,0}, {1,0,1}, {2,1,1}, {1,1,2}, {1,2,1}};
	int comb3[6][3] = {{0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}};
	int comb3_C[6][3] = {{0,1,2}, {2,1,0}, {0,2,1}, {1,2,0}, {2,0,1}, {1,0,2}};
	
	int e_a = instance.e_a, e_b = instance.e_b, e_c = instance.e_c;
	int S[3] = {instance.size_a, instance.size_b, instance.size_c};
	int C[3] = {instance.C_ab, instance.C_bc, instance.C_ca};
	int G = instance.g_abc;

	int counter[4][4][4] = {0};
	bool counted = false;

	vector<pair<double, int> > edges;
	for (const double &timestamp: E2times[e_a]) edges.push_back({timestamp, 1});
	for (const double &timestamp: E2times[e_b]) edges.push_back({timestamp, 2});
	for (const double &timestamp: E2times[e_c]) edges.push_back({timestamp, 3});
	sort(edges.begin(), edges.end());

	int start = 0;
	for (int end = 0; end < edges.size(); end++){
		while (edges[start].first + delta < edges[end].first){
			// decrement
			counter[edges[start].second][0][0] --;
			for (int i = 1; i <= 3; i++){
				counter[edges[start].second][i][0] -= counter[i][0][0];
			}
			start ++;
		}
		// increment
		for (int i = 1; i <= 3; i++){
			for (int j = 1; j <= 3; j++){
				counter[i][j][edges[end].second] += counter[i][j][0];
			}
			counter[i][edges[end].second][0] += counter[i][0][0];
		}
		counter[edges[end].second][0][0] ++;
	}

	for (int i = 0; i < 6; i++){
		int a = comb3[i][0], b = comb3[i][1], c = comb3[i][2];
		int a_C = comb3_C[i][0], b_C = comb3_C[i][1], c_C = comb3_C[i][2];
		int motif_index = get_motif_index(S[a], S[b], S[c], C[a_C], C[b_C], C[c_C], G, 1, 2, 3);
		motif2count[motif_index] += counter[a+1][b+1][c+1];
		if (counter[a+1][b+1][c+1]) counted = true;
	}
	return counted;
}

int main(int argc, char *argv[])
{
	clock_t start;

	// Configuration
	string dataset = argv[1];
	double delta = stod(argv[2]);

	// File path
	string graphFile = "../data/" + dataset + "/" + dataset + ".txt";
	string writeFile = "../results/" + dataset + "_" + to_string(delta) + "_dp.txt";
	
	// Read data
	start = clock();
	vector< vector<int> > node2hyperedge;
	vector< vector<int> > hyperedge2node;
	vector< unordered_set<int> > hyperedge2node_set;
	vector<double> hyperedge2time;
	read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set, hyperedge2time);

	int V = (int)node2hyperedge.size();
	int E = (int)hyperedge2node.size();
	cout << "# of nodes:\t\t" << V << endl;
	cout << "# of hyperedges:\t" << E << endl;
	cout << "Reading data done: "
	     << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "----------------------------------------------------------" << endl << endl;

	// Induced static hypergraph: preprocessing
	start = clock();
	unordered_map<string, int> hash2index;
	vector< vector<int> > V2E; V2E.resize(V);
	vector< vector<int> > E2V;
	vector< unordered_set<int> > E2V_set;
	vector< vector<double> > E2times;

	for (int i = 0; i < E; i++){
		string hash = hyperedge2hash(hyperedge2node[i]);
		if (hash2index.find(hash) == hash2index.end()){
			hash2index[hash] = (int)hash2index.size();
			vector<int> nodes;
			unordered_set<int> nodes_set;
			for (const int &v: hyperedge2node[i]){
				nodes.push_back(v);
				nodes_set.insert(v);
				V2E[v].push_back(hash2index[hash]);
			}
			E2V.push_back(nodes);
			E2V_set.push_back(nodes_set);
			E2times.push_back(vector<double>());
		}
		int index = hash2index[hash];
		E2times[index].push_back(hyperedge2time[i]);
	}
	int E_static = (int)E2V.size();

	for (int e = 0; e < E_static; e++)
		sort(E2times[e].begin(), E2times[e].end());

	cout << "# of static hyperedges:\t" << E_static << endl;

	cout << "Static hypergraph (1): preprocessing done:\t"
	     << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;

	// Induced static hypergraph: adjacency list
	start = clock();
	vector< vector<pair<int, int> > > E_adj;
	vector< unordered_map<int, int> > E_inter;
	vector<long long> upd_time(E_static, -1LL);
	E_adj.resize(E_static);
	E_inter.resize(E_static);
		
	for (int e_a = 0; e_a < E_static; e_a++){
		long long l_e_a = (long long)e_a;
		for (const int &v: E2V[e_a]){
			for (const int e_b: V2E[v]){
				if (e_a == e_b) continue;
				if ((upd_time[e_b] >> 31) ^ e_a){
					upd_time[e_b] = (l_e_a << 31) + (long long)E_adj[e_b].size();
					E_adj[e_b].push_back({e_a, 0});
				}
				E_adj[e_b][(int)(upd_time[e_b] & 0x7FFFFFFFLL)].second++;
			}
		}
	}

	for (int e_a = 0; e_a < E_static; e_a++){
		int deg_a = (int)E_adj[e_a].size();
		E_inter[e_a].rehash(deg_a);
		for (int i = 0; i < deg_a; i++){
			int e_b = E_adj[e_a][i].first;
			int C_ab = E_adj[e_a][i].second;
			E_inter[e_a].insert({e_b, C_ab});
		}
	}

	cout << "Static hypergraph (2): adjacency done:\t\t"
	     << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;

	// Induced static hypergraph: counting motifs
	start = clock();
	vector<int> intersection(E_static, 0);
	vector<long long> motif2count(96, 0);
	long long valid_cnt = 0, total_cnt = 0;
	
	int bin = (int)(E_static / 20);

	for (int e_a = 0; e_a < E_static; e_a++){
		if (e_a % bin == 0) cout << e_a << " / " << E_static << endl;

		long long l_e_a = (long long)e_a;
		int size_a = (int)E2V[e_a].size();
		int deg_a = (int)E_adj[e_a].size();

		for (int i = 0; i < deg_a; i++){
			int e_b = E_adj[e_a][i].first, C_ab = E_adj[e_a][i].second;
			int size_b = (int)E2V[e_b].size();
			int deg_b = (int)E_adj[e_b].size();

			const auto &nodes = E2V_set[e_b]; auto it_end = nodes.end(); int cnt = 0;
			for (const int &node: E2V[e_a]){
				if (nodes.find(node) != it_end) intersection[cnt++] = node;
			}

			for (int j = i+1; j < deg_a; j++){
				int e_c = E_adj[e_a][j].first, C_ca = E_adj[e_a][j].second;
				int size_c = (int)E2V[e_c].size();
				int deg_c = (int)E_adj[e_c].size();

				int C_bc = E_inter[e_b][e_c];
				if (C_bc){
					if (e_a < e_b){
						int g_abc = 0;
						const auto &nodes = E2V_set[e_c]; auto it_end = nodes.end();
						for (int k = 0; k < C_ab; k++){
							if (nodes.find(intersection[k]) != it_end) g_abc++;
						}
						hmotif H = hmotif(e_a, e_b, e_c, size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);
						bool valid = count_temporal_hmotif(H, motif2count, E2times, delta);
						if (valid) valid_cnt ++;
						total_cnt ++;
					}
				} else{
					int g_abc = 0;
					hmotif H = hmotif(e_a, e_b, e_c, size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc);
					bool valid = count_temporal_hmotif(H, motif2count, E2times, delta);
					if (valid) valid_cnt ++;
					total_cnt ++;
				}
			}
		}
	}

	cout << "Static hypergraph (3): counting done\t\t"
	     << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "----------------------------------------------------------" << endl << endl;

	cout << "Static h-motif: " << valid_cnt << " / " << total_cnt << endl;

	int comb2[6][3] = {{0,0,1}, {0,1,0}, {1,0,0}, {1,1,0}, {1,0,1}, {0,1,1}};
	int comb2_C[6][3] = {{0,1,1}, {1,1,0}, {1,0,1}, {2,1,1}, {1,1,2}, {1,2,1}};
	int comb3[6][3] = {{0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}};
	int comb3_C[6][3] = {{0,1,2}, {2,1,0}, {0,2,1}, {1,2,0}, {2,0,1}, {1,0,2}};
	
	for (int e_a = 0; e_a < E_static; e_a++){
		vector<pair<double, int> > _edges;
		for (const double &timestamp: E2times[e_a]) _edges.push_back({timestamp, 1});

		int counter[2][2][2] = {0};
		sort(_edges.begin(), _edges.end());

		int start = 0;
		for (int end = 0; end < _edges.size(); end++){
			while (_edges[start].first + delta < _edges[end].first){
				// decrement
				counter[1][0][0] --;
				counter[1][1][0] -= counter[1][0][0];
				start ++;
			}
			// increment
			counter[1][1][1] += counter[1][1][0];
			counter[1][1][0] += counter[1][0][0];
			counter[1][0][0] ++;
		}
		motif2count[95] += counter[1][1][1];

		for(const pair<int,int> &p_b: E_adj[e_a]){
			int e_b = p_b.first;
			if (e_a > e_b) continue;
			
			int S[2] = {(int)E2V[e_a].size(), (int)E2V[e_b].size()};
			int C[3] = {(int)E2V[e_a].size(), p_b.second, (int)E2V[e_b].size()};
			int G = p_b.second;

			int counter[3][3][3] = {0};

			vector<pair<double, int> > edges = _edges;
			for (const double &timestamp: E2times[e_b]) edges.push_back({timestamp, 2});
			sort(edges.begin(), edges.end());

			int start = 0;
			for (int end = 0; end < edges.size(); end++){
				while (edges[start].first + delta < edges[end].first){
					// decrement
					counter[edges[start].second][0][0] --;
					for (int i = 1; i <= 2; i++){
						counter[edges[start].second][i][0] -= counter[i][0][0];
					}
					start ++;
				}
				// increment
				for (int i = 1; i <= 2; i++){
					for (int j = 1; j <= 2; j++){
						counter[i][j][edges[end].second] += counter[i][j][0];
					}
					counter[i][edges[end].second][0] += counter[i][0][0];
				}
				counter[edges[end].second][0][0] ++;
			}
	
			for (int i = 0; i < 6; i++){
				int size_a = S[comb2[i][0]], size_b = S[comb2[i][1]], size_c = S[comb2[i][2]];
				int C_ab = C[comb2_C[i][0]], C_bc = C[comb2_C[i][1]], C_ca = C[comb2_C[i][2]];
				int motif_index = get_motif_index(size_a, size_b, size_c, C_ab, C_bc, C_ca, G, 1, 2, 3);
				motif2count[motif_index] += counter[comb2[i][0]+1][comb2[i][1]+1][comb2[i][2]+1];
			}
		}
	}

	double counting_runtime = (double)(clock() - start) / CLOCKS_PER_SEC;

	// File output
	start = clock();
	long long total_count = 0;
	ofstream resultFile(writeFile.c_str());

	resultFile << fixed << counting_runtime << endl;

	for (int i = 0; i < 96; i++){
		resultFile << i + 1 << "\t" << motif2count[i] << endl;
		total_count += motif2count[i];
	}
	resultFile.close();

	cout << "Counting temporal h-motifs done: "
	     << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "----------------------------------------------------------" << endl << endl;
	
	return 0;
}
