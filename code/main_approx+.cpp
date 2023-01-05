#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <numeric>
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
		int index_a, index_b, index_c;
		int size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc;
		hmotif(int index_a, int index_b, int index_c,
		       int size_a, int size_b, int size_c, int C_ab, int C_bc, int C_ca, int g_abc){
			this->index_a = index_a;
			this->index_b = index_b;
			this->index_c = index_c;
			this->size_a = size_a;
			this->size_b = size_b;
			this->size_c = size_c;
			this->C_ab = C_ab;
			this->C_bc = C_bc;
			this->C_ca = C_ca;
			this->g_abc = g_abc;
		}
};

inline long long convert_id(int x, int y)
{
	return x * (1LL << 31) + y;
}

string hyperedge2hash(vector<int>& nodes)
{
	// Generate unique hash value for each hyperedge.
	stringstream ss;
	sort(nodes.begin(), nodes.end());
	for (int i = 0; i < nodes.size(); i++){
		ss << nodes[i];
		if (i < nodes.size() - 1) ss << ",";
	}
	return ss.str();
}

void update_add(int e, double time_e,
		vector< vector<int> >& V2E, vector< vector<int> >& E2V, vector< vector<double> >& E2times,
		vector< vector< pair<int, int> > >& E_adj, vector< unordered_map<int, int> >& E_inter)
{
	// Add the hyperedge to the projected graph.

	// 1. Update V2E and E2times. 
	for (const int &v: E2V[e]){
		if (find(V2E[v].begin(), V2E[v].end(), e) == V2E[v].end()){
			V2E[v].push_back(e);
		}
	}
	E2times[e].push_back(time_e);
	
	// 2. Update E_adj and E_inter.
	vector<long long> upd_time(E2V.size(), -1LL);
	if (E2times[e].size() == 1){
		for (const int &v: E2V[e]){
			for (const int &_e: V2E[v]){
				if (e == _e) continue;
				if ((upd_time[_e] >> 31) ^ e){
					upd_time[_e] = ((long long)e << 31) + (int)E_adj[e].size();
					E_adj[e].push_back({_e, 0});
				}
				E_adj[e][(int)(upd_time[_e] & 0x7FFFFFFFLL)].second++;
			}
		}
		E_inter[e].rehash(E_adj[e].size());
		for (const pair<int, int>& p: E_adj[e]){
			int _e = p.first, C = p.second;
			E_adj[_e].push_back({e, C});
			E_inter[e].insert({_e, C});
			E_inter[_e].insert({e, C});
		}
	}
}

void update_delete(int e, int time_e, 
		   vector< vector<int> >& V2E, vector< vector<int> >& E2V, vector< vector<double> >& E2times,
		   vector< vector< pair<int, int> > >& E_adj, vector< unordered_map<int, int> >& E_inter)
{
	// Delete the hyperedge from the projected graph.

	E2times[e].erase(E2times[e].begin());
	if (E2times[e].size() == 0){
		for (const int &v: E2V[e]){
			V2E[v].erase(remove(V2E[v].begin(), V2E[v].end(), e), V2E[v].end());
		}
		for (const pair<int, int>& p: E_adj[e]){
			int _e = p.first;
			E_inter[_e].erase(e);
			for (int i = 0; i < E_adj[_e].size(); i++){
				if(E_adj[_e][i].first == e){
					E_adj[_e].erase(E_adj[_e].begin() + i);
					break;
				}
			}
		}
		E_adj[e].clear();
		E_inter[e].clear();
	}
}

int main(int argc, char *argv[])
{
	clock_t run_start;

	// Configuration
	string dataset = argv[1];
	double delta = stod(argv[2]);
	double _T = stod(argv[3]);
	int S = stoi(argv[4]);
	int SEED = 2022;

	// File path
	string graphFile = "../data/" + dataset + "/" + dataset;
	string writeFile = "../results/" + dataset + "_" + to_string(delta) + "_" + to_string(S) + "_" + to_string(_T);
	graphFile += ".txt";
	writeFile += "_approx+.txt";

	// Read data
	run_start = clock();
	vector< vector<int> > node2hyperedge;
	vector< vector<int> > hyperedge2node;
	vector< unordered_set<int> > hyperedge2node_set;
	vector<double> hyperedge2time;
	read_data(graphFile, node2hyperedge, hyperedge2node, hyperedge2node_set, hyperedge2time);
	
	int V = (int)node2hyperedge.size();
	int E = (int)hyperedge2node.size();
	cout << "# of nodes:\t\t\t" << V << endl;
	cout << "# of hyperedges:\t" << E << endl;
	cout << "Reading data done:\t"
	     << (double)(clock() - run_start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "-----------------------------------------------------------" << endl << endl;

	// Induced static hypergraph preprocessing
	run_start = clock();
	unordered_map<string, int> hash2index;
	vector<int> hyperedge2index;
	vector< vector<int> > E2V;
	vector< unordered_set<int> > E2V_set;

	for (int i = 0; i < E; i++){
		string hash = hyperedge2hash(hyperedge2node[i]);
		if (hash2index.find(hash) == hash2index.end()){
			hash2index[hash] = (int)hash2index.size();
			vector<int> nodes;
			unordered_set<int> nodes_set;
			for (const int &v: hyperedge2node[i]){
				nodes.push_back(v);
				nodes_set.insert(v);
			}
			E2V.push_back(nodes);
			E2V_set.push_back(nodes_set);
		}
		hyperedge2index.push_back(hash2index[hash]);
	}

	int E_static = (int)E2V.size();
	cout << "# of induced hyperedges:\t" << E_static << endl;

	cout << "Induced static hypergraph preprocessing done:\t"
	     << (double)(clock() - run_start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "-----------------------------------------------------------" << endl << endl;

	// Sorting temporal hyperedges
	run_start = clock();
	int x = 0;
	vector<int> hyperedge2order(E);
	
	iota(hyperedge2order.begin(), hyperedge2order.end(), x++);
	sort(hyperedge2order.begin(), hyperedge2order.end(), [&](long long i, long long j){
			if (hyperedge2time[i] == hyperedge2time[j]) return i < j;
			return hyperedge2time[i] < hyperedge2time[j];
			});

	double T_min = hyperedge2time[hyperedge2order[0]];
	double T_max = hyperedge2time[hyperedge2order[E-1]];
	double T = max(1.0, (T_max - T_min) * _T);

	// Valid time zone
	vector<pair<double, double> > valid_time;
	vector<long long> valid_acc, valid_acc_shifted = {0};
	vector<double> unique_timestamps;
	double t_a = hyperedge2time[hyperedge2order[0]] - T + 1, t_b = hyperedge2time[hyperedge2order[0]];
       	long long acc = 0;

	for (int i = 0; i < E; i++){
		int e_i = hyperedge2order[i];
		double time_i = hyperedge2time[e_i];
		unique_timestamps.push_back(time_i);

		if (time_i - T > t_b){
			valid_time.push_back({t_a, t_b});
			acc += (long long)(t_b - t_a + 1);
			valid_acc.push_back(acc);
			valid_acc_shifted.push_back(acc);
			t_a = time_i - T + 1;
		}
		t_b = time_i;
	}
	valid_time.push_back({t_a, t_b});
	acc += (long long)(t_b - t_a + 1);
	valid_acc.push_back(acc);
	valid_acc_shifted.push_back(acc);

	double valid_tot = (double)valid_acc[valid_acc.size()-1];
	cout << "Valid time range:\t" << valid_tot / (T_max - T_min + T) << endl;

	cout << "Sorting temporal hyperedges done:\t"
	     << (double)(clock() - run_start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "-----------------------------------------------------------" << endl << endl;

	// Counting temporal hypergraph motifs
	run_start = clock();
	vector<hmotif> instances;
	vector<long long> motif2count(96, 0);
	vector<int> intersection(E, 0);

	mt19937 gen(SEED);
	uniform_real_distribution<> dist(1, valid_acc[valid_acc.size()-1]+1);

	int bin = max(1, (int)(S / 20));

	for (int num_sample = 0; num_sample < S; num_sample++){
		if (num_sample % bin == 0){
			cout << num_sample << " / " << S << endl;
		}

		// Sample an interval
		long long s = (long long)dist(gen);

		int idx = lower_bound(valid_acc.begin(), valid_acc.end(), s) - valid_acc.begin();
		double T_start = valid_time[idx].first + (double)(s - valid_acc_shifted[idx] - 1);
		double T_end = T_start + T;
		int start_index = upper_bound(unique_timestamps.begin(), unique_timestamps.end(), T_start-1-delta) - unique_timestamps.begin();
	
		// Initialization
		int start = start_index;
		vector<long long> upd_time(E2V.size(), -1LL);
		vector< vector<int> > V2E; V2E.resize(V);
		vector< vector<double> > E2times; E2times.resize(E);
		vector< vector< pair<int, int> > > E_adj; E_adj.resize(E_static);
		vector< unordered_map<int, int> > E_inter; E_inter.resize(E_static);

		for (int end = start_index; end < E; end++){
			int e_a = hyperedge2order[end];
			int index_a = hyperedge2index[e_a];
			int size_a = (int)E2V[index_a].size();
			double time_a = hyperedge2time[e_a];

			if (time_a >= T_end) break;

			// Update the projected graph
			while (hyperedge2time[hyperedge2order[start]] + delta < time_a){
				int e_start = hyperedge2order[start];
				int index_start = hyperedge2index[e_start];
				double time_start = hyperedge2time[e_start];
				update_delete(index_start, time_start, V2E, E2V, E2times, E_adj, E_inter);
				start ++;
			}

			update_add(index_a, time_a, V2E, E2V, E2times, E_adj, E_inter);
		
			if ((T_start <= time_a) and (time_a < T_end)){
				// Count induced instances
				instances.clear();
				int deg_a = (int)E_adj[index_a].size();
				for (int i = 0; i < deg_a; i++){
					int index_b = E_adj[index_a][i].first;
					int C_ab = E_adj[index_a][i].second;
					int size_b = (int)E2V[index_b].size();
					int deg_b = (int)E_adj[index_b].size();

					long long ab_id = convert_id(index_a, index_b);
					const auto &nodes = E2V_set[index_b];
					auto it_end = nodes.end();
					int cnt = 0;
					for (const int &v: E2V[index_a]){
						if (nodes.find(v) != it_end) intersection[cnt++] = v;
					}

					for (int j = 0; j < deg_b; j++){
						int index_c = E_adj[index_b][j].first;
						int C_bc = E_adj[index_b][j].second;
						int size_c = (int)E2V[index_c].size();

						if (index_c == index_a or index_c == index_b) continue;
						upd_time[index_c] = ab_id;

						int C_ca = E_inter[index_a][index_c], g_abc = 0;
						const auto &nodes = E2V_set[index_c];
						auto it_end = nodes.end();
						for (int k = 0; k < C_ab; k++){
							if (nodes.find(intersection[k]) != it_end) g_abc++;
						}

						if (C_ca){
							if (index_b < index_c){
								instances.push_back(hmotif(index_a, index_b, index_c, size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc));
							}
						} else{
							instances.push_back(hmotif(index_a, index_b, index_c, size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc));
						}
					}

					for (int j = i+1; j < deg_a; j++){
						int index_c = E_adj[index_a][j].first;
						int C_ca = E_adj[index_a][j].second;
						int size_c = (int)E2V[index_c].size();

						if (upd_time[index_c] == ab_id) continue;
						if (index_c == index_a or index_c == index_b) continue;

						int C_bc = 0, g_abc = 0;
						instances.push_back(hmotif(index_a, index_b, index_c, size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc));
					}
				}
				
				// Count temporal h-motif instances
				for (const hmotif instance: instances){
					int index_a = instance.index_a;
					int index_b = instance.index_b;
					int index_c = instance.index_c;
					int size_a = instance.size_a;
					int size_b = instance.size_b;
					int size_c = instance.size_c;
					int C_ab = instance.C_ab;
					int C_bc = instance.C_bc;
					int C_ca = instance.C_ca;
					int g_abc = instance.g_abc;

					int motif_index_bc = get_motif_index(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc, 3, 1, 2);
					int motif_index_cb = get_motif_index(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc, 3, 2, 1);

					int bc_cnt = 0, cb_cnt = 0;
					int b_temp_cnt = (int)E2times[index_b].size(), c_temp_cnt = (int)E2times[index_c].size();
					for (int i = 0; i < b_temp_cnt; i++){
						for (int j = 0; j < c_temp_cnt; j++){
							if (E2times[index_b][i] < E2times[index_c][j]){
								bc_cnt += c_temp_cnt - j;
								break;
							}
						}
					}
					cb_cnt = b_temp_cnt * c_temp_cnt - bc_cnt;
					
					motif2count[motif_index_bc] += (double)bc_cnt;
					motif2count[motif_index_cb] += (double)cb_cnt;
				}

				// Temporal h-motifs 87 ~ 95
				for (int i = 0; i < deg_a; i++){
					int index_b = E_adj[index_a][i].first;
					int C_ab = E_adj[index_a][i].second;
					int size_b = (int)E2V[index_b].size();

					int motif_index_ab = get_motif_index(size_a, size_b, size_a, C_ab, C_ab, size_a, C_ab, 1, 2, 3);
					int motif_index_ba = get_motif_index(size_b, size_a, size_a, C_ab, size_a, C_ab, C_ab, 1, 2, 3);
					int motif_index_bb = get_motif_index(size_b, size_b, size_a, size_b, C_ab, C_ab, C_ab, 1, 2, 3);

					int ab_cnt = 0, ba_cnt = 0, bb_cnt;
					int a_temp_cnt = (int)E2times[index_a].size(), b_temp_cnt = (int)E2times[index_b].size();
					for (int j = 0; j < a_temp_cnt - 1; j++){
						for (int k = 0; k < b_temp_cnt; k++){
							if (E2times[index_a][j] < E2times[index_b][k]){
								ab_cnt += b_temp_cnt - k;
								break;
							}
						}
					}
					ba_cnt = (a_temp_cnt - 1) * b_temp_cnt - ab_cnt;
					bb_cnt = b_temp_cnt * (b_temp_cnt - 1) / 2;

					motif2count[motif_index_ab] += (double)ab_cnt;
					motif2count[motif_index_ba] += (double)ba_cnt;
					motif2count[motif_index_bb] += (double)bb_cnt;
				}

				// Temporal h-motif 96
				int a_temp_cnt = (int)E2times[index_a].size();
				if (a_temp_cnt >= 2) motif2count[95] += (double)((a_temp_cnt - 1) * (a_temp_cnt - 2) / 2);
			}
		}
	}

	vector<double> final_motif2count(96, 0);
	for (int i = 0; i < 96; i++){
		final_motif2count[i] = motif2count[i] * (valid_tot / (S * T));
	}

	double counting_runtime = (double)(clock() - run_start) / CLOCKS_PER_SEC;
	cout << "Counting temporal h-motifs done:\t"
	     << (double)(clock() - run_start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "-----------------------------------------------------------" << endl << endl;

	// File output
	run_start = clock();
	double total_count = 0;
	ofstream resultFile(writeFile.c_str());

	resultFile << fixed << counting_runtime << endl;

	for (int i = 0; i < 96; i++){
		resultFile << i + 1 << "\t" << final_motif2count[i] << endl;
		total_count += final_motif2count[i];
	}	
	resultFile << total_count << endl;
	resultFile.close();
	
	cout << "File output done:\t"
	     << (double)(clock() - run_start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "-----------------------------------------------------------" << endl << endl;

	return 0;
}
