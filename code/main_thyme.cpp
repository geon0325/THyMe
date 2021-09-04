#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>

#include "read_data.cpp"
#include "motif_id.cpp"

using namespace std;

inline long long convert_id(int x, int y)
{
	return x * (1LL << 31) + y;
}

void project_add(int e, double delta,
		 vector< vector<int> >& V2E, vector< vector<int> >& E2V, vector<double>& E2time,
		 vector< vector< pair<int, int> > >& E_adj, vector< unordered_map<int, int> >& E_inter)
{
	// Add a hyperedge to the projected graph
	vector<long long> upd_time(E2V.size(), -1LL);
	for (const int &v: E2V[e]){
		for (const int &_e: V2E[v]){
			if (e == _e) continue;
			double time_interval = E2time[e] - E2time[_e];
			if (time_interval < 0) break;
			if (time_interval > delta) continue;
			if (time_interval == 0 and e < _e) continue;
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

void project_delete(int e, vector< vector< pair<int, int> > >& E_adj, vector< unordered_map<int, int> >& E_inter)
{
	// Delete a hyperedge from the projected graph
	E_adj[e].clear();
	E_inter[e].clear();
}

int main(int argc, char *argv[])
{
	clock_t run_start;

	// Configuration
	string dataset = argv[1];
	double delta = stod(argv[2]);

	// File path
	string graphFile = "../data/" + dataset + "/" + dataset + ".txt";
	string writeFile = "../results/" + dataset + "_" + to_string(delta) + "_thyme.txt";

	// Read data
	run_start = clock();
	vector< vector<int> > V2E;
	vector< vector<int> > E2V;
	vector< unordered_set<int> > E2V_set;
	vector<double> E2time;
	read_data(graphFile, V2E, E2V, E2V_set, E2time);

	int V = (int)V2E.size();
	int E = (int)E2V.size();
	cout << "# of nodes:\t\t" << V << endl;
	cout << "# of hyperedges:\t" << E << endl;
	cout << "Reading data done: "
	     << (double)(clock() - run_start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "----------------------------------------------------------" << endl << endl;

	// Sort by timestamps
	run_start = clock();
	int x = 0;
	vector<int> E2index(E);
	iota(E2index.begin(), E2index.end(), x++);
	sort(E2index.begin(), E2index.end(), [&](long long i, long long j)
			{return E2time[i] + (double)i / V / 100 < E2time[j] + (double)j / V / 100;});

	for (int i = 0; i < V; i++){
		sort(V2E[i].begin(), V2E[i].end(), [&](long long x, long long y)
				{return E2time[x] + (double)x / V/ 100 < E2time[y] + (double)y / V / 100;});
	}

	cout << "Sorting done: "
	     << (double)(clock() - run_start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "----------------------------------------------------------" << endl << endl;
 
	// Projected graph initialization
	vector< vector< pair<int, int> > > E_adj;
	vector< unordered_map<int, int> > E_inter;
	E_adj.resize(E);
	E_inter.resize(E);

	// Temporal hypergraph motif counting
	run_start = clock();
	int start = 0;
	vector<long long> motif2count(96, 0);
	vector<long long> upd_time(E, -1LL);
	vector<int> intersection(E, 0);

	int bin = (int)(E / 20);

	for (int end = 0; end < E; end++){
		if (end % bin == 0) cout << end << " / " << E << endl;

		int e_a = E2index[end];
		int size_a = (int)E2V[e_a].size();
		double time_a = E2time[e_a];

		// Update the projected graph: add a hyperedge
		project_add(e_a, delta, V2E, E2V, E2time, E_adj, E_inter);

		// Update the projected graph: delete hyperedges
		while (E2time[start] + delta < time_a){
			project_delete(E2index[start++], E_adj, E_inter);
		}

		// Count temporal h-motifs
		int deg_a = (int)E_adj[e_a].size();
		for (int i = 0; i < deg_a; i++){
			int e_b = E_adj[e_a][i].first;
			int C_ab = E_adj[e_a][i].second;
			int deg_b = (int)E_adj[e_b].size();
			int size_b = (int)E2V[e_b].size();
			double time_b = E2time[e_b];

			long long ab_id = convert_id(e_a, e_b);

			const auto &nodes = E2V_set[e_b]; auto it_end = nodes.end(); int cnt = 0;
			for (const int &v: E2V[e_a]){
				if (nodes.find(v) != it_end) intersection[cnt++] = v;
			}

			for (int j = 0; j < deg_b; j++){
				int e_c = E_adj[e_b][j].first;
				int C_bc = E_adj[e_b][j].second;
				int size_c = (int)E2V[e_c].size();
				double time_c = E2time[e_c];

				if (e_c == e_a or e_c == e_b) continue;
				if (time_a - min(time_b, time_c) > delta) continue;
				upd_time[e_c] = ab_id;

				int C_ca = E_inter[e_a][e_c], g_abc = 0;
				const auto &nodes = E2V_set[e_c]; auto it_end = nodes.end(); 
				for (int k = 0; k < C_ab; k++){
					if (nodes.find(intersection[k]) != it_end) g_abc++;
				}

				int motif_index = get_motif_index(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc,
								  time_a, time_b, time_c);

				if (C_ca){
					if (e_b < e_c){
						motif2count[motif_index] ++;
					}
				} else {
					motif2count[motif_index] ++;
				}
			}

			for (int j = i+1; j < deg_a; j++){
				int e_c = E_adj[e_a][j].first;
				int C_ca = E_adj[e_a][j].second;
				int size_c = (int)E2V[e_c].size();
				double time_c = E2time[e_c];
				
				if (upd_time[e_c] == ab_id) continue;
				if (e_c == e_a or e_c == e_b) continue;
				if (time_a - min(time_b, time_c) > delta) continue;

				int C_bc = 0, g_abc = 0;
				int motif_index = get_motif_index(size_a, size_b, size_c, C_ab, C_bc, C_ca, g_abc,
								  time_a, time_b, time_c);
				motif2count[motif_index] ++;
			}
		}
	}

	double runtime = (double)(clock() - start) / CLOCKS_PER_SEC;
	cout << "Counting temporal h-motifs done: "
	     << runtime << " sec" << endl;
	cout << "----------------------------------------------------------" << endl << endl;

	// File output
	start = clock();
	ofstream resultFile(writeFile.c_str());
	long long total_count = 0;

	resultFile << fixed << runtime << endl;

	for (int i = 0; i < 96; i++){
		resultFile << i + 1 << "\t" << motif2count[i] << endl;
		total_count += motif2count[i];
	}
	resultFile.close();
	
	cout << "File output done: "
	     << (double)(clock() - start) / CLOCKS_PER_SEC << " sec" << endl;
	cout << "----------------------------------------------------------" << endl << endl;
 
	return 0;
}
