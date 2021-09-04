#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <algorithm>
#include <random>

using namespace std;

vector<string> split(string str, char delimiter)
{
	vector<string> internal;
	stringstream ss(str);
	string temp;
	while (getline(ss, temp, delimiter)){
		internal.push_back(temp);
	}
	return internal;
}

void read_data(string path,
	       vector< vector<int> >& node2hyperedge,
	       vector< vector<int> >& hyperedge2node,
	       vector< unordered_set<int> >& hyperedge2node_set,
	       vector<double>& hyperedge2time)
{
	ifstream graphFile(path.c_str());

	string line;
	vector<string> lines;
	while (getline(graphFile, line)){
		lines.push_back(line);
	}
	default_random_engine e(random_device{}());
	shuffle(lines.begin(), lines.end(), e);
	
	for (int i = 0; i < lines.size(); i++){
		string line = lines[i];
		vector<string> term = split(line, '\t');
		vector<string> nodes = split(term[0], ',');
		double timestamp = stod(term[1]);
		
		vector<int> tokens;
		unordered_set<int> tokens_set;
		for (int j = 0; j < nodes.size(); j++){
			tokens.push_back(stoi(nodes[j]));
			tokens_set.insert(stoi(nodes[j]));
			while (node2hyperedge.size() <= stoi(nodes[j])){
				node2hyperedge.push_back(vector<int>());
			}
			node2hyperedge[stoi(nodes[j])].push_back(i);
		}
		hyperedge2node.push_back(tokens);
		hyperedge2node_set.push_back(tokens_set);
		hyperedge2time.push_back(timestamp);
	}
}

void read_data_static(string path,
	       vector< vector<int> >& node2hyperedge,
	       vector< vector<int> >& hyperedge2node,
	       vector< unordered_set<int> >& hyperedge2node_set)
{
	ifstream graphFile(path.c_str());

	string line;
	vector<string> lines;
	while (getline(graphFile, line)){
		lines.push_back(line);
	}
	default_random_engine e(random_device{}());
	shuffle(lines.begin(), lines.end(), e);

	for (int i = 0; i < lines.size(); i++){
		string line = lines[i];
		vector<string> nodes = split(line, ',');
		
		vector<int> tokens;
		unordered_set<int> tokens_set;
		for (int j = 0; j < nodes.size(); j++){
			tokens.push_back(stoi(nodes[j]));
			tokens_set.insert(stoi(nodes[j]));
			while (node2hyperedge.size() <= stoi(nodes[j])){
				node2hyperedge.push_back(vector<int>());
			}
			node2hyperedge[stoi(nodes[j])].push_back(i);
		}
		hyperedge2node.push_back(tokens);
		hyperedge2node_set.push_back(tokens_set);
	}
}
