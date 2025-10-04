#pragma once
#include <iostream> 
#include <random>
#include <cmath>
#include <vector>
#include <queue>
#include <algorithm>
#include <set>
#include <iomanip>
#include <ctime>
#include <stack>
#include <fstream>
#include "graphtypes.h"
using namespace std;

void GenerateVertexDegreesPareto(vector<int>& degrees, int ver, double alpha, int xmin);
vector<vector<int>> Createadmatrix(int ver, const vector<int>& degrees);
vector<vector<int>> Createdostmatrix(vector<vector<int>> admatrix, vector<vector<int>> matrix, const int ver);
vector<vector<int>> Addmatrix(vector<vector<int>> matrix1, vector<vector<int>> matrix2, const int ver);
vector<vector<int>> Shimbel(vector<vector<int>> weightmatrix, vector<vector<int>> matrix, const int ver, const int parametr);
vector<int> PrintPath(vector<int> from, int finish);
vector<int> BFS(vector<vector<int>> weightmatrix, int ver, int start, int finish, vector<int>& path, int& iter_BFS);
vector<int> BellmanFord(vector<vector<int>>& weightmatrix, int ver, int start, int finish, vector<int>& path, int& iter_Bell);
void PrintWay(int ver, vector<int> dist, int from, int where_, vector<int> path);
void PrintMatrix(vector<vector<int>> matrix);
void menu();
void erase(vector<int>& vec, int index);
int generateParetoValue(double alpha, int xmin, mt19937& gen);
vector<vector<int>> GenerateWeightMatrix(int ver, const vector<vector<int>>& admatrix, double alpha, int xmin);
vector<vector<int>> GenerateCostMatrix(int ver, const vector<vector<int>>& admatrix, double alpha, int xmin);
vector<vector<int>> GenerateCapacityMatrix(int ver, const vector<vector<int>>& admatrix, double alpha, int xmin);
int MinPropusk(vector<vector<int>> bandwidth_matrix, int start, int finish, int ver, vector<int>& parent);
int Ford_Fulkerson(vector<vector<int>> bandwidth_matrix, int start, int finish, int ver);
int MinCost(vector<vector<int>>& cost_matrix, vector<vector<int>>& bandwidth_matrix,int s, int t, int ver, int need_flow);
vector<edge> Kruskal(vector<vector<int>> weightmatrix, int ver, int& mst_weight, int& iter_kruskal);
int determinant(int ver, vector<vector<int>> kirhgof_matrix);
vector<vector<int>> decodePrufer(vector<pair<int, int>> prufer);
vector<pair<int, int>> codePrufer(vector<edge> mst, int ver);
vector<pair<int, int>> MaxIndependentEdgeSet(vector<vector<int>>& matrix, int ver);
void printAllEdges(const vector<vector<int>>& matrix, int ver);
vector<vector<int>> Euler(int ver, vector<vector<int>> undirected_weight_matrix, vector<int> undirected_degrees);
void FindEuler(int ver, vector<vector<int>> eulergraph);
bool CanAddV(int v, vector<vector<int>>& gamiltongraph, vector<int>& path, int p);
bool FindGamilton(vector<vector<int>>& gamiltongraph, vector<int>& path, int p);
bool ExistGamiltonCycle(vector<vector<int>>& gamiltongraph);
vector<vector<int>> Gamilton(int ver, vector<vector<int>> undirected_weight_matrix);
void FindGamiltonKomivvv(vector<vector<int>>& gamiltongraph, vector<bool>& visited, vector<int>& path, int vert, int ver, int count, int cost, int& min_cost, vector<int>& min_path);
void Komivoyadger(vector<vector<int>>& gamiltongraph, int ver);
void PrintAllHamiltonCyclesWithWeights(vector<vector<int>>& gamiltongraph, vector<int>& path, int position, int& cycle_count, vector<vector<int>>& weight_matrix);