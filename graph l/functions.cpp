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
#include "functions.h"
using namespace std;

int trash1 = 1e9;
int trash2 = 2147483647; //���� ���?

void erase(vector<int>& vec, int index) {
    vec.erase(remove(vec.begin(), vec.end(), index), vec.end()); //��� ��������, ������ value, ������� � ������ �������
}

template<typename T>
const T& clamp(const T& value, const T& low, const T& high) {
    return (value < low) ? low : (high < value) ? high : value;
}

int generateParetoValue(double alpha, int xmin, mt19937& gen) {
    uniform_real_distribution<> dis(0.0, 1.0);
    double u = dis(gen);
    int value = static_cast<int>(xmin / pow(u, 1.0 / alpha));

    // ������������ �������� ���������� [1, 100]
    return clamp(value, 1, 100);
}

vector<vector<int>> GenerateWeightMatrix(int ver, const vector<vector<int>>& admatrix, double alpha, int xmin) {
    vector<vector<int>> weightmatrix(ver, vector<int>(ver, 0));
    static random_device rd;
    static mt19937 gen(rd());

    for (int i = 0; i < ver; i++) {
        for (int j = 0; j < ver; j++) {
            // ���������� �������� ������ ��� ������������ ����
            if (admatrix[i][j] == 1) {
                weightmatrix[i][j] = generateParetoValue(alpha, xmin, gen);
            }
        }
    }
    return weightmatrix;
}

vector<vector<int>> GenerateCostMatrix(int ver, const vector<vector<int>>& admatrix, double alpha, int xmin) {
    return GenerateWeightMatrix(ver, admatrix, alpha, xmin);
}

vector<vector<int>> GenerateCapacityMatrix(int ver, const vector<vector<int>>& admatrix, double alpha, int xmin) {
    return GenerateWeightMatrix(ver, admatrix, alpha, xmin);
}

vector<vector<int>> Createdostmatrix(vector<vector<int>> admatrix, vector<vector<int>> matrix, const int ver) {
    vector<vector<int>> new_matrix(ver, vector<int>(ver, 0));
    for (int i = 0; i < ver; i++) {
        for (int j = 0; j < ver; j++) {
            vector<int> vect;
            for (int k = 0; k < ver; k++) new_matrix[i][j] += (matrix[i][k] * admatrix[k][j]);
        }
    }
    return new_matrix;
}

vector<vector<int>> Addmatrix(vector<vector<int>> matrix1, vector<vector<int>> matrix2, const int ver) {
    vector<vector<int>> new_matrix(ver, vector<int>(ver, 0));
    for (int i = 0; i < ver; i++)
        for (int j = 0; j < ver; j++) new_matrix[i][j] = (matrix1[i][j] + matrix2[i][j]);
    return new_matrix;
}

void PrintMatrix(vector<vector<int>> matrix) {
    for (const auto& row : matrix) {
        for (int val : row) cout << setw(7) << val << setw(5);
        cout << "\n";
    }
    cout << "\n";
}

vector<int> PrintPath(vector<int> from, int finish) {
    vector<int> path;  // ������� ���� �� ��������� ������� �� finish
    // ��������������� ����, �������� "�����" �� ������� from
    // from[v] ������, �� ����� ������� �� ������ � ������� v ��� ������ �����
    for (int v = finish; v != -1; v = from[v]) {
        path.push_back(v);  // ��������� ������� ������� � ����
    }
    reverse(path.begin(), path.end());
    return path;
}

void PrintWay(int ver, vector<int> dist, int from, int where_, vector<int> path) {
    cout << "������ ����������: " << endl;
    for (int i = 0; i < ver; i++) {
        if (dist[i] != trash1 and dist[i] != trash2 and dist[i] != -trash1) cout << dist[i] << "\t";
        else cout << "-" << "\t";

    }
    if (dist[where_] != trash1 and dist[where_] != -trash1) {
        cout << "\n����� ��������: " << dist[where_] << ".";
        cout << "\n�������: ";
        for (int i = 0; i < path.size(); i++) {
            if (i == path.size() - 1) cout << path[i] << ".";
            else cout << path[i] << "->";
        }
    }
    else cout << "\n����� ������� ��������� ������" << endl;
}

//������������ ������� ���������� �� ������ ������
int determinant(int ver, vector<vector<int>> kirhgof_matrix) {
    int det = 0;
    vector<vector<int>> submatrix(ver, vector<int>(ver, 0));

    if (ver == 2) {
        return ((kirhgof_matrix[0][0] * kirhgof_matrix[1][1]) - (kirhgof_matrix[1][0] * kirhgof_matrix[0][1]));
    }
    else {
        for (int x = 0; x < ver; x++) {
            int subi = 0;
            for (int i = 1; i < ver; i++) {
                int subj = 0;
                for (int j = 0; j < ver; j++) {
                    if (j == x) {
                        continue;
                    }
                    submatrix[subi][subj] = kirhgof_matrix[i][j];
                    subj++;
                }
                subi++;
            }
            det = det + (pow(-1, x) * kirhgof_matrix[0][x] * determinant(ver - 1, submatrix));
        }
    }
    return det;
}


void GenerateVertexDegreesPareto(vector<int>& degrees, int ver, double alpha, int xmin) {
    vector<double> probability(ver); //������ ��� ������������, ��� ������������

    static random_device rd; // ��������� ���������
    static mt19937 gen(rd());

    double sum_prob = 0.0; //������������ ��������������� ����������� ��� ������������� ������
    for (int k = xmin; k <= ver; ++k) {
        double prob = pow(k, -(alpha + 1));  // ������ ��� ������� �������� �������
        probability[k - 1] = prob;
        sum_prob += prob; //��������� ��� ������������
    }

    for (int k = 0; k < ver; ++k) { //����������� �����������, ����� ����� ���� ����� 1
        probability[k] /= sum_prob;
    }

    // ������� ���������� ������������� � ����� ������������� (�����������)
    discrete_distribution<> distribution(probability.begin(), probability.end());

    bool flag = true; // ���������� ��������� �������
    while (flag) {
        for (int i = 0; i < ver; ++i) {
            degrees[i] = distribution(gen) + xmin; // ���������� �����, ��������������� ������������� ������
        }
        sort(degrees.begin(), degrees.end()); // ��������� �� �����������
        reverse(degrees.begin(), degrees.end()); // ��������� �� ��������

        // �������� �� ������������ ������������� �������� ������
        degrees[ver - 1] = 0;
        int pr_sum = 0;
        int sumVertexDegrees = 0;
        int max_deg = 0;
        for (int i = 0; i < ver; i++) {
            if (degrees[i] <= ver - i - 1) {
                pr_sum += 1;
                sumVertexDegrees += degrees[i];
            }
            // ��������� ���-�� ������ � ������������ ��������, ���� ����� ������� ���� => ����� ����
            else if (degrees[i] == ver - 1) max_deg += 1;
        }

        if (max_deg == 0) {
            degrees[0] = ver - 1;
            max_deg += 1;
        }

        if (pr_sum == ver && sumVertexDegrees <= (ver * ver - ver) / 2 && max_deg == 1) flag = false;
    }
    //����� ���� �������� ������
    for (int i = 0; i < ver; i++) cout << "������� ������� �" << i << " = " << degrees[i] << endl;
}

vector<vector<int>> Createadmatrix(int ver, const vector<int>& degrees) {
    //������� ������� ��������� - ������� �������� ���-�� ������*���-�� ������, ����������� ������
    vector<vector<int>> admatrix(ver, vector<int>(ver, 0));
    //������� ����� ������� �������� ������
    vector<int> copy_degrees(degrees);
    srand(time(0));
    for (int i = 0; i < ver; i++) {
        for (int d = 0; d < copy_degrees[i]; d++) {
            bool n_addded_edge = true;
            while (n_addded_edge) {
                double R = static_cast<double>(rand()) / RAND_MAX; //����� � ��������� [0;1)
                int j = static_cast<int>((ver)*R);
                //���������: ��� ����� ����� ��������� ��� �� �����������, ��� �� ������� �����, ��� ������� ������� ��� ������ ����
                if (admatrix[i][j] != 1 && i != j && copy_degrees[i] > 0 && j > i) {
                    admatrix[i][j] = 1;
                    n_addded_edge = false;
                }
            }
        }
    }
    return admatrix;
}

vector<vector<int>> Shimbel(vector<vector<int>> weightmatrix, vector<vector<int>> matrix, const int ver, const int parametr) {
    vector<vector<int>> new_matrix(ver, vector<int>(ver, 0));
    for (int i = 0; i < ver; i++) {
        for (int j = 0; j < ver; j++) {
            vector<int> vect;
            for (int k = 0; k < ver; k++) {
                if (matrix[i][k] == 0 or weightmatrix[k][j] == 0) vect.push_back(0);
                else vect.push_back(matrix[i][k] + weightmatrix[k][j]);
            }
            sort(vect.begin(), vect.end()); 
            if (vect.back() == vect.front() && vect.back() == 0) new_matrix[i][j] = 0;
            else {
                if (parametr == 1) new_matrix[i][j] = vect.back();
                else if (parametr == 2) {
                    erase(vect, 0);
                    new_matrix[i][j] = vect.front();
                }
            }
        }
    }
    return new_matrix;
}

//�������� ������ ������ � ������ BFS
vector<int> BFS(vector<vector<int>> weightmatrix, int ver, int start, int finish, vector<int>& path, int& iter_BFS) {
    vector<int> from(ver, -1); //����� ������������ ����: from[i] = �������, �� ������� ������ � i
    iter_BFS = 0;
    vector<int> distance(ver, trash1); //���������� �� ������, ���������� "�������������" (trash1), ����� ���������
    queue<int> q; //������� ��� ������ ������
    distance[start] = 0; //���������� �� ��������� ������� = 0
    q.push(start); //��������� ������� � �������

    while (!q.empty()) {
        int vert = q.front(); //��������� ������� ������� �� �������
        q.pop();
        for (int to = 0; to < ver; to++) { //��������� ��� ��������� �������, � ������� ����� ������� �� �������
            if (weightmatrix[vert][to] != 0) {  //���� ���� ����� �� vert � to 
                //���� ������ ����� �������� ���� �� ������� to � ��������� ���
                if (distance[to] > distance[vert] + weightmatrix[vert][to]) {
                    distance[to] = distance[vert] + weightmatrix[vert][to];  //��������� ���������� ����������
                    q.push(to); //��������� ������� � ������� ��� ����������� ������
                    from[to] = vert; //����������, ������ ������ � ������� to
                }
            }
            iter_BFS++;
        }
    }

    path = PrintPath(from, finish); //��������������� ���� �� start �� finish ����� ��������������� ������ from
    return distance; //���������� ������ ���������� �� start �� ���� ��������� ������
}

vector<int> BellmanFord(vector<vector<int>>& weightmatrix, int ver, int start, int finish, vector<int>& path, int& iter_Bell) {
    vector<int> dist(ver, trash1); // �������������� ���������� ��� "�������������"
    vector<int> from(ver, -1);  // ��� �������������� ����
    dist[start] = 0; // ���������� �� ��������� ������� 0
    iter_Bell = 0;
    // �������� ���� - ��������� ver-1 ��������
    for (int i = 0; i < ver - 1; ++i) 
    {
        bool updated = false;
        // �������� �� ���� ������ �����
        for (int u = 0; u < ver; ++u) {
            for (int v = 0; v < ver; ++v) {
                if (weightmatrix[u][v] != 0) { // ���� ����� ����������
                    if (dist[v] > dist[u] + weightmatrix[u][v]) {
                        dist[v] = dist[u] + weightmatrix[u][v];
                        from[v] = u;
                        updated = true;
                    }
                }
                iter_Bell++; // 1 �������� = 1 �������� �����
            }
        }
        // ���� �� ������� �������� �� ���� ���������, ������� ������
        if (!updated) break;
    }

    path = PrintPath(from, finish);
    return dist;
}

int MinPropusk(vector<vector<int>> bandwidth_matrix, int start, int finish, int ver, vector<int>& parent) {
    //parent ����� ������� ������������ ������� ��� ������ �������, ��������� ��� ���������� -1, ������� ��� �� ��������.
    fill(parent.begin(), parent.end(), -1);
    //��������� ������� �������� ��������� -2, ����� ����������, ��� ��� ��������� ����� ������.
    parent[start] = -2;
    //� ������� ����� ������� ����: {�������, ����������� ���������� ����������� �� ���}.
    queue<pair<int, int>> q;
    //��������� ��������� ������� � ������� � "�����������" ���������� ������������ (trash1).
    q.push({ start, trash1 });

    //�������� ����� � ������ (BFS)
    while (!q.empty()) {
        //��������� �� ������� ������� u � ������� ����������� ���������� ����������� cap.
        int u = q.front().first;
        int cap = q.front().second;
        q.pop();

        //��������� ���� ������� ������� u.
        for (int v = 0; v < ver; v++) {
            //���� ���������� ����� �� u � v (bandwidth_matrix[u][v] != 0), � ������� v ��� �� ���� �������� (parent[v] == -1)
            if (u != v && bandwidth_matrix[u][v] != 0 && parent[v] == -1) {
                // �������� ������� v ��� ����������, ���������� �� ������������ ������� (������ ������).
                parent[v] = u;
                int min_cap = min(cap, bandwidth_matrix[u][v]);
                //���� �� �������� ������� finish, ����� ���������� ����������� ���������� �����������.
                if (v == finish) return min_cap;
                //��������� ������� v � ������� � ����������� ���������� ������������.
                q.push({ v, min_cap });
            }
        }
    }
    return 0;
}

int Ford_Fulkerson(vector<vector<int>> bandwidth_matrix, int start, int finish, int ver) {
    // �������������� ������������ ������ (parent[i] = -1 ��������, ��� ������� i �� ���� ��������).
    vector<int> parent(ver, -1);
    int max_flow = 0; //��� ���������� ������������� ������
    int min_cap = 0;  //��� �������� ����������� ���������� ����������� ����

    // ���� ���������� ���� �� ��������� ������� � �������� ����� MinPropusk, ���������� ������ �����
    while (min_cap = MinPropusk(bandwidth_matrix, start, finish, ver, parent)) {
        max_flow += min_cap; //����������� ����� ����� �� �������� ����������� ���������� ����������� ��� ���������� ����

        //������������� ������� ���������� ������������:
        cout << "������������� ������� ���������� ������������:" << endl;
        for (int i = 0; i < ver; i++) {
            for (int j = 0; j < ver; j++) {
                cout << bandwidth_matrix[i][j] << "   ";
            }
            cout << endl;
        }
        cout << "-----------------" << endl;

        // ������� ������������� ���� � ���������� ������������
        cout << "������������� ���� � ���������� ������������: " << min_cap << endl;
        cout << "�������: ";

        int v_path = finish;  // ��� �������������� ���� ���������� ���������� v_path
        vector<int> path;
        while (v_path != start) {
            path.push_back(v_path);
            v_path = parent[v_path];
        }
        path.push_back(start);
        reverse(path.begin(), path.end());
        for (int p : path) {
            cout << p << "   ";
        }
        cout << endl;


        // ��������� ���������� ����������� �� ���������� ����
        int v = finish;
        while (v != start) {
            int u = parent[v];
            bandwidth_matrix[u][v] -= min_cap; //��������� ���������� ����������� �� ���������� ���� (�� ����������� ������)
            bandwidth_matrix[v][u] += min_cap; //����������� ���������� ����������� � �������� ����������� (�������� �����)
            v = u; //��������� � �������������� �������
        }
    }

    cout << "������������� ������� ���������� ������������:" << endl;
    for (int i = 0; i < ver; i++) {
        for (int j = 0; j < ver; j++) {
            cout << bandwidth_matrix[i][j] << "   ";
        }
        cout << endl;
    }
    cout << "-----------------" << endl;

    return max_flow; //���������� ������������ �����, ������� ������� �����
}



int MinCost(vector<vector<int>>& cost_matrix,
    vector<vector<int>>& bandwidth_matrix,
    int s, int t,
    int ver,
    int need_flow) {
    int total_cost = 0;  //��� �������� ����� ��������� ������
    int flow = 0; //��� �������� �������� ������
    vector<vector<int>> flow_matrix(ver, vector<int>(ver, 0));  //������� �������� ������ ����� ���������
    vector<vector<int>> residual_cap = bandwidth_matrix;  //������� ���������� ���������� ������������

    // ���� ������� ����� ������ �������, ���������� ������ ����
    while (flow < need_flow) {
        vector<vector<int>> residual_cost(ver, vector<int>(ver, 0));  //������� ���������� ���������� ����

        // ��������� residual_cost, ������� ���� � ������� � �������� �����������
        for (int i = 0; i < ver; ++i) {
            for (int j = 0; j < ver; ++j) {
                residual_cost[i][j] = cost_matrix[i][j];  // ��������� �� �������� �������
                if (flow_matrix[i][j] > 0) {  // ���� ���� ����� �� ����� �����
                    residual_cap[j][i] = flow_matrix[i][j];  // �������� �����
                    residual_cost[j][i] = -cost_matrix[i][j];  // �������� ��������� (�������������)
                }
            }
        }

        // ����� ���� ����������� ��������� (�������� ��������-�����)
        vector<int> dist(ver, trash1);  // ������ ���������� (���������)
        vector<int> parent(ver, -1);    // ������ ��������� ��� �������������� ����
        dist[s] = 0;  // ��������� �� ��������� ������� ����� 0

        // �������� ���� ��������-�����: ���� ���������� ����
        for (int i = 0; i < ver - 1; ++i) {
            for (int u = 0; u < ver; ++u) {
                for (int v = 0; v < ver; ++v) {
                    // ���� ���� �����, ������� �� ��������, � �� ����� ����� �������� ����
                    if (residual_cap[u][v] > 0 && dist[u] != trash1 &&
                        dist[v] > dist[u] + residual_cost[u][v]) {
                        dist[v] = dist[u] + residual_cost[u][v];  // ��������� ���������
                        parent[v] = u;  // ����������, ������ ������ � ������� v
                    }
                }
            }
        }

        // ���� �� ����� ���� �� �������� �������, ������, ������ ������ �� �����
        if (parent[t] == -1) break;

        // ���������� ����������� ���������� ����������� �� ����
        int delta = need_flow - flow;  // ����������� �� ���������� �����
        for (int v = t; v != s; v = parent[v]) {
            delta = min(delta, residual_cap[parent[v]][v]);  // ����������� �� ����������� ���������� ����������� �� ����
        }

        // ���������� ������ � ���������� ���������� ������������
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            flow_matrix[u][v] += delta;  // ����������� ����� �� ����� (u -> v)
            flow_matrix[v][u] -= delta;  // ��������� ����� � �������� ����������� (v -> u)
            total_cost += delta * cost_matrix[u][v];  // ����������� ����� ��������� ������
            residual_cap[u][v] -= delta;  // ��������� ���������� ����������� �� ����� (u -> v)
            residual_cap[v][u] += delta;  // ����������� ���������� ����������� � �������� ����������� (v -> u)
        }
        flow += delta;  // ����������� ����� �����
    }

    // ����� ������������� ������
    cout << "������������� ������:" << endl;
    for (int i = 0; i < ver; ++i) {
        for (int j = 0; j < ver; ++j) {
            if (flow_matrix[i][j] > 0) {  // ���� ���� ����� �� �����
                cout << i << "->" << j << ": " << flow_matrix[i][j]
                    << " (��������� " << cost_matrix[i][j] << ")" << endl;
            }
        }
    }

    // ���������� ����� ��������� ������������ ������
    return total_cost;
}



// ������� ��� ���������� ����� ���������, � �������� ����������� ������� k
int find_root(int k, vector<int> parent) {
    if (parent[k] == k) {  // ���� ������� �������� ������ ���������
        return k;
    }
    // ����� ���������� ������� ������ ���������
    return find_root(parent[k], parent);
}

// ������� ���������� ��� ���������, � ������� ������ ������� f_p � s_p
void union_(int f_p, int s_p, vector<int>& parent) {
    // ������� ����� ��������
    int x = find_root(f_p, parent);
    int y = find_root(s_p, parent);
    // ���������� ��� ���������, ������������ ������ ������ ��������� �� ������ �������
    parent[x] = y;
}

vector<edge> Kruskal(vector<vector<int>> weightmatrix, int ver, int& mst_weight, int& iter_kruskal) {
    vector<edge> edges;
    // ��������� ������ ����: ��� ������� ����� ��������� � ������ ���� <���, ������, ����>
    for (int i = 0; i < ver; i++) {
        for (int j = 0; j < ver; j++) {
            if (weightmatrix[i][j] != 0) {
                edges.push_back({ weightmatrix[i][j], i, j });
            }
        }
    }

    // ��������� ���� �� ���� (�� �����������)
    sort(edges.begin(), edges.end(), compare());

    cout << "\n��������������� ����: " << endl;
    for (const auto& e : edges) {
        cout << e.from << "-" << e.where_ << " ���: " << e.cost << endl;
    }

    vector<int> parent(ver);
    // �������������� ������������ ������, ��� ������ ������� - ��� ���� ������� (������ ������� ���������� �������� ����� ������)
    for (int i = 0; i < ver; i++) {
        parent[i] = i;
    }

    vector<edge> mst;
    // �������, ��� ����������� �������� ������ ������ ��������� ver-1 ����
    while (mst.size() != ver - 1) {
        // ���� ��������� ����� (� ����������� �����, ��� ��� �� ������������� ����)
        edge next_edge = edges.back();
        edges.pop_back();

        // ���� ����� ��������� ��� ������ ������� ����� � ������ �������
        int f_p = find_root(next_edge.from, parent);
        int s_p = find_root(next_edge.where_, parent);

        // ���� ����� ��������� �������� (�.�. ���������� ����� ����� �� ������� ����), �� ��������� ��� ����� � �������� ������
        if (f_p != s_p) {
            mst.push_back(next_edge);  // ��������� ����� � �������� ������
            mst_weight += next_edge.cost;  // ����������� ��������� ��� ��������� ������
            union_(f_p, s_p, parent);  // ���������� ���������
        }

        iter_kruskal++;  // ����������� ������� ���������� �������� (�� �����)
    }
    return mst;  // ���������� ����������� �������� ������
}


vector<pair<int, int>> codePrufer(vector<edge> mst, int ver) {
    // ������������� ������� ��������� ��� ������ (������������ ��������� ������)
    vector<vector<int>> mst_matrix(ver, vector<int>(ver, 0));

    // ��������� ������� ���������, ��� mst_matrix[i][j] - ��� ��� ����� ����� ��������� i � j
    for (int k = 0; k < mst.size(); k++) {
        int i = mst[k].from;
        int j = mst[k].where_;
        int w = mst[k].cost;
        mst_matrix[i][j] = w;
        mst_matrix[j][i] = w;
    }

    // ������������� ������� ��� �������� ���� �������
    vector<pair<int, int>> prufer;
    vector<int> degree_vert(ver);  // ������ ��� �������� �������� ������

    // ��������� ������� ���� ������
    for (int i = 0; i < ver; i++) {
        for (int j = 0; j < ver; j++) {
            if (mst_matrix[i][j] != 0) {
                degree_vert[i]++;
            }
        }
    }

    // ������� ���������� ����
    int r = 0;
    for (int i = 0; i < ver; i++) {
        r += degree_vert[i];
    }

    // ���� ���� ���� (r != 0), ���������� ��������� ����������� ������� � �������� 1
    while (r != 0) {
        for (int i = 0; i < ver; i++) {
            if (degree_vert[i] == 1) {  // ���� ������� ������� ����� 1, ��� �������� �� ��������
                for (int j = 0; j < ver; j++) {
                    if (mst_matrix[i][j] != 0) {
                        // ��������� ����� (��� � �������, � ������� ��� ���������)
                        prufer.push_back(make_pair(mst_matrix[i][j], j));
                        degree_vert[i]--;
                        degree_vert[j]--;
                        mst_matrix[i][j] = 0;
                        mst_matrix[j][i] = 0;
                        r -= 2;  // ��������� ���������� ���� �� 2 (�� ������� �������)
                    }
                }
            }
        }
    }

    return prufer;
}


vector<vector<int>> decodePrufer(vector<pair<int, int>> prufer) {
    // ���������� ������ � ������ (������ ���� ������� + 1)
    int ver = prufer.size() + 1;

    // ������������� ������� ��������� ��� ������������� ������
    vector<vector<int>> ans_matrix(ver, vector<int>(ver, 0));

    // ������ ��� ������������ �������� ������
    vector<int> vertexes(ver, 0);

    // ��������� ���������� ����, ��������� �� ������ ������� (�������)
    for (int i = 0; i < prufer.size(); i++) {
        vertexes[prufer[i].second] += 1;
    }

    int j = 0;
    int num = 0;

    // ��������������� ������ �� ���� �������
    for (int i = 0; i < ver - 1; i++) {
        for (j = 0; j < ver; j++) {
            if (i == ver - 2) {
                if (vertexes[j] == 0) {  // ������� ������� � �������� 0
                    vertexes[j] = -1;  // �������� ������� ��� ��������������
                    ans_matrix[j][prufer[i].second] = prufer[i].first;
                    ans_matrix[prufer[i].second][j] = prufer[i].first;
                    vertexes[prufer[i].second]--;  // ��������� ������� �������� �������
                    num = j;
                    break;
                }
            }
            else {
                if (vertexes[j] == 0 && num <= j) {
                    vertexes[j] = -1;
                    ans_matrix[j][prufer[i].second] = prufer[i].first;
                    ans_matrix[prufer[i].second][j] = prufer[i].first;
                    vertexes[prufer[i].second]--;
                    num = j;
                    break;
                }
            }
        }
    }
    return ans_matrix;  // ���������� ������� ��������� ���������������� ������
}

void printAllEdges(const vector<vector<int>>& matrix, int ver) {
    cout << "��� ���� �����:" << endl;
    for (int u = 0; u < ver; ++u) {
        for (int v = u + 1; v < ver; ++v) {
            if (matrix[u][v] != 0) {
                cout << u << " <-> " << v << " (���: " << matrix[u][v] << ")" << endl;
            }
        }
    }
}

vector<pair<int, int>> MaxIndependentEdgeSet(vector<vector<int>>& matrix, int ver) {
    vector<pair<int, int>> matching;  // ������ ��� �������� ��������� ����
    vector<bool> used(ver, false);  // ������ ��� ������������ �������������� ������

    // ���������� ��� �������
    for (int u = 0; u < ver; ++u) {
        if (!used[u]) {  // ���� ������� ��� �� ������������
            // ���� �����, ������� �� ���������� ������� u � ������� v
            for (int v = 0; v < ver; ++v) {
                if (matrix[u][v] != 0 && !used[v]) {  // ���� ���������� ����� � ������� v �� ������������
                    matching.push_back({ u, v });  // ��������� ����� � ���������
                    used[u] = true;  // �������� ������� u ��� ��������������
                    used[v] = true;  // �������� ������� v ��� ��������������
                    break;  // ��������� � ��������� ������� u
                }
            }
        }
    }

    return matching;  // ���������� ��������� ����
}

vector<vector<int>> Euler(int ver, vector<vector<int>> undirected_weight_matrix, vector<int> undirected_degrees) {
    vector<vector<int>> eulergraph(undirected_weight_matrix);
    int countevenvert = 0;
    bool flagEuler = false;
    cout << "\n������� �����: " << endl;
    PrintMatrix(eulergraph);
    //����� ���� �������� ������
    for (int i = 0; i < ver; i++) cout << "������� ������� �" << i << " = " << undirected_degrees[i] << endl;

    //���������, ��� ��� ������� ������ ����� ������ => ������� ����
    for (int i = 0; i < ver; i++) {
        if (undirected_degrees[i] % 2) countevenvert++;
    }
    if (countevenvert == ver) {
        cout << "���� �������� ���������\n" << endl;
        flagEuler = true;
    }
    //���� �� ��� ������� ������ ����� ������ => ���� �������������� ���� - ������� ���������
    else {
        cout << "\n���� �� �������� ���������. ������������ ����.\n" << endl;
        //����� ������� ��� ��������� ����� ���� ���� �� ������ ���������
        while (!flagEuler) {
            bool flag = 0;
            for (int i = 0; i < ver; i++) {
                for (int j = 0; j < ver; j++) {
                    //��������� ����� �����, ���� � ������� �������� ������� � ������ ����� ��� ���
                    if (undirected_degrees[i] % 2 != 0 and undirected_degrees[j] % 2 != 0 and eulergraph[i][j] == 0 and i != j) {
                        flag = true;
                        int r = rand() % (10 - (1) + 1) + 1;
                        eulergraph[i][j] = r;
                        eulergraph[j][i] = r;
                        undirected_degrees[i]++;
                        undirected_degrees[j]++;
                        cout << "�������� �����: " << i << "<->" << j << " � �����: " << r << endl;
                        break;
                    }
                }
                //���� �������� �����, �� ����� ��������� ������� ������ �� �������� � ����� �� ����� ���������� �������
                if (flag) break;
            }
            //���� �� ���������� �������� �����, �� ������� �������
            if (!flag) {
                for (int i = 0; i < ver; i++) {
                    for (int j = 0; j < ver; j++) {
                        if (undirected_degrees[i] % 2 != 0 and undirected_degrees[j] % 2 != 0 and eulergraph[i][j] != 0 and i != j) {
                            flag = true;
                            eulergraph[i][j] = 0;
                            eulergraph[j][i] = 0;
                            undirected_degrees[i]--;
                            undirected_degrees[j]--;
                            cout << "������� �����: " << i << "<->" << j << endl;
                            break;
                        }
                    }
                    //���� ������� �����, �� ����� ��������� ������� ������ �� �������� � ����� �� ����� �������� �������
                    if (flag) break;
                }
            }
            //��������� ���������� �� ������� ���� � ������� ��������� ������
            int countevenvert1 = 0;
            for (int i = 0; i < ver; i++) {
                if (undirected_degrees[i] % 2 == 0) countevenvert1++;
            }
            //���� ����������, �� ��������� ���� � ������� �� �����, ���� ��� - ���������� ���������/������� �����
            if (countevenvert1 == ver) flagEuler = true;
        }
        cout << "***********************" << endl;
        cout << "���������������� ������� ����." << endl;
        cout << "������� �����: " << endl;
        PrintMatrix(eulergraph);
        for (int i = 0; i < ver; i++) cout << "������� ������� �" << i << " = " << undirected_degrees[i] << endl;
        cout << "***********************" << endl;
    }
    return eulergraph;
}

void FindEuler(int ver, vector<vector<int>> eulergraph) {
    stack<int> vertexes;
    vector<int> cycle;
    int vert = 0;
    //���� ������ �������, �� ������� ���� �����
    for (vert = 0; vert < ver; vert++) {
        for (int j = 0; j < ver; j++) {
            if (eulergraph[vert][j] != 0) {
                break;
            }
        }
        if (vert < ver) {
            break;
        }
    }
    //��������� ��� ������� � ����
    vertexes.push(vert);
    //���� ���� �� ������
    while (!vertexes.empty()) {
        //����� ������� �������
        vert = vertexes.top();
        int i;
        //������� ������ ������� � ���
        for (i = 0; i < ver; i++) {
            if (eulergraph[vert][i] != 0) {
                break;
            }
        }
        //���� ������� �� �������, �� ��������� ��� ������� � ����, ������� �� �����
        if (i == ver) {
            cycle.push_back(vert);
            vertexes.pop();
        }
        //���� ����� �������, �� ��������� � ���� � ������� ����� ����� �������� ��������� � �������
        else {
            vertexes.push(i);
            eulergraph[vert][i] = 0;
            eulergraph[i][vert] = 0;
        }
    }
    cout << "\n������� ����:" << endl;
    for (int i = 0; i < cycle.size(); i++) {
        cout << cycle[i];
        if (i != cycle.size() - 1) {
            cout << " -> ";
        }
    }
}


//�������� ����������� ���������� ������� � ����
bool CanAddV(int v, vector<vector<int>>& gamiltongraph, vector<int>& path, int p) {
    if (gamiltongraph[path[p - 1]][v] == 0) return false; //����� ����� ��������� �������� path[p - 1] � �������� v � ������� ���������
    for (int i = 0; i < p; i++) {
        if (path[i] == v) return false;
    }
    return true;
}

bool FindGamilton(vector<vector<int>>& gamiltongraph, vector<int>& path, int p) {
    //���� ��� ������� �������� � ��������� ��������� � ������, �� ������������ true => ���� ����������� ���� � �����
    if (p == gamiltongraph.size()) return gamiltongraph[path[p - 1]][path[0]] != 0;

    //����� ���� ��������� ������� ��� ���������
    for (int v = 1; v < gamiltongraph.size(); v++) {
        //� ������� ���� ����������� ����� (���� ��������� ����� + �� ���� �������� ������� �������)
        if (CanAddV(v, gamiltongraph, path, p)) {
            //���� �������, �� ��������� ��� ������� � ����
            path[p] = v;
            //���� ����������� ���� ������, �� ���������� true
            if (FindGamilton(gamiltongraph, path, p + 1)) return true;
            //���� ����������� ���� �� ������, �� ������� ������� v �� ����
            path[p] = -1;
        }
    }
    return false;
}

// ������� ��� ������ ���� ������������� ������ � �� ��������� �����
void PrintAllHamiltonCyclesWithWeights(vector<vector<int>>& gamiltongraph, vector<int>& path, int position, int& cycle_count, vector<vector<int>>& weight_matrix) {
    int ver = gamiltongraph.size();
    if (position == ver) {
        // ���������, ���� �� ����� ����� ��������� � ������ ��������
        if (gamiltongraph[path[position - 1]][path[0]] != 0) {
            // ��������� ��������� ��� �����
            int total_weight = 0;
            for (int i = 0; i < ver - 1; i++) {
                total_weight += weight_matrix[path[i]][path[i + 1]];
            }
            total_weight += weight_matrix[path[ver - 1]][path[0]];  // �������� ����

            // ������� ���� � �����
            cout << "���� " << ++cycle_count << ": ";
            for (int i = 0; i < ver; i++) {
                cout << path[i] << (i < ver - 1 ? " -> " : "");
            }
            cout << " -> " << path[0] << " (��������� ���: " << total_weight << ")" << endl;
        }
        return;
    }

    // ���������� ��� ��������� ��������� �������
    for (int v = 0; v < ver; v++) {
        if (CanAddV(v, gamiltongraph, path, position)) {
            path[position] = v;
            PrintAllHamiltonCyclesWithWeights(gamiltongraph, path, position + 1, cycle_count, weight_matrix);
            path[position] = -1;  // ���������� ���������
        }
    }
}



//������� ��� �������� ������� ������������ ����� (�����������)
bool ExistGamiltonCycle(vector<vector<int>>& gamiltongraph) {
    vector<int> path(gamiltongraph.size(), -1);
    path[0] = 0;
    if (!FindGamilton(gamiltongraph, path, 1)) return false;
    return true;
}

vector<vector<int>> Gamilton(int ver, vector<vector<int>> undirected_weight_matrix) {
    vector<vector<int>> gamiltongraph(undirected_weight_matrix);

    // ������� ������� ���������
    cout << "\n������� �����: " << endl;
    PrintMatrix(gamiltongraph);

    // ��������� ������� ������������� ������
    vector<int> path(ver, -1);
    path[0] = 0;  // �������� � ������� 0
    int cycle_count = 0;

    cout << "\n����� ��������� ������������� ������." << endl;
    PrintAllHamiltonCyclesWithWeights(gamiltongraph, path, 1, cycle_count, undirected_weight_matrix);

    if (cycle_count > 0) {
        cout << "***********************" << endl;
        cout << "\n���� �������� �������������. ������� ������: " << cycle_count << endl;
        cout << "***********************" << endl;
        return gamiltongraph;
    }
    else {
        cout << "\n���� �� �������� �������������. ������������ ����." << endl;
        vector<int> vect(ver, -1);

        // ���������� ��������� ������������ ������
        for (int i = 0; i < ver; i++) {
            int r = rand() % ver;
            while (find(vect.begin(), vect.end(), r) != vect.end()) {
                r = rand() % ver;
            }
            vect[i] = r;
        }

        // ��������� ����������� �����
        for (int i = 0; i < ver; i++) {
            if (gamiltongraph[vect[i]][vect[(i + 1) % ver]] == 0) {
                int r = rand() % 10 + 1;
                gamiltongraph[vect[i]][vect[(i + 1) % ver]] = r;
                gamiltongraph[vect[(i + 1) % ver]][vect[i]] = r;
                cout << "��������� �����: " << vect[i] << " <-> " << vect[(i + 1) % ver]
                    << " � �����: " << r << endl;
            }
        }

        // ��������� ����� ����� �����������
        cycle_count = 0;
        fill(path.begin(), path.end(), -1);
        path[0] = 0;

        cout << "\n�������� ����� �����������." << endl;
        PrintAllHamiltonCyclesWithWeights(gamiltongraph, path, 1, cycle_count, undirected_weight_matrix);

        if (cycle_count > 0) {
            cout << "***********************" << endl;
            cout << "\n������ ���� �����������. ������� ������: " << cycle_count << endl;
            cout << "***********************" << endl;
        }
        else {
            cout << "\n�� ������� ������� ���� �������������." << endl;
        }
    }

    cout << "\n�������� ������� �����:" << endl;
    PrintMatrix(gamiltongraph);

    return gamiltongraph;
}

void FindGamiltonKomivvv(vector<vector<int>>& gamiltongraph, vector<bool>& visited, vector<int>& path, int vert, int ver, int count, int cost, int& min_cost, vector<int>& min_path) {
    path.push_back(vert);
    if (count == ver and gamiltongraph[vert][0]) {
        // ������� ������� ���� � �������
        for (int i = 0; i < path.size(); i++) cout << path[i] << " -> ";
        cout << path[0];
        cout << "\t��� �����: " << gamiltongraph[vert][0] + cost << endl;

        if (cost + gamiltongraph[vert][0] < min_cost) {
            min_cost = cost + gamiltongraph[vert][0];
            min_path = path;
        }
        path.pop_back();
        return;
    }
    for (int i = 0; i < ver; i++) {
        if (!visited[i] and gamiltongraph[vert][i] != 0) {
            visited[i] = true;
            FindGamiltonKomivvv(gamiltongraph, visited, path, i, ver, count + 1, cost + gamiltongraph[vert][i], min_cost, min_path);
            visited[i] = false;
        }
    }
    path.pop_back();
}


void Komivoyadger(vector<vector<int>>& gamiltongraph, int ver) {
    vector<bool> visited(ver); //������������ �������
    vector<int> path; //������� ����
    visited[0] = true;
    int min_cost = trash1; //����������� ��������� ����
    vector<int> min_path; //����������� ����

    cout << "***********************" << endl;
    cout << "��� ��������� ������������ �����:" << endl;

    FindGamiltonKomivvv(gamiltongraph, visited, path, 0, ver, 1, 0, min_cost, min_path);

    cout << "***********************" << endl;
    cout << "���� ������������: " << endl;
    for (int i = 0; i < min_path.size(); i++) cout << min_path[i] << " -> ";
    cout << min_path[0];
    cout << "\n���: " << min_cost << endl;
}


void menu() {

    printf("\n������� ����� �� 0 �� 19\n");
    printf("0) �����\n");
    printf("1) ������� ����\n");
    printf("2) ������� ������� ��������� ������\n");
    printf("3) ����� ��������\n");
    printf("4) ����������� ���������� �������� � �� ����������\n");
    printf("5) ����� ������ ����� ������� � ������\n");
    printf("6) ������� ������� �����\n");
    printf("7) �������� ��������-�����\n");
    printf("8) �������� �������� ������ 5 � 8\n");
    printf("9) ������� ���������� ������������\n");
    printf("10) ������� ���������\n");
    printf("11) ������������ ����� �� ��������� �����-����������\n");
    printf("12) ����� ����������� ���������\n");
    printf("13) ����� �������� �������� �� ��������\n");
    printf("14) ����������� �� ���� ����� �� ��������\n");
    printf("15) ��� ������� (���������� � ������������)\n");
    printf("16) ������������ ����������� ��������� �����\n");
    printf("17) ���������, �������� �� ���� ���������/�������������� ����, ���� �� �������.\n");
    printf("18) ���������, �������� �� ���� �������������/�������������� ����, ���� �� �����������.\n");
    printf("19) ������ ������ ������������ �� ������������� �����.\n");


}

