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
using namespace std;

int INF = 1e9;
//������� ��� ��������� �������� ������ � ����������� � �������������� ����
void GenerateVertexDegrees(vector<int>& degrees, int num_vertex) {
    float s = 2.0; //������ �� �������� �������� ������������ ��� ������� ������������� �������; ������ ���� >1 
    vector<double> probability(num_vertex); //������ ������������ ��������� ������ ��������� ������� ������ �� 1 �� (���� ���-�� ������-1)
    //������� ������ ������, ������������ ��� ��������� ��������� �����
    static random_device rd; //��������� ���������
    static mt19937 gen(rd());
    //��������� ������ ������������ � ������� ������� ����
    for (int k = 1; k <= num_vertex; ++k)
    {
        probability[k - 1] = 1.0 / pow(k, s); //k^(-s) = 1/(k^s) -1 ������ ��� ��� ���������� � probability ������ �������� � 0, � �� � 1 - ���� ��������� ������� �������
    }
    //������� ���������� ������������� � ����� ������������� (�����������, ����� ����� ���� ����� 1 � ������� ����������� ������)
    discrete_distribution<> distribution(probability.begin(), probability.end());
    //����� �������� distribution, ��� ����� ������������ ������ � ����������� ��������� ����� (������ gen ������ std::mt19937) 
    //��� ��������� ��������� ��������, ��� ����������� ������� �������� ��������������� ��� ����������� � �������������.
    bool flag = true;
    while (flag) {
        for (int i = 0; i < num_vertex; ++i) {
            degrees[i] = distribution(gen) + 1; //���������� ��������� ����� � �������������� ������������� ������������ (������ ����������� �� �������=���������� ����� ���������)
        }
        sort(degrees.begin(), degrees.end()); //��������� �� �����������
        reverse(degrees.begin(), degrees.end()); //��������� �� ��������
        degrees[num_vertex - 1] = 0;
        int pr_sum = 0;
        int sumVertexDegrees = 0;
        int max_deg = 0;
        for (int i = 0; i < num_vertex; i++) {
            if (degrees[i] <= num_vertex - i - 1) {
                pr_sum += 1;
                sumVertexDegrees += degrees[i];
            }
            //��������� ���-�� ������ � ���� ��������, ���� ����� ������� ���� => ����� ����
            else if (degrees[i] == num_vertex - 1) {
                max_deg += 1;
            }
        }
        if (max_deg == 0) {
            degrees[0] = num_vertex - 1;
            max_deg += 1;
        }

        if (pr_sum == num_vertex && sumVertexDegrees <= (num_vertex * num_vertex - num_vertex) / 2 && max_deg == 1) {
            flag = false;
        }
    }
}

//��������� ������� ��������� ���������������� �������� ������������� �����
//�������� ���-�� ������, ������ �������� ���� ������
//���������� ������� ���������
vector<vector<int>> Generate_adjacency_matrix(int num_vertex, const vector<int>& degrees) {
    //������� ������� ��������� - ������� �������� ���-�� ������*���-�� ������, ����������� ������
    vector<vector<int>> adjacency_matrix(num_vertex, vector<int>(num_vertex, 0));
    //������� ����� ������� �������� ������
    vector<int> copy_degrees(degrees);
    srand(time(0));
    for (int i = 0; i < num_vertex; i++) {
        for (int d = 0; d < copy_degrees[i]; d++) {
            bool n_addded_edge = true;
            while (n_addded_edge) {
                double R = static_cast<double>(rand()) / RAND_MAX; //����� � ��������� [0;1)
                int j = static_cast<int>((num_vertex)*R);
                //���������: ��� ����� ����� ��������� ��� �� �����������, ��� �� ������� �����, ��� ������� ������� ��� ������ ����
                if (adjacency_matrix[i][j] != 1 && i != j && copy_degrees[i] > 0 && j > i) {
                    adjacency_matrix[i][j] = 1;
                    n_addded_edge = false;
                }
            }
        }
    }
    return adjacency_matrix;
}

//�������� ��������
vector<vector<int>> Zhimbel(vector<vector<int>> weight_matrix, vector<vector<int>> matrix, const int num_vertex, const int parametr) {
    vector<vector<int>> new_matrix(num_vertex, vector<int>(num_vertex, 0));
    for (int i = 0; i < num_vertex; i++) {
        for (int j = 0; j < num_vertex; j++) {
            vector<int> vect;
            for (int k = 0; k < num_vertex; k++) {
                if (matrix[i][k] == 0 or weight_matrix[k][j] == 0) {
                    vect.push_back(0);
                }
                else {
                    vect.push_back(matrix[i][k] + weight_matrix[k][j]);
                }
            }
            sort(vect.begin(), vect.end()); //������������� �� �����������
            //���� ���� ������ �� ����� ==> �������� ����
            if (vect.back() == vect.front() && vect.back() == 0) {
                new_matrix[i][j] = 0;
            }
            else {
                if (parametr == 1) {
                    new_matrix[i][j] = vect.back(); //max
                }
                else if (parametr == 2) {
                    erase(vect, 0);
                    new_matrix[i][j] = vect.front(); //min
                }
            }
        }
    }
    return new_matrix;
}


//������������ ������� ������������
vector<vector<int>> Generate_reachability_matrix(vector<vector<int>> adjacency_matrix, vector<vector<int>> matrix, const int num_vertex) {
    vector<vector<int>> new_matrix(num_vertex, vector<int>(num_vertex, 0));
    //��������� ���������
    for (int i = 0; i < num_vertex; i++) {
        for (int j = 0; j < num_vertex; j++) {
            vector<int> vect;
            for (int k = 0; k < num_vertex; k++) {
                new_matrix[i][j] += (matrix[i][k] * adjacency_matrix[k][j]);
            }
        }
    }
    return new_matrix;
}

//�������� ������
vector<vector<int>> Addition_matrix(vector<vector<int>> matrix1, vector<vector<int>> matrix2, const int num_vertex) {
    vector<vector<int>> new_matrix(num_vertex, vector<int>(num_vertex, 0));
    for (int i = 0; i < num_vertex; i++) {
        for (int j = 0; j < num_vertex; j++) {
            new_matrix[i][j] = (matrix1[i][j] + matrix2[i][j]);
        }
    }
    return new_matrix;
}

//������� �������
void PrintMatrix(vector<vector<int>> matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            cout << setw(5) << val << setw(5);
        }
        cout << "\n";
    }
    cout << "\n";
}

//�������� �������
vector<int> GetPath(vector<int> from, int finish) {
    vector<int> path;
    //�������� � �����, ���� v!=-1, �� ����� ���������� ���� ����� v = from[v]
    for (int v = finish; v != -1; v = from[v]) {
        path.push_back(v);
    }
    reverse(path.begin(), path.end());
    return path;
}

//�������� ��������
vector<int> Dijkstra(vector<vector<int>> weight_matrix, int num_vertex, int start, int finish, vector<int>& path, int& iter_Dijkstra) {
    vector<int> from(num_vertex, -1);
    vector<int> distance(num_vertex, INF); //������ ����������
    distance[start] = 0; //�� ���� ������ ����� ��������� ������ ���������� ����, ��� ��������� ����
    set<pair<int, int>> unvisited_vertex; //������������� ��������� ����������� �������� � ������� ����������� (������� ������ ��������� ���� ������ ���������� �� �������, ������ - ����� �������)
    unvisited_vertex.insert(make_pair(distance[start], start));
    while (!unvisited_vertex.empty()) {
        int nearest = unvisited_vertex.begin()->second; //����� ��������� �������
        unvisited_vertex.erase(unvisited_vertex.begin());
        //�� ������ �������� �������� ���������� �� ����� �� ������
        for (int j = 0; j < num_vertex; j++) {
            iter_Dijkstra++;
            if (distance[j] > distance[nearest] + weight_matrix[nearest][j] && weight_matrix[nearest][j] != 0) {
                unvisited_vertex.erase(make_pair(distance[j], j)); //������ �������� ����� ��� ������� �������
                distance[j] = distance[nearest] + weight_matrix[nearest][j]; //��������� �����
                from[j] = nearest;
                unvisited_vertex.insert(make_pair(distance[j], j)); //���������� � ������������
            }
        }
    }
    path = GetPath(from, finish);
    return distance;
}

//�������� ������ � ������ BFS
vector<int> BFS(vector<vector<int>> weight_matrix, int num_vertex, int start, int finish, vector<int>& path, int& iter_BFS) {
    vector<int> from(num_vertex, -1);
    vector<int> distance(num_vertex, INF); //������ ���������� �� ������, ���������� �� ���� ������ ����, ����� ��������� �������
    queue<int> q; //������� ������, ��������� ���������
    distance[start] = 0;
    q.push(start);
    //���� ������� �� �����
    while (!q.empty()) {
        int vert = q.front(); //��������� ����� �������, ������� ����� ������������
        q.pop(); //������� ��� ������� �� �������
        for (int to = 0; to < num_vertex; to++) {
            iter_BFS++;
            //���������, ��� ������� ������ � ��������� � ��� ���������� �� ������� ������, ��� ���������� �� ����������+��� ����
            if (weight_matrix[vert][to] != 0 && distance[to] > distance[vert] + weight_matrix[vert][to]) {
                distance[to] = distance[vert] + weight_matrix[vert][to];
                q.push(to);
                from[to] = vert;
            }
        }
    }
    path = GetPath(from, finish);
    return distance;
}

//�������� ������ ������������� ���� 
vector<int> MaxWay(vector<vector<int>> weight_matrix, int num_vertex, int start, int finish, vector<int>& path) {
    vector<int> from(num_vertex, -1);
    vector<int> distance(num_vertex, INF); //������ ����������
    distance[start] = 0; //�� ���� ������ ����� ��������� ������ ���������� ����, ��� ��������� ����
    set<pair<int, int>> unvisited_vertex; //������������� ��������� ����������� �������� � ������� ����������� (������� ������ ��������� ���� ������ ���������� �� �������, ������ - ����� �������)
    unvisited_vertex.insert(make_pair(distance[start], start));
    while (!unvisited_vertex.empty()) {
        int nearest = unvisited_vertex.begin()->second; //����� ��������� �������
        unvisited_vertex.erase(unvisited_vertex.begin());
        //�� ������ �������� �������� ���������� �� ����� �� ������
        for (int j = 0; j < num_vertex; j++) {
            if (distance[j] > distance[nearest] + weight_matrix[nearest][j] * (-1) && weight_matrix[nearest][j] != 0) {
                unvisited_vertex.erase(make_pair(distance[j], j));
                distance[j] = distance[nearest] + weight_matrix[nearest][j] * (-1);
                from[j] = nearest;
                unvisited_vertex.insert(make_pair(distance[j], j));
            }
        }
    }
    path = GetPath(from, finish);
    for (int i = 0; i < num_vertex; i++) {
        distance[i] = distance[i] * (-1);
    }
    return distance;
}


int FF_bfs(vector<vector<int>> bandwidth_matrix, int start, int finish, int num_vertex, vector<int>& parent) {
    //��������� ������ ����������������, ��������� -1
    fill(parent.begin(), parent.end(), -1);
    //�������������� ��������� ������� = -2
    parent[start] = -2;
    //������� �������, �������� ����: 1�� ������� ���� - ����� �������, 2�� - ��� ��������������
    queue<pair<int, int>> q;
    //��������� � ������� ��������� �������
    q.push({ start,INF });

    //���� �� �����, ���� ������� �� �����
    while (!q.empty()) {
        //��������� ������� ���� ������� � ��� �������������� �� ������ ������
        int u = q.front().first;
        int cap = q.front().second;
        //������� �� ������� ���� ����
        q.pop();
        //�������� �� ���� �������� ������� � u
        for (int v = 0; v < num_vertex; v++) {
            //������� ����,� ������� �� ��� �� �������� (parent[v]==-1), ������� �� ����� u, � �������� ��� ����� ��������� (�������������� != 0)
            if (u != v && bandwidth_matrix[u][v] != 0 && parent[v] == -1) {
                //���������� u ��� �������� (���������������) v
                parent[v] = u;
                //��������� ����������� ��������������
                int min_cap = min(cap, bandwidth_matrix[u][v]);
                //���� ��������� �� �������� �������, �� ���������� ��������� ����������� ��������������
                if (v == finish) {
                    return min_cap;
                }
                //���� �� ��������� �� �������� �������, �� ��������� ������� v � ��� �������������� � ������ ������ � �������
                q.push({ v,min_cap });
            }
        }
    }
    //���� �� ����� ���� ����� ��������� � �������� ��������, �� ���������� ����
    return 0;
}

int Ford_Fulkerson(vector<vector<int>> bandwidth_matrix, int start, int finish, int num_vertex) {
    //�������������� ������������ ������ (parent[1]=0 ��������, ��� �� ���� 0 ����� ������� � 1�� = �������������� 1�� ���� - ������� ����) ��� ���������� ���� �� ��������� ������� � ��������
    vector<int> parent(num_vertex, -1);
    //�������� ������������� ������
    int max_flow = 0;
    //����� ��������� � ���������� �������� ��� �������������� ��� ������� ����
    int min_cap = 0;

    //���� �� �����, ���� ��� �������������� �� ������ ������� (���� ���� ���� �� ��������� ������� � ��������)
    //����� ����� ����� � ��� �������������� �������� ������� FF_bfs
    while (min_cap = FF_bfs(bandwidth_matrix, start, finish, num_vertex, parent)) {
        //��������� � ���� ����� ��� �������������� �������� ����
        max_flow += min_cap;
        //���������� � v �������� �������
        int v = finish;

        //���� �� ����� �� ��������� �������, ����� ���� �� �������� �� ���������
        while (v != start) {
            //���������� �������� ������� v
            int u = parent[v];
            //��������� �������������� �� u � v �� ��� �������������� �������� ����
            bandwidth_matrix[u][v] -= min_cap;
            //����������� (��� ��� � ��������������� �������) �������������� �� v � u �� ��� �������������� �������� ����
            bandwidth_matrix[v][u] += min_cap;
            //���������� � v ��� �������� (������������ �� ����� � ��������� �������)
            v = u;
        }
    }
    //���������� ���� �����
    return max_flow;
}


//���������������� �������� ��������
void modyDijkstra(vector<vector<int>>& weight_matrix, int num_vertex, int start, int finish, vector<int>& parent) {
    vector<int> distance(num_vertex, INF); //������ ����������
    distance[start] = 0; //�� ���� ������ ����� ��������� ������ ���������� ����, ��� ��������� ����
    set<pair<int, int>> unvisited_vertex; //������������� ��������� ����������� �������� � ������� ����������� (������� ������ ��������� ���� ������ ���������� �� �������, ������ - ����� �������)
    unvisited_vertex.insert(make_pair(distance[start], start));
    while (!unvisited_vertex.empty()) {
        int nearest = unvisited_vertex.begin()->second; //����� ��������� �������
        unvisited_vertex.erase(unvisited_vertex.begin());
        //�� ������ �������� �������� ���������� �� ����� �� ������
        for (int j = 0; j < num_vertex; j++) {
            if (distance[j] > distance[nearest] + weight_matrix[nearest][j] && weight_matrix[nearest][j] != 0) {
                unvisited_vertex.erase(make_pair(distance[j], j)); //������ �������� ����� ��� ������� �������
                distance[j] = distance[nearest] + weight_matrix[nearest][j]; //��������� �����
                parent[j] = nearest;
                unvisited_vertex.insert(make_pair(distance[j], j)); //���������� � ������������
            }
        }
    }
}

int mincost_maxflow(vector<vector<int>> cost_matrix, vector<vector<int>> bandwidth_matrix, int s, int t, int num_vertex, vector<int> parent, int need_flow) {
    int ans = 0;
    int flow = 0;
    vector<vector<int>> flow_matrix(num_vertex, vector<int>(num_vertex, 0));
    while (flow < need_flow) {
        vector<int> parent(num_vertex, -1);
        modyDijkstra(cost_matrix, num_vertex, s, t, parent);
        int delta = INF;
        for (int v = t; v != s; v = parent[v]) {
            delta = min(delta, bandwidth_matrix[parent[v]][v]); //����� ����� �� ����� ����
        }
        int razn = need_flow - flow;
        int flag = 0;
        if (delta > razn) {
            delta = razn;
            flag = 1;
        }
        flow += delta;

        for (int v = t; v != s; v = parent[v]) {
            bandwidth_matrix[parent[v]][v] -= delta;
            bandwidth_matrix[v][parent[v]] += delta; //��������� ���� � ��������������� �������
            cost_matrix[v][parent[v]] = -cost_matrix[parent[v]][v];
            flow_matrix[v][parent[v]] = delta;
            if (bandwidth_matrix[parent[v]][v] == 0) {
                cost_matrix[parent[v]][v] = INF;
            }
            ans += abs(cost_matrix[v][parent[v]] * flow_matrix[v][parent[v]]);
        }
    }
    for (int i = 0; i < t + 1; i++) {
        for (int j = 0; j < t + 1; j++) {
            if (flow_matrix[i][j] != 0) {
                std::cout << "�����, ���������� �� ����� " << j << "->" << i << " �����: " << flow_matrix[i][j] << ", ��������� �����: " << abs(cost_matrix[i][j]) << endl;
            }
        }
    }

    return ans;
}

void PrintWayInfo(int num_vertex, vector<int> dist, int from, int where_, vector<int> path) {
    cout << "������ ����������: " << endl;
    for (int i = 0; i < num_vertex; i++) {
        if (dist[i] != INF and dist[i] != -INF) {
            cout << dist[i] << "\t";
        }
        else {
            cout << "-" << "\t";
        }

    }
    if (dist[where_] != INF and dist[where_] != -INF) {
        cout << "\n����� ��������: " << dist[where_] << ".";
        cout << "\n�������: ";
        for (int i = 0; i < path.size(); i++) {
            if (i == path.size() - 1) {
                cout << path[i] << ".";
            }
            else {
                cout << path[i] << "->";
            }
        }
    }
    else {
        cout << "\n����� ������� ��������� ������" << endl;
    }
}



//��� �������� �� int 
int intinput()
{
    int val;
    while (true)
    {
        std::cin >> val;
        if (std::cin.peek() != '\n')
        {
            std::cin.clear();
            std::cin.ignore(std::cin.rdbuf()->in_avail());
            printf("\n�� ����� ����������� ��������. ��������� ����.\n");
        }
        else break;
    }
    return val;
}

void GetStartFinishVertex(int& start, int& finish, int num_vertex) {
    cout << "������� ����� ������� �� ������� ������ ��������� ������� �� 0 �� " << num_vertex - 1 << ":" << endl;
    start = intinput(); //�� ����� ������� ������ �������
    while (start < 0 || start > num_vertex - 1)
    {
        cout << "������������ ����. ������� ����� ������� �� ������� ������ ��������� ������� �� 0 �� " << num_vertex - 1 << ":" << endl;
        start = intinput();
    }

    cout << "������� ����� ������� � ������� ������ ��������� ������� �� 0 �� " << num_vertex - 1 << ":" << endl;
    finish = intinput(); //� ����� ������� ����� ���������
    while (finish < 0 || finish > num_vertex - 1)
    {
        cout << "������������ ����. ������� ����� ������� � ������� ������ ��������� ������� �� 0 �� " << num_vertex - 1 << ":" << endl;
        finish = intinput();
    }

}



//������������ �������
int determinant(int num_vertex, vector<vector<int>> kirhgof_matrix) {
    int det = 0;
    vector<vector<int>> submatrix(num_vertex, vector<int>(num_vertex, 0));

    if (num_vertex == 2) {
        return ((kirhgof_matrix[0][0] * kirhgof_matrix[1][1]) - (kirhgof_matrix[1][0] * kirhgof_matrix[0][1]));
    }
    else {
        for (int x = 0; x < num_vertex; x++) {
            int subi = 0;
            for (int i = 1; i < num_vertex; i++) {
                int subj = 0;
                for (int j = 0; j < num_vertex; j++) {
                    if (j == x) {
                        continue;
                    }
                    submatrix[subi][subj] = kirhgof_matrix[i][j];
                    subj++;
                }
                subi++;
            }
            det = det + (pow(-1, x) * kirhgof_matrix[0][x] * determinant(num_vertex - 1, submatrix));
        }
    }
    return det;
}

//��������
//��������� ��� �������� �����
struct edge {
    int cost;
    int from;
    int where_;
};

//����������, ����� ���������� � ������� ��������
struct compare {
    bool operator() (edge const& a, edge const& b) const {
        return a.cost > b.cost;
    }
};

//����������, ����� ���������� � ������� ��������
struct comparenum {
    bool operator() (edge const& a, edge const& b) const {
        return a.from < b.from;
    }
};

//����� ������ ��������� ������, � �������� ����������� �������
int find_root(int k, vector<int>parent) {
    if (parent[k] == k) {
        return k;
    }
    return find_root(parent[k], parent);
}

void union_(int f_p, int s_p, vector<int>& parent) {
    int x = find_root(f_p, parent);
    int y = find_root(s_p, parent);
    parent[x] = y;
}

vector<edge> Kruskal(vector<vector<int>> weight_matrix, int num_vertex, int& mst_weight, int& iter_kruskal) {
    vector<edge> edges;
    //��������� ������ edges, �������� ������ <���,������ �����, ���� �����>
    for (int i = 0; i < num_vertex; i++) {
        for (int j = 0; j < num_vertex; j++) {
            if (weight_matrix[i][j] != 0) {
                edges.push_back({ weight_matrix[i][j],i,j });
            }
        }
    }
    //��������� �� �������� ���� �����
    sort(edges.begin(), edges.end(), compare());

    vector<int> parent(num_vertex);
    for (int i = 0; i < num_vertex; i++) {
        parent[i] = i;
    }

    vector<edge> mst;
    while (mst.size() != num_vertex - 1) {
        edge next_edge = edges.back();
        edges.pop_back();

        //���� ������ ��������� ������ ������ ������� ����� � ������ ��������� ������ ������ ������� �����
        int f_p = find_root(next_edge.from, parent);
        int s_p = find_root(next_edge.where_, parent);

        //���� ����� �� �������, �� ���������� ����� �� ������� ����
        if (f_p != s_p) {
            mst.push_back(next_edge);
            mst_weight += next_edge.cost;
            union_(f_p, s_p, parent);
        }
        iter_kruskal++;
    }
    return mst;
}

//�������
void union_sets(int f_p, int s_p, vector<int>& parent, vector<int>& rank) {
    int x = find_root(f_p, parent);
    int y = find_root(s_p, parent);
    if (rank[x] < rank[y]) {
        parent[x] = y;
    }
    else if (rank[x] > rank[y]) {
        parent[y] = x;
    }
    else {
        parent[y] = x;
        rank[x]++;
    }
}

vector<edge> Boruvka(vector<vector<int>> weight_matrix, int num_vertex, int& mst_weight, int& iter_boruvka) {
    vector<edge> edges;
    vector<edge> mst;
    //��������� ������ edges, �������� ������ <���,������ �����, ���� �����>
    for (int i = 0; i < num_vertex; i++) {
        for (int j = 0; j < num_vertex; j++) {
            if (weight_matrix[i][j] != 0) {
                edges.push_back({ weight_matrix[i][j],i,j });
            }
        }
    }

    vector<int> parent(num_vertex);
    vector<int> rank(num_vertex);
    vector<int> cheapest(num_vertex);
    for (int i = 0; i < num_vertex; i++) {
        parent[i] = i;
        rank[i] = 0;
        cheapest[i] = 0;
    }
    int num_trees = num_vertex;

    while (num_trees > 1) {
        for (int i = 0; i < edges.size(); i++) {
            int f_p = find_root(edges[i].from, parent);
            int s_p = find_root(edges[i].where_, parent);
            if (f_p != s_p) {
                if (cheapest[f_p] == -1 || edges[cheapest[f_p]].cost > edges[i].cost) {
                    cheapest[f_p] = i;
                }
                if (cheapest[s_p] == -1 || edges[cheapest[s_p]].cost > edges[i].cost) {
                    cheapest[s_p] = i;
                }
            }
            iter_boruvka++;
        }
        for (int i = 0; i < num_vertex; i++) {
            if (cheapest[i] != -1) {
                int f_p = find_root(edges[cheapest[i]].from, parent);
                int s_p = find_root(edges[cheapest[i]].where_, parent);
                if (f_p != s_p) {
                    mst.push_back(edges[cheapest[i]]);
                    mst_weight += edges[cheapest[i]].cost;
                    union_sets(f_p, s_p, parent, rank);
                    num_trees--;
                }
                cheapest[i] = -1;
            }
            iter_boruvka++;
        }
    }
    return mst;
}

vector<pair<int, int>> codePrufer(vector<edge> mst, int num_vertex) {
    vector<vector<int>> mst_matrix(num_vertex, vector<int>(num_vertex, 0));
    for (int k = 0; k < mst.size(); k++) {
        int i = mst[k].from;
        int j = mst[k].where_;
        int w = mst[k].cost;
        mst_matrix[i][j] = w;
        mst_matrix[j][i] = w;
    }
    //����: ��� �����, ����� ��������
    vector<pair<int, int>> prufer;
    vector<int> degree_vert(num_vertex);
    for (int i = 0; i < num_vertex; i++) {
        for (int j = 0; j < num_vertex; j++) {
            if (mst_matrix[i][j] != 0) {
                degree_vert[i]++;
            }
        }
    }

    int r = 0;
    for (int i = 0; i < num_vertex; i++) {
        r += degree_vert[i];
    }

    while (r != 0) {
        for (int i = 0; i < num_vertex; i++) {
            if (degree_vert[i] == 1) {
                for (int j = 0; j < num_vertex; j++) {
                    if (mst_matrix[i][j] != 0) {
                        prufer.push_back(make_pair(mst_matrix[i][j], j));
                        degree_vert[i]--;
                        degree_vert[j]--;
                        mst_matrix[i][j] = 0;
                        mst_matrix[j][i] = 0;
                        r -= 2;
                    }
                }
            }
        }
    }

    return prufer;
}

vector<vector<int>> decodePrufer(vector<pair<int, int>> prufer) {
    int num_vertex = prufer.size() + 1;
    vector<vector<int>> ans_matrix(num_vertex, vector<int>(num_vertex, 0));
    vector<int> vertexes(num_vertex, 0);
    for (int i = 0; i < prufer.size(); i++) {
        vertexes[prufer[i].second] += 1;
    }

    int j = 0;
    int num = 0;
    for (int i = 0; i < num_vertex - 1; i++) {
        for (j = 0; j < num_vertex; j++) {
            if (i == num_vertex - 2) {
                if (vertexes[j] == 0) {
                    vertexes[j] = -1;
                    ans_matrix[j][prufer[i].second] = prufer[i].first;
                    ans_matrix[prufer[i].second][j] = prufer[i].first;
                    vertexes[prufer[i].second]--;
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
    return ans_matrix;
}

vector<vector<int>> MakeEulerGraph(int num_vertex, vector<vector<int>> undirected_weight_matrix, vector<int> undirected_degrees) {
    vector<vector<int>> eulergraph(undirected_weight_matrix);
    int countevenvert = 0;
    bool flagEuler = false;
    //������� ������� ���������
    cout << "\n������� �����: " << endl;
    PrintMatrix(eulergraph);
    //����� ���� �������� ������
    for (int i = 0; i < num_vertex; i++) {
        cout << "������� ������� �" << i << " = " << undirected_degrees[i] << endl;
    }
    //���������, ��� ��� ������� ������ ����� ������ => ������� ����
    for (int i = 0; i < num_vertex; i++) {
        if (undirected_degrees[i] % 2) {
            countevenvert++;
        }
    }
    if (countevenvert == num_vertex) {
        cout << "���� �������� ���������\n" << endl;
        flagEuler = true;
    }
    //���� �� ��� ������� ������ ����� ������ => ���� �������������� ���� - ������� ���������
    else {
        cout << "\n���� �� �������� ���������. ������������ ����.\n" << endl;
        //����� ������� ��� ��������� ����� ���� ���� �� ������ ���������
        while (!flagEuler) {
            bool flag = 0;
            for (int i = 0; i < num_vertex; i++) {
                for (int j = 0; j < num_vertex; j++) {
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
                if (flag) {
                    break;
                }
            }
            //���� �� ���������� �������� �����, �� ������� �������
            if (!flag) {
                for (int i = 0; i < num_vertex; i++) {
                    for (int j = 0; j < num_vertex; j++) {
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
                    if (flag) {
                        break;
                    }
                }
            }
            //��������� ���������� �� ������� ���� � ������� ��������� ������
            int countevenvert1 = 0;
            for (int i = 0; i < num_vertex; i++) {
                if (undirected_degrees[i] % 2 == 0) {
                    countevenvert1++;
                }
            }
            //���� ����������, �� ��������� ���� � ������� �� �����, ���� ��� - ���������� ���������/������� �����
            if (countevenvert1 == num_vertex) {
                flagEuler = true;
            }
        }
        cout << "������ ���� �������� ���������.\n" << endl;
        //������� ������� ���������
        cout << "������� �����: " << endl;
        PrintMatrix(eulergraph);
        //����� ���� �������� ������
        for (int i = 0; i < num_vertex; i++) {
            cout << "������� ������� �" << i << " = " << undirected_degrees[i] << endl;
        }
    }
    return eulergraph;
}

void FindEulerCycle(int num_vertex, vector<vector<int>> eulergraph) {
    stack<int> vertexes;
    vector<int> cycle;
    int vert = 0;
    //���� ������ �������, �� ������� ���� �����
    for (vert = 0; vert < num_vertex; vert++) {
        for (int j = 0; j < num_vertex; j++) {
            if (eulergraph[vert][j] != 0) {
                break;
            }
        }
        if (vert < num_vertex) {
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
        for (i = 0; i < num_vertex; i++) {
            if (eulergraph[vert][i] != 0) {
                break;
            }
        }
        //���� ������� �� �������, �� ��������� ��� ������� � ����, ������� �� �����
        if (i == num_vertex) {
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


//��������������� ������� ��� �������� ����������� ���������� ������� � ����
bool CanGoNext(int v, vector<vector<int>>& gamiltongraph, vector<int>& path, int position) {
    //���������, ���� �� ����� ����� ��������� �������� path[position - 1] � �������� v � ������� ���������
    //���� ���, �� ���������� false => ������� �� �� ��������� � ����
    if (gamiltongraph[path[position - 1]][v] == 0) {
        return false;
    }
    //����� ��������� �� ���� �� ������� ��� ��������
    //���� ����, �� ���������� false => ������� �� �� ��������� � ����
    for (int i = 0; i < position; i++) {
        if (path[i] == v) {
            return false;
        }
    }
    return true;
}


//������� ��� ������ ������������ ����� (�����������)
bool FindGamiltonCycle(vector<vector<int>>& gamiltongraph, vector<int>& path, int position) {
    //���� ��� ������� �������� � ��������� ��������� � ������, �� ������������ true => ���� ����������� ���� � �����
    if (position == gamiltongraph.size()) {
        return gamiltongraph[path[position - 1]][path[0]] != 0;
    }
    //����� ���� ��������� ������� ��� ���������
    for (int v = 1; v < gamiltongraph.size(); v++) {
        //� ������� ���� ����������� ����� (���� ��������� ����� + �� ���� �������� ������� �������)
        if (CanGoNext(v, gamiltongraph, path, position)) {
            //���� �������, �� ��������� ��� ������� � ����
            path[position] = v;
            //���� ����������� ���� ������, �� ���������� true
            if (FindGamiltonCycle(gamiltongraph, path, position + 1)) {
                return true;
            }
            //���� ����������� ���� �� ������, �� ������� ������� v �� ����
            path[position] = -1;
        }
    }
    return false;
}


//������� ��� �������� ������� ������������ ����� (�����������)
bool ExistGamiltonCycle(vector<vector<int>>& gamiltongraph) {
    vector<int> path(gamiltongraph.size(), -1);
    path[0] = 0;
    if (!FindGamiltonCycle(gamiltongraph, path, 1)) {
        return false;
    }
    return true;
}


//����������� ����: ��������, �����������
vector<vector<int>> MakeGamiltonGraph(int num_vertex, vector<vector<int>> undirected_weight_matrix) {
    vector<vector<int>> gamiltongraph(undirected_weight_matrix);
    //������� ������� ���������
    cout << "\n������� �����: " << endl;
    PrintMatrix(gamiltongraph);
    if (ExistGamiltonCycle(gamiltongraph)) {
        cout << "���� �������� �������������\n" << endl;
    }
    else {
        cout << "\n���� �� �������� �������������. ������������ ����." << endl;
        vector<int> vect(num_vertex, -1);
        //����� ��������� �����, ���� ���� �� ������ �������������
        for (int i = 0; i < num_vertex; i++) {
            int r = rand() % ((num_vertex - 1) - 0 + 1) + 0;
            //������� find ����� ��������������� ���� �� ������� � ������ �� �����, 
            //���� �� ������ �������, ������� ����� �������� ��� ������ r. 
            //���� �������� r �� ����� �������, �� find ������ ��������� �� ����� ������� => ������� ��� ��������
            if (find(vect.begin(), vect.end(), r) == vect.end()) {
                vect[i] = r;
            }
            else {
                i--;
            }
        }
        //�������� �� ���� ��������� ������� � ��������� ������� ����� ����� ��������������� ���������
        //���� ������ ����� ��� - ��������� ���������
        for (int i = 0; i < num_vertex; i++) {
            if (gamiltongraph[vect[i]][vect[(i + 1) % num_vertex]] == 0) {
                int r = rand() % (10 - (1) + 1) + 1;
                gamiltongraph[vect[i]][vect[(i + 1) % num_vertex]] = r;
                gamiltongraph[vect[(i + 1) % num_vertex]][vect[i]] = r;
                cout << "�������� �����: " << vect[i] << "<->" << vect[(i + 1) % num_vertex] << " � �����: " << r << endl;
            }
        }
        if (ExistGamiltonCycle(gamiltongraph)) {
            cout << "������ ���� �������� �������������." << endl;
            //������� ������� ���������
            cout << "������� �����: " << endl;
            PrintMatrix(gamiltongraph);
        }
        else {
            cout << "���� �� ������� �������������� � �� �� �������� �������������\n" << endl;
        }
    }
    return gamiltongraph;
}

void findGamiltonCycleforTSP(vector<vector<int>>& gamiltongraph, vector<bool>& visited, vector<int>& path, int vert, int num_vertex, int count, int cost, int& min_cost, vector<int>& min_path, ofstream& fout) {
    //��������� � ���� ������� �������
    path.push_back(vert);
    //���� ��� ������� ����� �������� � ��������� ������� ��������� � �������,
    //�� ���������� ���� � ��� ����� � ����
    if (count == num_vertex and gamiltongraph[vert][0]) {
        for (int i = 0; i < path.size(); i++) {
            fout << path[i] << " -> ";
            /* if (i != path.size() - 1) {
                 fout << path[i] << " -> ";
             }
             else {
                 fout << path[i];
             } */
        }
        fout << path[0];
        fout << "\t��� �����: " << gamiltongraph[vert][0] + cost << endl;
        //���� ��� ����� ������ ������������, �� ����������� ��������� � ���� ���������
        if (cost + gamiltongraph[vert][0] < min_cost) {
            min_cost = cost + gamiltongraph[vert][0];
            min_path = path;
        }
        path.pop_back();
        return;
    }
    for (int i = 0; i < num_vertex; i++) {
        if (!visited[i] and gamiltongraph[vert][i] != 0) {
            visited[i] = true;
            findGamiltonCycleforTSP(gamiltongraph, visited, path, i, num_vertex, count + 1, cost + gamiltongraph[vert][i], min_cost, min_path, fout);
            visited[i] = false;
        }
    }
    path.pop_back();
}

void TransportSalesmanProblem(vector<vector<int>>& gamiltongraph, int num_vertex, ofstream& fout) {
    vector<bool> visited(num_vertex); //������������ �������
    vector<int> path; //������� ����
    visited[0] = true;
    int min_cost = INF; //����������� ��������� ����
    vector<int> min_path; //����������� ����
    findGamiltonCycleforTSP(gamiltongraph, visited, path, 0, num_vertex, 1, 0, min_cost, min_path, fout);
    cout << "���� ������������: " << endl;
    for (int i = 0; i < min_path.size(); i++) {
        cout << min_path[i] << " -> ";
    }
    cout << min_path[0];
    cout << "\n���: " << min_cost << endl;
    fout.close();
}




//����� ����
void printmenu()
{
    printf("\n������� ����� �� 0 �� 18:\n");
    printf("0 - ��������� ������ ���������\n");
    printf("1 - ������������� ����� ����\n");
    printf("\n������������ �1\n");
    printf("2 - ������� ������� ��������� ������\n");
    printf("3 - ������� ������� �����\n");
    printf("4 - ����� ��������\n");
    printf("5 - ����� ���-�� ��������� �� ����� ������� �� ������\n");
    printf("\n������������ �2\n");
    printf("6 - �������� ��������\n");
    printf("7 - �������� ������ � ������ BFS\n");
    printf("8 - �������� ������ ������������� ����\n");
    printf("\n������������ �3\n");
    printf("9 - ������� ������� ���������� ������������\n");
    printf("10 - ������� ������� ���������\n");
    printf("11 - ����� ������������ ����� �� ��������� �����-����������\n");
    printf("12 - ����� �������� ����� ����������� ���������\n");
    printf("\n������������ �4\n");
    printf("13 - ������� ������� ��������\n");
    printf("14 - ����� ����� �������� ��������\n");
    printf("15 - ��������� ����������� �� ���� ����� �� ��������� ��������\n");
    printf("16 - ��������� ����������� �� ���� ����� �� ��������� �������\n");
    printf("17 - ������������ � ������������ ����� � ������� ���� �������\n");
    printf("\n������������ �5\n");
    printf("18 - ���������: �������� �� ���� ���������. ���� ���, �� ��������������. ��������� ������� ����.\n");
    printf("19 - ���������: �������� �� ���� �������������. ���� ���, �� ��������������. ��������� ������������ �����, �������� � ����. ������ ������������.\n\n");
}

//������� ����
void menu(bool& exitflag) {
    while (exitflag) {
        cout << "������� ���-�� ������ � �����:" << endl;
        int num_vertex = intinput(); //���-�� ������ - ������ ������������
        while (num_vertex < 2 || num_vertex > 100)
        {
            printf("\n���-�� ������ ����� ����� ���� �� 2 �� 100. ��������� ����.\n");
            num_vertex = intinput();
        }
        int num_edge = num_vertex - 1; //���-�� ����� � ������� ������������ ����� = ���-�� ������ - 1
        vector<int> degrees(num_vertex); //������ �������� ������  
        int res_sumVertexDegrees = 0;
        GenerateVertexDegrees(degrees, num_vertex); //���������� ������� ���� ������

        //����� ���� �������� ������
        for (int i = 0; i < num_vertex; i++) {
            cout << "������� ������� �" << i << " = " << degrees[i] << endl;
        }

        //������� ������� ���������
        vector<vector<int>> adjacency_matrix = Generate_adjacency_matrix(num_vertex, degrees);

        //������� ������� �����
        cout << "\n��������: 1 - ���� � �������������� ������, 2 - ���� � �������������� � �������������� ������:" << endl;
        int sign = intinput(); //���-�� ������ - ������ ������������
        while (sign < 1 || sign > 2)
        {
            cout << "������������ ����. ��������: 1 - ���� � �������������� ������, 2 - ���� � �������������� � �������������� ������:" << endl;
            sign = intinput();
        }
        srand(time(0));
        vector<vector<int>> weight_matrix(adjacency_matrix);
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (adjacency_matrix[i][j] != 0) {
                    if (sign == 1) {
                        //rand()%(end-start+1)+start
                        weight_matrix[i][j] = rand() % (10 - (1) + 1) + 1;
                    }
                    else if (sign == 2) {
                        //rand()%(end-start+1)+start
                        weight_matrix[i][j] = rand() % (10 - (-10) + 1) + (-10);
                    }
                }
            }
        }

        vector<vector<int>> ed_matrix(num_vertex, vector<int>(num_vertex, 0)); //��������� �������
        //��������� ��������� �������
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (i == j) {
                    ed_matrix[i][j] = 1;
                }
            }
        }

        //��������� ������� ������������
        vector<vector<int>> reachability_matrix(adjacency_matrix);
        vector<vector<int>> res_reachability_matrix(num_vertex, vector<int>(num_vertex, 0));
        res_reachability_matrix = Addition_matrix(ed_matrix, adjacency_matrix, num_vertex);
        for (int i = 0; i < num_vertex; i++) {
            reachability_matrix = Generate_reachability_matrix(adjacency_matrix, reachability_matrix, num_vertex);
            res_reachability_matrix = Addition_matrix(res_reachability_matrix, reachability_matrix, num_vertex);
        }

        //��������� ������� ���������� ������������ (������� �����)
        vector<vector<int>> bandwidth_matrix(num_vertex, vector<int>(num_vertex, 0));
        srand(time(0));
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (adjacency_matrix[i][j] != 0) {
                    //rand()%(end-start+1)+start
                    bandwidth_matrix[i][j] = rand() % (10 - (1) + 1) + 1;
                }
            }
        }

        //��������� ������� ��������� (������� � ����� �����)
        vector<vector<int>> cost_matrix(num_vertex, vector<int>(num_vertex, 0));
        srand(time(0) + 1);
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (adjacency_matrix[i][j] != 0) {
                    //rand()%(end-start+1)+start
                    cost_matrix[i][j] = rand() % (10 - (1) + 1) + 1;
                }
            }
        }


        //��������� ������� ��������� ������������������ �����
        vector<vector<int>> undirected_adjacency_matrix(num_vertex, vector<int>(num_vertex, 0));
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (adjacency_matrix[i][j] != 0) {
                    undirected_adjacency_matrix[j][i] = adjacency_matrix[i][j];
                    undirected_adjacency_matrix[i][j] = adjacency_matrix[i][j];
                }
            }
        }

        //��������� ������� �� ��������� ������ �� ������� ��������� ��� ������������������ �����
        vector<vector<int>> degree_matrix(num_vertex, vector<int>(num_vertex, 0));
        for (int i = 0; i < num_vertex; i++) {
            //int d = 0;
            for (int j = 0; j < num_vertex; j++) {
                if (undirected_adjacency_matrix[i][j] != 0) {
                    degree_matrix[i][i] += 1;
                }
            }

        }

        //��������� ������ �������� ������ ������������������ �����
        vector<int> undirected_degrees(num_vertex, 0);
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (degree_matrix[i][j] != 0) {
                    undirected_degrees[i] = degree_matrix[i][j];
                }
            }
        }

        //��������� ������� ��������
        vector<vector<int>> kirhgof_matrix(num_vertex, vector<int>(num_vertex, 0));
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                kirhgof_matrix[i][j] = degree_matrix[i][j] - undirected_adjacency_matrix[i][j];
            }
        }


        //��������� ������� ����� ������������������ �����
        vector<vector<int>> undirected_weight_matrix(num_vertex, vector<int>(num_vertex, 0));
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (cost_matrix[i][j] != 0) {
                    undirected_weight_matrix[j][i] = weight_matrix[i][j];
                    undirected_weight_matrix[i][j] = weight_matrix[i][j];
                }
            }
        }



        printmenu(); //������� ����
        int iter_Dijkstra = 0;
        int iter_BFS = 0;
        float iter_div = 0;
        ofstream fout;
        fout.open("GamiltonCycles.txt");
        int input = 0;
        bool flag = true;
        while (flag)
        {
            input = intinput();
            switch (input)
            {
            case 0: {
                printf("\n�� ������� ������� 0 - ���������� ������ ���������\n");
                flag = false;
                exitflag = false;
                break;
            }
            case 1: {
                printf("\n�� ������� ������� 1 - ������������� ����� ����\n");
                flag = false;
                break;
            }
            case 2: {
                printf("\n�� ������� ������� 2 - ������� ������� ��������� ������\n");
                PrintMatrix(adjacency_matrix);
                cout << "\n";
                printmenu();
                break;
            }
            case 3: {
                printf("\n�� ������� ������� 3 - ������� ������� �����\n");
                PrintMatrix(weight_matrix);
                cout << "\n";
                printmenu();
                break;
            }
            case 4: {
                printf("\n�� ������� ������� 4 - ����� ��������\n");
                if (num_vertex == 1) {
                    cout << "������ ������������ ����� �������� ��� ������� 1 �������" << endl;
                }
                else {
                    cout << "������� ����� �������� (����� �������� ������ ���� �� 1 �� " << num_vertex - 1 << " �����): " << endl;
                    int length = intinput();
                    while (length < 1 || length > num_vertex - 1)
                    {
                        cout << "����� �������� ������ ���� �� 1 �� " << num_vertex - 1 << " �����. ��������� ����." << endl;
                        length = intinput();
                    }

                    cout << "��������: 1 - ����� �������� ������������ �����; 2 - ����� �������� ����������� �����" << endl;
                    int parametr = intinput(); //���� 1 - ���� max, ���� 2 - ���� min
                    while (parametr < 1 || parametr > 2)
                    {
                        cout << "������������ ����. ��������: 1 - ����� �������� ������������ �����; 2 - ����� �������� ����������� �����" << endl;
                        parametr = intinput();
                    }
                    vector<vector<int>> zhimbel_matrix(weight_matrix);
                    for (int i = 0; i < length - 1; i++) {
                        zhimbel_matrix = Zhimbel(weight_matrix, zhimbel_matrix, num_vertex, parametr);
                    }
                    PrintMatrix(zhimbel_matrix);
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 5: {
                printf("\n�� ������� ������� 5 - ����� ���-�� ��������� �� ����� ������� �� ������\n");
                int from, where_;
                GetStartFinishVertex(from, where_, num_vertex);

                if (res_reachability_matrix[from][where_] != 0) {
                    cout << "\n���-�� ��������� �� ������� " << from << " � ������� " << where_ << ": " << res_reachability_matrix[from][where_] << endl;
                }
                else {
                    cout << "\n����� ������� ��������� ������" << endl;
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 6: {
                printf("\n�� ������� ������� 6 -  �������� ��������\n");
                int from, where_;
                GetStartFinishVertex(from, where_, num_vertex);

                vector<int> path;
                vector<int> dist = Dijkstra(weight_matrix, num_vertex, from, where_, path, iter_Dijkstra);
                cout << "\n���������� �������� ���������: " << iter_Dijkstra << endl;
                PrintWayInfo(num_vertex, dist, from, where_, path);
                cout << "\n";
                printmenu();
                break;
            }
            case 7: {
                printf("\n�� ������� ������� 7 -  �������� ������ � ������ BFS\n");
                int from, where_;
                GetStartFinishVertex(from, where_, num_vertex);

                vector<int> path;
                vector<int> dist = BFS(weight_matrix, num_vertex, from, where_, path, iter_BFS);
                cout << "\n���������� �������� ���������: " << iter_BFS << endl;
                iter_div = float(iter_BFS) / float(iter_Dijkstra);
                cout << "\n�������� �������� ������� ��������� ������ � ������ �: " << iter_div << " ���(�)." << endl;
                PrintWayInfo(num_vertex, dist, from, where_, path);
                cout << "\n";
                printmenu();
                break;
            }
            case 8: {
                printf("\n�� ������� ������� 8 -  �������� ������ ������������� ����\n");
                int from, where_;
                GetStartFinishVertex(from, where_, num_vertex);

                vector<int> path;
                vector<int> dist = MaxWay(weight_matrix, num_vertex, from, where_, path);
                PrintWayInfo(num_vertex, dist, from, where_, path);
                cout << "\n";
                printmenu();
                break;
            }
            case 9: {
                printf("\n�� ������� ������� 9 - ������� ������� ���������� ������������\n");
                PrintMatrix(bandwidth_matrix);
                cout << "\n";
                printmenu();
                break;
            }
            case 10: {
                printf("\n�� ������� ������� 10 - ������� ������� ���������\n");
                PrintMatrix(cost_matrix);
                cout << "\n";
                printmenu();
                break;
            }
            case 11: {
                printf("\n�� ������� ������� 11 - ����� ������������ ����� �� ��������� �����-����������\n");
                int from = 0;
                int where_ = num_vertex - 1;
                int max_flow = Ford_Fulkerson(bandwidth_matrix, from, where_, num_vertex);
                cout << "������������ �����: " << max_flow << endl;
                cout << "\n";
                printmenu();
                break;
            }
            case 12: {
                printf("\n�� ������� ������� 12 - ����� �������� ����� ����������� ���������\n");
                int from = 0;
                int where_ = num_vertex - 1;
                vector<int> path;
                int max_flow = Ford_Fulkerson(bandwidth_matrix, from, where_, num_vertex);
                int need_flow = (max_flow * 2) / 3;
                if (need_flow == 0) {
                    cout << "�������� ����� ����� 0, ����������� ��������� ������ ��������� ������" << endl;
                }
                else {
                    cout << "������������ �����: " << max_flow << endl;
                    cout << "�������� ����� [2/3*max]: " << need_flow << endl;
                    int cost = mincost_maxflow(cost_matrix, bandwidth_matrix, from, where_, num_vertex, path, need_flow);
                    cout << "����������� ��������� ������: " << cost << endl;
                }

                cout << "\n";
                printmenu();
                break;
            }
            case 13: {
                printf("\n�� ������� ������� 13 - ������� ������� ��������\n");
                PrintMatrix(kirhgof_matrix);
                cout << "\n";
                printmenu();
                break;
            }
            case 14: {
                printf("\n�� ������� ������� 14 - ����� ����� �������� ��������\n");
                vector<vector<int>> new_kirhgof_matrix(num_vertex - 1, vector<int>(num_vertex - 1, 0));
                for (int i = 0; i < num_vertex - 1; i++) {
                    for (int j = 0; j < num_vertex - 1; j++) {
                        new_kirhgof_matrix[i][j] = kirhgof_matrix[i + 1][j + 1];
                    }
                }
                if (num_vertex == 2) {
                    cout << "����� �������� ��������: " << 1 << "\n";
                }
                else {
                    cout << "����� �������� ��������: " << determinant(num_vertex - 1, new_kirhgof_matrix) << "\n";
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 15: {
                printf("\n�� ������� ������� 15 - ��������� ����������� �� ���� ����� �� ��������� ��������\n");
                int mst_weight = 0;
                int iter_kruskal = 0;
                vector<edge> kr = Kruskal(undirected_weight_matrix, num_vertex, mst_weight, iter_kruskal);
                cout << "������� ����� ������������������ �����: " << endl;
                PrintMatrix(undirected_weight_matrix);
                cout << "\n���������� �������� ���������: " << iter_kruskal << endl;
                cout << "����������� �� ���� �����: " << mst_weight << endl;
                for (int i = 0; i < kr.size(); i++) {
                    cout << kr[i].from << " <-> " << kr[i].where_ << " ���: " << "" << kr[i].cost << endl;
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 16: {
                printf("\n�� ������� ������� 16 - ��������� ����������� �� ���� ����� �� ��������� �������\n");
                int mst_weight = 0;
                int iter_boruvka = 0;
                vector<edge> kr = Boruvka(undirected_weight_matrix, num_vertex, mst_weight, iter_boruvka);
                cout << "������� ����� ������������������ �����: " << endl;
                PrintMatrix(undirected_weight_matrix);
                cout << "\n���������� �������� ���������: " << iter_boruvka << endl;
                cout << "����������� ��� ������: " << mst_weight << endl;
                for (int i = 0; i < kr.size(); i++) {
                    cout << kr[i].from << " <-> " << kr[i].where_ << " ���: " << "" << kr[i].cost << endl;
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 17: {
                printf("\n�� ������� ������� 17 - ������������ � ������������ ����� � ������� ���� �������\n");
                int mst_weight = 0;
                int iter_kruskal = 0;
                vector<edge> kr = Kruskal(undirected_weight_matrix, num_vertex, mst_weight, iter_kruskal);
                vector<pair<int, int>> pr = codePrufer(kr, num_vertex);
                cout << "��� �������: " << endl;
                for (int i = 0; i < pr.size(); i++) {
                    cout << "����� �������: " << pr[i].second << " ��� �����: " << pr[i].first << endl;
                }
                vector<vector<int>> dec_pr = decodePrufer(pr);
                cout << "\n�������������� ��� �������: " << endl;
                PrintMatrix(dec_pr);


                cout << "\n";
                printmenu();
                break;
            }
            case 18: {
                printf("\n�� ������� ������� 18 - ���������: �������� �� ���� ���������. ���� ���, �� ��������������. ��������� ������� ����.\n");
                vector<vector<int>> eulergraph = MakeEulerGraph(num_vertex, undirected_weight_matrix, undirected_degrees);
                if (num_vertex <= 2) {
                    cout << "���� �� ����� ���� ���������, ��� ��� ������� �� 2-� ������" << endl;
                }
                else {
                    FindEulerCycle(num_vertex, eulergraph);
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 19: {
                printf("\n�� ������� ������� 19 - ���������: �������� �� ���� �������������. ���� ���, �� ��������������. ��������� ������������ �����, �������� � ����. ������ ������������.\n");
                vector<vector<int>> gamiltongraph(num_vertex, vector<int>(num_vertex, 0));
                if (num_vertex <= 2) {
                    cout << "���� �� ����� ���� �������������, ��� ��� ������� �� 2-� ������" << endl;
                }
                else {
                    gamiltongraph = MakeGamiltonGraph(num_vertex, undirected_weight_matrix);
                    TransportSalesmanProblem(gamiltongraph, num_vertex, fout);
                }
                cout << "\n";
                printmenu();
                break;
            }

            default:
                printf("\n������������ ����. ���������� ��� ���.\n");
                break;
            }
        }
    }
}

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "����� ����������!" << endl;
    bool exitflag = true;
    menu(exitflag);
}

