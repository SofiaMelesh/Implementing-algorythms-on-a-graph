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
//Функция для генерации степеней вершин в соответсвии с распределением Зета
void GenerateVertexDegrees(vector<int>& degrees, int num_vertex) {
    float s = 2.0; //влияет на скорость убывания вероятностей для частоты встречаемости событий; должен быть >1 
    vector<double> probability(num_vertex); //вектор вероятностей появления каждой возможной степени вершин от 1 до (макс кол-во вершин-1)
    //создаем объект класса, используемый для получения случайных чисел
    static random_device rd; //начальное состояние
    static mt19937 gen(rd());
    //заполняем вектор вероятностей с помощью формулы Зета
    for (int k = 1; k <= num_vertex; ++k)
    {
        probability[k - 1] = 1.0 / pow(k, s); //k^(-s) = 1/(k^s) -1 делаем так как индексация в probability должна начаться с 0, а не с 1 - наим возможная степень вершины
    }
    //создаем дискретное распределение с этими вероятностями (нормализуем, чтобы сумма была равна 1 с помощью встроенного класса)
    discrete_distribution<> distribution(probability.begin(), probability.end());
    //После создания distribution, его можно использовать вместе с генератором случайных чисел (объект gen класса std::mt19937) 
    //для получения случайных значений, где вероятность каждого значения пропорциональна его вероятности в распределении.
    bool flag = true;
    while (flag) {
        for (int i = 0; i < num_vertex; ++i) {
            degrees[i] = distribution(gen) + 1; //генерируем случайное число с использованием распределения вероятностей (каждая вероятность по индексу=случайному числу смотрится)
        }
        sort(degrees.begin(), degrees.end()); //сортируем по возрастанию
        reverse(degrees.begin(), degrees.end()); //сортируем по убыванию
        degrees[num_vertex - 1] = 0;
        int pr_sum = 0;
        int sumVertexDegrees = 0;
        int max_deg = 0;
        for (int i = 0; i < num_vertex; i++) {
            if (degrees[i] <= num_vertex - i - 1) {
                pr_sum += 1;
                sumVertexDegrees += degrees[i];
            }
            //проверяем кол-во вершин с макс степенью, если такая вершина одна => исток один
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

//формируем матрицу смежности ориентированного связного ациклического графа
//передаем кол-во вершин, вектор степеней всех вершин
//возвращаем матрицу смежности
vector<vector<int>> Generate_adjacency_matrix(int num_vertex, const vector<int>& degrees) {
    //создаем матрицу смежности - матрица размером кол-во вершин*кол-во вершин, заполненная нулями
    vector<vector<int>> adjacency_matrix(num_vertex, vector<int>(num_vertex, 0));
    //создаем копию вектора степеней вершин
    vector<int> copy_degrees(degrees);
    srand(time(0));
    for (int i = 0; i < num_vertex; i++) {
        for (int d = 0; d < copy_degrees[i]; d++) {
            bool n_addded_edge = true;
            while (n_addded_edge) {
                double R = static_cast<double>(rand()) / RAND_MAX; //число в диапазоне [0;1)
                int j = static_cast<int>((num_vertex)*R);
                //проверяем: что связь между вершинами еще не установлена, что не создаем петлю, что степень вершины ещё больше нуля
                if (adjacency_matrix[i][j] != 1 && i != j && copy_degrees[i] > 0 && j > i) {
                    adjacency_matrix[i][j] = 1;
                    n_addded_edge = false;
                }
            }
        }
    }
    return adjacency_matrix;
}

//алгоритм Шимбелла
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
            sort(vect.begin(), vect.end()); //отсортировали по возрастанию
            //если весь вектор из нулей ==> выбираем ноль
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


//формирование матрицы достижимости
vector<vector<int>> Generate_reachability_matrix(vector<vector<int>> adjacency_matrix, vector<vector<int>> matrix, const int num_vertex) {
    vector<vector<int>> new_matrix(num_vertex, vector<int>(num_vertex, 0));
    //выполняем умножение
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

//сложение матриц
vector<vector<int>> Addition_matrix(vector<vector<int>> matrix1, vector<vector<int>> matrix2, const int num_vertex) {
    vector<vector<int>> new_matrix(num_vertex, vector<int>(num_vertex, 0));
    for (int i = 0; i < num_vertex; i++) {
        for (int j = 0; j < num_vertex; j++) {
            new_matrix[i][j] = (matrix1[i][j] + matrix2[i][j]);
        }
    }
    return new_matrix;
}

//вывести матрицу
void PrintMatrix(vector<vector<int>> matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            cout << setw(5) << val << setw(5);
        }
        cout << "\n";
    }
    cout << "\n";
}

//получить маршрут
vector<int> GetPath(vector<int> from, int finish) {
    vector<int> path;
    //начинаем с конца, если v!=-1, то после выполнения тела цикла v = from[v]
    for (int v = finish; v != -1; v = from[v]) {
        path.push_back(v);
    }
    reverse(path.begin(), path.end());
    return path;
}

//алгортим Дейкстры
vector<int> Dijkstra(vector<vector<int>> weight_matrix, int num_vertex, int start, int finish, vector<int>& path, int& iter_Dijkstra) {
    vector<int> from(num_vertex, -1);
    vector<int> distance(num_vertex, INF); //вектор расстояний
    distance[start] = 0; //до всех вершин кроме стартовой ставим расстояние беск, для стартовой ноль
    set<pair<int, int>> unvisited_vertex; //автоматически сортирует добавляемые элементы в порядке возрастания (поэтому первым элементом пары делаем расстояние до вершины, вторым - номер вершины)
    unvisited_vertex.insert(make_pair(distance[start], start));
    while (!unvisited_vertex.empty()) {
        int nearest = unvisited_vertex.begin()->second; //номер ближайшей вершины
        unvisited_vertex.erase(unvisited_vertex.begin());
        //на каждой итерации выясняем расстояние до одной из вершин
        for (int j = 0; j < num_vertex; j++) {
            iter_Dijkstra++;
            if (distance[j] > distance[nearest] + weight_matrix[nearest][j] && weight_matrix[nearest][j] != 0) {
                unvisited_vertex.erase(make_pair(distance[j], j)); //старое значение длины для вершины удаляем
                distance[j] = distance[nearest] + weight_matrix[nearest][j]; //формируем новое
                from[j] = nearest;
                unvisited_vertex.insert(make_pair(distance[j], j)); //записываем в непосещенные
            }
        }
    }
    path = GetPath(from, finish);
    return distance;
}

//алгоритм поиска в ширину BFS
vector<int> BFS(vector<vector<int>> weight_matrix, int num_vertex, int start, int finish, vector<int>& path, int& iter_BFS) {
    vector<int> from(num_vertex, -1);
    vector<int> distance(num_vertex, INF); //массив расстояний до вершин, изначально до всех вершин беск, кроме стартовой вершины
    queue<int> q; //очередь вершин, ожидающих обработку
    distance[start] = 0;
    q.push(start);
    //пока очередь не пуста
    while (!q.empty()) {
        int vert = q.front(); //извлекаем номер вершины, которую будет обрабатывать
        q.pop(); //достаем эту вершину из очереди
        for (int to = 0; to < num_vertex; to++) {
            iter_BFS++;
            //проверяем, что вершина смежна с выбранной и что расстояние до вершины больше, чем расстояние до предыдущей+вес пути
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

//алгортим поиска максимального пути 
vector<int> MaxWay(vector<vector<int>> weight_matrix, int num_vertex, int start, int finish, vector<int>& path) {
    vector<int> from(num_vertex, -1);
    vector<int> distance(num_vertex, INF); //вектор расстояний
    distance[start] = 0; //до всех вершин кроме стартовой ставим расстояние беск, для стартовой ноль
    set<pair<int, int>> unvisited_vertex; //автоматически сортирует добавляемые элементы в порядке возрастания (поэтому первым элементом пары делаем расстояние до вершины, вторым - номер вершины)
    unvisited_vertex.insert(make_pair(distance[start], start));
    while (!unvisited_vertex.empty()) {
        int nearest = unvisited_vertex.begin()->second; //номер ближайшей вершины
        unvisited_vertex.erase(unvisited_vertex.begin());
        //на каждой итерации выясняем расстояние до одной из вершин
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
    //Обновляем вектор предшественников, заполняем -1
    fill(parent.begin(), parent.end(), -1);
    //Предшественник стартовой вершины = -2
    parent[start] = -2;
    //Создаем очередь, хранящую пары: 1ый элемент пары - номер вершины, 2ой - мин пропускаемость
    queue<pair<int, int>> q;
    //Добавляем в очередь стартовую вершину
    q.push({ start,INF });

    //Идем по циклу, пока очередь не пуста
    while (!q.empty()) {
        //Сохраняем верхний узел очереди и мин пропускаемость на данный момент
        int u = q.front().first;
        int cap = q.front().second;
        //Удаляем из очереди этот узел
        q.pop();
        //Проходим по всем вершинам смежным с u
        for (int v = 0; v < num_vertex; v++) {
            //Находим узел,в который мы ещё не заходили (parent[v]==-1), который не равен u, к которому еще можно добраться (пропускаемость != 0)
            if (u != v && bandwidth_matrix[u][v] != 0 && parent[v] == -1) {
                //Записываем u как родителя (предшественника) v
                parent[v] = u;
                //Обновляем минимальную пропускаемость
                int min_cap = min(cap, bandwidth_matrix[u][v]);
                //Если добрались до конечной вершины, то возвращаем найденную минимальную пропускаемость
                if (v == finish) {
                    return min_cap;
                }
                //Если не добрались до конечной вершины, то добавляем вершину v и мин пропускаемость в данный момент в очередь
                q.push({ v,min_cap });
            }
        }
    }
    //Если не нашли пути между стартовой и конечной вершиной, то возвращаем ноль
    return 0;
}

int Ford_Fulkerson(vector<vector<int>> bandwidth_matrix, int start, int finish, int num_vertex) {
    //Инициализируем родительский вектор (parent[1]=0 означает, что из узла 0 можно попасть в 1ый = предшественник 1го узла - нулевой узел) для нахождения пути от стартовой вершины к конечной
    vector<int> parent(num_vertex, -1);
    //Значение максимального потока
    int max_flow = 0;
    //Будем сохранять в переменную значение мин пропускаемости для каждого пути
    int min_cap = 0;

    //Идем по циклу, пока мин пропускаемость не станет нулевой (пока есть путь от стартовой вершины к конечной)
    //Чтобы найти пусть и мин пропускаемость вызываем функцию FF_bfs
    while (min_cap = FF_bfs(bandwidth_matrix, start, finish, num_vertex, parent)) {
        //Добавляем в макс поток мин пропускаемость текущего пути
        max_flow += min_cap;
        //Записываем в v конечную вершину
        int v = finish;

        //Пока не дошли до стартовой вершины, будем идти от конечной до стартовой
        while (v != start) {
            //Записываем родителя вершины v
            int u = parent[v];
            //Уменьшаем пропускаемость из u в v на мин пропускаемость текущего пути
            bandwidth_matrix[u][v] -= min_cap;
            //Увеличиваем (так как в противоположную сторону) пропускаемость из v в u на мин пропускаемость текущего пути
            bandwidth_matrix[v][u] += min_cap;
            //Записываем в v его родителя (продвигаемся всё ближе к стартовой вершине)
            v = u;
        }
    }
    //Возвращаем макс поток
    return max_flow;
}


//модифицированный алгоритм Дейкстры
void modyDijkstra(vector<vector<int>>& weight_matrix, int num_vertex, int start, int finish, vector<int>& parent) {
    vector<int> distance(num_vertex, INF); //вектор расстояний
    distance[start] = 0; //до всех вершин кроме стартовой ставим расстояние беск, для стартовой ноль
    set<pair<int, int>> unvisited_vertex; //автоматически сортирует добавляемые элементы в порядке возрастания (поэтому первым элементом пары делаем расстояние до вершины, вторым - номер вершины)
    unvisited_vertex.insert(make_pair(distance[start], start));
    while (!unvisited_vertex.empty()) {
        int nearest = unvisited_vertex.begin()->second; //номер ближайшей вершины
        unvisited_vertex.erase(unvisited_vertex.begin());
        //на каждой итерации выясняем расстояние до одной из вершин
        for (int j = 0; j < num_vertex; j++) {
            if (distance[j] > distance[nearest] + weight_matrix[nearest][j] && weight_matrix[nearest][j] != 0) {
                unvisited_vertex.erase(make_pair(distance[j], j)); //старое значение длины для вершины удаляем
                distance[j] = distance[nearest] + weight_matrix[nearest][j]; //формируем новое
                parent[j] = nearest;
                unvisited_vertex.insert(make_pair(distance[j], j)); //записываем в непосещенные
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
            delta = min(delta, bandwidth_matrix[parent[v]][v]); //нашли поток по этому пути
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
            bandwidth_matrix[v][parent[v]] += delta; //добавляем дугу в противоположную сторону
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
                std::cout << "Поток, проходящий по ребру " << j << "->" << i << " равен: " << flow_matrix[i][j] << ", стоимость равна: " << abs(cost_matrix[i][j]) << endl;
            }
        }
    }

    return ans;
}

void PrintWayInfo(int num_vertex, vector<int> dist, int from, int where_, vector<int> path) {
    cout << "Вектор расстояний: " << endl;
    for (int i = 0; i < num_vertex; i++) {
        if (dist[i] != INF and dist[i] != -INF) {
            cout << dist[i] << "\t";
        }
        else {
            cout << "-" << "\t";
        }

    }
    if (dist[where_] != INF and dist[where_] != -INF) {
        cout << "\nДлина маршрута: " << dist[where_] << ".";
        cout << "\nМаршрут: ";
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
        cout << "\nТакой маршрут построить нельзя" << endl;
    }
}



//для проверки на int 
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
            printf("\nВы ввели некоректное значение. Повторите ввод.\n");
        }
        else break;
    }
    return val;
}

void GetStartFinishVertex(int& start, int& finish, int num_vertex) {
    cout << "Введите номер вершины из которой хотите построить маршрут от 0 до " << num_vertex - 1 << ":" << endl;
    start = intinput(); //из какой вершины строим маршрут
    while (start < 0 || start > num_vertex - 1)
    {
        cout << "Неправильный ввод. Введите номер вершины из которой хотите построить маршрут от 0 до " << num_vertex - 1 << ":" << endl;
        start = intinput();
    }

    cout << "Введите номер вершины в которую хотите построить маршрут от 0 до " << num_vertex - 1 << ":" << endl;
    finish = intinput(); //в какую вершину хотим добраться
    while (finish < 0 || finish > num_vertex - 1)
    {
        cout << "Неправильный ввод. Введите номер вершины в которую хотите построить маршрут от 0 до " << num_vertex - 1 << ":" << endl;
        finish = intinput();
    }

}



//определитель матрицы
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

//Краскала
//структура для хранения ребра
struct edge {
    int cost;
    int from;
    int where_;
};

//компаратор, чтобы расставить в порядке убывания
struct compare {
    bool operator() (edge const& a, edge const& b) const {
        return a.cost > b.cost;
    }
};

//компаратор, чтобы расставить в порядке убывания
struct comparenum {
    bool operator() (edge const& a, edge const& b) const {
        return a.from < b.from;
    }
};

//найти корень множества вершин, к которому принадлежит вершина
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
    //формируем вектор edges, хранящий тройку <вес,откуда ребро, куда ребро>
    for (int i = 0; i < num_vertex; i++) {
        for (int j = 0; j < num_vertex; j++) {
            if (weight_matrix[i][j] != 0) {
                edges.push_back({ weight_matrix[i][j],i,j });
            }
        }
    }
    //сортируем по убыванию веса ребер
    sort(edges.begin(), edges.end(), compare());

    vector<int> parent(num_vertex);
    for (int i = 0; i < num_vertex; i++) {
        parent[i] = i;
    }

    vector<edge> mst;
    while (mst.size() != num_vertex - 1) {
        edge next_edge = edges.back();
        edges.pop_back();

        //ищем корень множества вершин первой вершины ребра и корень множества вершин второй вершины ребра
        int f_p = find_root(next_edge.from, parent);
        int s_p = find_root(next_edge.where_, parent);

        //если корни не совпали, то добавление ребра не создаст цикл
        if (f_p != s_p) {
            mst.push_back(next_edge);
            mst_weight += next_edge.cost;
            union_(f_p, s_p, parent);
        }
        iter_kruskal++;
    }
    return mst;
}

//Борувка
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
    //формируем вектор edges, хранящий тройку <вес,откуда ребро, куда ребро>
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
    //пара: вес ребра, номер родителя
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
    //выводим матрицу стоимости
    cout << "\nМатрица весов: " << endl;
    PrintMatrix(eulergraph);
    //вывод всех степеней вершин
    for (int i = 0; i < num_vertex; i++) {
        cout << "Степень вершины №" << i << " = " << undirected_degrees[i] << endl;
    }
    //проверяем, что все степени вершин графа четные => эйлеров граф
    for (int i = 0; i < num_vertex; i++) {
        if (undirected_degrees[i] % 2) {
            countevenvert++;
        }
    }
    if (countevenvert == num_vertex) {
        cout << "Граф является Эйлеровым\n" << endl;
        flagEuler = true;
    }
    //если не все степени вершин графа четные => надо модифицировать граф - сделать эйлеровым
    else {
        cout << "\nГраф не является Эйлеровым. Модифицируем граф.\n" << endl;
        //будем удалять или добавлять ребра пока граф не станет эйлеровым
        while (!flagEuler) {
            bool flag = 0;
            for (int i = 0; i < num_vertex; i++) {
                for (int j = 0; j < num_vertex; j++) {
                    //добавляем новое ребро, если у вершины нечетная степень и такого ребра еще нет
                    if (undirected_degrees[i] % 2 != 0 and undirected_degrees[j] % 2 != 0 and eulergraph[i][j] == 0 and i != j) {
                        flag = true;
                        int r = rand() % (10 - (1) + 1) + 1;
                        eulergraph[i][j] = r;
                        eulergraph[j][i] = r;
                        undirected_degrees[i]++;
                        undirected_degrees[j]++;
                        cout << "Добавили ребро: " << i << "<->" << j << " с весом: " << r << endl;
                        break;
                    }
                }
                //если добавили ребро, то можно проверить степени вершин на четность и выйти из цикла добавления вершины
                if (flag) {
                    break;
                }
            }
            //если не получилось добавить ребро, то пробуем удалить
            if (!flag) {
                for (int i = 0; i < num_vertex; i++) {
                    for (int j = 0; j < num_vertex; j++) {
                        if (undirected_degrees[i] % 2 != 0 and undirected_degrees[j] % 2 != 0 and eulergraph[i][j] != 0 and i != j) {
                            flag = true;
                            eulergraph[i][j] = 0;
                            eulergraph[j][i] = 0;
                            undirected_degrees[i]--;
                            undirected_degrees[j]--;
                            cout << "Удалили ребро: " << i << "<->" << j << endl;
                            break;
                        }
                    }
                    //если удалили ребро, то можно проверить степени вершин на четность и выйти из цикла удаления вершины
                    if (flag) {
                        break;
                    }
                }
            }
            //проверяем получилось ли сделать граф с четными степенями вершин
            int countevenvert1 = 0;
            for (int i = 0; i < num_vertex; i++) {
                if (undirected_degrees[i] % 2 == 0) {
                    countevenvert1++;
                }
            }
            //если получилось, то поднимаем флаг и выходим из цикла, если нет - продолжаем добавлять/удалять ребро
            if (countevenvert1 == num_vertex) {
                flagEuler = true;
            }
        }
        cout << "Теперь граф является Эйлеровым.\n" << endl;
        //выводим матрицу стоимости
        cout << "Матрица весов: " << endl;
        PrintMatrix(eulergraph);
        //вывод всех степеней вершин
        for (int i = 0; i < num_vertex; i++) {
            cout << "Степень вершины №" << i << " = " << undirected_degrees[i] << endl;
        }
    }
    return eulergraph;
}

void FindEulerCycle(int num_vertex, vector<vector<int>> eulergraph) {
    stack<int> vertexes;
    vector<int> cycle;
    int vert = 0;
    //ищем первую вершину, из которой есть ребро
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
    //добавляем эту вершину в стек
    vertexes.push(vert);
    //пока стек не пустой
    while (!vertexes.empty()) {
        //берем верхнюю вершину
        vert = vertexes.top();
        int i;
        //находим первую смежную с ней
        for (i = 0; i < num_vertex; i++) {
            if (eulergraph[vert][i] != 0) {
                break;
            }
        }
        //если смежной не нашлось, то добавляем эту вершину в цикл, удаляем из стека
        if (i == num_vertex) {
            cycle.push_back(vert);
            vertexes.pop();
        }
        //если нашли смежную, то добавляем в стек и удаляем ребро между смежными вершинами в матрице
        else {
            vertexes.push(i);
            eulergraph[vert][i] = 0;
            eulergraph[i][vert] = 0;
        }
    }
    cout << "\nЭйлеров цикл:" << endl;
    for (int i = 0; i < cycle.size(); i++) {
        cout << cycle[i];
        if (i != cycle.size() - 1) {
            cout << " -> ";
        }
    }
}


//вспомогательная функция для проверки возможности добавления вершины в путь
bool CanGoNext(int v, vector<vector<int>>& gamiltongraph, vector<int>& path, int position) {
    //проверяем, есть ли ребро между последней вершиной path[position - 1] и вершиной v в матрице смежности
    //если нет, то возвращаем false => вершина не мб добавлена в путь
    if (gamiltongraph[path[position - 1]][v] == 0) {
        return false;
    }
    //затем проверяем не была ли вершина уже посещена
    //если была, то возвращаем false => вершина не мб добавлена в путь
    for (int i = 0; i < position; i++) {
        if (path[i] == v) {
            return false;
        }
    }
    return true;
}


//функция для поиска гамильтонова цикла (рекурсивный)
bool FindGamiltonCycle(vector<vector<int>>& gamiltongraph, vector<int>& path, int position) {
    //если все вершины посещены и последняя соединена с первой, то возвращается true => есть гамильтонов цикл в графе
    if (position == gamiltongraph.size()) {
        return gamiltongraph[path[position - 1]][path[0]] != 0;
    }
    //иначе ищем следующую вершину для посещения
    for (int v = 1; v < gamiltongraph.size(); v++) {
        //в которую есть возможность дойти (есть исходящее ребро + не была посещена смежная вершина)
        if (CanGoNext(v, gamiltongraph, path, position)) {
            //если найлена, то добавляем эту вершину в путь
            path[position] = v;
            //если гамильтонов цикл найден, то возвращаем true
            if (FindGamiltonCycle(gamiltongraph, path, position + 1)) {
                return true;
            }
            //если гамильтонов цикл не найден, то удаляем вершину v из пути
            path[position] = -1;
        }
    }
    return false;
}


//функция для проверки наличия гамильтонова цикла (рекурсивный)
bool ExistGamiltonCycle(vector<vector<int>>& gamiltongraph) {
    vector<int> path(gamiltongraph.size(), -1);
    path[0] = 0;
    if (!FindGamiltonCycle(gamiltongraph, path, 1)) {
        return false;
    }
    return true;
}


//гамильтонов граф: проверка, модификация
vector<vector<int>> MakeGamiltonGraph(int num_vertex, vector<vector<int>> undirected_weight_matrix) {
    vector<vector<int>> gamiltongraph(undirected_weight_matrix);
    //выводим матрицу стоимости
    cout << "\nМатрица весов: " << endl;
    PrintMatrix(gamiltongraph);
    if (ExistGamiltonCycle(gamiltongraph)) {
        cout << "Граф является гамильтоновым\n" << endl;
    }
    else {
        cout << "\nГраф не является гамильтоновым. Модифицируем граф." << endl;
        vector<int> vect(num_vertex, -1);
        //Будем добавлять ребра, пока граф не станет гамильтоновым
        for (int i = 0; i < num_vertex; i++) {
            int r = rand() % ((num_vertex - 1) - 0 + 1) + 0;
            //Функция find будет последовательно идти по массиву с начала до конца, 
            //пока не найдет элемент, который равен значению для поиска r. 
            //Если значение r не будет найдено, то find вернет указатель на конец массива => добавим это значение
            if (find(vect.begin(), vect.end(), r) == vect.end()) {
                vect[i] = r;
            }
            else {
                i--;
            }
        }
        //проходим по всем элементам вектора и проверяем наличие ребра между соответсвующими вершинами
        //если такого ребра нет - добавляем случайное
        for (int i = 0; i < num_vertex; i++) {
            if (gamiltongraph[vect[i]][vect[(i + 1) % num_vertex]] == 0) {
                int r = rand() % (10 - (1) + 1) + 1;
                gamiltongraph[vect[i]][vect[(i + 1) % num_vertex]] = r;
                gamiltongraph[vect[(i + 1) % num_vertex]][vect[i]] = r;
                cout << "Добавили ребро: " << vect[i] << "<->" << vect[(i + 1) % num_vertex] << " с весом: " << r << endl;
            }
        }
        if (ExistGamiltonCycle(gamiltongraph)) {
            cout << "Теперь граф является гамильтоновым." << endl;
            //выводим матрицу стоимости
            cout << "Матрица весов: " << endl;
            PrintMatrix(gamiltongraph);
        }
        else {
            cout << "Граф не удалось модифицировать и он не является гамильтоновым\n" << endl;
        }
    }
    return gamiltongraph;
}

void findGamiltonCycleforTSP(vector<vector<int>>& gamiltongraph, vector<bool>& visited, vector<int>& path, int vert, int num_vertex, int count, int cost, int& min_cost, vector<int>& min_path, ofstream& fout) {
    //добавляем в путь текущую вершину
    path.push_back(vert);
    //если все вершины графа посещены и последняя вершина соединена с нулевой,
    //то записываем путь и вес цикла в файл
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
        fout << "\tВес цикла: " << gamiltongraph[vert][0] + cost << endl;
        //если вес цикла меньше минимального, то минимальную стоимость и путь обновляем
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
    vector<bool> visited(num_vertex); //посещенность вершины
    vector<int> path; //текущий путь
    visited[0] = true;
    int min_cost = INF; //минимальная стоимость пути
    vector<int> min_path; //минимальный путь
    findGamiltonCycleforTSP(gamiltongraph, visited, path, 0, num_vertex, 1, 0, min_cost, min_path, fout);
    cout << "Цикл коммивояжера: " << endl;
    for (int i = 0; i < min_path.size(); i++) {
        cout << min_path[i] << " -> ";
    }
    cout << min_path[0];
    cout << "\nВес: " << min_cost << endl;
    fout.close();
}




//вывод меню
void printmenu()
{
    printf("\nВведите цифру от 0 до 18:\n");
    printf("0 - Завершить работу программы\n");
    printf("1 - Сгенерировать новый граф\n");
    printf("\nЛабораторная №1\n");
    printf("2 - Вывести матрицу смежности вершин\n");
    printf("3 - Вывести матрицу весов\n");
    printf("4 - Метод Шимбелла\n");
    printf("5 - Найти кол-во маршрутов от одной вершины до другой\n");
    printf("\nЛабораторная №2\n");
    printf("6 - Алгоритм Дейкстры\n");
    printf("7 - Алгоритм поиска в ширину BFS\n");
    printf("8 - Алгоритм поиска максимального пути\n");
    printf("\nЛабораторная №3\n");
    printf("9 - Вывести матрицу пропускных способностей\n");
    printf("10 - Вывести матрицу стоимости\n");
    printf("11 - Найти максимальный поток по алгоритму Форда-Фалкерсона\n");
    printf("12 - Найти заданный поток минимальной стоимости\n");
    printf("\nЛабораторная №4\n");
    printf("13 - Вывести Матрицу Кирхгофа\n");
    printf("14 - Найти число остовных деревьев\n");
    printf("15 - Построить минимальный по весу остов по алгоритму Краскала\n");
    printf("16 - Построить минимальный по весу остов по алгоритму Борувки\n");
    printf("17 - Закодировать и декодировать остов с помощью кода Прюфера\n");
    printf("\nЛабораторная №5\n");
    printf("18 - Проверить: является ли граф эйлеровым. Если нет, то модифицировать. Построить эйлеров цикл.\n");
    printf("19 - Проверить: является ли граф гамильтоновым. Если нет, то модифицировать. Построить гамильтоновы циклы, записать в файл. Задача коммивояжера.\n\n");
}

//функция меню
void menu(bool& exitflag) {
    while (exitflag) {
        cout << "Введите кол-во вершин в графе:" << endl;
        int num_vertex = intinput(); //кол-во вершин - задает пользователь
        while (num_vertex < 2 || num_vertex > 100)
        {
            printf("\nКол-во вершин графа может быть от 2 до 100. Повторите ввод.\n");
            num_vertex = intinput();
        }
        int num_edge = num_vertex - 1; //кол-во ребер в связном ациклическом графе = кол-во вершин - 1
        vector<int> degrees(num_vertex); //вектор степеней вершин  
        int res_sumVertexDegrees = 0;
        GenerateVertexDegrees(degrees, num_vertex); //генерируем степени всех вершин

        //вывод всех степеней вершин
        for (int i = 0; i < num_vertex; i++) {
            cout << "Степень вершины №" << i << " = " << degrees[i] << endl;
        }

        //создаем матрицу смежности
        vector<vector<int>> adjacency_matrix = Generate_adjacency_matrix(num_vertex, degrees);

        //создаем матрицу весов
        cout << "\nВыберите: 1 - граф с положительными весами, 2 - граф с положительными и отрицательными весами:" << endl;
        int sign = intinput(); //кол-во вершин - задает пользователь
        while (sign < 1 || sign > 2)
        {
            cout << "Неправильный ввод. Выберите: 1 - граф с положительными весами, 2 - граф с положительными и отрицательными весами:" << endl;
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

        vector<vector<int>> ed_matrix(num_vertex, vector<int>(num_vertex, 0)); //единичная матрица
        //формируем единичную матрицу
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (i == j) {
                    ed_matrix[i][j] = 1;
                }
            }
        }

        //формируем матрицу достижимости
        vector<vector<int>> reachability_matrix(adjacency_matrix);
        vector<vector<int>> res_reachability_matrix(num_vertex, vector<int>(num_vertex, 0));
        res_reachability_matrix = Addition_matrix(ed_matrix, adjacency_matrix, num_vertex);
        for (int i = 0; i < num_vertex; i++) {
            reachability_matrix = Generate_reachability_matrix(adjacency_matrix, reachability_matrix, num_vertex);
            res_reachability_matrix = Addition_matrix(res_reachability_matrix, reachability_matrix, num_vertex);
        }

        //формируем матрицу пропускных способностей (положит числа)
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

        //формируем матрицу стоимости (положит и отриц числа)
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


        //формируем матрицу смежности неориентированного графа
        vector<vector<int>> undirected_adjacency_matrix(num_vertex, vector<int>(num_vertex, 0));
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (adjacency_matrix[i][j] != 0) {
                    undirected_adjacency_matrix[j][i] = adjacency_matrix[i][j];
                    undirected_adjacency_matrix[i][j] = adjacency_matrix[i][j];
                }
            }
        }

        //формируем матрицу со степенями вершин на главной диагонали для неориентированного графа
        vector<vector<int>> degree_matrix(num_vertex, vector<int>(num_vertex, 0));
        for (int i = 0; i < num_vertex; i++) {
            //int d = 0;
            for (int j = 0; j < num_vertex; j++) {
                if (undirected_adjacency_matrix[i][j] != 0) {
                    degree_matrix[i][i] += 1;
                }
            }

        }

        //формируем вектор степеней вершин неориентированного графа
        vector<int> undirected_degrees(num_vertex, 0);
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (degree_matrix[i][j] != 0) {
                    undirected_degrees[i] = degree_matrix[i][j];
                }
            }
        }

        //формируем матрицу Кирхгофа
        vector<vector<int>> kirhgof_matrix(num_vertex, vector<int>(num_vertex, 0));
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                kirhgof_matrix[i][j] = degree_matrix[i][j] - undirected_adjacency_matrix[i][j];
            }
        }


        //формируем матрицу весов неориентированного графа
        vector<vector<int>> undirected_weight_matrix(num_vertex, vector<int>(num_vertex, 0));
        for (int i = 0; i < num_vertex; i++) {
            for (int j = 0; j < num_vertex; j++) {
                if (cost_matrix[i][j] != 0) {
                    undirected_weight_matrix[j][i] = weight_matrix[i][j];
                    undirected_weight_matrix[i][j] = weight_matrix[i][j];
                }
            }
        }



        printmenu(); //выводим меню
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
                printf("\nВы выбрали функцию 0 - завершение работы программы\n");
                flag = false;
                exitflag = false;
                break;
            }
            case 1: {
                printf("\nВы выбрали функцию 1 - Сгенерировать новый граф\n");
                flag = false;
                break;
            }
            case 2: {
                printf("\nВы выбрали функцию 2 - Вывести матрицу смежности вершин\n");
                PrintMatrix(adjacency_matrix);
                cout << "\n";
                printmenu();
                break;
            }
            case 3: {
                printf("\nВы выбрали функцию 3 - Вывести матрицу весов\n");
                PrintMatrix(weight_matrix);
                cout << "\n";
                printmenu();
                break;
            }
            case 4: {
                printf("\nВы выбрали функцию 4 - Метод Шимбелла\n");
                if (num_vertex == 1) {
                    cout << "Нельзя использовать метод Шимбелла при наличии 1 вершины" << endl;
                }
                else {
                    cout << "Введите длину маршрута (длина маршрута должна быть от 1 до " << num_vertex - 1 << " ребер): " << endl;
                    int length = intinput();
                    while (length < 1 || length > num_vertex - 1)
                    {
                        cout << "Длина маршрута должна быть от 1 до " << num_vertex - 1 << " ребер. Повторите ввод." << endl;
                        length = intinput();
                    }

                    cout << "Выберите: 1 - поиск маршрута максимальной длины; 2 - поиск маршрута минимальной длины" << endl;
                    int parametr = intinput(); //если 1 - ищем max, если 2 - ищем min
                    while (parametr < 1 || parametr > 2)
                    {
                        cout << "Неправильный ввод. Выберите: 1 - поиск маршрута максимальной длины; 2 - поиск маршрута минимальной длины" << endl;
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
                printf("\nВы выбрали функцию 5 - Найти кол-во маршрутов от одной вершины до другой\n");
                int from, where_;
                GetStartFinishVertex(from, where_, num_vertex);

                if (res_reachability_matrix[from][where_] != 0) {
                    cout << "\nКол-во маршрутов из вершины " << from << " в вершину " << where_ << ": " << res_reachability_matrix[from][where_] << endl;
                }
                else {
                    cout << "\nТакой маршрут построить нельзя" << endl;
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 6: {
                printf("\nВы выбрали функцию 6 -  Алгоритм Дейкстры\n");
                int from, where_;
                GetStartFinishVertex(from, where_, num_vertex);

                vector<int> path;
                vector<int> dist = Dijkstra(weight_matrix, num_vertex, from, where_, path, iter_Dijkstra);
                cout << "\nКоличество итераций алгоритма: " << iter_Dijkstra << endl;
                PrintWayInfo(num_vertex, dist, from, where_, path);
                cout << "\n";
                printmenu();
                break;
            }
            case 7: {
                printf("\nВы выбрали функцию 7 -  Алгоритм поиска в ширину BFS\n");
                int from, where_;
                GetStartFinishVertex(from, where_, num_vertex);

                vector<int> path;
                vector<int> dist = BFS(weight_matrix, num_vertex, from, where_, path, iter_BFS);
                cout << "\nКоличество итераций алгоритма: " << iter_BFS << endl;
                iter_div = float(iter_BFS) / float(iter_Dijkstra);
                cout << "\nАлгоритм Дейкстры быстрее алгоритма поиска в ширину в: " << iter_div << " раз(а)." << endl;
                PrintWayInfo(num_vertex, dist, from, where_, path);
                cout << "\n";
                printmenu();
                break;
            }
            case 8: {
                printf("\nВы выбрали функцию 8 -  Алгоритм поиска максимального пути\n");
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
                printf("\nВы выбрали функцию 9 - Вывести матрицу пропускных способностей\n");
                PrintMatrix(bandwidth_matrix);
                cout << "\n";
                printmenu();
                break;
            }
            case 10: {
                printf("\nВы выбрали функцию 10 - Вывести матрицу стоимости\n");
                PrintMatrix(cost_matrix);
                cout << "\n";
                printmenu();
                break;
            }
            case 11: {
                printf("\nВы выбрали функцию 11 - Найти максимальный поток по алгоритму Форда-Фалкерсона\n");
                int from = 0;
                int where_ = num_vertex - 1;
                int max_flow = Ford_Fulkerson(bandwidth_matrix, from, where_, num_vertex);
                cout << "Максимальный поток: " << max_flow << endl;
                cout << "\n";
                printmenu();
                break;
            }
            case 12: {
                printf("\nВы выбрали функцию 12 - Найти заданный поток минимальной стоимости\n");
                int from = 0;
                int where_ = num_vertex - 1;
                vector<int> path;
                int max_flow = Ford_Fulkerson(bandwidth_matrix, from, where_, num_vertex);
                int need_flow = (max_flow * 2) / 3;
                if (need_flow == 0) {
                    cout << "Заданный поток равен 0, минимальную стоимость потока вычислить нельзя" << endl;
                }
                else {
                    cout << "Максимальный поток: " << max_flow << endl;
                    cout << "Заданный поток [2/3*max]: " << need_flow << endl;
                    int cost = mincost_maxflow(cost_matrix, bandwidth_matrix, from, where_, num_vertex, path, need_flow);
                    cout << "Минимальная стоимость потока: " << cost << endl;
                }

                cout << "\n";
                printmenu();
                break;
            }
            case 13: {
                printf("\nВы выбрали функцию 13 - Вывести матрицу Кирхгофа\n");
                PrintMatrix(kirhgof_matrix);
                cout << "\n";
                printmenu();
                break;
            }
            case 14: {
                printf("\nВы выбрали функцию 14 - Найти число остовных деревьев\n");
                vector<vector<int>> new_kirhgof_matrix(num_vertex - 1, vector<int>(num_vertex - 1, 0));
                for (int i = 0; i < num_vertex - 1; i++) {
                    for (int j = 0; j < num_vertex - 1; j++) {
                        new_kirhgof_matrix[i][j] = kirhgof_matrix[i + 1][j + 1];
                    }
                }
                if (num_vertex == 2) {
                    cout << "Число остовных деревьев: " << 1 << "\n";
                }
                else {
                    cout << "Число остовных деревьев: " << determinant(num_vertex - 1, new_kirhgof_matrix) << "\n";
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 15: {
                printf("\nВы выбрали функцию 15 - Построить минимальный по весу остов по алгоритму Краскала\n");
                int mst_weight = 0;
                int iter_kruskal = 0;
                vector<edge> kr = Kruskal(undirected_weight_matrix, num_vertex, mst_weight, iter_kruskal);
                cout << "Матрица весов неориентированного графа: " << endl;
                PrintMatrix(undirected_weight_matrix);
                cout << "\nКоличество итераций алгоритма: " << iter_kruskal << endl;
                cout << "Минимальный по весу остов: " << mst_weight << endl;
                for (int i = 0; i < kr.size(); i++) {
                    cout << kr[i].from << " <-> " << kr[i].where_ << " вес: " << "" << kr[i].cost << endl;
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 16: {
                printf("\nВы выбрали функцию 16 - Построить минимальный по весу остов по алгоритму Борувки\n");
                int mst_weight = 0;
                int iter_boruvka = 0;
                vector<edge> kr = Boruvka(undirected_weight_matrix, num_vertex, mst_weight, iter_boruvka);
                cout << "Матрица весов неориентированного графа: " << endl;
                PrintMatrix(undirected_weight_matrix);
                cout << "\nКоличество итераций алгоритма: " << iter_boruvka << endl;
                cout << "Минимальный вес остова: " << mst_weight << endl;
                for (int i = 0; i < kr.size(); i++) {
                    cout << kr[i].from << " <-> " << kr[i].where_ << " вес: " << "" << kr[i].cost << endl;
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 17: {
                printf("\nВы выбрали функцию 17 - Закодировать и декодировать остов с помощью кода Прюфера\n");
                int mst_weight = 0;
                int iter_kruskal = 0;
                vector<edge> kr = Kruskal(undirected_weight_matrix, num_vertex, mst_weight, iter_kruskal);
                vector<pair<int, int>> pr = codePrufer(kr, num_vertex);
                cout << "Код Прюфера: " << endl;
                for (int i = 0; i < pr.size(); i++) {
                    cout << "номер вершины: " << pr[i].second << " вес ребра: " << pr[i].first << endl;
                }
                vector<vector<int>> dec_pr = decodePrufer(pr);
                cout << "\nДекодированный код Прюфера: " << endl;
                PrintMatrix(dec_pr);


                cout << "\n";
                printmenu();
                break;
            }
            case 18: {
                printf("\nВы выбрали функцию 18 - Проверить: является ли граф эйлеровым. Если нет, то модифицировать. Построить эйлеров цикл.\n");
                vector<vector<int>> eulergraph = MakeEulerGraph(num_vertex, undirected_weight_matrix, undirected_degrees);
                if (num_vertex <= 2) {
                    cout << "Граф не может быть эйлеровым, так как состоит их 2-х вершин" << endl;
                }
                else {
                    FindEulerCycle(num_vertex, eulergraph);
                }
                cout << "\n";
                printmenu();
                break;
            }
            case 19: {
                printf("\nВы выбрали функцию 19 - Проверить: является ли граф гамильтоновым. Если нет, то модифицировать. Построить гамильтоновы циклы, записать в файл. Задача коммивояжера.\n");
                vector<vector<int>> gamiltongraph(num_vertex, vector<int>(num_vertex, 0));
                if (num_vertex <= 2) {
                    cout << "Граф не может быть гамильтоновым, так как состоит их 2-х вершин" << endl;
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
                printf("\nНеправильный ввод. Попробуйте ещё раз.\n");
                break;
            }
        }
    }
}

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Добро пожаловать!" << endl;
    bool exitflag = true;
    menu(exitflag);
}

