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
int trash2 = 2147483647; //макс инт?

void erase(vector<int>& vec, int index) {
    vec.erase(remove(vec.begin(), vec.end(), index), vec.end()); //все элементы, равные value, начиная с начала вектора
}

template<typename T>
const T& clamp(const T& value, const T& low, const T& high) {
    return (value < low) ? low : (high < value) ? high : value;
}

int generateParetoValue(double alpha, int xmin, mt19937& gen) {
    uniform_real_distribution<> dis(0.0, 1.0);
    double u = dis(gen);
    int value = static_cast<int>(xmin / pow(u, 1.0 / alpha));

    // Ограничиваем значение диапазоном [1, 100]
    return clamp(value, 1, 100);
}

vector<vector<int>> GenerateWeightMatrix(int ver, const vector<vector<int>>& admatrix, double alpha, int xmin) {
    vector<vector<int>> weightmatrix(ver, vector<int>(ver, 0));
    static random_device rd;
    static mt19937 gen(rd());

    for (int i = 0; i < ver; i++) {
        for (int j = 0; j < ver; j++) {
            // Генерируем значения ТОЛЬКО для существующих рёбер
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
    vector<int> path;  // Будущий путь от начальной вершины до finish
    // Восстанавливаем путь, двигаясь "назад" по вектору from
    // from[v] хранит, из какой вершины мы пришли в вершину v при обходе графа
    for (int v = finish; v != -1; v = from[v]) {
        path.push_back(v);  // Добавляем текущую вершину в путь
    }
    reverse(path.begin(), path.end());
    return path;
}

void PrintWay(int ver, vector<int> dist, int from, int where_, vector<int> path) {
    cout << "Вектор расстояний: " << endl;
    for (int i = 0; i < ver; i++) {
        if (dist[i] != trash1 and dist[i] != trash2 and dist[i] != -trash1) cout << dist[i] << "\t";
        else cout << "-" << "\t";

    }
    if (dist[where_] != trash1 and dist[where_] != -trash1) {
        cout << "\nДлина маршрута: " << dist[where_] << ".";
        cout << "\nМаршрут: ";
        for (int i = 0; i < path.size(); i++) {
            if (i == path.size() - 1) cout << path[i] << ".";
            else cout << path[i] << "->";
        }
    }
    else cout << "\nТакой маршрут построить нельзя" << endl;
}

//определитель матрицы разложение по первой строке
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
    vector<double> probability(ver); //вектор для вероятностей, для нормализации

    static random_device rd; // начальное состояние
    static mt19937 gen(rd());

    double sum_prob = 0.0; //рассчитываем нормализованные вероятности для распределения Парето
    for (int k = xmin; k <= ver; ++k) {
        double prob = pow(k, -(alpha + 1));  // Расчёт для каждого значения степени
        probability[k - 1] = prob;
        sum_prob += prob; //суммируем для нормализации
    }

    for (int k = 0; k < ver; ++k) { //нормализуем вероятности, чтобы сумма была равна 1
        probability[k] /= sum_prob;
    }

    // создаем дискретное распределение с этими вероятностями (нормализуем)
    discrete_distribution<> distribution(probability.begin(), probability.end());

    bool flag = true; // генерируем случайные степени
    while (flag) {
        for (int i = 0; i < ver; ++i) {
            degrees[i] = distribution(gen) + xmin; // Генерируем число, соответствующее распределению Парето
        }
        sort(degrees.begin(), degrees.end()); // сортируем по возрастанию
        reverse(degrees.begin(), degrees.end()); // сортируем по убыванию

        // проверка на допустимость распределения степеней вершин
        degrees[ver - 1] = 0;
        int pr_sum = 0;
        int sumVertexDegrees = 0;
        int max_deg = 0;
        for (int i = 0; i < ver; i++) {
            if (degrees[i] <= ver - i - 1) {
                pr_sum += 1;
                sumVertexDegrees += degrees[i];
            }
            // проверяем кол-во вершин с максимальной степенью, если такая вершина одна => исток один
            else if (degrees[i] == ver - 1) max_deg += 1;
        }

        if (max_deg == 0) {
            degrees[0] = ver - 1;
            max_deg += 1;
        }

        if (pr_sum == ver && sumVertexDegrees <= (ver * ver - ver) / 2 && max_deg == 1) flag = false;
    }
    //вывод всех степеней вершин
    for (int i = 0; i < ver; i++) cout << "Степень вершины №" << i << " = " << degrees[i] << endl;
}

vector<vector<int>> Createadmatrix(int ver, const vector<int>& degrees) {
    //создаем матрицу смежности - матрица размером кол-во вершин*кол-во вершин, заполненная нулями
    vector<vector<int>> admatrix(ver, vector<int>(ver, 0));
    //создаем копию вектора степеней вершин
    vector<int> copy_degrees(degrees);
    srand(time(0));
    for (int i = 0; i < ver; i++) {
        for (int d = 0; d < copy_degrees[i]; d++) {
            bool n_addded_edge = true;
            while (n_addded_edge) {
                double R = static_cast<double>(rand()) / RAND_MAX; //число в диапазоне [0;1)
                int j = static_cast<int>((ver)*R);
                //проверяем: что связь между вершинами еще не установлена, что не создаем петлю, что степень вершины ещё больше нуля
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

//алгоритм обхода вершин в ширину BFS
vector<int> BFS(vector<vector<int>> weightmatrix, int ver, int start, int finish, vector<int>& path, int& iter_BFS) {
    vector<int> from(ver, -1); //чтобы восстановить путь: from[i] = вершина, из которой пришли в i
    iter_BFS = 0;
    vector<int> distance(ver, trash1); //расстояния до вершин, изначально "бесконечность" (trash1), кроме стартовой
    queue<int> q; //очередь для обхода вершин
    distance[start] = 0; //расстояние до стартовой вершины = 0
    q.push(start); //стартовую вершину в очередь

    while (!q.empty()) {
        int vert = q.front(); //извлекаем текущую вершину из очереди
        q.pop();
        for (int to = 0; to < ver; to++) { //проверяем все возможные вершины, в которые можно попасть из текущей
            if (weightmatrix[vert][to] != 0) {  //если есть ребро из vert в to 
                //если найден более короткий путь до вершины to — обновляем его
                if (distance[to] > distance[vert] + weightmatrix[vert][to]) {
                    distance[to] = distance[vert] + weightmatrix[vert][to];  //обновляем кратчайшее расстояние
                    q.push(to); //добавляем вершину в очередь для дальнейшего обхода
                    from[to] = vert; //запоминаем, откуда пришли в вершину to
                }
            }
            iter_BFS++;
        }
    }

    path = PrintPath(from, finish); //восстанавливаем путь от start до finish через вспомогательный вектор from
    return distance; //возвращаем вектор расстояний от start до всех остальных вершин
}

vector<int> BellmanFord(vector<vector<int>>& weightmatrix, int ver, int start, int finish, vector<int>& path, int& iter_Bell) {
    vector<int> dist(ver, trash1); // Инициализируем расстояния как "бесконечность"
    vector<int> from(ver, -1);  // Для восстановления пути
    dist[start] = 0; // Расстояние до стартовой вершины 0
    iter_Bell = 0;
    // Основной цикл - выполняем ver-1 итераций
    for (int i = 0; i < ver - 1; ++i) 
    {
        bool updated = false;
        // Проходим по всем ребрам графа
        for (int u = 0; u < ver; ++u) {
            for (int v = 0; v < ver; ++v) {
                if (weightmatrix[u][v] != 0) { // Если ребро существует
                    if (dist[v] > dist[u] + weightmatrix[u][v]) {
                        dist[v] = dist[u] + weightmatrix[u][v];
                        from[v] = u;
                        updated = true;
                    }
                }
                iter_Bell++; // 1 итерация = 1 проверка ребра
            }
        }
        // Если на текущей итерации не было улучшений, выходим раньше
        if (!updated) break;
    }

    path = PrintPath(from, finish);
    return dist;
}

int MinPropusk(vector<vector<int>> bandwidth_matrix, int start, int finish, int ver, vector<int>& parent) {
    //parent будет хранить родительскую вершину для каждой вершины, заполняем его значениями -1, вершины еще не посещены.
    fill(parent.begin(), parent.end(), -1);
    //стартовую вершину помечаем значением -2, чтобы обозначить, что это начальная точка поиска.
    parent[start] = -2;
    //в очереди будем хранить пары: {вершина, минимальная пропускная способность до нее}.
    queue<pair<int, int>> q;
    //добавляем начальную вершину в очередь с "бесконечной" пропускной способностью (trash1).
    q.push({ start, trash1 });

    //проводим поиск в ширину (BFS)
    while (!q.empty()) {
        //извлекаем из очереди вершину u и текущую минимальную пропускную способность cap.
        int u = q.front().first;
        int cap = q.front().second;
        q.pop();

        //проверяем всех соседей вершины u.
        for (int v = 0; v < ver; v++) {
            //если существует ребро от u к v (bandwidth_matrix[u][v] != 0), и вершина v еще не была посещена (parent[v] == -1)
            if (u != v && bandwidth_matrix[u][v] != 0 && parent[v] == -1) {
                // Помечаем вершину v как посещенную, запоминаем ее родительскую вершину (откуда пришли).
                parent[v] = u;
                int min_cap = min(cap, bandwidth_matrix[u][v]);
                //если мы достигли вершины finish, сразу возвращаем минимальную пропускную способность.
                if (v == finish) return min_cap;
                //вставляем вершину v в очередь с обновленной пропускной способностью.
                q.push({ v, min_cap });
            }
        }
    }
    return 0;
}

int Ford_Fulkerson(vector<vector<int>> bandwidth_matrix, int start, int finish, int ver) {
    // Инициализируем родительский вектор (parent[i] = -1 означает, что вершина i не была посещена).
    vector<int> parent(ver, -1);
    int max_flow = 0; //для накопления максимального потока
    int min_cap = 0;  //для хранения минимальной пропускной способности пути

    // Пока существует путь от стартовой вершины к конечной через MinPropusk, продолжаем искать поток
    while (min_cap = MinPropusk(bandwidth_matrix, start, finish, ver, parent)) {
        max_flow += min_cap; //увеличиваем общий поток на величину минимальной пропускной способности для найденного пути

        //промежуточная матрица пропускных способностей:
        cout << "Промежуточная матрица пропускных способностей:" << endl;
        for (int i = 0; i < ver; i++) {
            for (int j = 0; j < ver; j++) {
                cout << bandwidth_matrix[i][j] << "   ";
            }
            cout << endl;
        }
        cout << "-----------------" << endl;

        // Выводим промежуточный путь с пропускной способностью
        cout << "Промежуточный путь с пропускной способностью: " << min_cap << endl;
        cout << "Маршрут: ";

        int v_path = finish;  // Для восстановления пути используем переменную v_path
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


        // Обновляем пропускные способности по найденному пути
        int v = finish;
        while (v != start) {
            int u = parent[v];
            bandwidth_matrix[u][v] -= min_cap; //уменьшаем пропускную способность по найденному пути (по направлению потока)
            bandwidth_matrix[v][u] += min_cap; //увеличиваем пропускную способность в обратном направлении (обратный поток)
            v = u; //переходим к предшествующей вершине
        }
    }

    cout << "Промежуточная матрица пропускных способностей:" << endl;
    for (int i = 0; i < ver; i++) {
        for (int j = 0; j < ver; j++) {
            cout << bandwidth_matrix[i][j] << "   ";
        }
        cout << endl;
    }
    cout << "-----------------" << endl;

    return max_flow; //возвращаем максимальный поток, который удалось найти
}



int MinCost(vector<vector<int>>& cost_matrix,
    vector<vector<int>>& bandwidth_matrix,
    int s, int t,
    int ver,
    int need_flow) {
    int total_cost = 0;  //для хранения общей стоимости потока
    int flow = 0; //для хранения текущего потока
    vector<vector<int>> flow_matrix(ver, vector<int>(ver, 0));  //матрица текущего потока между вершинами
    vector<vector<int>> residual_cap = bandwidth_matrix;  //матрица остаточных пропускных способностей

    // Пока текущий поток меньше нужного, продолжаем искать пути
    while (flow < need_flow) {
        vector<vector<int>> residual_cost(ver, vector<int>(ver, 0));  //матрица остаточных стоимостей рёбер

        // Заполняем residual_cost, включая рёбра с потоком в обратном направлении
        for (int i = 0; i < ver; ++i) {
            for (int j = 0; j < ver; ++j) {
                residual_cost[i][j] = cost_matrix[i][j];  // Стоимость из исходной матрицы
                if (flow_matrix[i][j] > 0) {  // Если есть поток по этому ребру
                    residual_cap[j][i] = flow_matrix[i][j];  // Обратный поток
                    residual_cost[j][i] = -cost_matrix[i][j];  // Обратная стоимость (отрицательная)
                }
            }
        }

        // Поиск пути минимальной стоимости (алгоритм Беллмана-Форда)
        vector<int> dist(ver, trash1);  // Массив расстояний (стоимости)
        vector<int> parent(ver, -1);    // Массив родителей для восстановления пути
        dist[s] = 0;  // Стоимость до стартовой вершины равна 0

        // Основной цикл Беллмана-Форда: ищем кратчайшие пути
        for (int i = 0; i < ver - 1; ++i) {
            for (int u = 0; u < ver; ++u) {
                for (int v = 0; v < ver; ++v) {
                    // Если есть ребро, которое не насыщено, и мы нашли более короткий путь
                    if (residual_cap[u][v] > 0 && dist[u] != trash1 &&
                        dist[v] > dist[u] + residual_cost[u][v]) {
                        dist[v] = dist[u] + residual_cost[u][v];  // Обновляем стоимость
                        parent[v] = u;  // Запоминаем, откуда пришли в вершину v
                    }
                }
            }
        }

        // Если не нашли путь до конечной вершины, значит, потока больше не будет
        if (parent[t] == -1) break;

        // Нахождение минимальной пропускной способности по пути
        int delta = need_flow - flow;  // Ограничение на оставшийся поток
        for (int v = t; v != s; v = parent[v]) {
            delta = min(delta, residual_cap[parent[v]][v]);  // Ограничение по минимальной пропускной способности на пути
        }

        // Обновление потока и остаточных пропускных способностей
        for (int v = t; v != s; v = parent[v]) {
            int u = parent[v];
            flow_matrix[u][v] += delta;  // Увеличиваем поток по ребру (u -> v)
            flow_matrix[v][u] -= delta;  // Уменьшаем поток в обратном направлении (v -> u)
            total_cost += delta * cost_matrix[u][v];  // Увеличиваем общую стоимость потока
            residual_cap[u][v] -= delta;  // Уменьшаем пропускную способность по ребру (u -> v)
            residual_cap[v][u] += delta;  // Увеличиваем пропускную способность в обратном направлении (v -> u)
        }
        flow += delta;  // Увеличиваем общий поток
    }

    // Вывод распределения потока
    cout << "Распределение потока:" << endl;
    for (int i = 0; i < ver; ++i) {
        for (int j = 0; j < ver; ++j) {
            if (flow_matrix[i][j] > 0) {  // Если есть поток по ребру
                cout << i << "->" << j << ": " << flow_matrix[i][j]
                    << " (стоимость " << cost_matrix[i][j] << ")" << endl;
            }
        }
    }

    // Возвращаем общую стоимость минимального потока
    return total_cost;
}



// Функция для нахождения корня множества, к которому принадлежит вершина k
int find_root(int k, vector<int> parent) {
    if (parent[k] == k) {  // Если вершина является корнем множества
        return k;
    }
    // Иначе рекурсивно находим корень множества
    return find_root(parent[k], parent);
}

// Функция объединяет два множества, в которые входят вершины f_p и s_p
void union_(int f_p, int s_p, vector<int>& parent) {
    // Находим корни множеств
    int x = find_root(f_p, parent);
    int y = find_root(s_p, parent);
    // Объединяем два множества, устанавливая корень одного множества на корень другого
    parent[x] = y;
}

vector<edge> Kruskal(vector<vector<int>> weightmatrix, int ver, int& mst_weight, int& iter_kruskal) {
    vector<edge> edges;
    // Формируем список рёбер: для каждого ребра добавляем в список пару <вес, откуда, куда>
    for (int i = 0; i < ver; i++) {
        for (int j = 0; j < ver; j++) {
            if (weightmatrix[i][j] != 0) {
                edges.push_back({ weightmatrix[i][j], i, j });
            }
        }
    }

    // Сортируем рёбра по весу (по возрастанию)
    sort(edges.begin(), edges.end(), compare());

    cout << "\nОтсортированные рёбра: " << endl;
    for (const auto& e : edges) {
        cout << e.from << "-" << e.where_ << " Вес: " << e.cost << endl;
    }

    vector<int> parent(ver);
    // Инициализируем родительский вектор, где каждый элемент - это сама вершина (каждая вершина изначально является своим корнем)
    for (int i = 0; i < ver; i++) {
        parent[i] = i;
    }

    vector<edge> mst;
    // Покажем, что минимальное остовное дерево должно содержать ver-1 рёбер
    while (mst.size() != ver - 1) {
        // Берём следующее ребро (с минимальным весом, так как мы отсортировали рёбра)
        edge next_edge = edges.back();
        edges.pop_back();

        // Ищем корни множества для первой вершины ребра и второй вершины
        int f_p = find_root(next_edge.from, parent);
        int s_p = find_root(next_edge.where_, parent);

        // Если корни различных множеств (т.е. добавление этого ребра не создаст цикл), то добавляем это ребро в остовное дерево
        if (f_p != s_p) {
            mst.push_back(next_edge);  // Добавляем ребро в остовное дерево
            mst_weight += next_edge.cost;  // Увеличиваем суммарный вес остовного дерева
            union_(f_p, s_p, parent);  // Объединяем множества
        }

        iter_kruskal++;  // Увеличиваем счётчик количества итераций (по рёбрам)
    }
    return mst;  // Возвращаем минимальное остовное дерево
}


vector<pair<int, int>> codePrufer(vector<edge> mst, int ver) {
    // Инициализация матрицы смежности для дерева (минимального остовного дерева)
    vector<vector<int>> mst_matrix(ver, vector<int>(ver, 0));

    // Заполняем матрицу смежности, где mst_matrix[i][j] - это вес ребра между вершинами i и j
    for (int k = 0; k < mst.size(); k++) {
        int i = mst[k].from;
        int j = mst[k].where_;
        int w = mst[k].cost;
        mst_matrix[i][j] = w;
        mst_matrix[j][i] = w;
    }

    // Инициализация вектора для хранения кода Прюфера
    vector<pair<int, int>> prufer;
    vector<int> degree_vert(ver);  // Массив для хранения степеней вершин

    // Вычисляем степени всех вершин
    for (int i = 0; i < ver; i++) {
        for (int j = 0; j < ver; j++) {
            if (mst_matrix[i][j] != 0) {
                degree_vert[i]++;
            }
        }
    }

    // Считаем количество рёбер
    int r = 0;
    for (int i = 0; i < ver; i++) {
        r += degree_vert[i];
    }

    // Пока есть рёбра (r != 0), продолжаем извлекать минимальные вершины с степенью 1
    while (r != 0) {
        for (int i = 0; i < ver; i++) {
            if (degree_vert[i] == 1) {  // Если степень вершины равна 1, это кандидат на удаление
                for (int j = 0; j < ver; j++) {
                    if (mst_matrix[i][j] != 0) {
                        // Добавляем ребро (вес и вершину, с которой оно соединено)
                        prufer.push_back(make_pair(mst_matrix[i][j], j));
                        degree_vert[i]--;
                        degree_vert[j]--;
                        mst_matrix[i][j] = 0;
                        mst_matrix[j][i] = 0;
                        r -= 2;  // Уменьшаем количество рёбер на 2 (мы удалили вершину)
                    }
                }
            }
        }
    }

    return prufer;
}


vector<vector<int>> decodePrufer(vector<pair<int, int>> prufer) {
    // Количество вершин в дереве (размер кода Прюфера + 1)
    int ver = prufer.size() + 1;

    // Инициализация матрицы смежности для возвращаемого дерева
    vector<vector<int>> ans_matrix(ver, vector<int>(ver, 0));

    // Массив для отслеживания степеней вершин
    vector<int> vertexes(ver, 0);

    // Вычисляем количество рёбер, исходящих из каждой вершины (степени)
    for (int i = 0; i < prufer.size(); i++) {
        vertexes[prufer[i].second] += 1;
    }

    int j = 0;
    int num = 0;

    // Восстанавливаем дерево по коду Прюфера
    for (int i = 0; i < ver - 1; i++) {
        for (j = 0; j < ver; j++) {
            if (i == ver - 2) {
                if (vertexes[j] == 0) {  // Находим вершину с степенью 0
                    vertexes[j] = -1;  // Помечаем вершину как использованную
                    ans_matrix[j][prufer[i].second] = prufer[i].first;
                    ans_matrix[prufer[i].second][j] = prufer[i].first;
                    vertexes[prufer[i].second]--;  // Уменьшаем степень соседней вершины
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
    return ans_matrix;  // Возвращаем матрицу смежности восстановленного дерева
}

void printAllEdges(const vector<vector<int>>& matrix, int ver) {
    cout << "Все рёбра графа:" << endl;
    for (int u = 0; u < ver; ++u) {
        for (int v = u + 1; v < ver; ++v) {
            if (matrix[u][v] != 0) {
                cout << u << " <-> " << v << " (Вес: " << matrix[u][v] << ")" << endl;
            }
        }
    }
}

vector<pair<int, int>> MaxIndependentEdgeSet(vector<vector<int>>& matrix, int ver) {
    vector<pair<int, int>> matching;  // Вектор для хранения выбранных рёбер
    vector<bool> used(ver, false);  // Массив для отслеживания использованных вершин

    // Перебираем все вершины
    for (int u = 0; u < ver; ++u) {
        if (!used[u]) {  // Если вершина еще не использована
            // Ищем ребро, которое не использует вершину u и вершину v
            for (int v = 0; v < ver; ++v) {
                if (matrix[u][v] != 0 && !used[v]) {  // Если существует ребро и вершина v не использована
                    matching.push_back({ u, v });  // Добавляем ребро в множество
                    used[u] = true;  // Отмечаем вершину u как использованную
                    used[v] = true;  // Отмечаем вершину v как использованную
                    break;  // Переходим к следующей вершине u
                }
            }
        }
    }

    return matching;  // Возвращаем множество рёбер
}

vector<vector<int>> Euler(int ver, vector<vector<int>> undirected_weight_matrix, vector<int> undirected_degrees) {
    vector<vector<int>> eulergraph(undirected_weight_matrix);
    int countevenvert = 0;
    bool flagEuler = false;
    cout << "\nМатрица весов: " << endl;
    PrintMatrix(eulergraph);
    //вывод всех степеней вершин
    for (int i = 0; i < ver; i++) cout << "Степень вершины №" << i << " = " << undirected_degrees[i] << endl;

    //проверяем, что все степени вершин графа четные => эйлеров граф
    for (int i = 0; i < ver; i++) {
        if (undirected_degrees[i] % 2) countevenvert++;
    }
    if (countevenvert == ver) {
        cout << "Граф является Эйлеровым\n" << endl;
        flagEuler = true;
    }
    //если не все степени вершин графа четные => надо модифицировать граф - сделать эйлеровым
    else {
        cout << "\nГраф не является Эйлеровым. Модифицируем граф.\n" << endl;
        //будем удалять или добавлять ребра пока граф не станет эйлеровым
        while (!flagEuler) {
            bool flag = 0;
            for (int i = 0; i < ver; i++) {
                for (int j = 0; j < ver; j++) {
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
                if (flag) break;
            }
            //если не получилось добавить ребро, то пробуем удалить
            if (!flag) {
                for (int i = 0; i < ver; i++) {
                    for (int j = 0; j < ver; j++) {
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
                    if (flag) break;
                }
            }
            //проверяем получилось ли сделать граф с четными степенями вершин
            int countevenvert1 = 0;
            for (int i = 0; i < ver; i++) {
                if (undirected_degrees[i] % 2 == 0) countevenvert1++;
            }
            //если получилось, то поднимаем флаг и выходим из цикла, если нет - продолжаем добавлять/удалять ребро
            if (countevenvert1 == ver) flagEuler = true;
        }
        cout << "***********************" << endl;
        cout << "Модифицированный эйлеров граф." << endl;
        cout << "Матрица весов: " << endl;
        PrintMatrix(eulergraph);
        for (int i = 0; i < ver; i++) cout << "Степень вершины №" << i << " = " << undirected_degrees[i] << endl;
        cout << "***********************" << endl;
    }
    return eulergraph;
}

void FindEuler(int ver, vector<vector<int>> eulergraph) {
    stack<int> vertexes;
    vector<int> cycle;
    int vert = 0;
    //ищем первую вершину, из которой есть ребро
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
    //добавляем эту вершину в стек
    vertexes.push(vert);
    //пока стек не пустой
    while (!vertexes.empty()) {
        //берем верхнюю вершину
        vert = vertexes.top();
        int i;
        //находим первую смежную с ней
        for (i = 0; i < ver; i++) {
            if (eulergraph[vert][i] != 0) {
                break;
            }
        }
        //если смежной не нашлось, то добавляем эту вершину в цикл, удаляем из стека
        if (i == ver) {
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


//проверка возможности добавления вершины в путь
bool CanAddV(int v, vector<vector<int>>& gamiltongraph, vector<int>& path, int p) {
    if (gamiltongraph[path[p - 1]][v] == 0) return false; //ребро между последней вершиной path[p - 1] и вершиной v в матрице смежности
    for (int i = 0; i < p; i++) {
        if (path[i] == v) return false;
    }
    return true;
}

bool FindGamilton(vector<vector<int>>& gamiltongraph, vector<int>& path, int p) {
    //если все вершины посещены и последняя соединена с первой, то возвращается true => есть гамильтонов цикл в графе
    if (p == gamiltongraph.size()) return gamiltongraph[path[p - 1]][path[0]] != 0;

    //иначе ищем следующую вершину для посещения
    for (int v = 1; v < gamiltongraph.size(); v++) {
        //в которую есть возможность дойти (есть исходящее ребро + не была посещена смежная вершина)
        if (CanAddV(v, gamiltongraph, path, p)) {
            //если найлена, то добавляем эту вершину в путь
            path[p] = v;
            //если гамильтонов цикл найден, то возвращаем true
            if (FindGamilton(gamiltongraph, path, p + 1)) return true;
            //если гамильтонов цикл не найден, то удаляем вершину v из пути
            path[p] = -1;
        }
    }
    return false;
}

// Функция для вывода всех гамильтоновых циклов с их суммарным весом
void PrintAllHamiltonCyclesWithWeights(vector<vector<int>>& gamiltongraph, vector<int>& path, int position, int& cycle_count, vector<vector<int>>& weight_matrix) {
    int ver = gamiltongraph.size();
    if (position == ver) {
        // Проверяем, есть ли ребро между последней и первой вершиной
        if (gamiltongraph[path[position - 1]][path[0]] != 0) {
            // Вычисляем суммарный вес цикла
            int total_weight = 0;
            for (int i = 0; i < ver - 1; i++) {
                total_weight += weight_matrix[path[i]][path[i + 1]];
            }
            total_weight += weight_matrix[path[ver - 1]][path[0]];  // Замыкаем цикл

            // Выводим цикл с весом
            cout << "Цикл " << ++cycle_count << ": ";
            for (int i = 0; i < ver; i++) {
                cout << path[i] << (i < ver - 1 ? " -> " : "");
            }
            cout << " -> " << path[0] << " (Суммарный вес: " << total_weight << ")" << endl;
        }
        return;
    }

    // Перебираем все возможные следующие вершины
    for (int v = 0; v < ver; v++) {
        if (CanAddV(v, gamiltongraph, path, position)) {
            path[position] = v;
            PrintAllHamiltonCyclesWithWeights(gamiltongraph, path, position + 1, cycle_count, weight_matrix);
            path[position] = -1;  // Откатываем изменения
        }
    }
}



//функция для проверки наличия гамильтонова цикла (рекурсивный)
bool ExistGamiltonCycle(vector<vector<int>>& gamiltongraph) {
    vector<int> path(gamiltongraph.size(), -1);
    path[0] = 0;
    if (!FindGamilton(gamiltongraph, path, 1)) return false;
    return true;
}

vector<vector<int>> Gamilton(int ver, vector<vector<int>> undirected_weight_matrix) {
    vector<vector<int>> gamiltongraph(undirected_weight_matrix);

    // Выводим матрицу стоимости
    cout << "\nМатрица весов: " << endl;
    PrintMatrix(gamiltongraph);

    // Проверяем наличие гамильтоновых циклов
    vector<int> path(ver, -1);
    path[0] = 0;  // Начинаем с вершины 0
    int cycle_count = 0;

    cout << "\nВывод найденных гамильтоновых циклов." << endl;
    PrintAllHamiltonCyclesWithWeights(gamiltongraph, path, 1, cycle_count, undirected_weight_matrix);

    if (cycle_count > 0) {
        cout << "***********************" << endl;
        cout << "\nГраф является гамильтоновым. Найдено циклов: " << cycle_count << endl;
        cout << "***********************" << endl;
        return gamiltongraph;
    }
    else {
        cout << "\nГраф не является гамильтоновым. Модифицируем граф." << endl;
        vector<int> vect(ver, -1);

        // Генерируем случайную перестановку вершин
        for (int i = 0; i < ver; i++) {
            int r = rand() % ver;
            while (find(vect.begin(), vect.end(), r) != vect.end()) {
                r = rand() % ver;
            }
            vect[i] = r;
        }

        // Добавляем недостающие ребра
        for (int i = 0; i < ver; i++) {
            if (gamiltongraph[vect[i]][vect[(i + 1) % ver]] == 0) {
                int r = rand() % 10 + 1;
                gamiltongraph[vect[i]][vect[(i + 1) % ver]] = r;
                gamiltongraph[vect[(i + 1) % ver]][vect[i]] = r;
                cout << "Добавлено ребро: " << vect[i] << " <-> " << vect[(i + 1) % ver]
                    << " с весом: " << r << endl;
            }
        }

        // Проверяем снова после модификации
        cycle_count = 0;
        fill(path.begin(), path.end(), -1);
        path[0] = 0;

        cout << "\nПроверка после модификации." << endl;
        PrintAllHamiltonCyclesWithWeights(gamiltongraph, path, 1, cycle_count, undirected_weight_matrix);

        if (cycle_count > 0) {
            cout << "***********************" << endl;
            cout << "\nТеперь граф гамильтонов. Найдено циклов: " << cycle_count << endl;
            cout << "***********************" << endl;
        }
        else {
            cout << "\nНе удалось сделать граф гамильтоновым." << endl;
        }
    }

    cout << "\nИтоговая матрица весов:" << endl;
    PrintMatrix(gamiltongraph);

    return gamiltongraph;
}

void FindGamiltonKomivvv(vector<vector<int>>& gamiltongraph, vector<bool>& visited, vector<int>& path, int vert, int ver, int count, int cost, int& min_cost, vector<int>& min_path) {
    path.push_back(vert);
    if (count == ver and gamiltongraph[vert][0]) {
        // Выводим текущий цикл в консоль
        for (int i = 0; i < path.size(); i++) cout << path[i] << " -> ";
        cout << path[0];
        cout << "\tВес цикла: " << gamiltongraph[vert][0] + cost << endl;

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
    vector<bool> visited(ver); //посещенность вершины
    vector<int> path; //текущий путь
    visited[0] = true;
    int min_cost = trash1; //минимальная стоимость пути
    vector<int> min_path; //минимальный путь

    cout << "***********************" << endl;
    cout << "Все возможные гамильтоновы циклы:" << endl;

    FindGamiltonKomivvv(gamiltongraph, visited, path, 0, ver, 1, 0, min_cost, min_path);

    cout << "***********************" << endl;
    cout << "Цикл коммивояжера: " << endl;
    for (int i = 0; i < min_path.size(); i++) cout << min_path[i] << " -> ";
    cout << min_path[0];
    cout << "\nВес: " << min_cost << endl;
}


void menu() {

    printf("\nВведите число от 0 до 19\n");
    printf("0) Выход\n");
    printf("1) Создать граф\n");
    printf("2) Вывести матрицу смежности вершин\n");
    printf("3) Метод Шимбелла\n");
    printf("4) Возможность построения маршрута и их количество\n");
    printf("5) Обход вершин графа поиском в ширину\n");
    printf("6) Вывести матрицу весов\n");
    printf("7) Алгоритм Беллмана-Форда\n");
    printf("8) Сравнить скорости работы 5 и 8\n");
    printf("9) Матрица пропускных способностей\n");
    printf("10) Матрица стоимости\n");
    printf("11) Максимальный поток по алгоритму Форда-Фалкерсона\n");
    printf("12) Поток минимальной стоимости\n");
    printf("13) Число остовных деревьев по Кирхгофу\n");
    printf("14) Минимальный по весу остов по Краскалу\n");
    printf("15) Код Прюфера (кодировать и декодировать)\n");
    printf("16) Максимальное независимое множество ребер\n");
    printf("17) Проверить, является ли граф эйлеровым/модифицировать граф, если не эйлеров.\n");
    printf("18) Проверить, является ли граф гамильтоновым/модифицировать граф, если не гамильтонов.\n");
    printf("19) Решить задачу коммивояжера на гамильтоновом графе.\n");


}

