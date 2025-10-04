#include "functions.h"
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

double alpha = 1; //������� ������ �������������� ����� ����������, �� ���� ������� � ������� ��������

//?�� ?������������ ����� �������� ����� [1,3] (����, ��������, ����� ������)
//����� ���������� "��������" �������������.
//��� ������� - ������������� ���������� ����� ��������������, � ����������� ��������� � ������� �������� � ���������� � ������ ��������
//��� ������ - ������������� ���������� ����� �����������, � ������� ������ �������������� ����� ����������


int xmin = 1;
//�������� ������� ������
//xmin ��������� ������� � ������� ��������� �� ��� ������ � �� ������ �� ����������:(

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "����� ����������!" << endl;
    cout << "����������:" << endl;
    cout << "����������� ������� ������������ ���� � ������������ � �������������� ������ � �������� ������������� ����������� ������" << endl;
    cout << "alpha = " << alpha << "   xmin = " << xmin << endl << endl;


    if (xmin != 1) return 0;
    if (alpha <= 0) return 0;

    bool continuef = true;
    while (continuef) {
        cout << "������� ������ ���������� ������ � ����� (�� 2 �� 100):" << endl;
        int ver;
        bool validInput = false;
        while (!validInput) {
            cin >> ver;
            if (cin.fail() || ver > 100 || ver < 2) {
                cin.clear(); //������� ����� ������
                cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                cout << "������ �����! ������� ����� �� 2 �� 100 ������������!" << endl;
            }
            else {
                validInput = true;
            }
        }
        vector<int> degrees(ver); //������ �������� ������ 

        GenerateVertexDegreesPareto(degrees, ver, alpha, xmin); //���������� ������� ���� ������

        vector<vector<int>> admatrix = Createadmatrix(ver, degrees);

        cout << "1 - ���� � �������������� ������, 2 - ���� � �������������� � �������������� ������" << endl;
        int sign;
        bool validInput1 = false;
        while (!validInput1) {
            cin >> sign;
            if (cin.fail() || sign > 2 || sign < 1) {
                cin.clear(); //������� ����� ������
                cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                cout << "������ �����! ������� ���� 1, ���� 2!" << endl;
            }
            else {
                validInput1 = true;
            }
        }

        vector<vector<int>> ed_matrix(ver, vector<int>(ver, 0)); //��������� �������
        //��������� ��������� �������
        for (int i = 0; i < ver; i++) {
            for (int j = 0; j < ver; j++) {
                if (i == j) {
                    ed_matrix[i][j] = 1;
                }
            }
        }
        //��������� ������� ������������
        vector<vector<int>> dostmatrix(admatrix);
        vector<vector<int>> res_dostmatrix(ver, vector<int>(ver, 0));
        res_dostmatrix = Addmatrix(ed_matrix, admatrix, ver);
        for (int i = 0; i < ver; i++) {
            dostmatrix = Createdostmatrix(admatrix, dostmatrix, ver);
            res_dostmatrix = Addmatrix(res_dostmatrix, dostmatrix, ver);
        }

        // ��������� ������ � �������������� ������
        vector<vector<int>> weightmatrix = GenerateWeightMatrix(ver, admatrix, alpha, xmin);
        if (sign == 2) {
            // ��� ������������� ����� ����������� ����� ��������
            random_device rd;
            mt19937 gen(rd());
            uniform_int_distribution<> dis(0, 1); // 50% chance to negate

            for (int i = 0; i < ver; i++) {
                for (int j = 0; j < ver; j++) {
                    if (weightmatrix[i][j] != 0 && dis(gen) == 1) {
                        weightmatrix[i][j] *= -1;
                    }
                }
            }
        }



        vector<vector<int>> cost_matrix = GenerateCostMatrix(ver, admatrix, alpha, xmin);
        vector<vector<int>> bandwidth_matrix = GenerateCapacityMatrix(ver, admatrix, alpha, xmin);

        //��������� ������� ��������� ������������������ �����
        vector<vector<int>> undirected_admatrix(ver, vector<int>(ver, 0));
        for (int i = 0; i < ver; i++) {
            for (int j = 0; j < ver; j++) {
                if (admatrix[i][j] != 0) {
                    undirected_admatrix[j][i] = admatrix[i][j];
                    undirected_admatrix[i][j] = admatrix[i][j];
                }
            }
        }

        //��������� ������� �� ��������� ������ �� ������� ��������� ��� ������������������ �����
        vector<vector<int>> degree_matrix(ver, vector<int>(ver, 0));
        for (int i = 0; i < ver; i++) {
            //int d = 0;
            for (int j = 0; j < ver; j++) {
                if (undirected_admatrix[i][j] != 0) {
                    degree_matrix[i][i] += 1;
                }
            }

        }

        //��������� ������ �������� ������ ������������������ �����
        vector<int> undirected_degrees(ver, 0);
        for (int i = 0; i < ver; i++) {
            for (int j = 0; j < ver; j++) {
                if (degree_matrix[i][j] != 0) {
                    undirected_degrees[i] = degree_matrix[i][j];
                }
            }
        }

        //��������� ������� ��������
        vector<vector<int>> kirhgof_matrix(ver, vector<int>(ver, 0));
        for (int i = 0; i < ver; i++) {
            for (int j = 0; j < ver; j++) {
                kirhgof_matrix[i][j] = degree_matrix[i][j] - undirected_admatrix[i][j];
            }
        }


        //��������� ������� ����� ������������������ �����
        vector<vector<int>> undirected_weightmatrix(ver, vector<int>(ver, 0));
        for (int i = 0; i < ver; i++) {
            for (int j = 0; j < ver; j++) {
                if (cost_matrix[i][j] != 0) {
                    undirected_weightmatrix[j][i] = weightmatrix[i][j];
                    undirected_weightmatrix[i][j] = weightmatrix[i][j];
                }
            }
        }

        menu();

        int iter_Bell = 0;
        int iter_BFS = 0;
        float iter_div = 0;
        ofstream fout;
        int input = 0;
        bool flag = true;
        vector<vector<int>> gamiltongraph(ver, vector<int>(ver, 0));
        int max_flow = -1;

        while (flag)
        {
            bool validInput2 = false;
            while (!validInput2) {
                cin >> input;
                if (cin.fail() || input > 19 || input < 0) {
                    cin.clear(); //������� ����� ������
                    cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                    cout << "������ �����! ����� �� 0 �� 19 ������������!" << endl;
                }
                else {
                    validInput2 = true;
                }
                switch (input)
                {
                case 0: {
                    flag = false;
                    continuef = false;
                    break;
                }
                case 1: {
                    flag = false;
                    break;
                }
                case 2: {
                    cout << "***********************" << endl;
                    PrintMatrix(admatrix);
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;
                }
                case 3: {
                    cout << "������� ����� �������� (����� �������� ������ ���� �� 1 �� " << ver - 1 << " �����): " << endl;
                    int length;
                    bool validInput3 = false;
                    while (!validInput3) {
                        cin >> length;
                        if (cin.fail() || length > (ver - 1) || length < 1) {
                            cin.clear(); //������� ����� ������
                            cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                            cout << "������ �����! ������� ����� �� 1 �� (���-�� ������ - 1) ������������!" << endl;
                        }
                        else {
                            validInput3 = true;
                        }
                    }

                    cout << "��������: 1 - ����� �������� ������������ �����; 2 - ����� �������� ����������� �����" << endl;
                    int parametr; //���� 1 - ���� max, ���� 2 - ���� min
                    bool validInput4 = false;
                    while (!validInput4) {
                        cin >> parametr;
                        if (cin.fail() || parametr > 2 || parametr < 1) {
                            cin.clear(); //������� ����� ������
                            cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                            cout << "������ �����! ������� ���� 1, ���� 2!" << endl;
                        }
                        else {
                            validInput4 = true;
                        }
                    }
                    vector<vector<int>> Shimbel_matrix(weightmatrix);
                    for (int i = 0; i < length - 1; i++) {
                        Shimbel_matrix = Shimbel(weightmatrix, Shimbel_matrix, ver, parametr);
                    }
                    cout << "***********************" << endl;
                    PrintMatrix(Shimbel_matrix);
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;

                }
                case 4: {
                    int start, finish;
                    cout << "������� ����� ������� �� ������� ������ ��������� ������� �� 0 �� " << ver - 1 << ":" << endl;
                    bool validInput3 = false;
                    while (!validInput3) {
                        cin >> start;
                        if (cin.fail() || start > (ver - 1) || start < 0) {
                            cin.clear(); //������� ����� ������
                            cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                            cout << "������ �����! ������� ����� �� 0 �� (���-�� ������ - 1) ������������!" << endl;
                        }
                        else validInput3 = true;
                    }
                    cout << "������� ����� ������� � ������� ������ ��������� ������� �� 0 �� " << ver - 1 << ":" << endl;
                    bool validInput4 = false;
                    while (!validInput4) {
                        cin >> finish;
                        if (cin.fail() || finish > (ver - 1) || finish < 0) {
                            cin.clear(); //������� ����� ������
                            cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                            cout << "������ �����! ������� ����� �� 0 �� (���-�� ������ - 1) ������������!" << endl;
                        }
                        else validInput4 = true;

                    }
                    cout << "***********************" << endl;
                    if (res_dostmatrix[start][finish] != 0) cout << "���-�� ��������� �� ������� " << start << " � ������� " << finish << ": " << res_dostmatrix[start][finish] << endl;
                    else cout << "����� ������� ��������� ������" << endl;
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;
                }
                case 5: {

                    int start, finish;
                    cout << "������� ����� ������� �� ������� ������ ��������� ������� �� 0 �� " << ver - 1 << ":" << endl;
                    bool validInput3 = false;
                    while (!validInput3) {
                        cin >> start;
                        if (cin.fail() || start > (ver - 1) || start < 0) {
                            cin.clear(); //������� ����� ������
                            cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                            cout << "������ �����! ������� ����� �� 0 �� (���-�� ������ - 1) ������������!" << endl;
                        }
                        else validInput3 = true;
                    }
                    finish = ver - 1;
                    vector<int> path;
                    vector<int> dist = BFS(weightmatrix, ver, start, finish, path, iter_BFS);
                    cout << "***********************" << endl;
                    PrintWay(ver, dist, start, finish, path);
                    cout << "\n";
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;

                }
                case 6: {
                    cout << "***********************" << endl;
                    PrintMatrix(weightmatrix);
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;

                }
                case 7: {
                    //������ ���� �����
                    int start, finish;
                    cout << "������� ����� ������� �� ������� ������ ��������� ������� �� 0 �� " << ver - 1 << ":" << endl;
                    bool validInput3 = false;
                    while (!validInput3) {
                        cin >> start;
                        if (cin.fail() || start > (ver - 1) || start < 0) {
                            cin.clear(); //������� ����� ������
                            cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                            cout << "������ �����! ������� ����� �� 0 �� (���-�� ������ - 1) ������������!" << endl;
                        }
                        else validInput3 = true;
                    }
                    cout << "������� ����� ������� � ������� ������ ��������� ������� �� 0 �� " << ver - 1 << ":" << endl;
                    bool validInput4 = false;
                    while (!validInput4) {
                        cin >> finish;
                        if (cin.fail() || finish > (ver - 1) || finish < 0) {
                            cin.clear(); //������� ����� ������
                            cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                            cout << "������ �����! ������� ����� �� 0 �� (���-�� ������ - 1) ������������!" << endl;
                        }
                        else validInput4 = true;

                    }
                    vector<int> path;
                    vector<int> dist = BellmanFord(weightmatrix, ver, start, finish, path, iter_BFS);
                    cout << "***********************" << endl;
                    PrintWay(ver, dist, start, finish, path);
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;
                }
                case 8: {
                    //��������� 
                    int start = 0;
                    int finish = ver - 1;

                    vector<int> path1;
                    vector<int> path2;
                    vector<int> dist1 = BFS(weightmatrix, ver, start, finish, path1, iter_BFS);
                    vector<int> dist2 = BellmanFord(weightmatrix, ver, start, finish, path2, iter_Bell);
                    cout << "\n���������� �������� ��������� ����� � � ������: " << iter_BFS << endl;
                    cout << "���������� �������� ��������� ��������-�����: " << iter_Bell << endl;
                    cout << "***********************" << endl;
                    if (iter_BFS > iter_Bell) cout << "�������� ��������-����� ������� � " << (double)iter_BFS / iter_Bell << " ����." << endl;
                    if (iter_BFS < iter_Bell) cout << "����� ������ ������� � ������ ������� � " << (double)iter_Bell / iter_BFS << " ����." << endl;
                    if (iter_BFS == iter_Bell) cout << "��������� ��������� �� ��������." << endl;
                    cout << "***********************" << endl;
                    menu();
                    break;
                }
                case 9: {
                    cout << "***********************" << endl;
                    PrintMatrix(bandwidth_matrix);
                    cout << "\n";
                    cout << "***********************" << endl;
                    menu();
                    break;
                }
                case 10: {
                    cout << "***********************" << endl;
                    PrintMatrix(cost_matrix);
                    cout << "\n";
                    cout << "***********************" << endl;
                    menu();
                    break;
                }
                case 11: {
                    int from = 0;
                    int where_ = ver - 1;
                    max_flow = Ford_Fulkerson(bandwidth_matrix, from, where_, ver);
                    cout << "***********************" << endl;
                    cout << "������������ �����: " << max_flow << endl;
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;
                }
                case 12: {

                    int start = 0;
                    int finish = ver - 1;
                    vector<int> path;
                    if (max_flow == -1) max_flow = Ford_Fulkerson(bandwidth_matrix, start, finish, ver);
                    int need_flow = (max_flow * 2) / 3;
                    if (max_flow == 1) need_flow = 1;
                    cout << "�������� ������ [2/3*max]: " << need_flow << endl;
                    int cost = MinCost(cost_matrix, bandwidth_matrix, start, finish, ver, need_flow);
                    cout << "***********************" << endl;
                    cout << "����������� ��������� ������: " << cost << endl;
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;
                
                }
                case 13: {
                    cout << "�������� ����������������� ���� �� ����������������." << endl;
                    cout << "������� ���������: " << endl;
                    PrintMatrix(undirected_admatrix);

                    cout << "������� ��������: " << endl;
                    PrintMatrix(kirhgof_matrix);
                    cout << "\n";

                    vector<vector<int>> kirhgof_matrix2(ver - 1, vector<int>(ver - 1, 0));
                    for (int i = 0; i < ver - 1; i++) {
                        for (int j = 0; j < ver - 1; j++) kirhgof_matrix2[i][j] = kirhgof_matrix[i + 1][j + 1];
                    }
                    cout << "***********************" << endl;
                    if (ver == 2) cout << "����� �������� �������� �� ��������: 1" << "\n";
                    else cout << "����� �������� �������� �� ��������: " << determinant(ver - 1, kirhgof_matrix2) << "\n";
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;
                }
                case 14: {
                    int mst_weight = 0;
                    int iter_kruskal = 0;
                    vector<edge> kr = Kruskal(undirected_weightmatrix, ver, mst_weight, iter_kruskal);
                    cout << "***********************" << endl;
                    cout << "����������� �� ���� �����: " << mst_weight << endl;
                    cout << "***********************" << endl;
                    for (int i = 0; i < kr.size(); i++) {
                        cout << kr[i].from << "-" << kr[i].where_ << " ���: " << "" << kr[i].cost << endl;
                    }
                    cout << "\n";
                    menu();
                    break;

                }
                case 15: {
                    //����� � �����
                    int mst_weight = 0;
                    int iter_kruskal = 0;
                    vector<edge> kr = Kruskal(undirected_weightmatrix, ver, mst_weight, iter_kruskal);
                    vector<pair<int, int>> pr = codePrufer(kr, ver);
                    cout << "�����: " << mst_weight << endl;
                    for (int i = 0; i < kr.size(); i++) {
                        cout << kr[i].from << "-" << kr[i].where_ << " ���: " << "" << kr[i].cost << endl;
                    }
                    cout << "\n";
                    cout << "***********************" << endl;
                    cout << "��� �������: " << endl;
                    for (int i = 0; i < pr.size(); i++) cout << "����� �������: " << pr[i].second << " ��� �����: " << pr[i].first << endl;
                    cout << "***********************" << endl;
                    vector<vector<int>> dec_pr = decodePrufer(pr);
                    cout << "***********************" << endl;
                    cout << "\n�������������� ��� �������: " << endl;
                    PrintMatrix(dec_pr);
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;
                }
                case 16: {
                    //���� ����� ����� �����
                    cout << "��������: 1 - �������� ����, 2 - ����� (�� ��������� ��������): ";
                    int choice;
                    bool validInput9 = false;
                    while (!validInput9) {
                        cin >> choice;
                        if (cin.fail() || choice > 2 || choice < 1) {
                            cin.clear(); //������� ����� ������
                            cin.ignore(numeric_limits<streamsize>::max(), '\n'); //������� ������ �����
                            cout << "������ �����! ������� ���� 1, ���� 2!" << endl;
                        }
                        else validInput9 = true;

                    }

                    vector<vector<int>> target_matrix;
                    if (choice == 1) {
                        target_matrix = undirected_weightmatrix;
                        cout << "�������� ����:" << endl;
                        PrintMatrix(target_matrix);
                        
                    }
                    else {
                        int mst_weight = 0;
                        int iter_kruskal = 0;
                        vector<edge> kr = Kruskal(undirected_weightmatrix, ver, mst_weight, iter_kruskal);
                        target_matrix = vector<vector<int>>(ver, vector<int>(ver, 0));
                        for (const auto& e : kr) {
                            target_matrix[e.from][e.where_] = e.cost;
                            target_matrix[e.where_][e.from] = e.cost;
                        }
                        cout << "�����:" << endl;
                        PrintMatrix(target_matrix);

                    }
                    printAllEdges(target_matrix, ver);
                    vector<pair<int, int>> matching = MaxIndependentEdgeSet(target_matrix, ver);
                    cout << "������������ ����������� ��������� ����� (�������������):" << endl;
                    for (const auto& edge : matching) {
                        cout << edge.first << " <-> " << edge.second << endl;
                    }
                    cout << "***********************" << endl;
                    cout << "������ �������������: " << matching.size() << endl;
                    cout << "***********************" << endl;
                    cout << "\n";
                    menu();
                    break;
                }
                case 17: {
                    //���������, �������� �� ���� ���������

                    if (ver <= 2) {
                        cout << "***********************" << endl;
                        cout << "���� �������� 2 �������. ���� �� �������� ��������� � ��� ������ ��������������." << endl;
                        cout << "***********************" << endl;
                    }
                    else {
                        vector<vector<int>> eulergraph = Euler(ver, undirected_weightmatrix, undirected_degrees);
                        FindEuler(ver, eulergraph);
                    }

                    cout << "\n";
                    menu();
                    break;
                }
                case 18: {
                    //���������, �������� �� ���� �������������. ���� ���, �� �������������� ���� (����������, ��� ��������)
                    gamiltongraph.assign(ver, vector<int>(ver, 0));
                    if (ver <= 2) {
                        cout << "***********************" << endl;
                        cout << "���� �������� 2 �������. ���� �� �������� ������������� � ��� ������ ��������������." << endl;
                        cout << "***********************" << endl;
                    }
                    else gamiltongraph = Gamilton(ver, undirected_weightmatrix);
                    cout << "\n";
                    menu();
                    break;
                }
                case 19: {
                    vector<vector<int>> gamiltongraph(ver, vector<int>(ver, 0));
                    if (ver <= 2) {
                        cout << "���� �� ����� ���� �������������, ��� ��� ������� �� 2-� ������" << endl;
                    }
                    else {
                        gamiltongraph = Gamilton(ver, undirected_weightmatrix);
                        Komivoyadger(gamiltongraph, ver);
                    }
                    cout << "\n";
                    menu();
                    break;
                }

                }
            }
        }
    }
    cout << "�� ����� ������!" << endl;
    return 0;
}

