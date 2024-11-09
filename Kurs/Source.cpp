#include <iostream>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <cmath>
#include<Windows.h>
using namespace std;
//Функция создания симплека n+1
vector<vector<double>> CreateSimplex(int n) {
	vector<vector<double>> simplex(n + 1, vector<double>(n));

	for (int i = 0; i < n + 1; ++i) {
		for (int j = 0; j < n; ++j) {
			simplex[i][j] = -5 + static_cast<double>(rand()) / (RAND_MAX / 10); // Случайные значения от -5 до 5
		}
	}
	return simplex;

}
//Функция выводы симлекса
void printSimplex(const vector<vector<double>>& simplex) {
	for (size_t i = 0; i < simplex.size(); ++i) {
		cout << "Vertex " << i + 1 << ": ";
		for (size_t j = 0; j < simplex[i].size(); ++j) {
			cout << simplex[i][j] << " ";
		}
		cout << endl;
	}
}

double function_value(vector<double> Simplex) {
	return pow(Simplex[0],2) - 3*Simplex[0] + 2* Simplex[1] - 20;
}

vector<vector<double>>Method_Neldera_Mida(vector<vector<double>> Simplex, double alpha, double betha, double gamma, int n) {
	double Sr = 0;
	double F = 0;
	double S = 10;
	double exp = 0.00000001;
	int vertices = Simplex.size();
	int dimension = Simplex[0].size();

	//Сортировка точек симплкса в соответствии со значениями функций (от большей к меньшей)
	sort(Simplex.begin(), Simplex.end(), [](const vector<double>& a, const vector<double>& b) {return function_value(a) > function_value(b); });
	//Точки с max значением функции, следующим за max, и min
	vector<double> Xh = Simplex[0];
	vector<double> Xg = Simplex[1];
	vector<double> Xl = Simplex[Simplex.size() - 1];
	//Вектора нужные для вычислений 
	vector<double> center_of_mass(dimension, 0.0);
	vector<double> Xr(dimension, 0.0);
	vector<double> Xe(dimension, 0.0);
	vector<double> Xc(dimension, 0.0);
	//Значение функций в этих точках
	while (S > exp) {

		double FuncXh = function_value(Xh);
		double FuncXg = function_value(Xg);
		double FuncXl = function_value(Xl);
	
		//Вычисление центра тяжести 
		for (int i = 1; i < vertices;i++) {
			for (int j = 0; j < dimension; j++) {
				center_of_mass[j] += Simplex[i][j];
			}
		}
		for (int j = 0; j < dimension; j++) {
			center_of_mass[j] /= (vertices - 1);
		}
		//Вычисление функции в точке центра массы и точки Xr
		double Func_of_center = function_value(center_of_mass);
		double FuncXr = function_value(Xr);
		//Операция отражения точки Xh относительно центра масс
		for (int j = 0; j < dimension; j++) {
			Xr[j] = (1 - alpha) * center_of_mass[j] - alpha * Xh[j];
		}
		//Блок условий
		if (FuncXr < FuncXl) {
			for (int j = 0; j < dimension; j++) {
				Xe[j] = gamma * Xr[j] + (1 - gamma) * center_of_mass[j];
			}
			double FuncXe = function_value(Xe);
			if (FuncXe < FuncXl) {
				for (int j = 0; j < dimension; j++) {
					Simplex[0][j] = Xe[j];
				}
				// Переходи к проверке условия
				for (int i = 0; i < dimension;i++) {
					F += function_value(Simplex[i]);
				}
				F = F / (n + 1);
				for (int i = 0; i < dimension; i++) {
					Sr += function_value(Simplex[i]) - F;
				}
				Sr = Sr * Sr / (n + 1);
				S = sqrt(S);
			}
			else {
				for (int j = 0; j < dimension; j++) {
					Simplex[0][j] = Xr[j];
				}
				//Переход к проверке условия
				for (int i = 0; i < dimension;i++) {
					F += function_value(Simplex[i]);
				}
				F = F / (n + 1);
				for (int i = 0; i < dimension; i++) {
					Sr += function_value(Simplex[i]) - F;
				}
				Sr = Sr * Sr / (n + 1);
				S = sqrt(S);
			}
		}
		else {
			if (FuncXr > FuncXg) {
				for (int j = 0; j < dimension; j++) {
					Simplex[0][j] = Xr[j];
				}
				//Переход к проверке условия
				for (int i = 0; i < dimension;i++) {
					F += function_value(Simplex[i]);
				}
				F = F / (n + 1);
				for (int i = 0; i < dimension; i++) {
					S += function_value(Simplex[i]) - F;
				}
				Sr = Sr * Sr / (n + 1);
				S = sqrt(Sr);
			}
			else {
				;
			}
		}
		if (FuncXr < FuncXh) {
			for (int j = 0; j < dimension; j++) {
				Xr[j] = Simplex[0][j];
			}
			double prom = FuncXh;
			FuncXh = FuncXr;
			FuncXr = prom;
		}
		if (FuncXr >= FuncXh) {
			for (int j = 0; j < dimension; j++) {
				Xc[j] = betha*Xh[j] + (1 - betha)* center_of_mass[j];
			}
		}
		double FuncXc = function_value(Xc);
		if (FuncXc < FuncXh) {
			for (int j = 0; j < dimension; j++) {
				Xh[j] = Xc[j];
			}
		}
		else {
			for (int i = 0; i < vertices - 1; i++) {
				for (int j = 0; j < dimension; j++) {
					Simplex[i][j] = Simplex[i][j] + 0.5 * (Simplex[i][j] - Xl[j]);
				}
			}
		}
		//Проверка условия выхода
		for (int i = 0; i < dimension;i++) {
			F += function_value(Simplex[i]);
		}
		F = F / (n + 1);
		for (int i = 0; i < dimension; i++) {
			Sr += pow(function_value(Simplex[i]) - F,2);
		}
		Sr = Sr  / (n + 1);
		S = sqrt(Sr);
		//Очистка точек: центр масс, Xr, Xc, Xe
		
		Sr = 0;
		F = 0;
		
		
	}

	return Simplex;
}



int main() {
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);
	srand(static_cast<unsigned int>(time(0)));
	int n = 2;
	int alpha = 1;
	int betha = 0.5;
	int gamma = 2.0;
	vector<double> min;
	vector<vector<double>> simplex = CreateSimplex(n);
	printSimplex(simplex);
	simplex = Method_Neldera_Mida(simplex, alpha, betha, gamma, n);
	cout << "------------" << endl;
	sort(simplex.begin(), simplex.end(), [](const vector<double>& a, const vector<double>& b) {return function_value(a) > function_value(b); });
	printSimplex(simplex);
	cout << "Координаты точки минимума" << endl;
	min = simplex[simplex.size() - 1];
	for (int j = 0; j < n; j++) {
		cout << min[j] << endl;
	}

}