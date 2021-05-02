#include <io.h>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <time.h>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace chrono;


void Submatrix(vector<float>& matrixOrigin, vector<float>& matrix11, vector<float>& matrix12, vector<float>& matrix21, vector<float>& matrix22, int Rows, int Cols)
{
	// int Rows, int Cols �κ������ ������ 

	for (int i = 0; i < Rows; i++)
	{
		int temp = i * (Cols * 2);
		int temp2 = i * (Cols);
		for (int j = 0; j < Cols; j++)
		{
			int gidx = temp + j;
			int gidx2 = temp2 + j;
			matrix11[gidx2] = matrixOrigin[gidx];									//�� ������ // 11
			matrix12[gidx2] = matrixOrigin[gidx + Cols];							//�� ������ // 12  //
			matrix21[gidx2] = matrixOrigin[Rows * Cols * 2 + gidx];					//�� �ϴ���� // 21
			matrix22[gidx2] = matrixOrigin[Rows * Cols * 2 + gidx + Cols];			//�� �ϴ���� // 22  //
		}
	}

}

int main() {

	vector<float> test =
	{   1, 1, 2, 2, 1, 1, 2, 2,
		1, 1, 2, 2, 1, 1, 2, 2,
		1, 1, 2, 2, 1, 1, 2, 2,
		3, 3, 4, 4, 3, 3, 4, 4
	};
		vector<float> test1(8);
		vector<float> test2(8);
		vector<float> test3(8);
		vector<float> test4(8);


		Submatrix(test, test1, test2, test3, test4, 2, 4);


		cout << test2[0];
		cout << test2[1];
		cout << test2[2];
		cout << test2[3];
		cout <<endl;


	return 0;

}
