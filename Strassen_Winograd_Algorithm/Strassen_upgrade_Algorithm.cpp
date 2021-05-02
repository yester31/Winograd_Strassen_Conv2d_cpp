#include <stdio.h>
#include <stdlib.h>
#include <vector> //vector�� �������� ���̺귯��
#include <math.h> //log,floor ���
#include <time.h> //clock()
#include <io.h>
#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>

using namespace std;
using namespace chrono;

/***************************************************************************
	Strassen Algorithm

	��ó :: http://yimoyimo.tk/Strassen/
*****************************************************************************/

// ���̳ʸ� ������ �о�� ���� 2���� �迭�� �����Ű�� �Լ�.
// parameter : ���� ���� �����ų vector �ּ� , ���ϸ�, , �������, �������� 
void MatrixInit(vector<int>& matrix, int Row, int Col)
{
	int idx = 1;

	for (int i = 0; i < Row; i++)
	{
		int temp = i * Col;
		for (int j = 0; j < Col; j++)
		{
			int gidx = temp + j;
			matrix[gidx] = idx;
			idx++;
		}
	}
}
void MatrixInitone(vector<int>& matrix, int Row, int Col)
{
	for (int i = 0; i < Row; i++)
	{
		int temp = i * Col;
		for (int j = 0; j < Col; j++)
		{
			int gidx = temp + j;
			matrix[gidx] = 1;
		}
	}
}


// ��� A�� B�� ���Ͽ� C�� �����Ű�� �Լ�(A��� B��� C����� ���� ũ��)
void MatrixSum(vector<int>& matrixA, vector<int>& matrixB, vector<int>& matrixC, int Rows, int Cols)
{
	for (int i = 0; i < Rows; i++)
	{
		int temp = i * Cols;
		for (int j = 0; j < Cols; j++)
		{
			int gidx = temp + j;
			matrixC[gidx] = matrixA[gidx] + matrixB[gidx];
		}
	}
}

// ��� A�� B�� ���� C�� �����Ű�� �Լ�(A��� B��� C����� ���� ũ��)
void MatrixSub(vector<int>& matrixA, vector<int>& matrixB, vector<int>& matrixC, int Rows, int Cols)
{
	for (int i = 0; i < Rows; i++)
	{
		int temp = i * Cols;
		for (int j = 0; j < Cols; j++) 
		{
			int gidx = temp + j;
			matrixC[gidx] = matrixA[gidx] - matrixB[gidx];
		}
	}


}

// ��� A�� B�� ���Ͽ� C�� �����Ű�� �Լ�(A��� B��� C����� ���� ũ��)
// parameter : ��� A, B, C 

void MatrixMul(vector<int> Prev_matrix, vector<int> Post_matrix, vector<int>& Y_output, int prev_rows, int prev_colsAndPost_rows, int post_cols)
{

	for (int i = 0; i < prev_rows; ++i) {
		int temp1 = i * prev_colsAndPost_rows;
		int temp2 = i * post_cols;
		for (int j = 0; j < post_cols; ++j)
		{
			int sum = 0;
			for (int k = 0; k < prev_colsAndPost_rows; ++k)
			{
				sum += Prev_matrix[temp1 + k] * Post_matrix[k * post_cols + j];
			}
			Y_output[temp2 + j] = (int)sum;
		}
	}
}



// �Ӱ谪 ���ϴ� �Լ�
// parameter : ����� ���� �� ���� ���� (ex:1024)
// return : �Ӱ谪 ��ȯ
int getThreshold(int n)
{
	int th;
	double k = floor(log(n) / log(2) - 6);
	th = (int)floor(n / pow(2.0, k)) + 1;
	return th;
}

int getThreshold2(int n)
{
	int th;
	double k = floor(log(n) / log(2) - 4);// lgn -4
	th = (int)floor(n / pow(2.0, k)) + 1;
	return th;
}

// 4���� �κ���ķ� ������ �Լ�
// parameter : ���� ���, ������ ��� ���� 4��
void Submatrix(vector<int>& matrixOrigin, vector<int>& matrix11, vector<int>& matrix12, vector<int>& matrix21, vector<int>& matrix22, int Rows, int Cols)
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
			matrix11[gidx2] = matrixOrigin[gidx];									//�� ������
			matrix12[gidx2] = matrixOrigin[gidx + Cols];							//�� ������
			matrix21[gidx2] = matrixOrigin[Rows * Cols * 2 + gidx];				//�� �ϴ����
			matrix22[gidx2] = matrixOrigin[Rows * Cols * 2 + gidx + Cols];		//�� �ϴ����
		}
	}

}

// 4���� �κ���ĵ��� ����� ���ִ� �Լ�
// parameter : ��ģ ����� ������ ��� , �κ���� 11, 12, 21, 22
void Mergematrix(vector<int>& matrixOrigin, vector<int>& matrix11, vector<int>& matrix12, vector<int>& matrix21, vector<int>& matrix22, int Rows, int Cols)
{
	for (int i = 0; i < Rows; i++)
	{
		int temp = i * (Cols * 2);
		int temp2 = i * Cols;
		for (int j = 0; j < Cols; j++)
		{
			int gidx = temp + j;
			int gidx2 = temp2 + j;
			matrixOrigin[gidx] = matrix11[gidx2];										//�� ������
			matrixOrigin[gidx + Cols] = matrix12[gidx2];								//�� ������
			matrixOrigin[Rows * Cols * 2 + gidx] = matrix21[gidx2];         	    //�� �ϴ����
			matrixOrigin[Rows * Cols * 2 + gidx + Cols] = matrix22[gidx2];	    //�� �ϴ����
		}
	}
}

void ElementMul(vector<int>& Y_output, vector<int> Prev_matrix, vector<int> Post_matrix) {
	for (int i = 0; i < 16; i++)
		Y_output[i] = Prev_matrix[i] * Post_matrix[i];

}


// ��Ʈ�� �˰��� �Լ�
void Strassen(vector<int>& matrixU, vector<int>& matrixV, vector<int>& matrixM, int U_Row, int U_Col, int V_Row, int V_Col)
{

	if (1 == U_Row % 2 || 1 == U_Col % 2 || 1 == V_Row % 2 || 1 == V_Col % 2)
	{
		//cout << "U_Row : " << U_Row << " , " << "U_Col : " << U_Col << " , " << "V_Row : " << V_Row << " , " << "V_Col : " << V_Col << " )" << endl;
		//cout << "Can't divide the matrix anymore" << endl;
		//ElementMul(matrixC, matrixA, matrixB);
		MatrixMul(matrixU, matrixV, matrixM, U_Row, U_Col, V_Col);
		return;
	}else {
		if (V_Col <= getThreshold(V_Col)) {
			MatrixMul(matrixU, matrixV, matrixM, U_Row, U_Col, V_Col);
			return;
		}
		else {
		
			//cout << " ss" << endl;
			int newU_Row = U_Row / 2;					//4����� �ϱ� ����
			int newU_Col = U_Col / 2;
			int newV_Row = V_Row / 2;
			int newV_Col = V_Col / 2;

			//a11~a22 �κ����, b11~b22 �κ���� 
			vector<int> a11(newU_Row * newU_Col), a12(newU_Row * newU_Col), a21(newU_Row * newU_Col), a22(newU_Row * newU_Col);
			vector<int> b11(newV_Row * newV_Col), b12(newV_Row * newV_Col), b21(newV_Row * newV_Col), b22(newV_Row * newV_Col);

			//�κ���ĵ��� �������� m1~m7 ������
			vector<int>  m1(newU_Row * newV_Col), m2(newU_Row * newV_Col), m3(newU_Row * newV_Col), m4(newU_Row * newV_Col), m5(newU_Row * newV_Col), m6(newU_Row * newV_Col), m7(newU_Row * newV_Col);

			//a11~b22 �� ���������� �ӽ÷� ������ �׸�
			vector<int>  tempA(newU_Row * newU_Col), tempB(newV_Row * newV_Col);

			vector<int>  tempAc(newU_Row * newV_Col), tempBc(newU_Row * newV_Col);

			// m1~m7 ���� ����� C�� ���ϱ� ���� ���� �� ���
			vector<int>  c11(newU_Row * newV_Col), c12(newU_Row * newV_Col), c21(newU_Row * newV_Col), c22(newU_Row * newV_Col);


			//A�� �κ���� 4��, B�� �κ���� 4�� ����
			Submatrix(matrixU, a11, a12, a21, a22, newU_Row, newU_Col);
			Submatrix(matrixV, b11, b12, b21, b22, newV_Row, newV_Col);

			MatrixSum(a11, a22, tempA, newU_Row, newU_Col);				// a11+a22
			MatrixSum(b11, b22, tempB, newV_Row, newV_Col);		       // b11+b22
			Strassen(tempA, tempB, m1, newU_Row, newU_Col, newV_Row, newV_Col);    // m1=(a11+a11)(b11+b22)

			MatrixSum(a21, a22, tempA, newU_Row, newU_Col);            // a21+a22
			Strassen(tempA, b11, m2, newU_Row, newU_Col, newV_Row, newV_Col);      // m2=(a21+a22)b11

			MatrixSub(b12, b22, tempB, newV_Row, newV_Col);            // b12-b22
			Strassen(a11, tempB, m3, newU_Row, newU_Col, newV_Row, newV_Col);      // m3=a11(b12-b22)

			MatrixSub(b21, b11, tempB, newV_Row, newV_Col);            // b21-b11
			Strassen(a22, tempB, m4, newU_Row, newU_Col, newV_Row, newV_Col);      // m4=a22(b21-11)

			MatrixSum(a11, a12, tempA, newU_Row, newU_Col);            //  a11+a12
			Strassen(tempA, b22, m5, newU_Row, newU_Col, newV_Row, newV_Col); 	   // m5=(a11+a12)b22

			MatrixSub(a21, a11, tempA, newU_Row, newU_Col);            // a21-a11
			MatrixSum(b11, b12, tempB, newV_Row, newV_Col);            // b11+b12
			Strassen(tempA, tempB, m6, newU_Row, newU_Col, newV_Row, newV_Col);    // m6=(a21-a11)(b11+b12)

			MatrixSub(a12, a22, tempA, newU_Row, newU_Col);            // a12-a22
			MatrixSum(b21, b22, tempB, newV_Row, newV_Col);            // b21+b22
			Strassen(tempA, tempB, m7, newU_Row, newU_Col, newV_Row, newV_Col);    // m7 = (a12 - a22)(a12 - a22)


			// ������ ���� m1~m7 �����  c11 ~ c22 �� �����.
			MatrixSum(m1, m4, tempAc, newU_Row, newV_Col); //m1 + m4
			MatrixSum(tempAc, m7, tempBc, newU_Row, newV_Col); //m1 + m4 + m7
			MatrixSub(tempBc, m5, c11, newU_Row, newV_Col); //c11 = m1 + m4 - m5 + m7

			MatrixSum(m3, m5, c12, newU_Row, newV_Col); //c12 = m3 + m5

			MatrixSum(m2, m4, c21, newU_Row, newV_Col); //c21 = m2 + m4

			MatrixSum(m1, m3, tempAc, newU_Row, newV_Col); //m1 + m3
			MatrixSum(tempAc, m6, tempBc, newU_Row, newV_Col); //m1 + m3 + m6
			MatrixSub(tempBc, m2, c22, newU_Row, newV_Col); //c22 = m1 + m3 - m2 + m6

			//�� ����
			Mergematrix(matrixM, c11, c12, c21, c22, newU_Row, newV_Col);
		}
	}

	/*

	if (1 == U_Row % 2 || 1 == U_Col % 2 || 1 == V_Row % 2 || 1 == V_Col % 2)
	{

		//if (1 == U_Row % 2 || 1 == U_Col % 2 || 1 == V_Row % 2 || 1 == V_Col % 2)
		//cout << "U_Row : " << U_Row << " , " << "U_Col : " << U_Col << " , " << "V_Row : " << V_Row << " , " << "V_Col : " << V_Col << " )" << endl;
		//cout << "Can't divide the matrix anymore" << endl;
		//ElementMul(matrixC, matrixA, matrixB);
		MatrixMul(matrixU, matrixV, matrixM, U_Row, U_Col, V_Col);
		return;
	}
	else
	{
		//cout << " ss" << endl;
		int newU_Row = U_Row / 2;					//4����� �ϱ� ����
		int newU_Col = U_Col / 2;
		int newV_Row = V_Row / 2;
		int newV_Col = V_Col / 2;

		//a11~a22 �κ����, b11~b22 �κ���� 
		vector<int> a11(newU_Row * newU_Col), a12(newU_Row * newU_Col), a21(newU_Row * newU_Col), a22(newU_Row * newU_Col);
		vector<int> b11(newV_Row * newV_Col), b12(newV_Row * newV_Col), b21(newV_Row * newV_Col), b22(newV_Row * newV_Col);

		//�κ���ĵ��� �������� m1~m7 ������
		vector<int>  m1(newU_Row * newV_Col), m2(newU_Row * newV_Col), m3(newU_Row * newV_Col), m4(newU_Row * newV_Col), m5(newU_Row * newV_Col), m6(newU_Row * newV_Col), m7(newU_Row * newV_Col);

		//a11~b22 �� ���������� �ӽ÷� ������ �׸�
		vector<int>  tempA(newU_Row * newU_Col), tempB(newV_Row * newV_Col);

		vector<int>  tempAc(newU_Row * newV_Col), tempBc(newU_Row * newV_Col);

		// m1~m7 ���� ����� C�� ���ϱ� ���� ���� �� ���
		vector<int>  c11(newU_Row * newV_Col), c12(newU_Row * newV_Col), c21(newU_Row * newV_Col), c22(newU_Row * newV_Col);


		//A�� �κ���� 4��, B�� �κ���� 4�� ����
		Submatrix(matrixU, a11, a12, a21, a22, newU_Row, newU_Col);
		Submatrix(matrixV, b11, b12, b21, b22, newV_Row, newV_Col);

		MatrixSum(a11, a22, tempA, newU_Row, newU_Col);				// a11+a22
		MatrixSum(b11, b22, tempB, newV_Row, newV_Col);		       // b11+b22
		Strassen(tempA, tempB, m1, newU_Row, newU_Col, newV_Row, newV_Col);    // m1=(a11+a11)(b11+b22)

		MatrixSum(a21, a22, tempA, newU_Row, newU_Col);            // a21+a22
		Strassen(tempA, b11, m2, newU_Row, newU_Col, newV_Row, newV_Col);      // m2=(a21+a22)b11

		MatrixSub(b12, b22, tempB, newV_Row, newV_Col);            // b12-b22
		Strassen(a11, tempB, m3, newU_Row, newU_Col, newV_Row, newV_Col);      // m3=a11(b12-b22)

		MatrixSub(b21, b11, tempB, newV_Row, newV_Col);            // b21-b11
		Strassen(a22, tempB, m4, newU_Row, newU_Col, newV_Row, newV_Col);      // m4=a22(b21-11)

		MatrixSum(a11, a12, tempA, newU_Row, newU_Col);            //  a11+a12
		Strassen(tempA, b22, m5, newU_Row, newU_Col, newV_Row, newV_Col); 	   // m5=(a11+a12)b22

		MatrixSub(a21, a11, tempA, newU_Row, newU_Col);            // a21-a11
		MatrixSum(b11, b12, tempB, newV_Row, newV_Col);            // b11+b12
		Strassen(tempA, tempB, m6, newU_Row, newU_Col, newV_Row, newV_Col);    // m6=(a21-a11)(b11+b12)

		MatrixSub(a12, a22, tempA, newU_Row, newU_Col);            // a12-a22
		MatrixSum(b21, b22, tempB, newV_Row, newV_Col);            // b21+b22
		Strassen(tempA, tempB, m7, newU_Row, newU_Col, newV_Row, newV_Col);    // m7 = (a12 - a22)(a12 - a22)


		// ������ ���� m1~m7 �����  c11 ~ c22 �� �����.
		MatrixSum(m1, m4, tempAc, newU_Row, newV_Col); //m1 + m4
		MatrixSum(tempAc, m7, tempBc, newU_Row, newV_Col); //m1 + m4 + m7
		MatrixSub(tempBc, m5, c11, newU_Row, newV_Col); //c11 = m1 + m4 - m5 + m7

		MatrixSum(m3, m5, c12, newU_Row, newV_Col); //c12 = m3 + m5

		MatrixSum(m2, m4, c21, newU_Row, newV_Col); //c21 = m2 + m4

		MatrixSum(m1, m3, tempAc, newU_Row, newV_Col); //m1 + m3
		MatrixSum(tempAc, m6, tempBc, newU_Row, newV_Col); //m1 + m3 + m6
		MatrixSub(tempBc, m2, c22, newU_Row, newV_Col); //c22 = m1 + m3 - m2 + m6

		//�� ����
		Mergematrix(matrixM, c11, c12, c21, c22, newU_Row, newV_Col);
	}*/

}


void FileOutput(vector<int>& matrix, int row, int col)
{
	
	for (int i = 0; i < row; i++)
	{
		int temp = col * i;
		for (int j = 0; j < col; j++)
		{
			int gidx = temp + j;
			cout << setw(5) << matrix[gidx];
		}
		cout << endl;
	}
	cout << endl;

}


int main()
{
	int U_Row = 512;
	int U_Col = 512;
	int V_Row = 512;
	int V_Col = 512;

	cout << "U matrix ( " << U_Row << ", " << U_Col << " ) " << endl;
	cout << "V matrix ( " << V_Row << ", " << V_Col << " ) " << endl;

	vector <int> U(U_Row * U_Col); // 
	vector <int> V(V_Row * V_Col); // 

	vector <int> M(U_Row * V_Col); //
	vector <int> Mg(U_Row * V_Col); //

	MatrixInit(U, U_Row, U_Col);
	//MatrixInit(V, V_Row, V_Col);
	//MatrixInitone(U, U_Row, U_Col);
	MatrixInitone(V, V_Row, V_Col);

	//FileOutput(U, U_Row, U_Col);
	//FileOutput(V, V_Row, V_Col);
	
	long long start_usec3 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	MatrixMul(U, V, Mg, U_Row, U_Col, V_Col);
	long long end_usec3 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	int frame_sec3 = int(end_usec3 - start_usec3);
	cout << frame_sec3 << "u sec (MatrixMul) v1" << endl;
	//FileOutput(Mg, U_Row, V_Col);
	
	long long start_usec2 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	Strassen(U, V, M, U_Row, U_Col, V_Row, V_Col);
	long long end_usec2 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	int frame_sec2 = int(end_usec2 - start_usec2);
	cout << frame_sec2 << "u sec (Strassen) v1" << endl;
	
	//FileOutput(M, U_Row, V_Col);
	//cout << endl << endl;
	//FileOutput(D1);

	return 0;
}