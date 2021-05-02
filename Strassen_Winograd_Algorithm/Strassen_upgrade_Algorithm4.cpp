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
void MatrixInit(vector<float>& matrix, int Row, int Col)
{
	float idx = 1.f;

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
void MatrixInitone(vector<float>& matrix, int Row, int Col)
{
	for (int i = 0; i < Row; i++)
	{
		int temp = i * Col;
		for (int j = 0; j < Col; j++)
		{
			int gidx = temp + j;
			matrix[gidx] = 1.f;
		}
	}
}


// ��� A�� B�� ���Ͽ� C�� �����Ű�� �Լ�(A��� B��� C����� ���� ũ��)
void MatrixSum(vector<float>& matrixA, vector<float>& matrixB, vector<float>& matrixC, int Rows, int Cols)
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
void MatrixSub(vector<float>& matrixA, vector<float>& matrixB, vector<float>& matrixC, int Rows, int Cols)
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

void MatrixMul(vector<float> Prev_matrix, vector<float> Post_matrix, vector<float>& Y_output, int prev_rows, int prev_colsAndPost_rows, int post_cols)
{

	for (int i = 0; i < prev_rows; ++i) {
		int temp1 = i * prev_colsAndPost_rows;
		int temp2 = i * post_cols;
		for (int j = 0; j < post_cols; ++j)
		{
			float sum = 0;
			for (int k = 0; k < prev_colsAndPost_rows; ++k)
			{
				sum += Prev_matrix[temp1 + k] * Post_matrix[k * post_cols + j];
			}
			Y_output[temp2 + j] = sum;
		}
	}
}

void MatrixMulElement2X2(vector<float> matrix_U, vector<float> matrix_V, vector<float>& matrix_M, int U_rows, int UcolsVrows, int V_cols)
{

	for (int i = 0; i < U_rows / 2; ++i) { // ��

		//int temp1 = i * UcolsVrows ;
		//int temp2 = i * V_cols ;

		for (int j = 0; j < V_cols / 2; ++j) // �� 
		{
			int m1 = 0;
			int	m2 = 0;
			int	m3 = 0;
			int	m4 = 0;

			for (int k = 0; k < UcolsVrows / 2; ++k)
			{
				m1 += matrix_U[i * UcolsVrows * 2 + k * 2] * matrix_V[k * V_cols * 2 + j * 2];
				m2 += matrix_U[i * UcolsVrows * 2 + k * 2 + 1] * matrix_V[k * V_cols * 2 + j * 2 + 1];
				m3 += matrix_U[i * UcolsVrows * 2 + k * 2 + UcolsVrows] * matrix_V[k * V_cols * 2 + j * 2 + V_cols];
				m4 += matrix_U[i * UcolsVrows * 2 + k * 2 + 1 + UcolsVrows] * matrix_V[k * V_cols * 2 + j * 2 + 1 + V_cols];
			}
			matrix_M[i * V_cols * 2 + j * 2] = m1;
			matrix_M[i * V_cols * 2 + j * 2 + 1] = m2;
			matrix_M[i * V_cols * 2 + j * 2 + V_cols] = m3;
			matrix_M[i * V_cols * 2 + j * 2 + V_cols + 1] = m4;
		}
	}
}

void MatrixMulElement4X4(vector<float> matrix_U, vector<float> matrix_V, vector<float>& matrix_M, int U_rows, int UcolsVrows, int V_cols)
{

	for (int i = 0; i < U_rows / 4; ++i) { // ��

		//int temp1 = i * UcolsVrows ;
		//int temp2 = i * V_cols ;

		for (int j = 0; j < V_cols / 4; ++j) // �� 
		{
			int m1 = 0;int	m2 = 0;int	m3 = 0;int	m4 = 0;int m5 = 0;int	m6 = 0;int	m7 = 0;int	m8 = 0;
			int m9 = 0;int	m10 = 0;int	m11 = 0;int	m12 = 0;int m13 = 0;int	m14 = 0;int	m15 = 0;int	m16 = 0;

			for (int k = 0; k < UcolsVrows / 4; ++k)
			{
				int midx1 = i * UcolsVrows * 4 + k * 4;
				int midx2 = k * V_cols * 4 + j * 4;
				m1 += matrix_U[midx1] * matrix_V[midx2];
				m2 += matrix_U[midx1 + 1] * matrix_V[midx2 + 1];
				m3 += matrix_U[midx1 + 2] * matrix_V[midx2 + 2];
				m4 += matrix_U[midx1 + 3] * matrix_V[midx2 + 3];

				int midx3 = midx1 + UcolsVrows;
				int midx4 = midx2 + V_cols;
				m5 += matrix_U[midx3] * matrix_V[midx4];
				m6 += matrix_U[midx3 + 1] * matrix_V[midx4 + 1];
				m7 += matrix_U[midx3 + 2] * matrix_V[midx4 + 2];
				m8 += matrix_U[midx3 + 3] * matrix_V[midx4 + 3];

				int midx5 = midx3 + UcolsVrows;
				int midx6 = midx4 + V_cols;
				m9 += matrix_U[midx5] * matrix_V[midx6];
				m10 += matrix_U[midx5 + 1] * matrix_V[midx6 + 1];
				m11 += matrix_U[midx5 + 2] * matrix_V[midx6 + 2];
				m12 += matrix_U[midx5 + 3] * matrix_V[midx6 + 3];

				int midx7 = midx5 + UcolsVrows;
				int midx8 = midx6 + V_cols;
				m13 += matrix_U[midx7] * matrix_V[midx8];
				m14 += matrix_U[midx7 + 1] * matrix_V[midx8 + 1];
				m15 += matrix_U[midx7 + 2] * matrix_V[midx8 + 2];
				m16 += matrix_U[midx7 + 3] * matrix_V[midx8 + 3];

			}
			int idx1 = i * V_cols * 4 + j * 4;
			matrix_M[idx1] = m1;
			matrix_M[idx1 + 1] = m2;
			matrix_M[idx1 + 2] = m3;
			matrix_M[idx1 + 3] = m4;

			int idx2 = idx1 + V_cols;
			matrix_M[idx2] = m5;
			matrix_M[idx2 + 1] = m6;
			matrix_M[idx2 + 2] = m7;
			matrix_M[idx2 + 3] = m8;

			int idx3 = idx2 + V_cols;
			matrix_M[idx3] = m9;
			matrix_M[idx3 + 1] = m10;
			matrix_M[idx3 + 2] = m11;
			matrix_M[idx3 + 3] = m12;

			int idx4 = idx3 + V_cols;
			matrix_M[idx4] = m13;
			matrix_M[idx4 + 1] = m14;
			matrix_M[idx4 + 2] = m15;
			matrix_M[idx4 + 3] = m16;
		}
	}
}

// �Ӱ谪 ���ϴ� �Լ�
int getThreshold(int n)
{
	int th;
	double k = floor(log(n) / log(2) - 4);
	th = (int)floor(n / pow(2.0, k)) + 1;
	return th;
}

// 4���� �κ���ķ� ������ �Լ�
// parameter : ���� ���, ������ ��� ���� 4��
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
			matrix11[gidx2] = matrixOrigin[gidx];									//�� ������
			matrix12[gidx2] = matrixOrigin[gidx + Cols];							//�� ������
			matrix21[gidx2] = matrixOrigin[Rows * Cols * 2 + gidx];					//�� �ϴ����
			matrix22[gidx2] = matrixOrigin[Rows * Cols * 2 + gidx + Cols];			//�� �ϴ����
		}
	}

}

// 4���� �κ���ĵ��� ����� ���ִ� �Լ�
// parameter : ��ģ ����� ������ ��� , �κ���� 11, 12, 21, 22
void Mergematrix(vector<float>& matrixOrigin, vector<float>& matrix11, vector<float>& matrix12, vector<float>& matrix21, vector<float>& matrix22, int Rows, int Cols)
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
			matrixOrigin[Rows * Cols * 2 + gidx] = matrix21[gidx2];         			//�� �ϴ����
			matrixOrigin[Rows * Cols * 2 + gidx + Cols] = matrix22[gidx2];				//�� �ϴ����
		}
	}
}

// ��Ʈ�� �˰��� �Լ�
void Strassen(vector<float>& matrixU, vector<float>& matrixV, vector<float>& matrixM, int U_Row, int U_Col, int V_Row, int V_Col)
{

	if (1 == (U_Row / 4) % 2 || 1 == (U_Col / 4) % 2 || 1 == (V_Row / 4) % 2 || 1 == (V_Col / 4) % 2)
	{
		//cout << "U_Row : " << U_Row << " , " << "U_Col : " << U_Col << " , " << "V_Row : " << V_Row << " , " << "V_Col : " << V_Col << " )" << endl;
		//cout << "Can't divide the matrix anymore" << endl;
		MatrixMulElement4X4(matrixU, matrixV, matrixM, U_Row, U_Col, V_Col);
		return;
	}
	else {
		//int th = V_Col / 4;
		if (V_Col <= getThreshold(V_Col)) {

			MatrixMulElement4X4(matrixU, matrixV, matrixM, U_Row, U_Col, V_Col);

			return;
		}
		else {

			//cout << " ss" << endl;
			int newU_Row = U_Row / 2;					//4����� �ϱ� ����
			int newU_Col = U_Col / 2;
			int newV_Row = V_Row / 2;
			int newV_Col = V_Col / 2;

			//a11~a22 �κ����, b11~b22 �κ���� 
			vector<float> a11(newU_Row * newU_Col), a12(newU_Row * newU_Col), a21(newU_Row * newU_Col), a22(newU_Row * newU_Col);
			vector<float> b11(newV_Row * newV_Col), b12(newV_Row * newV_Col), b21(newV_Row * newV_Col), b22(newV_Row * newV_Col);

			//�κ���ĵ��� �������� m1~m7 ������
			vector<float>  m1(newU_Row * newV_Col), m2(newU_Row * newV_Col), m3(newU_Row * newV_Col), m4(newU_Row * newV_Col), m5(newU_Row * newV_Col), m6(newU_Row * newV_Col), m7(newU_Row * newV_Col);

			//a11~b22 �� ���������� �ӽ÷� ������ �׸�
			vector<float>  tempA(newU_Row * newU_Col), tempB(newV_Row * newV_Col);

			vector<float>  tempAc(newU_Row * newV_Col), tempBc(newU_Row * newV_Col);

			// m1~m7 ���� ����� C�� ���ϱ� ���� ���� �� ���
			vector<float>  c11(newU_Row * newV_Col), c12(newU_Row * newV_Col), c21(newU_Row * newV_Col), c22(newU_Row * newV_Col);


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
}

void FileOutput(vector<float>& matrix, int row, int col)
{

	for (int i = 0; i < row; i++)
	{
		int temp = col * i;
		for (int j = 0; j < col; j++)
		{
			int gidx = temp + j;
			cout << setw(10) << matrix[gidx];
		}
		cout << endl;
	}
	cout << endl;

}


int main()
{
	int OC = 16;
	int IC = 8;
	int N = 8;

	cout << "U matrix ( " << OC << ", " << IC << " ) " << endl;
	cout << "V matrix ( " << IC << ", " << N << " ) " << endl;

	vector<float> U(OC * IC); // 
	vector<float> V(IC * N); // 

	vector<float> M(OC * N); //
	vector<float> Mg(OC * N); //

	//MatrixInit(U, OC, IC);
	MatrixInit(V, IC, N);
	MatrixInitone(U, OC, IC);
	//MatrixInitone(V, IC, N);

	FileOutput(U, OC, IC);
	FileOutput(V, IC, N);


	//MatrixMulElement4X4(U, V, Mg, OC, IC, N);
	//FileOutput(Mg, OC, N);


	/*
	long long start_usec3 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	MatrixMul(U, V, Mg, OC, IC, N);
	long long end_usec3 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	int frame_sec3 = int(end_usec3 - start_usec3);
	cout << frame_sec3 << "u sec (MatrixMul) v1" << endl;
	//FileOutput(Mg, OC, N);
	*/


	long long start_usec2 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	Strassen(U, V, M, OC, IC, IC, N);
	long long end_usec2 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	int frame_sec2 = int(end_usec2 - start_usec2);
	cout << frame_sec2 << "u sec (Strassen) v1" << endl;
	FileOutput(M, OC, N);

	//cout << endl << endl;
	//FileOutput(D1);

	return 0;
}