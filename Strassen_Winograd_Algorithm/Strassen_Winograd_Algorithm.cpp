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

/***************************************************************************
	Strassen-Winograd algotirhm
*****************************************************************************/


void filterTransform2(vector<float>& Output_U, vector<float> filter, int out_ch, int in_ch) {

	int temp1u = in_ch * 4;
	int temp1f = in_ch * 3 * 3;
	for (int oc_idx = 0; oc_idx < out_ch; oc_idx++)
	{
		int temp2u = oc_idx * temp1u * 4;
		int temp2f = oc_idx * temp1f;
		for (int ic_idx = 0; ic_idx < in_ch; ic_idx++)
		{
			int u_idx = ic_idx * 4 + temp2u;
			int f_idx = ic_idx * 3 * 3 + temp2f;

			float* g1 = &filter[f_idx];
			float* g2 = &filter[f_idx + 1];
			float* g3 = &filter[f_idx + 2];
			float* g4 = &filter[f_idx + 3];
			float* g5 = &filter[f_idx + 4];
			float* g6 = &filter[f_idx + 5];
			float* g7 = &filter[f_idx + 6];
			float* g8 = &filter[f_idx + 7];
			float* g9 = &filter[f_idx + 8];

			Output_U[u_idx] = *g1;
			Output_U[u_idx + 1] = (*g1 + *g2 + *g3) / 2.f;
			Output_U[u_idx + 2] = (*g1 - *g2 + *g3) / 2.f;
			Output_U[u_idx + 3] = *g3;

			float temp1 = *g1 + *g4 + *g7;
			float temp2 = *g2 + *g5 + *g8;
			float temp3 = *g3 + *g6 + *g9;

			Output_U[u_idx + temp1u] = (temp1) / 2.f;
			Output_U[u_idx + 1 + temp1u] = (temp1 + temp2 + temp3) / 4.f;
			Output_U[u_idx + 2 + temp1u] = (temp1 - temp2 + temp3) / 4.f;
			Output_U[u_idx + 3 + temp1u] = (temp3) / 2.f;

			float temp4 = *g1 - *g4 + *g7;
			float temp5 = *g2 - *g5 + *g8;
			float temp6 = *g3 - *g6 + *g9;

			Output_U[u_idx + temp1u * 2] = (temp4) / 2.f;
			Output_U[u_idx + 1 + temp1u * 2] = (temp4 + temp5 + temp6) / 4.f;
			Output_U[u_idx + 2 + temp1u * 2] = (temp4 - temp5 + temp6) / 4.f;
			Output_U[u_idx + 3 + temp1u * 2] = (temp6) / 2.f;

			Output_U[u_idx + temp1u * 3] = *g7;
			Output_U[u_idx + 1 + temp1u * 3] = (*g7 + *g8 + *g9) / 2.f;
			Output_U[u_idx + 2 + temp1u * 3] = (*g7 - *g8 + *g9) / 2.f;
			Output_U[u_idx + 3 + temp1u * 3] = *g9;


		}
	}
}

void inputTransform2(vector<float>& V, vector<float>input, int input_n, int input_c, int input_h, int input_w) {

	int output_H = input_h - 2;
	int	output_W = input_w - 2;
	int area_o = output_H * output_W;
	//int tile_N = (input_h - 2) / 2 * (input_w - 2) / 2;
	// int temp_v = 4 * area_o;// 16 * (input_h - 2)/2 * (input_w - 2)/2
	// int tile_N = area_o / 4;
	int tile_h = (input_h - 2) / 2;
	int tile_w = (input_w - 2) / 2;
	int temp_i = input_c * input_h * input_w;
	int temp_n = input_n * 4;
	for (int n = 0; n < input_n; n++) // N
	{
		int temp_i2 = n * temp_i;

		for (int inch = 0; inch < input_c; inch++) // IC
		{
			int temp_i3 = inch * input_h * input_w + temp_i2;

			int tilecount = 0;
			for (int row = 0; row < output_H; row += 2)  // T_N
			{
				tilecount = tile_w * row / 2;

				int row_idx1 = row * input_w;
				int row_idx2 = row_idx1 + input_w;
				int row_idx3 = row_idx2 + input_w;
				int row_idx4 = row_idx3 + input_w;

				for (int col = 0; col < output_W; col += 2)
				{
					int gidx = tilecount * input_c * input_n * 16 + inch * input_n * 16 + n * 4;
					tilecount ++;

					int t_idx1 = temp_i3 + col + row_idx1;
					int t_idx2 = temp_i3 + col + row_idx2;
					int t_idx3 = temp_i3 + col + row_idx3;
					int t_idx4 = temp_i3 + col + row_idx4;

					float* d1 = &input[t_idx1];
					float* d2 = &input[t_idx1 + 1];
					float* d3 = &input[t_idx1 + 2];
					float* d4 = &input[t_idx1 + 3];

					float* d5 = &input[t_idx2];
					float* d6 = &input[t_idx2 + 1];
					float* d7 = &input[t_idx2 + 2];
					float* d8 = &input[t_idx2 + 3];

					float* d9 = &input[t_idx3];
					float* d10 = &input[t_idx3 + 1];
					float* d11 = &input[t_idx3 + 2];
					float* d12 = &input[t_idx3 + 3];

					float* d13 = &input[t_idx4];
					float* d14 = &input[t_idx4 + 1];
					float* d15 = &input[t_idx4 + 2];
					float* d16 = &input[t_idx4 + 3];

					float dd1 = *d11 - (*d3);
					float dd2 = *d2 - (*d10);
					float dd3 = *d7 + (*d11);
					float dd4 = *d6 + (*d10);
					float dd5 = *d7 - (*d11);
					float dd6 = *d10 - (*d6);
					float dd7 = *d15 - (*d7);
					float dd8 = *d6 - (*d14);

					V[gidx] = *d1 - *d9 + dd1;
					V[gidx + 1] = dd2 - dd1;
					V[gidx + 2] = -dd1 - dd2;
					V[gidx + 3] = dd2 - *d4 + *d12;

					int temp_v4 = gidx + temp_n;
					V[temp_v4] = *d5 + *d9 - dd3;
					V[temp_v4 + 1] = dd4 + dd3;
					V[temp_v4 + 2] = dd3 - dd4;
					V[temp_v4 + 3] = dd4 - *d8 - *d12;

					int temp_v5 = temp_v4 + temp_n;
					V[temp_v5] = *d9 - *d5 + dd5;
					V[temp_v5 + 1] = dd6 - dd5;
					V[temp_v5 + 2] = -(dd6 + dd5);
					V[temp_v5 + 3] = dd6 + *d8 - *d12;

					int temp_v6 = temp_v5 + temp_n;
					V[temp_v6] = *d5 - *d13 + dd7;
					V[temp_v6 + 1] = dd8 - dd7;
					V[temp_v6 + 2] = -dd7 - dd8;
					V[temp_v6 + 3] = dd8 - *d8 + *d16;

					/*
					cout << setw(5) << V[gidx] << " " << setw(5) << V[gidx + 1] << " " << setw(5) << V[gidx + 2] << " " << setw(5) << V[gidx + 3];
					cout << setw(5) << V[temp_v4] << " " << setw(5) << V[temp_v4 + 1] << " " << setw(5) << V[temp_v4 + 2] << " " << setw(5) << V[temp_v4 + 3];
					cout << setw(5) << V[temp_v5] << " " << setw(5) << V[temp_v5 + 1] << " " << setw(5) << V[temp_v5 + 2] << " " << setw(5) << V[temp_v5 + 3];
					cout << setw(5) << V[temp_v6] << " " << setw(5) << V[temp_v6 + 1] << " " << setw(5) << V[temp_v6 + 2] << " " << setw(5) << V[temp_v6 + 3];
					cout << endl;
					*/
				}
			}
		}
	}
}


void valueCheck(vector<float> &valueCheckInput, int input_n, int input_c, int input_h, int input_w, int offset = 0) {
	if (offset == 1) { input_n = 1; }

	int temp1 = input_w * input_h * input_c;
	for (int ⁠n_idx = 0; ⁠n_idx < input_n; ⁠n_idx++)
	{
		int temp2 = ⁠n_idx * temp1;
		for (int ⁠c_idx = 0; ⁠c_idx < input_c; ⁠c_idx++)
		{
			int temp3 = ⁠c_idx * input_w * input_h + temp2;
			for (int ⁠h_idx = 0; ⁠h_idx < input_h; ⁠h_idx++)
			{
				int temp4 = ⁠h_idx * input_w + temp3;
				for (int w_idx = 0; w_idx < input_w; w_idx++)
				{
					int g_idx = w_idx + temp4;
					cout << setw(5) << valueCheckInput[g_idx] << " ";
				}cout << endl;
			}cout << endl; cout << endl;
		}
	}
}


void valueCheckV(vector<float>& valueCheckInput, int tile_n, int input_c, int input_n, int input_r, int offset = 0) {
	if (offset == 1) { input_n = 1; }

	int temp1 = input_c * input_n * input_r;
	for (int ⁠n_idx = 0; ⁠n_idx < tile_n; ⁠n_idx++)
	{
		int temp2 = ⁠n_idx * temp1;
		for (int ⁠c_idx = 0; ⁠c_idx < input_c; ⁠c_idx++)
		{
			int temp3 = ⁠c_idx * input_n * input_r + temp2;
			for (int ⁠h_idx = 0; ⁠h_idx < 4; ⁠h_idx++)
			{
				int temp4 = ⁠h_idx * input_n * 4 + temp3;
				for (int w_idx = 0; w_idx < input_n * 4; w_idx++)
				{
					int g_idx = w_idx + temp4;
					cout << setw(5) << valueCheckInput[g_idx] << " ";
				}cout << endl;
			}cout << endl; cout << endl;
		}
	}
}

// 행렬 A와 B를 더하여 C에 저장시키는 함수(A행렬 B행렬 C행렬은 같은 크기)
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

// 행렬 A와 B를 빼서 C에 저장시키는 함수(A행렬 B행렬 C행렬은 같은 크기)
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

// 행렬 A와 B를 곱하여 C에 저장시키는 함수(A행렬 B행렬 C행렬은 같은 크기)
// parameter : 행렬 A, B, C 

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

	for (int i = 0; i < U_rows / 2; ++i) { // 행

		//int temp1 = i * UcolsVrows ;
		//int temp2 = i * V_cols ;

		for (int j = 0; j < V_cols / 2; ++j) // 열 
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

	for (int i = 0; i < U_rows / 4; ++i) { // 행

		//int temp1 = i * UcolsVrows ;
		//int temp2 = i * V_cols ;

		for (int j = 0; j < V_cols / 4; ++j) // 열 
		{
			int m1 = 0; int	m2 = 0; int	m3 = 0; int	m4 = 0; int m5 = 0; int	m6 = 0; int	m7 = 0; int	m8 = 0;
			int m9 = 0; int	m10 = 0; int	m11 = 0; int	m12 = 0; int m13 = 0; int	m14 = 0; int	m15 = 0; int	m16 = 0;

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

// 임계값 구하는 함수
int getThreshold(int n)
{
	int th;
	double k = floor(log(n) / log(2) - 4);
	th = (int)floor(n / pow(2.0, k)) + 1;
	return th;
}

// 4개의 부분행렬로 나누는 함수
// parameter : 나눌 행렬, 저장할 행렬 공간 4개
void Submatrix(vector<float>& matrixOrigin, vector<float>& matrix11, vector<float>& matrix12, vector<float>& matrix21, vector<float>& matrix22, int Rows, int Cols)
{
	// int Rows, int Cols 부분행렬의 사이즈 

	for (int i = 0; i < Rows; i++)
	{
		int temp = i * (Cols * 2);
		int temp2 = i * (Cols);
		for (int j = 0; j < Cols; j++)
		{
			int gidx = temp + j;
			int gidx2 = temp2 + j;
			matrix11[gidx2] = matrixOrigin[gidx];									//좌 상단행렬
			matrix12[gidx2] = matrixOrigin[gidx + Cols];							//우 상단행렬
			matrix21[gidx2] = matrixOrigin[Rows * Cols * 2 + gidx];					//좌 하단행렬
			matrix22[gidx2] = matrixOrigin[Rows * Cols * 2 + gidx + Cols];			//우 하단행렬
		}
	}

}

// 4개의 부분행렬들을 재결합 해주는 함수
// parameter : 합친 결과를 저장할 행렬 , 부분행렬 11, 12, 21, 22
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
			matrixOrigin[gidx] = matrix11[gidx2];										//좌 상단행렬
			matrixOrigin[gidx + Cols] = matrix12[gidx2];								//우 상단행렬
			matrixOrigin[Rows * Cols * 2 + gidx] = matrix21[gidx2];         			//좌 하단행렬
			matrixOrigin[Rows * Cols * 2 + gidx + Cols] = matrix22[gidx2];				//우 하단행렬
		}
	}
}

// 쉬트라센 알고리즘 함수
void Strassen(vector<float>& matrixU, vector<float>& matrixV, vector<float>& matrixM, int U_Row, int U_Col, int V_Row, int V_Col)
{

	if (1 == (U_Row/4) % 2 || 1 == (U_Col/4) % 2 || 1 == (V_Row/4) % 2 || 1 == (V_Col/4) % 2)
	{
		//cout << "U_Row : " << U_Row << " , " << "U_Col : " << U_Col << " , " << "V_Row : " << V_Row << " , " << "V_Col : " << V_Col << " )" << endl;
		//cout << "Can't divide the matrix anymore" << endl;
		MatrixMulElement4X4(matrixU, matrixV, matrixM, U_Row, U_Col, V_Col);
		return;
	}
	else {
		//int th = V_Col / 4;
		if (V_Col/4 <= getThreshold(V_Col/4)) {

			MatrixMulElement4X4(matrixU, matrixV, matrixM, U_Row, U_Col, V_Col);

			return;
		}
		else {

			//cout << " ss" << endl;
			int newU_Row = U_Row / 2;					//4등분을 하기 위해
			int newU_Col = U_Col / 2;
			int newV_Row = V_Row / 2;
			int newV_Col = V_Col / 2;

			//a11~a22 부분행렬, b11~b22 부분행렬 
			vector<float> a11(newU_Row * newU_Col), a12(newU_Row * newU_Col), a21(newU_Row * newU_Col), a22(newU_Row * newU_Col);
			vector<float> b11(newV_Row * newV_Col), b12(newV_Row * newV_Col), b21(newV_Row * newV_Col), b22(newV_Row * newV_Col);

			//부분행렬들의 연산결과를 m1~m7 에저장
			vector<float>  m1(newU_Row * newV_Col), m2(newU_Row * newV_Col), m3(newU_Row * newV_Col), m4(newU_Row * newV_Col), m5(newU_Row * newV_Col), m6(newU_Row * newV_Col), m7(newU_Row * newV_Col);

			//a11~b22 의 연산결과들을 임시로 저장할 그릇
			vector<float>  tempA(newU_Row * newU_Col), tempB(newV_Row * newV_Col);

			vector<float>  tempAc(newU_Row * newV_Col), tempBc(newU_Row * newV_Col);

			// m1~m7 연산 결과로 C를 구하기 위해 저장 할 행렬
			vector<float>  c11(newU_Row * newV_Col), c12(newU_Row * newV_Col), c21(newU_Row * newV_Col), c22(newU_Row * newV_Col);


			//A의 부분행렬 4개, B의 부분행렬 4개 생성
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


			// 위에서 계산된 m1~m7 결과로  c11 ~ c22 를 만든다.
			MatrixSum(m1, m4, tempAc, newU_Row, newV_Col); //m1 + m4
			MatrixSum(tempAc, m7, tempBc, newU_Row, newV_Col); //m1 + m4 + m7
			MatrixSub(tempBc, m5, c11, newU_Row, newV_Col); //c11 = m1 + m4 - m5 + m7

			MatrixSum(m3, m5, c12, newU_Row, newV_Col); //c12 = m3 + m5

			MatrixSum(m2, m4, c21, newU_Row, newV_Col); //c21 = m2 + m4

			MatrixSum(m1, m3, tempAc, newU_Row, newV_Col); //m1 + m3
			MatrixSum(tempAc, m6, tempBc, newU_Row, newV_Col); //m1 + m3 + m6
			MatrixSub(tempBc, m2, c22, newU_Row, newV_Col); //c22 = m1 + m3 - m2 + m6

			//재 병합
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
	cout << "Winograd Convolutions filterTransform function! \n\n";

	int input_c = 8;
	int output_c = 8;

	// h[O_Ch][ln_Ch][KernelSize_h][KernelSize_w] 
	// h[4][3][3][3]
	// 임시 filter 값 
	vector<float> h(output_c*input_c * 3 * 3);
	float count = 1.f;
	int temp1 = input_c * 3 * 3;
	for (int ⁠n_idx = 0; ⁠n_idx < output_c; ⁠n_idx++)
	{
		int temp2 = ⁠n_idx * temp1;
		for (int ⁠c_idx = 0; ⁠c_idx < input_c; ⁠c_idx++)
		{
			int temp3 = ⁠c_idx * 3 * 3 + temp2;
			for (int ⁠h_idx = 0; ⁠h_idx < 3; ⁠h_idx++)
			{
				int temp4 = ⁠h_idx * 3 + temp3;
				for (int w_idx = 0; w_idx < 3; w_idx++)
				{
					int g_idx = w_idx + temp4;
					h[g_idx] = count;
					count++;
				}
			}
		}
	}

	// filater 값 체크
	//valueCheck(h, output_c, input_c, 3, 3);

	// U[O_Ch][ln_Ch][4][4] 
	vector<float> U(output_c * input_c * 4 * 4);



	///////////////////////////////
	// U 값 체크 [OC][IC][u]
	////////////////////////////////

	//valueCheck(U, output_c, input_c, 4, 4);

	int input_n = 16;
	int input_H = 16;
	int input_W = 16;
	int area_i = input_H * input_W;

	// d[N][IC][input_H][input_W] 
	// 임시 input 값 
	vector<float> input(input_n * input_c * area_i);
	 count = 1.f;
	int tempi = input_c * input_H * input_W;
	for (int ⁠n_idx = 0; ⁠n_idx < input_n; ⁠n_idx++)
	{
		int tempi2 = ⁠n_idx * tempi;
		for (int ⁠c_idx = 0; ⁠c_idx < input_c; ⁠c_idx++)
		{
			int tempi3 = ⁠c_idx * input_H * input_W + tempi2;
			for (int ⁠h_idx = 0; ⁠h_idx < input_H; ⁠h_idx++)
			{
				int tempi4 = ⁠h_idx * input_W + tempi3;
				for (int w_idx = 0; w_idx < input_W; w_idx++)
				{
					int g_idx = w_idx + tempi4;
					input[g_idx] = count;
					count++;
				}
			}
		}
	}

	// filater 값 체크
	//valueCheck(input, input_n, input_c, input_H, input_W);

	// output[N][out_ch][4][4] 
	int output_H = input_H - 2;
	int output_W = input_W - 2;
	int area_o = output_H * output_W;
	int tile_n = area_o / 4;
	vector<float> V(input_n * input_c * tile_n * 16);



	long long start_usec = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();

	// filterTransform 수행
	filterTransform2(U, h, output_c, input_c);
	// inputTransform 수행
	inputTransform2(V, input, input_n, input_c, input_H, input_W);

	//cout << endl;	cout << endl;

	//valueCheckV(V, tile_n, input_c, input_n, 16);


	///////////////////////////////
	// V 값 체크 [N][IC][T][v]
	////////////////////////////////

	vector<float> M; //
	vector<float> subM(input_n * output_c * 16); //
	vector<float> convOutput(input_n * output_c * area_o); //
	//for () // tile 개수 만큼
	
	
	for(int i = 0 ; i < tile_n; i++){
		
		vector<float> subV; //
		int tile_idx = input_n * input_c * 16 * i;

		for (int j = 0; j < input_n * input_c * 16; j++) {
			subV.push_back(V[j + tile_idx]);
		}

	Strassen(U, subV, subM, output_c * 4, input_c * 4, input_c * 4, input_n * 4);
	
		for (int j = 0; j < input_n * output_c * 16; j++) {
			M.push_back(subM[j]);
		}
	}

	//valueCheck(M, tile_n, output_c, input_n, 16);






	int temp_o = output_c * output_H * output_W;
	int temp_i = input_c * input_H * input_W;

	int temp_u = input_c * 16;

	for (int n = 0; n < input_n; n++)
	{
		int temp_i2 = n * temp_i;
		int temp_o2 = n * temp_o;

		int tilecount = 0;
		for (int row = 0; row < output_H; row += 2)
		{

			tilecount = output_H * row / 4;

			int row_idxo1 = row * output_W;
			int row_idxo2 = row_idxo1 + output_W;

			for (int col = 0; col < output_W; col += 2)
			{

				for (int outch = 0; outch < output_c; outch++)
				{
					int ot_idx1 = outch * output_H * output_W + temp_o2 + col + row_idxo1;
					int ot_idx3 = outch * output_H * output_W + temp_o2 + col + row_idxo2;

					int u_idx = tilecount * output_c * input_n * 16 + outch * input_n * 16 + n * 4;

						// U . V
						float m1 = M[u_idx];
						float m2 = M[u_idx + 1];
						float m3 = M[u_idx + 2];
						float m4 = M[u_idx + 3];

						int u_idx2 = u_idx + input_n * 4;
						float m5 = M[u_idx2];
						float m6 = M[u_idx2 + 1];
						float m7 = M[u_idx2 + 2];
						float m8 = M[u_idx2 + 3];

						int u_idx3 = u_idx2 + input_n * 4;
						float m9 = M[u_idx3];
						float m10 = M[u_idx3 + 1];
						float m11 = M[u_idx3 + 2];
						float m12 = M[u_idx3 + 3];

						int u_idx4 = u_idx3 + input_n * 4;
						float m13 = M[u_idx4];
						float m14 = M[u_idx4 + 1];
						float m15 = M[u_idx4 + 2];
						float m16 = M[u_idx4 + 3];

						// output transfom
						float sub_y1 = m2 + m6 + m10;
						float sub_y2 = m3 + m7 + m11;
						float sub_y3 = m6 - m10 - m14;
						float sub_y4 = m7 - m11 - m15;
					

					convOutput[ot_idx1] = m1 + m5 + m9 + sub_y1 + sub_y2;
					convOutput[ot_idx1 + 1] = sub_y1 - sub_y2 - m4 - m8 - m12;
					convOutput[ot_idx3] = m5 - m9 - m13 + sub_y3 + sub_y4;
					convOutput[ot_idx3 + 1] = sub_y3 - sub_y4 - m8 + m12 + m16;
				}
				tilecount++;
			}
		}
	}

	long long end_usec = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	int frame_sec = int(end_usec - start_usec);

	//valueCheck(convOutput, input_n, output_c, output_H, output_W);

	cout << "====================================="<< endl;
	cout << frame_sec << "u sec (Strassen - Winograd Convolution)" << endl;

}

