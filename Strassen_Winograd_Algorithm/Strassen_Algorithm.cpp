#include <stdio.h>
#include <stdlib.h>
#include <vector> //vector를 쓰기위한 라이브러리
#include <math.h> //log,floor 사용
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

	출처 :: http://yimoyimo.tk/Strassen/
*****************************************************************************/

#define SIZE 1024

// 바이너리 파일을 읽어와 벡터 2차원 배열에 저장시키는 함수.
// parameter : 읽은 값을 저장시킬 vector 주소 , 파일명, , 행사이즈, 열사이즈 

void MatrixInit(vector<vector<int>> & matrix, int Row, int Col)
{
	int idx = 1;
	for (int i = 0; i < Row; i++)
	{
		for (int j = 0; j < Col; j++)
		{
			matrix[i][j] = idx;
			//idx++;
		}
	}
}

void MatrixInit(vector<int> &matrix, int Row, int Col)
{
	//int idx = 1;
	
	for (int i = 0; i < Row; i++)
	{
		int temp = i * Col;
		for (int j = 0; j < Col; j++)
		{
			int gidx = temp + j; 
			matrix[gidx] = 1;
			//idx++;
		}
	}
}


// 행렬 A와 B를 더하여 C에 저장시키는 함수(A행렬 B행렬 C행렬은 같은 크기)
// parameter :행렬 A,B,C 
// return : 없음
void MatrixSum(vector< vector<int> >& matrixA, vector< vector<int> >& matrixB, vector< vector<int> >& matrixC)
{
	for (int i = 0; i < (int)matrixA.size(); i++) // matrixA.size() = 행 길이
	{
		for (int j = 0; j < (int)matrixA[i].size(); j++) // matrixA.size() = 열 길이
		{
			matrixC[i][j] = matrixA[i][j] + matrixB[i][j];
		}
	}

}

void MatrixSum(vector<int> & matrixA,  vector<int>& matrixB, vector<int> & matrixC)
{
	int side = sqrt(matrixA.size());

	for (int i = 0; i < side; i++) // matrixA.size() = 행 길이
	{
		int temp = i * side;
		for (int j = 0; j < side; j++) // matrixA.size() = 열 길이
		{
			int gidx = temp + j;
			matrixC[gidx] = matrixA[gidx] + matrixB[gidx];
		}
	}
}

// 행렬 A와 B를 빼서 C에 저장시키는 함수(A행렬 B행렬 C행렬은 같은 크기)
// parameter : 행렬 A,B,C 
// return : 없음
void MatrixSub(vector<vector<int>>& matrixA, vector<vector<int>>& matrixB, vector< vector<int> >& matrixC)
{
	for (int i = 0; i < (int)matrixA.size(); i++) // matrixA.size() = 행 길이
	{
		for (int j = 0; j < (int)matrixA[i].size(); j++) // matrixA.size() = 열 길이
		{
			matrixC[i][j] = matrixA[i][j] - matrixB[i][j];
		}
	}
}

void MatrixSub(vector<int>& matrixA, vector<int>& matrixB,  vector<int>& matrixC)
{
	int side = sqrt(matrixA.size());

	for (int i = 0; i < side; i++) // matrixA.size() = 행 길이
	{
		int temp = i * side;
		for (int j = 0; j < side; j++) // matrixA.size() = 열 길이
		{
			int gidx = temp + j;
			matrixC[gidx] = matrixA[gidx] - matrixB[gidx];
		}
	}


}

// 행렬 A와 B를 곱하여 C에 저장시키는 함수(A행렬 B행렬 C행렬은 같은 크기)
// parameter : 행렬 A, B, C 
// return : 없음
void MatrixMul(vector< vector<int> >& matrixA, vector< vector<int> >& matrixB, vector< vector<int> >& matrixC)
{
	for (int i = 0; i < (int)matrixA.size(); i++) // matrixA.size() = 행 길이
	{
		for (int j = 0; j < (int)matrixA[i].size(); j++) // matrixA.size() = 열 길이
		{
			for (int k = 0; k < (int)matrixA[i].size(); k++)
			{
				matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
			}
		}
	}
}


void MatrixMul( vector<int> Prev_matrix, vector<int> Post_matrix, vector<int>& Y_output, int prev_rows, int prev_colsAndPost_rows, int post_cols)
{
	int newNum = sqrt(Prev_matrix.size());  //부분행렬의 사이즈 

	prev_rows = newNum;
	prev_colsAndPost_rows = newNum;
	post_cols = newNum;
	for (int i = 0; i < prev_rows; ++i) {
		int temp1 = i * prev_colsAndPost_rows;
		int temp2 = i * post_cols;
		for (int j = 0; j < post_cols; ++j)
		{
			int sum = 0;
			for (int k = 0; k < prev_colsAndPost_rows; ++k)
			{
				sum += Prev_matrix[temp1 + k] * Post_matrix[k* post_cols + j];
			}
			Y_output[temp2 + j] = (int)sum;
		}
	}
}



// 임계값 구하는 함수
// parameter : 행렬의 가로 및 세로 길이 (ex:1024)
// return : 임계값 반환
int getThreshold(int n)
{
	int th;
	double k = floor(log(n) / log(2) - 6);
	th = (int)floor(n / pow(2.0, k)) + 1;
	return th;
}

// 4개의 부분행렬로 나누는 함수
// parameter : 나눌 행렬, 저장할 행렬 공간 4개
// return : 없음
void Submatrix(vector< vector<int> >& matrixOrigin, vector< vector<int> >& matrix11, vector< vector<int> >& matrix12,
	vector< vector<int> >& matrix21, vector< vector<int> >& matrix22)
{
	int newNum = matrix11.size();  //부분행렬의 사이즈 

	for (int i = 0; i < newNum; i++)
	{
		for (int j = 0; j < newNum; j++)
		{
			matrix11[i][j] = matrixOrigin[i][j];						//좌 상단행렬
			matrix12[i][j] = matrixOrigin[i][j + newNum];				//우 상단행렬
			matrix21[i][j] = matrixOrigin[i + newNum][j];				//좌 하단행렬
			matrix22[i][j] = matrixOrigin[i + newNum][j + newNum];		//우 하단행렬
		}
	}

}

void Submatrix( vector<int> & matrixOrigin,  vector<int> & matrix11,  vector<int> & matrix12, vector<int> & matrix21,  vector<int> & matrix22)
{
	int newNum = sqrt(matrix11.size());  //부분행렬의 사이즈 

	for (int i = 0; i < newNum; i++)
	{
		int temp = i * (newNum*2);
		int temp2 = i * (newNum );
		for (int j = 0; j < newNum; j++)
		{
			int gidx = temp + j;
			int gidx2 = temp2 + j;
			matrix11[gidx2] = matrixOrigin[gidx];									//좌 상단행렬
			matrix12[gidx2] = matrixOrigin[gidx + newNum];							//우 상단행렬
			matrix21[gidx2] = matrixOrigin[newNum * newNum * 2 + gidx];				//좌 하단행렬
			matrix22[gidx2] = matrixOrigin[newNum * newNum * 2 + gidx + newNum];		//우 하단행렬
		}
	}

}

// 4개의 부분행렬들을 재결합 해주는 함수
// parameter : 합친 결과를 저장할 행렬 , 부분행렬 11, 12, 21, 22
// return : 없음
void Mergematrix(vector< vector<int> >& matrixOrigin, vector< vector<int> >& matrix11, vector< vector<int> >& matrix12,
	vector< vector<int> >& matrix21, vector< vector<int> >& matrix22)
{
	int newNum = matrix11.size();  //부분행렬의 사이즈

	for (int i = 0; i < newNum; i++)
	{
		for (int j = 0; j < newNum; j++)
		{
			matrixOrigin[i][j] = matrix11[i][j];						//좌 상단행렬
			matrixOrigin[i][j + newNum] = matrix12[i][j];				//우 상단행렬
			matrixOrigin[i + newNum][j] = matrix21[i][j];         	    //좌 하단행렬
			matrixOrigin[i + newNum][j + newNum] = matrix22[i][j];	    //우 하단행렬
		}
	}
}

void Mergematrix( vector<int> & matrixOrigin,  vector<int> & matrix11,  vector<int> & matrix12, vector<int> & matrix21,  vector<int> & matrix22)
{
	int newNum = sqrt(matrix11.size());  //부분행렬의 사이즈 

	for (int i = 0; i < newNum; i++)
	{
		int temp = i * (newNum * 2);
		int temp2 = i * newNum;
		for (int j = 0; j < newNum; j++)
		{
			int gidx = temp + j;
			int gidx2 = temp2 + j;
			matrixOrigin[gidx] = matrix11[gidx2];										//좌 상단행렬
			matrixOrigin[gidx + newNum] = matrix12[gidx2];								//우 상단행렬
			matrixOrigin[newNum * newNum * 2 + gidx] = matrix21[gidx2];         	    //좌 하단행렬
			matrixOrigin[newNum * newNum * 2 + gidx + newNum] = matrix22[gidx2];	    //우 하단행렬
		}
	}
}

void ElementMul(vector<int>& Y_output, vector<int> Prev_matrix, vector<int> Post_matrix){
	for(int i = 0 ; i < 16 ; i ++)
	Y_output[i] = Prev_matrix[i] * Post_matrix[i];

}


// 쉬트라센 알고리즘 함수
// parameter : 행렬의 사이즈(ex: 1024x1024 -> 1024입력) , 행렬 A, B, C
// return : 없음
void Strassen(int n, vector< vector<int> >& matrixA, vector< vector<int> >& matrixB, vector< vector<int> >& matrixC)
{
	
	if (n <= getThreshold(n))
	{
		//cout << " getThreshold" << endl;
		//ElementMul(matrixC, matrixA, matrixB);
		MatrixMul(matrixA, matrixB, matrixC);
		return;
	}
	else
	{
		//cout << " ss" << endl;
		int newRow = n / 2;					//4등분을 하기 위해
		vector<int> newCol(newRow, 0);

		//a11~a22 부분행렬, b11~b22 부분행렬 
		vector < vector<int> > a11(newRow, newCol), a12(newRow, newCol), a21(newRow, newCol), a22(newRow, newCol);
		vector < vector<int> > b11(newRow, newCol), b12(newRow, newCol), b21(newRow, newCol), b22(newRow, newCol);

		//부분행렬들의 연산결과를 m1~m7 에저장
		vector<vector<int>> m1(newRow, newCol), m2(newRow, newCol), m3(newRow, newCol), m4(newRow, newCol)
			, m5(newRow, newCol), m6(newRow, newCol), m7(newRow, newCol);

		//a11~b22 의 연산결과들을 임시로 저장할 그릇
		vector < vector<int> > tempA(newRow, newCol), tempB(newRow, newCol);

		// m1~m7 연산 결과로 C를 구하기 위해 저장 할 행렬
		vector < vector<int> > c11(newRow, newCol), c12(newRow, newCol), c21(newRow, newCol), c22(newRow, newCol);


		//A의 부분행렬 4개, B의 부분행렬 4개 생성
		Submatrix(matrixA, a11, a12, a21, a22);
		Submatrix(matrixB, b11, b12, b21, b22);


		MatrixSum(a11, a22, tempA);		       // a11+a22
		MatrixSum(b11, b22, tempB);		       // b11+b22
		Strassen(newRow, tempA, tempB, m1);    // m1=(a11+a11)(b11+b22)

		MatrixSum(a21, a22, tempA);            // a21+a22
		Strassen(newRow, tempA, b11, m2);      // m2=(a21+a22)b11

		MatrixSub(b12, b22, tempB);            // b12-b22
		Strassen(newRow, a11, tempB, m3);      // m3=a11(b12-b22)

		MatrixSub(b21, b11, tempB);            // b21-b11
		Strassen(newRow, a22, tempB, m4);      // m4=a22(b21-11)

		MatrixSum(a11, a12, tempA);            //  a11+a12
		Strassen(newRow, tempA, b22, m5); 	   // m5=(a11+a12)b22

		MatrixSub(a21, a11, tempA);            // a21-a11
		MatrixSum(b11, b12, tempB);            // b11+b12
		Strassen(newRow, tempA, tempB, m6);    // m6=(a21-a11)(b11+b12)

		MatrixSub(a12, a22, tempA);            // a12-a22
		MatrixSum(b21, b22, tempB);            // b21+b22
		Strassen(newRow, tempA, tempB, m7);    // m7 = (a12 - a22)(a12 - a22)


		// 위에서 계산된 m1~m7 결과로  c11 ~ c22 를 만든다.
		MatrixSum(m1, m4, tempA); //m1 + m4
		MatrixSum(tempA, m7, tempB); //m1 + m4 + m7
		MatrixSub(tempB, m5, c11); //c11 = m1 + m4 - m5 + m7

		MatrixSum(m3, m5, c12); //c12 = m3 + m5

		MatrixSum(m2, m4, c21); //c21 = m2 + m4

		MatrixSum(m1, m3, tempA); //m1 + m3
		MatrixSum(tempA, m6, tempB); //m1 + m3 + m6
		MatrixSub(tempB, m2, c22); //c22 = m1 + m3 - m2 + m6

		//재 병합
		Mergematrix(matrixC, c11, c12, c21, c22);
	}
}

void Strassen(int n,  vector<int> & matrixA,  vector<int> & matrixB,  vector<int> & matrixC)
{

	if (n <= getThreshold(n))
	{
		//cout << " getThreshold" << endl;
		//ElementMul(matrixC, matrixA, matrixB);
		MatrixMul(matrixA, matrixB, matrixC, n,n,n);
		return;
	}
	else
	{
		//cout << " ss" << endl;
		int newRow = n / 2;					//4등분을 하기 위해

		//a11~a22 부분행렬, b11~b22 부분행렬 
		vector<int> a11(newRow * newRow), a12(newRow * newRow), a21(newRow * newRow), a22(newRow * newRow);
		vector<int> b11(newRow * newRow), b12(newRow * newRow), b21(newRow * newRow), b22(newRow * newRow);

		//부분행렬들의 연산결과를 m1~m7 에저장
		vector<int>  m1(newRow * newRow), m2(newRow * newRow), m3(newRow * newRow), m4(newRow * newRow), m5(newRow * newRow), m6(newRow * newRow), m7(newRow * newRow);

		//a11~b22 의 연산결과들을 임시로 저장할 그릇
		vector<int>  tempA(newRow * newRow), tempB(newRow * newRow);

		// m1~m7 연산 결과로 C를 구하기 위해 저장 할 행렬
		vector<int>  c11(newRow * newRow), c12(newRow * newRow), c21(newRow * newRow), c22(newRow * newRow);


		//A의 부분행렬 4개, B의 부분행렬 4개 생성
		Submatrix(matrixA, a11, a12, a21, a22);
		Submatrix(matrixB, b11, b12, b21, b22);


		MatrixSum(a11, a22, tempA);		       // a11+a22
		MatrixSum(b11, b22, tempB);		       // b11+b22
		Strassen(newRow, tempA, tempB, m1);    // m1=(a11+a11)(b11+b22)

		MatrixSum(a21, a22, tempA);            // a21+a22
		Strassen(newRow, tempA, b11, m2);      // m2=(a21+a22)b11

		MatrixSub(b12, b22, tempB);            // b12-b22
		Strassen(newRow, a11, tempB, m3);      // m3=a11(b12-b22)

		MatrixSub(b21, b11, tempB);            // b21-b11
		Strassen(newRow, a22, tempB, m4);      // m4=a22(b21-11)

		MatrixSum(a11, a12, tempA);            //  a11+a12
		Strassen(newRow, tempA, b22, m5); 	   // m5=(a11+a12)b22

		MatrixSub(a21, a11, tempA);            // a21-a11
		MatrixSum(b11, b12, tempB);            // b11+b12
		Strassen(newRow, tempA, tempB, m6);    // m6=(a21-a11)(b11+b12)

		MatrixSub(a12, a22, tempA);            // a12-a22
		MatrixSum(b21, b22, tempB);            // b21+b22
		Strassen(newRow, tempA, tempB, m7);    // m7 = (a12 - a22)(a12 - a22)


		// 위에서 계산된 m1~m7 결과로  c11 ~ c22 를 만든다.
		MatrixSum(m1, m4, tempA); //m1 + m4
		MatrixSum(tempA, m7, tempB); //m1 + m4 + m7
		MatrixSub(tempB, m5, c11); //c11 = m1 + m4 - m5 + m7

		MatrixSum(m3, m5, c12); //c12 = m3 + m5

		MatrixSum(m2, m4, c21); //c21 = m2 + m4

		MatrixSum(m1, m3, tempA); //m1 + m3
		MatrixSum(tempA, m6, tempB); //m1 + m3 + m6
		MatrixSub(tempB, m2, c22); //c22 = m1 + m3 - m2 + m6

		//재 병합
		Mergematrix(matrixC, c11, c12, c21, c22);
	}
}


void FileOutput(vector<vector<int> >& matrix)
{

	for (int i = 0; i < (int)matrix.size(); i++)
	{
		for (int j = 0; j < (int)matrix.size(); j++)
		{
			printf("%d  ", matrix[i][j]);
		}
		printf("\n");
	}

}

void FileOutput(vector<int>& matrix)
{
	int side = sqrt(matrix.size());
	for (int i = 0; i < side; i++)
	{
		int temp = side * i;
		for (int j = 0; j < side; j++)
		{
			int gidx = temp + j;
			printf("%d  ", matrix[gidx]);
		}
		printf("\n");
	}

}


int main()
{
	int Row = SIZE;
	int Col = SIZE;


	// ex) vector<vector<int> >  arr(6, vector<int>(5, 0));
	// 설명: vector< vector<int> > 형 벡터 6개(가로 6줄)를 할당 한다는 뜻
	//       vector<int>(5,0) 은 모든 가로줄은 5개짜리 0으로 초기화 된 익명의 int 형 벡터 배열을 생성해 초기값을 넣는다.

	vector < vector<int> > A(Row, vector<int>(Col, 0)); // A[Row][Col] 인 동적 배열
	vector < vector<int> > B(Row, vector<int>(Col, 0)); // B[Row][Col] 인 동적 배열
	vector < vector<int> > C(Row, vector<int>(Col, 0)); // C[Row][Col] 인 동적 배열
	vector < vector<int> > D(Row, vector<int>(Col, 0)); // D[Row][Col] 인 동적 배열

	vector <int> A1(Row*Col); // 
	vector <int> B1(Row*Col); // 
	vector <int> D1(Row*Col); //

	MatrixInit(A, Row, Col);
	MatrixInit(B, Row, Col);

	MatrixInit(A1, Row, Col);
	MatrixInit(B1, Row, Col);

	/*
	long long start_usec = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	MatrixMul(A, B, D);
	long long end_usec = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	int frame_sec = int(end_usec - start_usec);
	cout << frame_sec << "u sec (MatrixMul) v2" << endl;
	*/
	
	long long start_usec3 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	MatrixMul( A1, B1, D1, Row, Row, Row);
	long long end_usec3 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	int frame_sec3 = int(end_usec3 - start_usec3);
	cout << frame_sec3 << "u sec (MatrixMul) v1" << endl;
	
	/*
	long long start_usec4 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	Strassen(SIZE, A, B, D);
	long long end_usec4 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	int frame_sec4 = int(end_usec4 - start_usec4);
	cout << frame_sec4 << "u sec (Strassen) v2" << endl;
	//FileOutput(D);
	*/

	long long start_usec2 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	Strassen(SIZE, A1, B1, D1);
	long long end_usec2 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	int frame_sec2 = int(end_usec2 - start_usec2);
	cout << frame_sec2 << "u sec (Strassen) v1" << endl;

	//FileOutput(D1);

	
	//cout << endl;FileOutput(C);cout << endl;



	return 0;
}