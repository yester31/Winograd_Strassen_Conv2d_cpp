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


void filterTransform(vector<float>& Output_U, vector<float> filter, int output_c, int input_c) {

	int temp1u = input_c * 16;
	int temp1f = input_c * 9;
	for (int oc_idx = 0; oc_idx < output_c; oc_idx++)
	{
		int temp2u = oc_idx * temp1u;
		int temp2f = oc_idx * temp1f;
		for (int ⁠c_idx = 0; ⁠c_idx < input_c; ⁠c_idx++)
		{
			int u_idx = ⁠c_idx * 16 + temp2u;
			int f_idx = ⁠c_idx * 9 + temp2f;

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

			Output_U[u_idx + 4] = (temp1) / 2.f;
			Output_U[u_idx + 5] = (temp1 + temp2 + temp3) / 4.f;
			Output_U[u_idx + 6] = (temp1 - temp2 + temp3) / 4.f;
			Output_U[u_idx + 7] = (temp3) / 2.f;

			float temp4 = *g1 - *g4 + *g7;
			float temp5 = *g2 - *g5 + *g8;
			float temp6 = *g3 - *g6 + *g9;

			Output_U[u_idx + 8] = (temp4) / 2.f;
			Output_U[u_idx + 9] = (temp4 + temp5 + temp6) / 4.f;
			Output_U[u_idx + 10] = (temp4 - temp5 + temp6) / 4.f;
			Output_U[u_idx + 11] = (temp6) / 2.f;

			Output_U[u_idx + 12] = *g7;
			Output_U[u_idx + 13] = (*g7 + *g8 + *g9) / 2.f;
			Output_U[u_idx + 14] = (*g7 - *g8 + *g9) / 2.f;
			Output_U[u_idx + 15] = *g9;

		}
	}
}

void zeroPadding(vector<float>& zeroPaddingOutput, vector<float>& zeroPaddingInput, int input_n, int input_c, int input_h, int input_w, int leftPadingSize, int rightPadingSize, int topPadingSize, int bottomPadingSize) {

	int temp1 = input_w * input_h * input_c;
	int temp1o = (input_h + topPadingSize + bottomPadingSize) * (input_w + leftPadingSize + rightPadingSize) * input_c;
	for (int ⁠n_idx = 0; ⁠n_idx < input_n; ⁠n_idx++)
	{
		int temp2 = ⁠n_idx * temp1;
		int temp2o = ⁠n_idx * temp1o;
		for (int ⁠c_idx = 0; ⁠c_idx < input_c; ⁠c_idx++)
		{
			int temp3 = ⁠c_idx * input_w * input_h + temp2;
			int temp3o = ⁠c_idx * (input_w + leftPadingSize + rightPadingSize) * (input_h + topPadingSize + bottomPadingSize) + temp2o;
			for (int ⁠h_idx = 0; ⁠h_idx < input_h; ⁠h_idx++)
			{
				int temp4 = ⁠h_idx * input_w + temp3;
				int temp4o = (⁠h_idx + topPadingSize) * (input_w + leftPadingSize + rightPadingSize) + leftPadingSize + temp3o;

				for (int w_idx = 0; w_idx < input_w; w_idx++)
				{
					int ⁠g_idx = w_idx + temp4;
					int g_idx_Output = w_idx + temp4o;
					zeroPaddingOutput[g_idx_Output] = zeroPaddingInput[⁠g_idx];
				}
			}
		}
	}
}

void winogradConv2d(vector<float>& convOutput, vector<float>& convInput, vector<float> kernel, int input_n, int input_c, int input_h, int input_w, int output_c, int leftPadingSize, int rightPadingSize, int topPadingSize, int bottomPadingSize) {

	int input_h_z = input_h + topPadingSize + bottomPadingSize;
	int input_w_z = input_w + leftPadingSize + rightPadingSize;

	//1. zeropadding 수행
	vector<float> convInputWithZP(input_n * input_c * input_h_z * input_w_z);
	zeroPadding(convInputWithZP, convInput, input_n, input_c, input_h, input_w, leftPadingSize, rightPadingSize, topPadingSize, bottomPadingSize);
	//cout << "zeropadding" << endl;
	//valueCheck(convInputWithZP, input_n, input_c, input_h_z, input_w_z);

	vector<float> U(output_c * input_c * 16);
	//2. filterTransform 수행
	filterTransform(U, kernel, output_c, input_c);

	//3. input transform 및 element wise multiplication, output transform수행
	int output_h_z = input_h_z - 2;
	int output_w_z = input_w_z - 2;

	int temp_o = output_c * output_h_z * output_w_z;
	int temp_i = input_c * input_h_z * input_w_z;
	int temp_u = input_c * 16;

	if (input_w_z % 2 == 0 && input_h_z % 2 == 0) {
		//cout << "가로 짝, 세로 짝" << endl;
		for (int n = 0; n < input_n; n++)
		{
			int temp_i2 = n * temp_i;
			int temp_o2 = n * temp_o;

			for (int row = 0; row < output_h_z; row += 2)
			{
				int row_idx1 = row * input_w_z;
				int row_idx2 = row_idx1 + input_w_z;
				int row_idx3 = row_idx2 + input_w_z;
				int row_idx4 = row_idx3 + input_w_z;

				int row_idxo1 = row * output_w_z;
				int row_idxo2 = row_idxo1 + output_w_z;

				for (int col = 0; col < output_w_z; col += 2)
				{
					for (int outch = 0; outch < output_c; outch++)
					{
						int temp_u2 = outch * temp_u;
						int ot_idx1 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo1;
						int ot_idx3 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo2;

						float y1 = 0;
						float y2 = 0;
						float y3 = 0;
						float y4 = 0;

						for (int inch = 0; inch < input_c; inch++)
						{
							int temp_ic = inch * input_h_z * input_w_z;
							int u_idx = inch * 16 + temp_u2; // U idex

							int t_idx1 = temp_ic + temp_i2 + row_idx1 + col;
							int t_idx2 = temp_ic + temp_i2 + row_idx2 + col;
							int t_idx3 = temp_ic + temp_i2 + row_idx3 + col;
							int t_idx4 = temp_ic + temp_i2 + row_idx4 + col;

							float* d1 = &convInputWithZP[t_idx1];
							float* d2 = &convInputWithZP[t_idx1 + 1];
							float* d3 = &convInputWithZP[t_idx1 + 2];
							float* d4 = &convInputWithZP[t_idx1 + 3];

							float* d5 = &convInputWithZP[t_idx2];
							float* d6 = &convInputWithZP[t_idx2 + 1];
							float* d7 = &convInputWithZP[t_idx2 + 2];
							float* d8 = &convInputWithZP[t_idx2 + 3];

							float* d9 = &convInputWithZP[t_idx3];
							float* d10 = &convInputWithZP[t_idx3 + 1];
							float* d11 = &convInputWithZP[t_idx3 + 2];
							float* d12 = &convInputWithZP[t_idx3 + 3];

							float* d13 = &convInputWithZP[t_idx4];
							float* d14 = &convInputWithZP[t_idx4 + 1];
							float* d15 = &convInputWithZP[t_idx4 + 2];
							float* d16 = &convInputWithZP[t_idx4 + 3];

							float dd1 = *d11 - (*d3);
							float dd2 = *d2 - (*d10);
							float dd3 = *d7 + (*d11);
							float dd4 = *d6 + (*d10);
							float dd5 = *d7 - (*d11);
							float dd6 = *d10 - (*d6);
							float dd7 = *d15 - (*d7);
							float dd8 = *d6 - (*d14);

							float v1 = *d1 - *d9 + dd1;
							float v2 = dd2 - dd1;//
							float v3 = -dd1 - dd2;//
							float v4 = dd2 - *d4 + *d12;

							float v5 = *d5 + *d9 - dd3;
							float v6 = dd4 + dd3;
							float v7 = dd3 - dd4;
							float v8 = dd4 - *d8 - *d12;

							float v9 = *d9 - *d5 + dd5;
							float v10 = dd6 - dd5;
							float v11 = -(dd6 + dd5);
							float v12 = dd6 + *d8 - *d12;

							float v13 = *d5 - *d13 + dd7;
							float v14 = dd8 - dd7;
							float v15 = -dd7 - dd8;
							float v16 = dd8 - *d8 + *d16;

							/*
							cout << setw(5) << v1 << " " << setw(5) << v2 << " " << setw(5) << v3 << " " << setw(5) << v4 ;
							cout << setw(5) << v5 << " " << setw(5) << v6 << " " << setw(5) << v7 << " " << setw(5) << v8 ;
							cout << setw(5) << v9 << " " << setw(5) << v10 << " " << setw(5) << v11 << " " << setw(5) << v12 ;
							cout << setw(5) << v13 << " " << setw(5) << v14 << " " << setw(5) << v15 << " " << setw(5) << v16;
							cout << endl;
							*/



							// U . V
							float m1 = v1 * U[u_idx];
							float m2 = v2 * U[u_idx + 1];
							float m3 = v3 * U[u_idx + 2];
							float m4 = v4 * U[u_idx + 3];
							float m5 = v5 * U[u_idx + 4];
							float m6 = v6 * U[u_idx + 5];
							float m7 = v7 * U[u_idx + 6];
							float m8 = v8 * U[u_idx + 7];
							float m9 = v9 * U[u_idx + 8];
							float m10 = v10 * U[u_idx + 9];
							float m11 = v11 * U[u_idx + 10];
							float m12 = v12 * U[u_idx + 11];
							float m13 = v13 * U[u_idx + 12];
							float m14 = v14 * U[u_idx + 13];
							float m15 = v15 * U[u_idx + 14];
							float m16 = v16 * U[u_idx + 15];

							// output transfom
							float sub_y1 = m2 + m6 + m10;
							float sub_y2 = m3 + m7 + m11;
							float sub_y3 = m6 - m10 - m14;
							float sub_y4 = m7 - m11 - m15;

							y1 += m1 + m5 + m9 + sub_y1 + sub_y2;
							y2 += sub_y1 - sub_y2 - m4 - m8 - m12;
							y3 += m5 - m9 - m13 + sub_y3 + sub_y4;
							y4 += sub_y3 - sub_y4 - m8 + m12 + m16;
						}

						convOutput[ot_idx1] = y1;
						convOutput[ot_idx1 + 1] = y2;
						convOutput[ot_idx3] = y3;
						convOutput[ot_idx3 + 1] = y4;
					}
				}
			}
		}

	}
	else if (input_w_z % 2 == 1 && input_h_z % 2 == 0) {
		//cout << "가로 홀, 세로 짝" << endl;

		for (int n = 0; n < input_n; n++)
		{
			int temp_i2 = n * temp_i;
			int temp_o2 = n * temp_o;

			for (int row = 0; row < output_h_z; row += 2)
			{
				int row_idx1 = row * input_w_z;
				int row_idx2 = row_idx1 + input_w_z;
				int row_idx3 = row_idx2 + input_w_z;
				int row_idx4 = row_idx3 + input_w_z;

				int row_idxo1 = row * output_w_z;
				int row_idxo2 = row_idxo1 + output_w_z;

				int col;
				for (col = 0; col < output_w_z - 1; col += 2)
				{
					for (int outch = 0; outch < output_c; outch++)
					{
						int temp_u2 = outch * temp_u;
						int ot_idx1 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo1;
						int ot_idx3 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo2;

						float y1 = 0;
						float y2 = 0;
						float y3 = 0;
						float y4 = 0;

						for (int inch = 0; inch < input_c; inch++)
						{
							int temp_ic = inch * input_h_z * input_w_z;
							int u_idx = inch * 16 + temp_u2; // U idex

							int t_idx1 = temp_ic + temp_i2 + row_idx1 + col;
							int t_idx2 = temp_ic + temp_i2 + row_idx2 + col;
							int t_idx3 = temp_ic + temp_i2 + row_idx3 + col;
							int t_idx4 = temp_ic + temp_i2 + row_idx4 + col;

							float* d1 = &convInputWithZP[t_idx1];
							float* d2 = &convInputWithZP[t_idx1 + 1];
							float* d3 = &convInputWithZP[t_idx1 + 2];
							float* d4 = &convInputWithZP[t_idx1 + 3];

							float* d5 = &convInputWithZP[t_idx2];
							float* d6 = &convInputWithZP[t_idx2 + 1];
							float* d7 = &convInputWithZP[t_idx2 + 2];
							float* d8 = &convInputWithZP[t_idx2 + 3];

							float* d9 = &convInputWithZP[t_idx3];
							float* d10 = &convInputWithZP[t_idx3 + 1];
							float* d11 = &convInputWithZP[t_idx3 + 2];
							float* d12 = &convInputWithZP[t_idx3 + 3];

							float* d13 = &convInputWithZP[t_idx4];
							float* d14 = &convInputWithZP[t_idx4 + 1];
							float* d15 = &convInputWithZP[t_idx4 + 2];
							float* d16 = &convInputWithZP[t_idx4 + 3];

							float dd1 = *d11 - (*d3);
							float dd2 = *d2 - (*d10);
							float dd3 = *d7 + (*d11);
							float dd4 = *d6 + (*d10);
							float dd5 = *d7 - (*d11);
							float dd6 = *d10 - (*d6);
							float dd7 = *d15 - (*d7);
							float dd8 = *d6 - (*d14);

							float v1 = *d1 - *d9 + dd1;
							float v2 = dd2 - dd1;//
							float v3 = -dd1 - dd2;//
							float v4 = dd2 - *d4 + *d12;

							float v5 = *d5 + *d9 - dd3;
							float v6 = dd4 + dd3;
							float v7 = dd3 - dd4;
							float v8 = dd4 - *d8 - *d12;

							float v9 = *d9 - *d5 + dd5;
							float v10 = dd6 - dd5;
							float v11 = -(dd6 + dd5);
							float v12 = dd6 + *d8 - *d12;

							float v13 = *d5 - *d13 + dd7;
							float v14 = dd8 - dd7;
							float v15 = -dd7 - dd8;
							float v16 = dd8 - *d8 + *d16;

							// U . V
							float m1 = v1 * U[u_idx];
							float m2 = v2 * U[u_idx + 1];
							float m3 = v3 * U[u_idx + 2];
							float m4 = v4 * U[u_idx + 3];
							float m5 = v5 * U[u_idx + 4];
							float m6 = v6 * U[u_idx + 5];
							float m7 = v7 * U[u_idx + 6];
							float m8 = v8 * U[u_idx + 7];
							float m9 = v9 * U[u_idx + 8];
							float m10 = v10 * U[u_idx + 9];
							float m11 = v11 * U[u_idx + 10];
							float m12 = v12 * U[u_idx + 11];
							float m13 = v13 * U[u_idx + 12];
							float m14 = v14 * U[u_idx + 13];
							float m15 = v15 * U[u_idx + 14];
							float m16 = v16 * U[u_idx + 15];

							// output transfom
							float sub_y1 = m2 + m6 + m10;
							float sub_y2 = m3 + m7 + m11;
							float sub_y3 = m6 - m10 - m14;
							float sub_y4 = m7 - m11 - m15;

							y1 += m1 + m5 + m9 + sub_y1 + sub_y2;
							y2 += sub_y1 - sub_y2 - m4 - m8 - m12;
							y3 += m5 - m9 - m13 + sub_y3 + sub_y4;
							y4 += sub_y3 - sub_y4 - m8 + m12 + m16;
						}

						convOutput[ot_idx1] = y1;
						convOutput[ot_idx1 + 1] = y2;
						convOutput[ot_idx3] = y3;
						convOutput[ot_idx3 + 1] = y4;
					}
				}

				for (int outch = 0; outch < output_c; outch++)
				{
					int temp_u2 = outch * temp_u;
					int ot_idx1 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo1;
					int ot_idx3 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo2;

					float y1 = 0;
					float y3 = 0;

					for (int inch = 0; inch < input_c; inch++)
					{
						int temp_ic = inch * input_h_z * input_w_z;
						int u_idx = inch * 16 + temp_u2; // U idex

						int t_idx1 = temp_ic + temp_i2 + row_idx1 + col;
						int t_idx2 = temp_ic + temp_i2 + row_idx2 + col;
						int t_idx3 = temp_ic + temp_i2 + row_idx3 + col;
						int t_idx4 = temp_ic + temp_i2 + row_idx4 + col;

						float* d1 = &convInputWithZP[t_idx1];
						float* d2 = &convInputWithZP[t_idx1 + 1];
						float* d3 = &convInputWithZP[t_idx1 + 2];

						float* d5 = &convInputWithZP[t_idx2];
						float* d6 = &convInputWithZP[t_idx2 + 1];
						float* d7 = &convInputWithZP[t_idx2 + 2];

						float* d9 = &convInputWithZP[t_idx3];
						float* d10 = &convInputWithZP[t_idx3 + 1];
						float* d11 = &convInputWithZP[t_idx3 + 2];

						float* d13 = &convInputWithZP[t_idx4];
						float* d14 = &convInputWithZP[t_idx4 + 1];
						float* d15 = &convInputWithZP[t_idx4 + 2];

						float dd1 = *d11 - (*d3);
						float dd2 = *d2 - (*d10);
						float dd3 = *d7 + (*d11);
						float dd4 = *d6 + (*d10);
						float dd5 = *d7 - (*d11);
						float dd6 = *d10 - (*d6);
						float dd7 = *d15 - (*d7);
						float dd8 = *d6 - (*d14);

						float v1 = *d1 - *d9 + dd1;
						float v2 = dd2 - dd1;//
						float v3 = -dd1 - dd2;//

						float v5 = *d5 + *d9 - dd3;
						float v6 = dd4 + dd3;
						float v7 = dd3 - dd4;

						float v9 = *d9 - *d5 + dd5;
						float v10 = dd6 - dd5;
						float v11 = -(dd6 + dd5);

						float v13 = *d5 - *d13 + dd7;
						float v14 = dd8 - dd7;
						float v15 = -dd7 - dd8;

						// U . V
						float m1 = v1 * U[u_idx];
						float m2 = v2 * U[u_idx + 1];
						float m3 = v3 * U[u_idx + 2];
						float m5 = v5 * U[u_idx + 4];
						float m6 = v6 * U[u_idx + 5];
						float m7 = v7 * U[u_idx + 6];
						float m9 = v9 * U[u_idx + 8];
						float m10 = v10 * U[u_idx + 9];
						float m11 = v11 * U[u_idx + 10];
						float m13 = v13 * U[u_idx + 12];
						float m14 = v14 * U[u_idx + 13];
						float m15 = v15 * U[u_idx + 14];

						// output transfom
						float sub_y1 = m2 + m6 + m10;
						float sub_y2 = m3 + m7 + m11;
						float sub_y3 = m6 - m10 - m14;
						float sub_y4 = m7 - m11 - m15;

						y1 += m1 + m5 + m9 + sub_y1 + sub_y2;
						y3 += m5 - m9 - m13 + sub_y3 + sub_y4;
					}
					convOutput[ot_idx1] = y1;
					convOutput[ot_idx3] = y3;
				}
			}
		}
	}
	else if (input_w_z % 2 == 0 && input_h_z % 2 == 1) {

		//cout << "가로 짝, 세로 홀" << endl;

		for (int n = 0; n < input_n; n++)
		{
			int temp_i2 = n * temp_i;
			int temp_o2 = n * temp_o;

			int row;
			for (row = 0; row < output_h_z - 1; row += 2)
			{
				int row_idx1 = row * input_w_z;
				int row_idx2 = row_idx1 + input_w_z;
				int row_idx3 = row_idx2 + input_w_z;
				int row_idx4 = row_idx3 + input_w_z;

				int row_idxo1 = row * output_w_z;
				int row_idxo2 = row_idxo1 + output_w_z;

				for (int col = 0; col < output_w_z; col += 2)
				{
					for (int outch = 0; outch < output_c; outch++)
					{
						int temp_u2 = outch * temp_u;
						int ot_idx1 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo1;
						int ot_idx3 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo2;

						float y1 = 0;
						float y2 = 0;
						float y3 = 0;
						float y4 = 0;

						for (int inch = 0; inch < input_c; inch++)
						{
							int temp_ic = inch * input_h_z * input_w_z;
							int u_idx = inch * 16 + temp_u2; // U idex

							int t_idx1 = temp_ic + temp_i2 + row_idx1 + col;
							int t_idx2 = temp_ic + temp_i2 + row_idx2 + col;
							int t_idx3 = temp_ic + temp_i2 + row_idx3 + col;
							int t_idx4 = temp_ic + temp_i2 + row_idx4 + col;

							float* d1 = &convInputWithZP[t_idx1];
							float* d2 = &convInputWithZP[t_idx1 + 1];
							float* d3 = &convInputWithZP[t_idx1 + 2];
							float* d4 = &convInputWithZP[t_idx1 + 3];

							float* d5 = &convInputWithZP[t_idx2];
							float* d6 = &convInputWithZP[t_idx2 + 1];
							float* d7 = &convInputWithZP[t_idx2 + 2];
							float* d8 = &convInputWithZP[t_idx2 + 3];

							float* d9 = &convInputWithZP[t_idx3];
							float* d10 = &convInputWithZP[t_idx3 + 1];
							float* d11 = &convInputWithZP[t_idx3 + 2];
							float* d12 = &convInputWithZP[t_idx3 + 3];

							float* d13 = &convInputWithZP[t_idx4];
							float* d14 = &convInputWithZP[t_idx4 + 1];
							float* d15 = &convInputWithZP[t_idx4 + 2];
							float* d16 = &convInputWithZP[t_idx4 + 3];

							float dd1 = *d11 - (*d3);
							float dd2 = *d2 - (*d10);
							float dd3 = *d7 + (*d11);
							float dd4 = *d6 + (*d10);
							float dd5 = *d7 - (*d11);
							float dd6 = *d10 - (*d6);
							float dd7 = *d15 - (*d7);
							float dd8 = *d6 - (*d14);

							float v1 = *d1 - *d9 + dd1;
							float v2 = dd2 - dd1;//
							float v3 = -dd1 - dd2;//
							float v4 = dd2 - *d4 + *d12;

							float v5 = *d5 + *d9 - dd3;
							float v6 = dd4 + dd3;
							float v7 = dd3 - dd4;
							float v8 = dd4 - *d8 - *d12;

							float v9 = *d9 - *d5 + dd5;
							float v10 = dd6 - dd5;
							float v11 = -(dd6 + dd5);
							float v12 = dd6 + *d8 - *d12;

							float v13 = *d5 - *d13 + dd7;
							float v14 = dd8 - dd7;
							float v15 = -dd7 - dd8;
							float v16 = dd8 - *d8 + *d16;

							// U . V
							float m1 = v1 * U[u_idx];
							float m2 = v2 * U[u_idx + 1];
							float m3 = v3 * U[u_idx + 2];
							float m4 = v4 * U[u_idx + 3];
							float m5 = v5 * U[u_idx + 4];
							float m6 = v6 * U[u_idx + 5];
							float m7 = v7 * U[u_idx + 6];
							float m8 = v8 * U[u_idx + 7];
							float m9 = v9 * U[u_idx + 8];
							float m10 = v10 * U[u_idx + 9];
							float m11 = v11 * U[u_idx + 10];
							float m12 = v12 * U[u_idx + 11];
							float m13 = v13 * U[u_idx + 12];
							float m14 = v14 * U[u_idx + 13];
							float m15 = v15 * U[u_idx + 14];
							float m16 = v16 * U[u_idx + 15];

							// output transfom
							float sub_y1 = m2 + m6 + m10;
							float sub_y2 = m3 + m7 + m11;
							float sub_y3 = m6 - m10 - m14;
							float sub_y4 = m7 - m11 - m15;

							y1 += m1 + m5 + m9 + sub_y1 + sub_y2;
							y2 += sub_y1 - sub_y2 - m4 - m8 - m12;
							y3 += m5 - m9 - m13 + sub_y3 + sub_y4;
							y4 += sub_y3 - sub_y4 - m8 + m12 + m16;
						}

						convOutput[ot_idx1] = y1;
						convOutput[ot_idx1 + 1] = y2;
						convOutput[ot_idx3] = y3;
						convOutput[ot_idx3 + 1] = y4;
					}
				}
			}

			int row_idx1 = row * input_w_z;
			int row_idx2 = row_idx1 + input_w_z;
			int row_idx3 = row_idx2 + input_w_z;
			int row_idx4 = row_idx3 + input_w_z;

			int row_idxo1 = row * output_w_z;
			int row_idxo2 = row_idxo1 + output_w_z;

			for (int col = 0; col < output_w_z; col += 2)
			{
				for (int outch = 0; outch < output_c; outch++)
				{
					int temp_u2 = outch * temp_u;
					int ot_idx1 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo1;
					int ot_idx3 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo2;

					float y1 = 0;
					float y2 = 0;

					for (int inch = 0; inch < input_c; inch++)
					{
						int temp_ic = inch * input_h_z * input_w_z;
						int u_idx = inch * 16 + temp_u2; // U idex

						int t_idx1 = temp_ic + temp_i2 + row_idx1 + col;
						int t_idx2 = temp_ic + temp_i2 + row_idx2 + col;
						int t_idx3 = temp_ic + temp_i2 + row_idx3 + col;
						int t_idx4 = temp_ic + temp_i2 + row_idx4 + col;

						float* d1 = &convInputWithZP[t_idx1];
						float* d2 = &convInputWithZP[t_idx1 + 1];
						float* d3 = &convInputWithZP[t_idx1 + 2];
						float* d4 = &convInputWithZP[t_idx1 + 3];

						float* d5 = &convInputWithZP[t_idx2];
						float* d6 = &convInputWithZP[t_idx2 + 1];
						float* d7 = &convInputWithZP[t_idx2 + 2];
						float* d8 = &convInputWithZP[t_idx2 + 3];

						float* d9 = &convInputWithZP[t_idx3];
						float* d10 = &convInputWithZP[t_idx3 + 1];
						float* d11 = &convInputWithZP[t_idx3 + 2];
						float* d12 = &convInputWithZP[t_idx3 + 3];

						float dd1 = *d11 - (*d3);
						float dd2 = *d2 - (*d10);
						float dd3 = *d7 + (*d11);
						float dd4 = *d6 + (*d10);
						float dd5 = *d7 - (*d11);
						float dd6 = *d10 - (*d6);

						float v1 = *d1 - *d9 + dd1;
						float v2 = dd2 - dd1;//
						float v3 = -dd1 - dd2;//
						float v4 = dd2 - *d4 + *d12;

						float v5 = *d5 + *d9 - dd3;
						float v6 = dd4 + dd3;
						float v7 = dd3 - dd4;
						float v8 = dd4 - *d8 - *d12;

						float v9 = *d9 - *d5 + dd5;
						float v10 = dd6 - dd5;
						float v11 = -(dd6 + dd5);
						float v12 = dd6 + *d8 - *d12;

						// U . V
						float m1 = v1 * U[u_idx];
						float m2 = v2 * U[u_idx + 1];
						float m3 = v3 * U[u_idx + 2];
						float m4 = v4 * U[u_idx + 3];
						float m5 = v5 * U[u_idx + 4];
						float m6 = v6 * U[u_idx + 5];
						float m7 = v7 * U[u_idx + 6];
						float m8 = v8 * U[u_idx + 7];
						float m9 = v9 * U[u_idx + 8];
						float m10 = v10 * U[u_idx + 9];
						float m11 = v11 * U[u_idx + 10];
						float m12 = v12 * U[u_idx + 11];

						// output transfom
						float sub_y1 = m2 + m6 + m10;
						float sub_y2 = m3 + m7 + m11;

						y1 += m1 + m5 + m9 + sub_y1 + sub_y2;
						y2 += sub_y1 - sub_y2 - m4 - m8 - m12;
					}

					convOutput[ot_idx1] = y1;
					convOutput[ot_idx1 + 1] = y2;
				}
			}
		}
	}
	else {
		//cout << "가로 홀, 세로 홀" << endl;

		for (int n = 0; n < input_n; n++)
		{
			int temp_i2 = n * temp_i;
			int temp_o2 = n * temp_o;

			int row;
			for (row = 0; row < output_h_z - 1; row += 2)
			{
				int row_idx1 = row * input_w_z;
				int row_idx2 = row_idx1 + input_w_z;
				int row_idx3 = row_idx2 + input_w_z;
				int row_idx4 = row_idx3 + input_w_z;

				int row_idxo1 = row * output_w_z;
				int row_idxo2 = row_idxo1 + output_w_z;

				int col;
				for (col = 0; col < output_w_z - 1; col += 2)
				{

					for (int outch = 0; outch < output_c; outch++)
					{
						int temp_u2 = outch * temp_u;
						int ot_idx1 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo1;
						int ot_idx3 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo2;

						float y1 = 0;
						float y2 = 0;
						float y3 = 0;
						float y4 = 0;

						for (int inch = 0; inch < input_c; inch++)
						{
							int temp_ic = inch * input_h_z * input_w_z;
							int u_idx = inch * 16 + temp_u2; // U idex

							int t_idx1 = temp_ic + temp_i2 + row_idx1 + col;
							int t_idx2 = temp_ic + temp_i2 + row_idx2 + col;
							int t_idx3 = temp_ic + temp_i2 + row_idx3 + col;
							int t_idx4 = temp_ic + temp_i2 + row_idx4 + col;

							float* d1 = &convInputWithZP[t_idx1];
							float* d2 = &convInputWithZP[t_idx1 + 1];
							float* d3 = &convInputWithZP[t_idx1 + 2];
							float* d4 = &convInputWithZP[t_idx1 + 3];

							float* d5 = &convInputWithZP[t_idx2];
							float* d6 = &convInputWithZP[t_idx2 + 1];
							float* d7 = &convInputWithZP[t_idx2 + 2];
							float* d8 = &convInputWithZP[t_idx2 + 3];

							float* d9 = &convInputWithZP[t_idx3];
							float* d10 = &convInputWithZP[t_idx3 + 1];
							float* d11 = &convInputWithZP[t_idx3 + 2];
							float* d12 = &convInputWithZP[t_idx3 + 3];

							float* d13 = &convInputWithZP[t_idx4];
							float* d14 = &convInputWithZP[t_idx4 + 1];
							float* d15 = &convInputWithZP[t_idx4 + 2];
							float* d16 = &convInputWithZP[t_idx4 + 3];

							float dd1 = *d11 - (*d3);
							float dd2 = *d2 - (*d10);
							float dd3 = *d7 + (*d11);
							float dd4 = *d6 + (*d10);
							float dd5 = *d7 - (*d11);
							float dd6 = *d10 - (*d6);
							float dd7 = *d15 - (*d7);
							float dd8 = *d6 - (*d14);

							float v1 = *d1 - *d9 + dd1;
							float v2 = dd2 - dd1;//
							float v3 = -dd1 - dd2;//
							float v4 = dd2 - *d4 + *d12;

							float v5 = *d5 + *d9 - dd3;
							float v6 = dd4 + dd3;
							float v7 = dd3 - dd4;
							float v8 = dd4 - *d8 - *d12;

							float v9 = *d9 - *d5 + dd5;
							float v10 = dd6 - dd5;
							float v11 = -(dd6 + dd5);
							float v12 = dd6 + *d8 - *d12;

							float v13 = *d5 - *d13 + dd7;
							float v14 = dd8 - dd7;
							float v15 = -dd7 - dd8;
							float v16 = dd8 - *d8 + *d16;

							// U . V
							float m1 = v1 * U[u_idx];
							float m2 = v2 * U[u_idx + 1];
							float m3 = v3 * U[u_idx + 2];
							float m4 = v4 * U[u_idx + 3];
							float m5 = v5 * U[u_idx + 4];
							float m6 = v6 * U[u_idx + 5];
							float m7 = v7 * U[u_idx + 6];
							float m8 = v8 * U[u_idx + 7];
							float m9 = v9 * U[u_idx + 8];
							float m10 = v10 * U[u_idx + 9];
							float m11 = v11 * U[u_idx + 10];
							float m12 = v12 * U[u_idx + 11];
							float m13 = v13 * U[u_idx + 12];
							float m14 = v14 * U[u_idx + 13];
							float m15 = v15 * U[u_idx + 14];
							float m16 = v16 * U[u_idx + 15];

							// output transfom
							float sub_y1 = m2 + m6 + m10;
							float sub_y2 = m3 + m7 + m11;
							float sub_y3 = m6 - m10 - m14;
							float sub_y4 = m7 - m11 - m15;

							y1 += m1 + m5 + m9 + sub_y1 + sub_y2;
							y2 += sub_y1 - sub_y2 - m4 - m8 - m12;
							y3 += m5 - m9 - m13 + sub_y3 + sub_y4;
							y4 += sub_y3 - sub_y4 - m8 + m12 + m16;
						}

						convOutput[ot_idx1] = y1;
						convOutput[ot_idx1 + 1] = y2;
						convOutput[ot_idx3] = y3;
						convOutput[ot_idx3 + 1] = y4;
					}
				}

				for (int outch = 0; outch < output_c; outch++)
				{
					int temp_u2 = outch * temp_u;
					int ot_idx1 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo1;
					int ot_idx3 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo2;

					float y1 = 0;
					float y3 = 0;

					for (int inch = 0; inch < input_c; inch++)
					{
						int temp_ic = inch * input_h_z * input_w_z;
						int u_idx = inch * 16 + temp_u2; // U idex

						int t_idx1 = temp_ic + temp_i2 + row_idx1 + col;
						int t_idx2 = temp_ic + temp_i2 + row_idx2 + col;
						int t_idx3 = temp_ic + temp_i2 + row_idx3 + col;
						int t_idx4 = temp_ic + temp_i2 + row_idx4 + col;

						float* d1 = &convInputWithZP[t_idx1];
						float* d2 = &convInputWithZP[t_idx1 + 1];
						float* d3 = &convInputWithZP[t_idx1 + 2];

						float* d5 = &convInputWithZP[t_idx2];
						float* d6 = &convInputWithZP[t_idx2 + 1];
						float* d7 = &convInputWithZP[t_idx2 + 2];

						float* d9 = &convInputWithZP[t_idx3];
						float* d10 = &convInputWithZP[t_idx3 + 1];
						float* d11 = &convInputWithZP[t_idx3 + 2];

						float* d13 = &convInputWithZP[t_idx4];
						float* d14 = &convInputWithZP[t_idx4 + 1];
						float* d15 = &convInputWithZP[t_idx4 + 2];

						float dd1 = *d11 - (*d3);
						float dd2 = *d2 - (*d10);
						float dd3 = *d7 + (*d11);
						float dd4 = *d6 + (*d10);
						float dd5 = *d7 - (*d11);
						float dd6 = *d10 - (*d6);
						float dd7 = *d15 - (*d7);
						float dd8 = *d6 - (*d14);

						float v1 = *d1 - *d9 + dd1;
						float v2 = dd2 - dd1;//
						float v3 = -dd1 - dd2;//

						float v5 = *d5 + *d9 - dd3;
						float v6 = dd4 + dd3;
						float v7 = dd3 - dd4;

						float v9 = *d9 - *d5 + dd5;
						float v10 = dd6 - dd5;
						float v11 = -(dd6 + dd5);

						float v13 = *d5 - *d13 + dd7;
						float v14 = dd8 - dd7;
						float v15 = -dd7 - dd8;

						// U . V
						float m1 = v1 * U[u_idx];
						float m2 = v2 * U[u_idx + 1];
						float m3 = v3 * U[u_idx + 2];
						float m5 = v5 * U[u_idx + 4];
						float m6 = v6 * U[u_idx + 5];
						float m7 = v7 * U[u_idx + 6];
						float m9 = v9 * U[u_idx + 8];
						float m10 = v10 * U[u_idx + 9];
						float m11 = v11 * U[u_idx + 10];
						float m13 = v13 * U[u_idx + 12];
						float m14 = v14 * U[u_idx + 13];
						float m15 = v15 * U[u_idx + 14];

						// output transfom
						float sub_y1 = m2 + m6 + m10;
						float sub_y2 = m3 + m7 + m11;
						float sub_y3 = m6 - m10 - m14;
						float sub_y4 = m7 - m11 - m15;

						y1 += m1 + m5 + m9 + sub_y1 + sub_y2;
						y3 += m5 - m9 - m13 + sub_y3 + sub_y4;
					}
					convOutput[ot_idx1] = y1;
					convOutput[ot_idx3] = y3;
				}
			}

			int row_idx1 = row * input_w_z;
			int row_idx2 = row_idx1 + input_w_z;
			int row_idx3 = row_idx2 + input_w_z;
			int row_idx4 = row_idx3 + input_w_z;

			int row_idxo1 = row * output_w_z;
			int row_idxo2 = row_idxo1 + output_w_z;

			int col;
			for (col = 0; col < output_w_z - 1; col += 2)
			{
				for (int outch = 0; outch < output_c; outch++)
				{
					int temp_u2 = outch * temp_u;
					int ot_idx1 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo1;
					int ot_idx3 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo2;

					float y1 = 0;
					float y2 = 0;

					for (int inch = 0; inch < input_c; inch++)
					{
						int temp_ic = inch * input_h_z * input_w_z;
						int u_idx = inch * 16 + temp_u2; // U idex

						int t_idx1 = temp_ic + temp_i2 + row_idx1 + col;
						int t_idx2 = temp_ic + temp_i2 + row_idx2 + col;
						int t_idx3 = temp_ic + temp_i2 + row_idx3 + col;
						int t_idx4 = temp_ic + temp_i2 + row_idx4 + col;

						float* d1 = &convInputWithZP[t_idx1];
						float* d2 = &convInputWithZP[t_idx1 + 1];
						float* d3 = &convInputWithZP[t_idx1 + 2];
						float* d4 = &convInputWithZP[t_idx1 + 3];

						float* d5 = &convInputWithZP[t_idx2];
						float* d6 = &convInputWithZP[t_idx2 + 1];
						float* d7 = &convInputWithZP[t_idx2 + 2];
						float* d8 = &convInputWithZP[t_idx2 + 3];

						float* d9 = &convInputWithZP[t_idx3];
						float* d10 = &convInputWithZP[t_idx3 + 1];
						float* d11 = &convInputWithZP[t_idx3 + 2];
						float* d12 = &convInputWithZP[t_idx3 + 3];

						float dd1 = *d11 - (*d3);
						float dd2 = *d2 - (*d10);
						float dd3 = *d7 + (*d11);
						float dd4 = *d6 + (*d10);
						float dd5 = *d7 - (*d11);
						float dd6 = *d10 - (*d6);

						float v1 = *d1 - *d9 + dd1;
						float v2 = dd2 - dd1;//
						float v3 = -dd1 - dd2;//
						float v4 = dd2 - *d4 + *d12;

						float v5 = *d5 + *d9 - dd3;
						float v6 = dd4 + dd3;
						float v7 = dd3 - dd4;
						float v8 = dd4 - *d8 - *d12;

						float v9 = *d9 - *d5 + dd5;
						float v10 = dd6 - dd5;
						float v11 = -(dd6 + dd5);
						float v12 = dd6 + *d8 - *d12;

						// U . V
						float m1 = v1 * U[u_idx];
						float m2 = v2 * U[u_idx + 1];
						float m3 = v3 * U[u_idx + 2];
						float m4 = v4 * U[u_idx + 3];
						float m5 = v5 * U[u_idx + 4];
						float m6 = v6 * U[u_idx + 5];
						float m7 = v7 * U[u_idx + 6];
						float m8 = v8 * U[u_idx + 7];
						float m9 = v9 * U[u_idx + 8];
						float m10 = v10 * U[u_idx + 9];
						float m11 = v11 * U[u_idx + 10];
						float m12 = v12 * U[u_idx + 11];

						// output transfom
						float sub_y1 = m2 + m6 + m10;
						float sub_y2 = m3 + m7 + m11;

						y1 += m1 + m5 + m9 + sub_y1 + sub_y2;
						y2 += sub_y1 - sub_y2 - m4 - m8 - m12;
					}

					convOutput[ot_idx1] = y1;
					convOutput[ot_idx1 + 1] = y2;
				}
			}

			for (int outch = 0; outch < output_c; outch++)
			{
				int temp_u2 = outch * temp_u;
				int ot_idx1 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo1;
				int ot_idx3 = outch * output_h_z * output_w_z + temp_o2 + col + row_idxo2;

				float y1 = 0;

				for (int inch = 0; inch < input_c; inch++)
				{
					int temp_ic = inch * input_h_z * input_w_z;
					int u_idx = inch * 16 + temp_u2; // U idex

					int t_idx1 = temp_ic + temp_i2 + row_idx1 + col;
					int t_idx2 = temp_ic + temp_i2 + row_idx2 + col;
					int t_idx3 = temp_ic + temp_i2 + row_idx3 + col;

					float* d1 = &convInputWithZP[t_idx1];
					float* d2 = &convInputWithZP[t_idx1 + 1];
					float* d3 = &convInputWithZP[t_idx1 + 2];

					float* d5 = &convInputWithZP[t_idx2];
					float* d6 = &convInputWithZP[t_idx2 + 1];
					float* d7 = &convInputWithZP[t_idx2 + 2];

					float* d9 = &convInputWithZP[t_idx3];
					float* d10 = &convInputWithZP[t_idx3 + 1];
					float* d11 = &convInputWithZP[t_idx3 + 2];

					float dd1 = *d11 - (*d3);
					float dd2 = *d2 - (*d10);
					float dd3 = *d7 + (*d11);
					float dd4 = *d6 + (*d10);
					float dd5 = *d7 - (*d11);
					float dd6 = *d10 - (*d6);

					float v1 = *d1 - *d9 + dd1;
					float v2 = dd2 - dd1;//
					float v3 = -dd1 - dd2;//

					float v5 = *d5 + *d9 - dd3;
					float v6 = dd4 + dd3;
					float v7 = dd3 - dd4;

					float v9 = *d9 - *d5 + dd5;
					float v10 = dd6 - dd5;
					float v11 = -(dd6 + dd5);

					// U . V
					float m1 = v1 * U[u_idx];
					float m2 = v2 * U[u_idx + 1];
					float m3 = v3 * U[u_idx + 2];
					float m5 = v5 * U[u_idx + 4];
					float m6 = v6 * U[u_idx + 5];
					float m7 = v7 * U[u_idx + 6];
					float m9 = v9 * U[u_idx + 8];
					float m10 = v10 * U[u_idx + 9];
					float m11 = v11 * U[u_idx + 10];

					// output transfom
					float sub_y1 = m2 + m6 + m10;
					float sub_y2 = m3 + m7 + m11;

					y1 += m1 + m5 + m9 + sub_y1 + sub_y2;
				}
				convOutput[ot_idx1] = y1;
			}
		}
	}
}

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

void inputTransform3(vector<float>& V, vector<float>input, int input_n, int input_c, int input_h, int input_w) {

	int output_H = input_h - 2;
	int	output_W = input_w - 2;
	int area_o = output_H * output_W;
	//	int tile_N = (input_h - 2) / 2 * (input_w - 2) / 2;
	int tile_N = area_o / 4;
	int tile_wn = output_W / 2;

	int temp_i = input_c * input_h * input_w;
	int temp_v = 4 * area_o;// 16 * (input_h - 2)/2 * (input_w - 2)/2
	int tempic = input_n * tile_N * 16;
	for (int n = 0; n < input_n; n++) // N
	{
		int temp_i2 = n * temp_i;
		int tempic3 = n * tile_N * 16;

		for (int inch = 0; inch < input_c; inch++) // IC
		{
			int temp_i3 = inch * input_h * input_w + temp_i2;
			int tempic2 = tempic * inch;

			int tilecount = 0;
			for (int row = 0; row < output_H; row += 2)  // T_N
			{
				tilecount = tile_wn * (row / 2);

				int row_idx1 = row * input_w;
				int row_idx2 = row_idx1 + input_w;
				int row_idx3 = row_idx2 + input_w;
				int row_idx4 = row_idx3 + input_w;

				for (int col = 0; col < output_W; col += 2)
				{
					int temp_v3 = tempic2 + tempic3 + (tilecount + col / 2) * 4;

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

					V[temp_v3] = *d1 - *d9 + dd1;
					V[temp_v3 + 1] = dd2 - dd1;
					V[temp_v3 + 2] = -dd1 - dd2;
					V[temp_v3 + 3] = dd2 - *d4 + *d12;

					int temp_v4 = temp_v3 + tile_N * 4;
					V[temp_v4] = *d5 + *d9 - dd3;
					V[temp_v4 + 1] = dd4 + dd3;
					V[temp_v4 + 2] = dd3 - dd4;
					V[temp_v4 + 3] = dd4 - *d8 - *d12;

					int temp_v5 = temp_v3 + tile_N * 4 + tile_N * 4;
					V[temp_v5] = *d9 - *d5 + dd5;
					V[temp_v5 + 1] = dd6 - dd5;
					V[temp_v5 + 2] = -(dd6 + dd5);
					V[temp_v5 + 3] = dd6 + *d8 - *d12;

					int temp_v6 = temp_v3 + tile_N * 4 + tile_N * 4 + tile_N * 4;
					V[temp_v6] = *d5 - *d13 + dd7;
					V[temp_v6 + 1] = dd8 - dd7;
					V[temp_v6 + 2] = -dd7 - dd8;
					V[temp_v6 + 3] = dd8 - *d8 + *d16;

					/*
					cout << setw(5) << V[temp_v3] << " " << setw(5) << V[temp_v3 + 1] << " " << setw(5) << V[temp_v3 + 2] << " " << setw(5) << V[temp_v3 + 3];
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

void valueCheck(vector<float>& valueCheckInput, int input_n, int input_c, int input_h, int input_w, int offset = 0) {
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

void MatrixMulElement4X4n(vector<float> matrix_U, vector<float> matrix_V, vector<float>& matrix_M, int U_rows, int UcolsVrows, int V_cols, int tile_n)
{
	int Output_c = U_rows / 4;
	int Input_c = UcolsVrows / 4;
	int Input_n = V_cols / (4 * tile_n);
	//int tile_n = tile_n; // tile 총 개수 
	int Output_W = (int)sqrt(tile_n) * 2;
	int temp_o = Output_c * Output_W * Output_W;

	for (int oc = 0; oc < Output_c; ++oc) { // OC

		int tempo1 =  oc * Output_W * Output_W ;
		//int temp2 = i * V_cols ;

		for (int n = 0; n < Input_n; ++n) // N
		{
			int temp_o2 = n * temp_o;

			for (int row = 0; row < Output_W; row += 2) //TILE
			{
				int newrow = row / 2;
				int row_idxo1 = row * Output_W;
				int row_idxo2 = row_idxo1 + Output_W;
				int row_idxo3 = newrow * sqrt(tile_n);

				for (int col = 0; col < Output_W; col += 2)//TILE
				{
					int m1 = 0; int	m2 = 0; int	m3 = 0; int	m4 = 0; int m5 = 0; int	m6 = 0; int	m7 = 0; int	m8 = 0;
					int m9 = 0; int	m10 = 0; int m11 = 0; int	m12 = 0; int m13 = 0; int	m14 = 0; int	m15 = 0; int	m16 = 0;

					for (int ic = 0; ic < Input_c; ++ic) // IC
					{
						int midx1 = oc * Input_c * 16 + ic * 4;

						int midx2 = ic * Input_n * tile_n * 16 + n * tile_n * 16 + (col/2 + row_idxo3)*4;

						m1 += matrix_U[midx1] * matrix_V[midx2];
						m2 += matrix_U[midx1 + 1] * matrix_V[midx2 + 1];
						m3 += matrix_U[midx1 + 2] * matrix_V[midx2 + 2];
						m4 += matrix_U[midx1 + 3] * matrix_V[midx2 + 3];

						int midx3 = midx1 + Input_c * 4;
						int midx4 = midx2 + tile_n * 4;
						m5 += matrix_U[midx3] * matrix_V[midx4];
						m6 += matrix_U[midx3 + 1] * matrix_V[midx4 + 1];
						m7 += matrix_U[midx3 + 2] * matrix_V[midx4 + 2];
						m8 += matrix_U[midx3 + 3] * matrix_V[midx4 + 3];

						int midx5 = midx3 + Input_c * 4;
						int midx6 = midx4 + tile_n * 4;
						m9 += matrix_U[midx5] * matrix_V[midx6];
						m10 += matrix_U[midx5 + 1] * matrix_V[midx6 + 1];
						m11 += matrix_U[midx5 + 2] * matrix_V[midx6 + 2];
						m12 += matrix_U[midx5 + 3] * matrix_V[midx6 + 3];

						int midx7 = midx5 + Input_c * 4;
						int midx8 = midx6 + tile_n * 4;
						m13 += matrix_U[midx7] * matrix_V[midx8];
						m14 += matrix_U[midx7 + 1] * matrix_V[midx8 + 1];
						m15 += matrix_U[midx7 + 2] * matrix_V[midx8 + 2];
						m16 += matrix_U[midx7 + 3] * matrix_V[midx8 + 3];

					}

					int ot_idx1 = temp_o2 + tempo1 + col + row_idxo1;
					int ot_idx2 = temp_o2 + tempo1 + col + row_idxo2;

					// output transfom
					float sub_y1 = m2 + m6 + m10;
					float sub_y2 = m3 + m7 + m11;
					float sub_y3 = m6 - m10 - m14;
					float sub_y4 = m7 - m11 - m15;

					matrix_M[ot_idx1] = m1 + m5 + m9 + sub_y1 + sub_y2;
					matrix_M[ot_idx1 + 1] = sub_y1 - sub_y2 - m4 - m8 - m12;
					matrix_M[ot_idx2] = m5 - m9 - m13 + sub_y3 + sub_y4;
					matrix_M[ot_idx2 + 1] = sub_y3 - sub_y4 - m8 + m12 + m16;

				}
			}
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
			matrix11[gidx2] = matrixOrigin[gidx];									//좌 상단행렬 // 11
			matrix12[gidx2] = matrixOrigin[gidx + Cols];							//우 상단행렬 // 12  //
			matrix21[gidx2] = matrixOrigin[Rows * Cols * 2 + gidx];					//좌 하단행렬 // 21
			matrix22[gidx2] = matrixOrigin[Rows * Cols * 2 + gidx + Cols];			//우 하단행렬 // 22  //
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
			matrixOrigin[gidx + Cols] = matrix12[gidx2];								//우 상단행렬 ////
			matrixOrigin[Rows * Cols * 2 + gidx] = matrix21[gidx2];         			//좌 하단행렬
			matrixOrigin[Rows * Cols * 2 + gidx + Cols] = matrix22[gidx2];				//우 하단행렬 ////
		}
	}
}
void Mergematrix2(vector<float>& matrixOrigin, vector<float>& matrix11, vector<float>& matrix12, vector<float>& matrix21, vector<float>& matrix22, int Rows, int Cols)
{
	for (int i = 0; i < Rows; i++)
	{
		int temp = i * Cols ;

		for (int j = 0; j < Cols; j++)
		{
			int gidx = temp + j;

			matrixOrigin[gidx] = matrix11[gidx];										//좌 상단행렬 11
			matrixOrigin[Rows * Cols  + gidx] = matrix12[gidx];						//우 상단행렬 12
			matrixOrigin[Rows * 2 * Cols  + gidx] = matrix21[gidx];         			//좌 하단행렬 21
			matrixOrigin[Rows * 3 * Cols  + gidx] = matrix22[gidx];					//우 하단행렬 22
		}
	}
}
// 쉬트라센 알고리즘 함수
void Strassen(vector<float>& matrixU, vector<float>& matrixV, vector<float>& matrixM, int U_Row, int U_Col, int V_Row, int V_Col, int tile_n )
{

	if (1 == (U_Row / 4) % 2 || 1 == (U_Col / 4) % 2 || 1 == (V_Row / 4) % 2 || 1 == (V_Col / 4) % 2)
	{
		MatrixMulElement4X4n(matrixU, matrixV, matrixM, U_Row, U_Col, V_Col, tile_n);
		return;
	}
	else {
		int th = V_Row/4;
		if (th <= getThreshold(th)) {
			MatrixMulElement4X4n(matrixU, matrixV, matrixM, U_Row, U_Col, V_Col, tile_n);
			return;
		}
		else {
			int newU_Row = U_Row / 2;					//4등분을 하기 위해
			int newU_Col = U_Col / 2;
			int newV_Row = V_Row / 2;
			int newV_Col = V_Col / 2;

			//a11~a22 부분행렬, b11~b22 부분행렬 
			vector<float> a11(newU_Row * newU_Col), a12(newU_Row * newU_Col), a21(newU_Row * newU_Col), a22(newU_Row * newU_Col);
			vector<float> b11(newV_Row * newV_Col), b12(newV_Row * newV_Col), b21(newV_Row * newV_Col), b22(newV_Row * newV_Col);

			//부분행렬들의 연산결과를 m1~m7 에저장
			vector<float>  m1(newU_Row/2 * newV_Col / 2), m2(newU_Row / 2 * newV_Col / 2), m3(newU_Row / 2 * newV_Col / 2), m4(newU_Row / 2 * newV_Col / 2), m5(newU_Row / 2 * newV_Col / 2), m6(newU_Row / 2 * newV_Col / 2), m7(newU_Row / 2 * newV_Col / 2);

			//a11~b22 의 연산결과들을 임시로 저장할 그릇
			vector<float>  tempA(newU_Row * newU_Col), tempB(newV_Row * newV_Col);

			vector<float>  tempAc(newU_Row / 2 * newV_Col / 2), tempBc(newU_Row / 2 * newV_Col / 2);

			// m1~m7 연산 결과로 C를 구하기 위해 저장 할 행렬
			vector<float>  c11(newU_Row / 2 * newV_Col / 2), c12(newU_Row / 2 * newV_Col / 2), c21(newU_Row / 2 * newV_Col / 2), c22(newU_Row / 2 * newV_Col / 2);


			//A의 부분행렬 4개, B의 부분행렬 4개 생성
			Submatrix(matrixU, a11, a12, a21, a22, newU_Row, newU_Col);
			Submatrix(matrixV, b11, b12, b21, b22, newV_Row, newV_Col);

			MatrixSum(a11, a22, tempA, newU_Row, newU_Col);				// a11+a22
			MatrixSum(b11, b22, tempB, newV_Row, newV_Col);		       // b11+b22
			Strassen(tempA, tempB, m1, newU_Row, newU_Col, newV_Row, newV_Col, tile_n);    // m1=(a11+a11)(b11+b22)

			MatrixSum(a21, a22, tempA, newU_Row, newU_Col);            // a21+a22
			Strassen(tempA, b11, m2, newU_Row, newU_Col, newV_Row, newV_Col, tile_n);      // m2=(a21+a22)b11

			MatrixSub(b12, b22, tempB, newV_Row, newV_Col);            // b12-b22
			Strassen(a11, tempB, m3, newU_Row, newU_Col, newV_Row, newV_Col, tile_n);      // m3=a11(b12-b22)//////

			MatrixSub(b21, b11, tempB, newV_Row, newV_Col);            // b21-b11
			Strassen(a22, tempB, m4, newU_Row, newU_Col, newV_Row, newV_Col, tile_n);      // m4=a22(b21-11)

			MatrixSum(a11, a12, tempA, newU_Row, newU_Col);            //  a11+a12
			Strassen(tempA, b22, m5, newU_Row, newU_Col, newV_Row, newV_Col, tile_n); 	   // m5=(a11+a12)b22////////

			MatrixSub(a21, a11, tempA, newU_Row, newU_Col);            // a21-a11
			MatrixSum(b11, b12, tempB, newV_Row, newV_Col);            // b11+b12
			Strassen(tempA, tempB, m6, newU_Row, newU_Col, newV_Row, newV_Col, tile_n);    // m6=(a21-a11)(b11+b12)

			MatrixSub(a12, a22, tempA, newU_Row, newU_Col);            // a12-a22
			MatrixSum(b21, b22, tempB, newV_Row, newV_Col);            // b21+b22
			Strassen(tempA, tempB, m7, newU_Row, newU_Col, newV_Row, newV_Col, tile_n);    // m7 = (a12 - a22)(a12 - a22)


			// 위에서 계산된 m1~m7 결과로  c11 ~ c22 를 만든다.
			MatrixSum(m1, m4, tempAc, newU_Row/2, newV_Col/2); //m1 + m4
			MatrixSum(tempAc, m7, tempBc, newU_Row/2, newV_Col/2); //m1 + m4 + m7
			MatrixSub(tempBc, m5, c11, newU_Row/2, newV_Col/2); //c11 = m1 + m4 - m5 + m7

			MatrixSum(m3, m5, c12, newU_Row/2, newV_Col/2); //c12 = m3 + m5////////////////////////////////////////

			MatrixSum(m2, m4, c21, newU_Row/2, newV_Col/2); //c21 = m2 + m4

			MatrixSum(m1, m3, tempAc, newU_Row/2, newV_Col/2); //m1 + m3
			MatrixSum(tempAc, m6, tempBc, newU_Row/2, newV_Col/2); //m1 + m3 + m6
			MatrixSub(tempBc, m2, c22, newU_Row/2, newV_Col/2); //c22 = m1 + m3 - m2 + m6///////////////////////////////

			//재 병합
			Mergematrix(matrixM, c11, c12, c21, c22, newU_Row/2, newV_Col/2);
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
	cout << "Strassen-Winograd Convolutions filterTransform function! \n\n";

	int output_c = 4;
	int input_c = 32;
	int input_n = 8;
	int input_H = 6;
	int input_W = 6;

	// h[O_Ch][ln_Ch][KernelSize_h][KernelSize_w] 
	// h[4][3][3][3]
	// 임시 filter 값 
	vector<float> h(output_c * input_c * 3 * 3);
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
					//count++;
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

	//////////////////////////////////////////////////////////////////////////////////////////////////
	// filterTransform 수행
	filterTransform2(U, h, output_c, input_c);
	// inputTransform 수행
	inputTransform3(V, input, input_n, input_c, input_H, input_W);
	vector<float> M(input_n * output_c * output_H * output_W); //??

	// Strassen - Winograd 수행
	Strassen(U, V, M, output_c * 4, input_c * 4, input_c * 4, input_n * 4 * tile_n, tile_n);
	//////////////////////////////////////////////////////////////////////////////////////////////////

	long long end_usec = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();	
	int frame_sec = int(end_usec - start_usec);

	valueCheck(M, input_n, output_c, output_H, output_W, 1);
	//valueCheck(convOutput, input_n, output_c, output_H, output_W);

	

	cout << "======================================================" << endl;
	cout << frame_sec << "u sec (Strassen - Winograd Convolution)" << endl;

	int leftPadingSize = 0;
	int rightPadingSize = 0;
	int topPadingSize = 0;
	int bottomPadingSize = 0;

	int output_h_z = (input_H + topPadingSize + bottomPadingSize) - 2;
	int output_w_z = (input_W + leftPadingSize + rightPadingSize) - 2;
	// output[input_n][output_c][output_h_z][output_w_z] 
	vector<float> output(input_n * output_c * output_h_z * output_w_z);

	long long start_usec2 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();

	winogradConv2d(output, input, h, input_n, input_c, input_H, input_W, output_c, leftPadingSize, rightPadingSize, topPadingSize, bottomPadingSize);

	long long end_usec2 = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
	int frame_sec2 = int(end_usec2 - start_usec2);

	// output 값 체크
	cout << "======================================================" << endl;
	cout << "===== Winograd Convolution ===== \n";
	cout << frame_sec2 << "u sec ( Winograd Convolution)" << endl;

	valueCheck(output, input_n, output_c, output_h_z, output_w_z, 1);

	








}

