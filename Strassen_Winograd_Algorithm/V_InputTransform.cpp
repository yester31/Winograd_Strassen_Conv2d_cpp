#include <io.h>
#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <iomanip>
#include <time.h>


using namespace std;

/***************************************************************************
	Input Transform 함수,  V = Bt d B
*****************************************************************************/

void inputTransform(vector<float> &V, vector<float>input, int input_n, int input_c, int input_h, int input_w) {

	int output_H = input_h - 2;
	int	output_W = input_w - 2;
	int area_o = output_H * output_W;

	int temp_i = input_c * input_h * input_w;
	int temp_v = input_c * 4 * area_o;
	for (int n = 0; n < input_n; n++) // N
	{
		int temp_i2 = n * temp_i;
		int temp_v2 = n * temp_v;

		for (int inch = 0; inch < input_c; inch++) // IC
		{
			int temp_i3 = inch * input_h * input_w + temp_i2;
			int temp_v3 = inch * 4 * area_o + temp_v2;

			for (int row = 0; row < output_H; row += 2)  // T_N
			{
				int row_idx1 = row * input_w;
				int row_idx2 = row_idx1 + input_w;
				int row_idx3 = row_idx2 + input_w;
				int row_idx4 = row_idx3 + input_w;

				int row_idxo1 = row * output_W;
				int row_idxo2 = row_idxo1 + output_W;

				for (int col = 0; col < output_W; col += 2)
				{

					//int u_idx = inch * 16 + temp_u2; // U idex
					int t_idx1 = temp_i3 + col + row_idx1;
					int t_idx2 = temp_i3 + col + row_idx2;
					int t_idx3 = temp_i3 + col + row_idx3;
					int t_idx4 = temp_i3 + col + row_idx4;

					float *d1 = &input[t_idx1];
					float *d2 = &input[t_idx1 + 1];
					float *d3 = &input[t_idx1 + 2];
					float *d4 = &input[t_idx1 + 3];

					float *d5 = &input[t_idx2];
					float *d6 = &input[t_idx2 + 1];
					float *d7 = &input[t_idx2 + 2];
					float *d8 = &input[t_idx2 + 3];

					float *d9 = &input[t_idx3];
					float *d10 = &input[t_idx3 + 1];
					float *d11 = &input[t_idx3 + 2];
					float *d12 = &input[t_idx3 + 3];

					float *d13 = &input[t_idx4];
					float *d14 = &input[t_idx4 + 1];
					float *d15 = &input[t_idx4 + 2];
					float *d16 = &input[t_idx4 + 3];

					/*
					cout << setw(5) << *d1 << " " << setw(5) << *d2 << " " << setw(5) << *d3 << " " << setw(5) << *d4 << endl;
					cout << setw(5) << *d5 << " " << setw(5) << *d6 << " " << setw(5) << *d7 << " " << setw(5) << *d8 << endl;
					cout << setw(5) << *d9 << " " << setw(5) << *d10 << " " << setw(5) << *d11 << " " << setw(5) << *d12 << endl;
					cout << setw(5) << *d13 << " " << setw(5) << *d14 << " " << setw(5) << *d15 << " " << setw(5) << *d16 << endl;
					cout << endl;
					*/
					
					float dd1 = *d11 - (*d3);
					float dd2 = *d2 - (*d10);
					float dd3 = *d7 + (*d11);
					float dd4 = *d6 + (*d10);
					float dd5 = *d7 - (*d11);
					float dd6 = *d10 - (*d6);
					float dd7 = *d15 - (*d7);
					float dd8 = *d6 - (*d14);

					V[temp_v3] = *d1 - *d9 + dd1; temp_v3++;
					V[temp_v3] = dd2 - dd1; temp_v3++;
					V[temp_v3] = -dd1 - dd2; temp_v3++;
					V[temp_v3] = dd2 - *d4 + *d12; temp_v3++;

					V[temp_v3] = *d5 + *d9 - dd3; temp_v3++;
					V[temp_v3] = dd4 + dd3; temp_v3++;
					V[temp_v3] = dd3 - dd4; temp_v3++;
					V[temp_v3] = dd4 - *d8 - *d12; temp_v3++;

					V[temp_v3] = *d9 - *d5 + dd5; temp_v3++;
					V[temp_v3] = dd6 - dd5; temp_v3++;
					V[temp_v3] = -(dd6 + dd5); temp_v3++;
					V[temp_v3] = dd6 + *d8 - *d12; temp_v3++;

					V[temp_v3] = *d5 - *d13 + dd7; temp_v3++;
					V[temp_v3] = dd8 - dd7; temp_v3++;
					V[temp_v3] = -dd7 - dd8; temp_v3++;
					V[temp_v3] = dd8 - *d8 + *d16; temp_v3++;
				


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

int main()
{
	cout << "Winograd Convolutions filterTransform function! \n\n";

	int batch_size = 5;
	int FeatureCh = 3;
	int input_H = 6;
	int input_W = 6;
	int area_i = input_H * input_W;
	//int out_ch = 4;

	// d[N][IC][input_H][input_W] 
	// d[1][3][10][10]
	// 임시 input 값 

	vector<float> input(batch_size * FeatureCh * area_i);
	float count = 1.f;
	int temp1 = FeatureCh * input_H * input_W;
	for (int ⁠n_idx = 0; ⁠n_idx < batch_size; ⁠n_idx++)
	{
		int temp2 = ⁠n_idx * temp1;
		for (int ⁠c_idx = 0; ⁠c_idx < FeatureCh; ⁠c_idx++)
		{
			int temp3 = ⁠c_idx * input_H * input_W + temp2;
			for (int ⁠h_idx = 0; ⁠h_idx < input_H; ⁠h_idx++)
			{
				int temp4 = ⁠h_idx * input_W + temp3;
				for (int w_idx = 0; w_idx < input_W; w_idx++)
				{
					int g_idx = w_idx + temp4;
					input[g_idx] = count;
					count++;
				}
			}
		}
	}


	// filater 값 체크
	valueCheck(input, batch_size, FeatureCh, input_H, input_W);

	// output[N][out_ch][4][4] 
	int output_H = input_H - 2;
	int output_W = input_W - 2;
	int area_o = output_H * output_W;

	vector<float> V(batch_size * FeatureCh * 4 * area_o);


	inputTransform(V, input, batch_size, FeatureCh, input_H, input_W);


	// V 값 체크
	valueCheck(V, FeatureCh, batch_size,  area_o/4, 16);

	//cout << " = ======= = = = == = ==" << endl;

	//valueCheck(V, batch_size, FeatureCh,  area_o / 4, 16);


	return 0;
}