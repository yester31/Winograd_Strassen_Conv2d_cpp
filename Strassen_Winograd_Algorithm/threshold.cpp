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


void getThreshold(int n)
{
	int th;
	double k = floor(log(n) / log(2) - 4);// lgn -4
	th = (int)floor(n / pow(2.0, k)) + 1;
	cout << "n :: " << n << "  k ::  " << k << "  th :: " << th << endl;
}


int main()
{
	getThreshold(4);
	getThreshold(16);
	getThreshold(32);
	getThreshold(64);
	getThreshold(128);



	return 0;
}