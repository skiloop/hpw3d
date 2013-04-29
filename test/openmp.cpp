#include <iostream>
#include <omp.h>
#include "../src/datastruct.h"

using namespace std;

int main(){
	int thread_count=5;
	data1d<double> data(20);
#pragma omp parallel for  num_threads(thread_count)
	for(int i=0;i<data.n;i++){
		printf("%f\n",data.p[i]);
	}
	return 0;
}


