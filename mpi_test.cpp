// mpi_test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>  
#include <vector>
#include <string>
#include "childnode.h"
#include "mainnode.h"
#include "sais.h"
using namespace std;

//#define DEBUG_HEAP

#define MPIDEBUG
//#define RESULT

//#define LOSERTREE
int main(int argc, char* argv[])
{
#ifdef MPIDEBUG
	int myid, numproces;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numproces);
	//cout << "sendData " << sizeof(SendData<unsigned char>) << "," << sizeof(SendData<uint_type>) << endl;
	//cout << "recvData " << sizeof(RecvData<unsigned char>) << "," << sizeof(RecvData<uint_type>) << endl;
	string path("D://test//process");
	if (myid == MAINNODEID) {
		MainNode mainNode(numproces - 1);
		mainNode.compute();
		cout << "mainNode is over, communicate traffic is "<<(double)(mainNode.commColumn)/1024/1024 <<" M, communicate time is "<<(double)(mainNode.commTime)/CLOCKS_PER_SEC << endl;
	}
	else {
		ChildNode childNode(myid);
		childNode.compute(path);
		cout << "child " <<myid<<" is over"<< endl;
	}
	MPI_Finalize();
#endif // !MPIDEBUG

#ifdef RESULT
	string fileName1("D://test//process1.txt");
	string fileName2("D://test//process2.txt");
	int64 n1 = MyIO<char>::getFileLenCh(fileName1);
	int64 n2 = MyIO<char>::getFileLenCh(fileName2);
	int64 n = n1 + n2 + 3;
	unsigned char *data = new unsigned char[n];
	MyIO<unsigned char>::read(data, fileName1, n1, std::ios_base::in | std::ios_base::binary, 0);
	data[n1] = 1;
	MyIO<unsigned char>::read(data+n1+1, fileName2, n2, std::ios_base::in | std::ios_base::binary, 0);
	data[n - 1] = 0;
	data[n - 2] = 2;
	uint_type *SA = new uint_type[n];
	SA_IS(data,SA,n,256,sizeof(char),0);

	/*cout << "data is :";
	for (int64 i = 0; i < n; i++) {
		cout << (int32)data[i] << ",";
	}
	cout << endl;
	cout << "SA is :";
	for (int64 i = 0; i < n; i++) {
	cout << SA[i] << ",";
	}
	cout << endl;*/

	fileName1 = "D://test//result1.txt";
	fileName2 = "D://test//result2.txt";

	uint_type *GRANK = new uint_type[n];
	n1 = MyIO<char>::getFileLenCh(fileName1);
	n2 = MyIO<char>::getFileLenCh(fileName2);
	//SACP[0] = n - 1;
	//int64 *GRANK2 = new int64[n2];
	MyIO<uint_type>::read(GRANK, fileName1, n1/sizeof(uint_type), std::ios_base::in | std::ios_base::binary, 0);

	/*cout << "GRank1 is :";
	for (int64 i = 0; i < n1; i++) {
		cout << GRANK[i] << ",";
	}
	cout << endl;*/

	MyIO<uint_type>::read(GRANK+ n1 / sizeof(uint_type), fileName2, n2/ sizeof(uint_type), std::ios_base::in | std::ios_base::binary, 0);

	/*cout << "GRank2 is :";
	for (int64 i = 0; i < n2; i++) {
		cout << GRANK2[i] << ",";
	}
	cout << endl;

	for (int64 i = 0,j=n1+1; i < n2; i++) {
		GRANK[j++]= GRANK2[i];
	}*/

	GRANK[n - 1] = 0;
	int64 *sacmp = new int64[n];
	/*cout << "GRANK is :";
	for (int64 i = 0; i < n; i++) {
		cout << GRANK[i] << ",";
		sacmp[i] = 0;
	}
	cout << endl;*/

	
	for (int64 i = 0; i < n; i++) {
		sacmp[GRANK[i]] = i;
	}

	/*cout << "sacmp is :";
	for (int64 i = 0; i < n; i++) {
		cout << sacmp[i] << ",";
	}
	cout << endl;*/


	for (int64 i = 0; i < n; i++) {
		if (SA[i] != sacmp[i])cout << "i is "<<i<<" error " <<"SA is "<< SA[i] <<" sacmp is "<<sacmp[i]<< endl;
	}
	cout << "success.." << endl;
	cin.get();
#endif
	
#ifdef LOSERTREE

#endif // LOSERTREE







	//system("pause");  
	return 0;


}
