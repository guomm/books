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

//#define MPIDEBUG
#define RESULT

int main(int argc, char* argv[])
{
#ifdef MPIDEBUG
	int myid, numproces;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numproces);

	string path("D://test//process");
	if (myid == MAINNODEID) {
		MainNode mainNode(numproces - 1);
		mainNode.compute();
		cout << "mainNode is over" << endl;
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
	int64 *SA = new int64[n];
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

	int64 *GRANK = new int64[n];
	n1 = MyIO<char>::getFileLenCh(fileName1);
	n2 = MyIO<char>::getFileLenCh(fileName2);
	//SACP[0] = n - 1;
	//int64 *GRANK2 = new int64[n2];
	MyIO<int64>::read(GRANK, fileName1, n1/8, std::ios_base::in | std::ios_base::binary, 0);

	/*cout << "GRank1 is :";
	for (int64 i = 0; i < n1; i++) {
		cout << GRANK[i] << ",";
	}
	cout << endl;*/

	MyIO<int64>::read(GRANK+ n1 / 8, fileName2, n2/8, std::ios_base::in | std::ios_base::binary, 0);

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
		if (SA[i] != sacmp[i])cout << "i is "<<i<<" error" <<"SA is "<< SA[i] <<" sacmp is "<<sacmp[i]<< endl;
	}
	cout << "success.." << endl;
	cin.get();
#endif
	
	//system("pause");  
	return 0;

#ifdef DEBUG_HEAP
	MinHeap hp;
	CommunicateData d1(1, 'a');
	hp.Insert(d1);
	//int max;
	std::cout << "目前堆中最大值为：" << hp.Get_Min().data << std::endl;
	hp.DeleteMin();
	std::cout << "删除的堆中最大的元素为：" << hp.Get_Min().data << std::endl;
	//cout << "删除后，现在堆中最大的元素为：" << hp.Get_Max() << endl;
	std::cout << "现在堆的大小为：" << hp.Get_Size() << std::endl;
	std::cout << "向堆中插入新元素" << std::endl;
//	CommunicateData<8> d1(1, 2);
	CommunicateData d2(1, 'c');
	CommunicateData d3(1, 'b');
	CommunicateData d4(1, 'f');
	CommunicateData d5(2, 'e');
	CommunicateData d6(2, 'a');
	CommunicateData d7(2, 'c');
	hp.Insert(d1);
	hp.Insert(d2);
	hp.Insert(d3);
	hp.Insert(d4);
	hp.Insert(d5);
	hp.Insert(d6);
	hp.Insert(d7);
	cout << (d1 < d2) << endl;
	cout << (d1 < d6) << endl;
	cout << (d2 > d1) << endl;
	cout << (d6 > d2) << endl;
	std::cout << "插入后现在堆的大小为：" << hp.Get_Size() << endl;

	std::cout << "现在由小到大输出堆中元素" << endl;
	//hp.printHeap();
	//cout << ".................\n";
	do
	{
		CommunicateData temp = hp.DeleteMin();

		temp.printA();
		//hp.printHeap();
		//cout << ".................\n";
	} while (hp.Get_Size()>0);

	//std::cout << endl;
#endif
	//cout << sizeof(CommunicateData) << endl;
	/*string path("D://test//process");
	checkGSA(path, 3);
	std::cin.get();*/
}
