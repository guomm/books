// mpi_test.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>  
#include <vector>
#include <string>
#include "childnode.h"
#include "mainnode.h"
using namespace std;

//#define DEBUG_HEAP

int main(int argc, char* argv[])
{
	int myid, numproces;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numproces);

	string path("D://test//process");
	//pasis(path, myid, numproces);
	//checkGSA(path, 3);
	//cout << myid << endl;
	//cout << numproces << endl;
	if (myid == MAINNODEID) {
		MainNode mainNode(numproces - 1);
	//	cout << "....MainNode..." << endl;
		mainNode.compute();
	}
	else {
		ChildNode childNode(myid);
		//cout << "....childNode..." << endl;
		childNode.compute(path);
	}
	MPI_Finalize();
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
