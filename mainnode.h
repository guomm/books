#ifndef  __MAINNODE
#define  __MAINNODE
#include <iostream>
#include <time.h>
#include "mycommon.h"
//#include "mystack.h"

//#define tget(j,i) ( (recvBuf[j+i/8]&mask[(i)%8]) ? 1 : 0 )
//#define isGget(i) ( (recvBuf[(i)/8]&mask[(i)%8]) ? 1 : 0 )

//#define SADEBUG
#define SALDEBUG
//#define SASDEBUG
//#define RENAMEDEBUG

#define LEVELL -1
//#define PQQ
using namespace std;

class MainNode {
private:
	int32 numOfChild;
	MPI_Status status;
	int64 *sizeOfChild;

	uint_type *recvBuf;
	RecvData<unsigned char> *recvData0;
	RecvData<uint_type> *recvData1;
	uint_type *sendBuf;;
	//RecvData *recvData;
	MarkData *markData;

	int32 *readCur;

	LoserTree<RecvData<unsigned char>> loserTree0;
	VictoryTree<unsigned char> victoryTree0;
	LoserTree<RecvData<uint_type>> loserTree1;
	VictoryTree<uint_type> victoryTree1;
	LoserTree<uint_type> rename;
	//MarkData *markData;
public:
	int64 commColumn = 0;
	clock_t commTime, sTime;
	MainNode(int32 num_child_node) {
		numOfChild = num_child_node;
		//	uint_type recvLen = 2 * commSize * sizeof(uint_type) + commSize + commSize * sizeof(uint_type);
		int32 sumLen = commSize * numOfChild;
		recvBuf = new uint_type[sumLen];
		sendBuf = new uint_type[sumLen];
		recvData0 = new RecvData<unsigned char>[sumLen];
		recvData1 = new RecvData<uint_type>[sumLen];
		readCur = new int32[sumLen];
		sizeOfChild = new int64[numOfChild];
		markData = new MarkData[numOfChild];
		for (int32 i = 0; i < numOfChild; i++) {
			sizeOfChild[i] = 0;
			MarkData temp;
			markData[i] = temp;
		}
		for (int32 i = 0; i < sumLen; i++) {
			readCur[i] = 0;
			sendBuf[i] = 0;
			recvBuf[i] = 0;
			RecvData<unsigned char> tt;
			recvData0[i] = tt;
			RecvData<uint_type> hh;
			recvData1[i] = hh;
		}

		RecvData<unsigned char> maxData0(0, MAXU, 1, (std::numeric_limits<unsigned char>::max)());
		RecvData<uint_type> maxData1(0, MAXU, 1, (std::numeric_limits<uint_type>::max)());
		RecvData<unsigned char> minData0(0, MINU, 1, (std::numeric_limits<unsigned char>::min)());
		RecvData<uint_type> minData1(0, MINU, 1, (std::numeric_limits<uint_type>::min)());
		//MinData.setMember(0, MINU, 0, (std::numeric_limits<T>::min)());

		loserTree0.setMAXMIN(minData0, maxData0);
		loserTree1.setMAXMIN(minData1, maxData1);
		rename.setMAXMIN((std::numeric_limits<uint_type>::min)(), (std::numeric_limits<uint_type>::max)());

		loserTree0.initTree(numOfChild);
		victoryTree0.initTree(numOfChild);
		loserTree1.initTree(numOfChild);
		victoryTree1.initTree(numOfChild);
		rename.initTree(numOfChild);
	}

	void compute() {
		computeGSA(0);
	}

	void reName(int32 level, int64 &rank) {
#ifdef RENAMEDEBUG
		cout << "主节点开始重命名: numOfChild is " << numOfChild << endl;
#endif
		uint_type  offset = 0;
		rename.initTree(numOfChild);

		uint_type i, j = 0;
		for (i = 0; i < numOfChild; i++) {

			offset = i * commSize;
			readCur[i] = offset;
			sTime = clock();
			MPI_Recv((char *)(recvBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			commTime += clock() - sTime;
			commColumn += commSize * sizeof(uint_type);

#ifdef RENAMEDEBUG
			cout << "重命名：从 " << i << "接收到的数据：" << endl;
			for (j = 0; j < commSize; j++) {
				rnData[offset + j].printA();
			}
#endif
		}

		for (i = 0; i < numOfChild; i++) {
			rename.setVal(i, recvBuf[readCur[i]]);
		}
		rename.CreateLoserTree();

		uint_type minVal = 0;
		while (!rename.isFinish()) {
			i = rename.getMinIndx();
			offset = i*commSize;
			if (recvBuf[readCur[i]] == 0) {
				sTime = clock();
				//发送数据
				MPI_Send(((char *)(sendBuf + offset)), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				commColumn += commSize * sizeof(uint_type);
				commTime += clock() - sTime;
				rename.setPosFinsih(i);
			}
			else {
				if (minVal != recvBuf[readCur[i]]) {
					rank++;
					minVal = recvBuf[readCur[i]];
				}

				sendBuf[readCur[i]++] = rank;
				if (readCur[i] == (offset + commSize)) {
					//发送数据
					sTime = clock();
					MPI_Send((char *)(sendBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					commColumn += commSize * sizeof(uint_type);
#ifdef RENAMEDEBUG
					std::cout << "重命名:向 " << i << "发送的数据：" << endl;
					for (j = 0; j < commSize; j++) {
						cout << sendBuf[offset + j] << ",";
					}
					cout << endl;
#endif

					MPI_Recv((char *)(recvBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
					commColumn += commSize * sizeof(uint_type);
					commTime += clock() - sTime;
					readCur[i] = offset;

#ifdef RENAMEDEBUG
					std::cout << "重命名：从 " << i << "接收到的数据：" << endl;
					for (j = 0; j < commSize; j++) {
						rnData[offset + j].printA();
					}
#endif
				}
				//rnPQ.push(rnData[readCur[i]]);
				rename.setAndFresh(i, recvBuf[readCur[i]]);
			}
		}

	}

	void reNameNoCycle(int32 level, int64 &rank) {
#ifdef RENAMEDEBUG
		cout << "主节点开始重命名: numOfChild is " << numOfChild << endl;
#endif
		uint_type  offset = 0;
		rename.initTree(numOfChild);

		uint_type i, j = 0;
		for (i = 0; i < numOfChild; i++) {

			offset = i * commSize;
			readCur[i] = offset;
			sTime = clock();
			MPI_Recv((char *)(recvBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			commTime += clock() - sTime;
			commColumn += commSize * sizeof(uint_type);

#ifdef RENAMEDEBUG
			cout << "重命名：从 " << i << "接收到的数据：" << endl;
			for (j = 0; j < commSize; j++) {
				rnData[offset + j].printA();
			}
#endif
		}

		for (i = 0; i < numOfChild; i++) {
			rename.setVal(i, recvBuf[readCur[i]]);
		}
		rename.CreateLoserTree();

		//uint_type minVal = 0;
		while (!rename.isFinish()) {
			i = rename.getMinIndx();
			offset = i*commSize;
			if (recvBuf[readCur[i]] == 0) {
				sTime = clock();
				//发送数据
				MPI_Send(((char *)(sendBuf + offset)), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				commColumn += commSize * sizeof(uint_type);
				commTime += clock() - sTime;
				rename.setPosFinsih(i);
			}
			else {
				rank++;
				sendBuf[readCur[i]++] = rank;
				if (readCur[i] == (offset + commSize)) {
					//发送数据
					sTime = clock();
					MPI_Send((char *)(sendBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					commColumn += commSize * sizeof(uint_type);
#ifdef RENAMEDEBUG
					std::cout << "重命名:向 " << i << "发送的数据：" << endl;
					for (j = 0; j < commSize; j++) {
						cout << sendBuf[offset + j] << ",";
					}
					cout << endl;
#endif

					MPI_Recv((char *)(recvBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
					commColumn += commSize * sizeof(uint_type);
					commTime += clock() - sTime;
					readCur[i] = offset;

#ifdef RENAMEDEBUG
					std::cout << "重命名：从 " << i << "接收到的数据：" << endl;
					for (j = 0; j < commSize; j++) {
						rnData[offset + j].printA();
					}
#endif
				}
				//rnPQ.push(rnData[readCur[i]]);
				rename.setAndFresh(i, recvBuf[readCur[i]]);
			}
		}

	}
	void computeGSA(int32 level) {

		uint_type n = 0;
		//接受每个子节点的大小
		int32 i, j = 0;
		for (i = 0; i < numOfChild; i++) {
			MPI_Recv(&(sizeOfChild[i]), 1, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD, &status);
			n += sizeOfChild[i];
			MarkData temp;
			markData[i] = temp;
		}
#ifdef	SADEBUG
		for (i = 0; i < numOfChild; i++) {
			std::cout << "the size of " << i << "is " << sizeOfChild[i] << endl;
		}

		std::cout << "n is " << n << endl;
#endif

		//开始计算SA
		uint_type gPos = numOfChild;

		if (level == 0)computeGSAlTree0(level, gPos);
		else computeGSAlTree1(level, gPos);

		gPos = n + 1;

		if (level == 0)computeGSAsTree0(level, gPos);
		else computeGSAsTree1(level, gPos);

		bool mark = 0;
		//int32 i = 0;
		for (i = 0; i < numOfChild; i++) {
			if (markData[i].isCycle == 1)mark = 1;
		}
		if (mark) {
			for (i = 0; i < numOfChild; i++)markData[i].isCycle = 1;
		}
		for (i = 0; i < numOfChild; i++) {
			//cout << "i  " << i << "--------------" << endl;
			MPI_Send((char *)(markData + i), 1 * sizeof(MarkData), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
		}
		//bool isCycle = 0;
		int64 rank = 0;
		if (mark) {
			reName(level, rank);
			computeGSA(level + 1);
		}
		else {
			reNameNoCycle(level, rank);
		}
		
		//rename
		/*reName(level, rank);

		if (isCycle) {
			for (i = 0; i < numOfChild; i++)MPI_Send(&rank, 1, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD);
			computeGSA(level + 1);
		}
		else {
			rank = 0;
			for (i = 0; i < numOfChild; i++)MPI_Send(&rank, 1, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD);
		}*/

		gPos = numOfChild;

		if (level == 0)computeGSAlTree0(level, gPos);
		else computeGSAlTree1(level, gPos);

		gPos = n + 1;

		if (level == 0)computeGSAsTree0(level, gPos);
		else computeGSAsTree1(level, gPos);
	}

	void computeGSAlTree0(int32 level, uint_type gPos) {
		uint_type i, j, offset = 0;
		//cout << "kkkk" << endl;
		loserTree0.initTree(numOfChild);
		for (i = 0; i < numOfChild; i++) {
			readCur[i] = i * commSize;
			sTime = clock();
			MPI_Recv((char*)(recvData0 + readCur[i]), commSize * sizeof(RecvData<unsigned char>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			commColumn += commSize * sizeof(RecvData<unsigned char>);
			commTime += clock() - sTime;
			/*SendData<unsigned char> aa[commSize];
			MPI_Recv((char*)aa, commSize * sizeof(SendData<unsigned char>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			cout << i << " 接收到的数据是：" << endl;
			for (j = 0; j < commSize; j++) {
			aa[j].printA();
			}*/
		}

#ifdef SALDEBUG
		if (level > LEVELL) {

			for (i = 0; i < numOfChild; i++) {
				cout << "level is " << level << " computeSAL get data from " << i << " ：" << endl;
				for (j = 0; j < commSize; j++) {
					recvData0[i * commSize + j].printA();
				}
			}
		}

#endif

		for (i = 0; i < numOfChild; i++) {
			//mndata0[i].setMember(recvData0[readCur[i]], i);
			//PQL0.push(mndata0[i]);
			loserTree0.setVal(i, recvData0[readCur[i]]);
		}

		//MNData<unsigned char> temp, minVal;
		//int32 minPos = 0;
		loserTree0.CreateLoserTree();

		RecvData<unsigned char> minVal;
		while (!loserTree0.isFinish()) {


			i = loserTree0.getMinIndx();
			//cout << "当前的最小值是：" <<i<< endl;
			//recvData0[readCur[i]].printA();

			offset = i * commSize;
			if (recvData0[readCur[i]].prePos == commSize) {
				sTime = clock();
				MPI_Send((char *)(sendBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				commColumn += commSize * sizeof(uint_type);
				commTime += clock() - sTime;
				loserTree0.setPosFinsih(i);
#ifdef SALDEBUG
				if (level > LEVELL) {
					cout << "level is " << level << " computeSAL last send data to  " << i << " ：" << endl;
					for (j = 0; j < commSize; j++) {
						cout << sendBuf[i * commSize + j] << ",";
					}
					cout << endl;
				}

#endif
			}
			else {
				if (minVal != recvData0[readCur[i]]) {
					gPos++;
					minVal.setMember(recvData0[readCur[i]]);
				}

				j = recvData0[readCur[i]].prePos;
				//cout << "hhh  " <<j<< endl;
				if (j != -1)recvData0[offset + j].suffixGIndex = gPos;
				sendBuf[readCur[i]++] = gPos;
				//PQL0.pop();
				//cout << ",,," << endl;
				if (readCur[i] == (offset + commSize)) {
					//发送数据
					sTime = clock();
					MPI_Send((char *)(sendBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					commColumn += commSize * sizeof(uint_type);

#ifdef SALDEBUG
					if (level > LEVELL) {
						cout << "level is " << level << " computeSAL send data to  " << i << " ：" << endl;
						for (j = 0; j < commSize; j++) {
							cout << sendBuf[i * commSize + j] << ",";
						}
						cout << endl;
					}

#endif
					readCur[i] = offset;
					//接收数据
					MPI_Recv((char*)(recvData0 + readCur[i]), commSize * sizeof(RecvData<unsigned char>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
					commColumn += commSize * sizeof(RecvData<unsigned char>);
					commTime += clock() - sTime;

#ifdef SALDEBUG
					if (level > LEVELL) {
						cout << "level is " << level << " computeSAL get data from " << i << " ：" << endl;
						for (int64 t = 0; t < numOfChild; t++) {
							for (j = 0; j < commSize; j++) {
								recvData0[t * commSize + j].printA();
							}
						}
					}

#endif

				}
				loserTree0.setAndFresh(i, recvData0[readCur[i]]);
				//cout << "i is " << i << endl;
				//mndata0[i].setMember(recvData0[readCur[i]], i);
				//mndata0[i].printA();
				//MNData<unsigned char> temp(recvData0[readCur[i]], i);
				//PQL0.push(mndata0[i]);
			}
		}

		//cout << "computeGSAlTree0 is over" << endl;
	}
	void computeGSAlTree1(int32 level, uint_type gPos) {
		uint_type i, j, offset = 0;
		//cout << "kkkk" << endl;
		loserTree1.initTree(numOfChild);
		for (i = 0; i < numOfChild; i++) {
			readCur[i] = i * commSize;
			sTime = clock();
			MPI_Recv((char*)(recvData1 + readCur[i]), commSize * sizeof(RecvData<uint_type>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			commColumn += commSize * sizeof(RecvData<uint_type>);
			commTime += clock() - sTime;
			/*SendData<unsigned char> aa[commSize];
			MPI_Recv((char*)aa, commSize * sizeof(SendData<unsigned char>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			cout << i << " 接收到的数据是：" << endl;
			for (j = 0; j < commSize; j++) {
			aa[j].printA();
			}*/
		}

#ifdef SALDEBUG
		if (level > LEVELL) {

			for (i = 0; i < numOfChild; i++) {
				cout << "level is " << level << " computeSAL get data from " << i << " ：" << endl;
				for (j = 0; j < commSize; j++) {
					recvData1[i * commSize + j].printA();
				}
			}
		}

#endif

		for (i = 0; i < numOfChild; i++) {
			//mndata0[i].setMember(recvData0[readCur[i]], i);
			//PQL0.push(mndata0[i]);
			loserTree1.setVal(i, recvData1[readCur[i]]);
		}

		//MNData<unsigned char> temp, minVal;
		//int32 minPos = 0;
		loserTree1.CreateLoserTree();

		RecvData<uint_type> minVal;
		while (!loserTree1.isFinish()) {


			i = loserTree1.getMinIndx();
			//cout << "当前的最小值是：" <<i<< endl;
			//recvData0[readCur[i]].printA();

			offset = i * commSize;
			if (recvData1[readCur[i]].prePos == commSize) {
				sTime = clock();
				MPI_Send((char *)(sendBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				commColumn += commSize * sizeof(uint_type);
				commTime += clock() - sTime;
				loserTree1.setPosFinsih(i);
#ifdef SALDEBUG
				if (level > LEVELL) {
					cout << "level is " << level << " computeSAL last send data to  " << i << " ：" << endl;
					for (j = 0; j < commSize; j++) {
						cout << sendBuf[i * commSize + j] << ",";
					}
					cout << endl;
				}

#endif
			}
			else {
				if (minVal != recvData1[readCur[i]]) {
					gPos++;
					minVal.setMember(recvData1[readCur[i]]);
				}

				j = recvData1[readCur[i]].prePos;
				//cout << "hhh  " <<j<< endl;
				if (j != -1)recvData1[offset + j].suffixGIndex = gPos;
				sendBuf[readCur[i]++] = gPos;
				//PQL0.pop();
				//cout << ",,," << endl;
				if (readCur[i] == (offset + commSize)) {
					//发送数据
					sTime = clock();
					MPI_Send((char *)(sendBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					commColumn += commSize * sizeof(uint_type);

#ifdef SALDEBUG
					if (level > LEVELL) {
						cout << "level is " << level << " computeSAL send data to  " << i << " ：" << endl;
						for (j = 0; j < commSize; j++) {
							cout << sendBuf[i * commSize + j] << ",";
						}
						cout << endl;
					}

#endif
					readCur[i] = offset;
					//接收数据
					MPI_Recv((char*)(recvData1 + readCur[i]), commSize * sizeof(RecvData<uint_type>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
					commColumn += commSize * sizeof(RecvData<uint_type>);
					commTime += clock() - sTime;

#ifdef SALDEBUG
					if (level > LEVELL) {
						cout << "level is " << level << " computeSAL get data from " << i << " ：" << endl;
						for (int64 t = 0; t < numOfChild; t++) {
							for (j = 0; j < commSize; j++) {
								recvData1[t * commSize + j].printA();
							}
						}
					}

#endif

				}
				loserTree1.setAndFresh(i, recvData1[readCur[i]]);
				//cout << "i is " << i << endl;
				//mndata0[i].setMember(recvData0[readCur[i]], i);
				//mndata0[i].printA();
				//MNData<unsigned char> temp(recvData0[readCur[i]], i);
				//PQL0.push(mndata0[i]);
			}
		}

		//cout << "computeGSAlTree0 is over" << endl;
	}
	void computeGSAsTree0(int32 level, uint_type gPos) {
		uint_type i, j, offset = 0;
		//	cout << "computeGSAsTree0 is start" << endl;
		victoryTree0.initTree(numOfChild);
		for (i = 0; i < numOfChild; i++) {
			readCur[i] = i * commSize;
			sTime = clock();
			MPI_Recv((char*)(recvData0 + readCur[i]), commSize * sizeof(RecvData<unsigned char>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			commColumn += commSize * sizeof(RecvData<unsigned char>);
			commTime += clock() - sTime;
			/*SendData<unsigned char> aa[commSize];
			MPI_Recv((char*)aa, commSize * sizeof(SendData<unsigned char>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			cout << i << " 接收到的数据是：" << endl;
			for (j = 0; j < commSize; j++) {
			aa[j].printA();
			}*/
		}

#ifdef SASDEBUG
		if (level > LEVELL) {

			for (i = 0; i < numOfChild; i++) {
				cout << "level is " << level << " computeSAS get data from " << i << " ：" << endl;
				for (j = 0; j < commSize; j++) {
					recvData0[i * commSize + j].printA();
				}
			}
		}

#endif

		for (i = 0; i < numOfChild; i++) {
			//mndata0[i].setMember(recvData0[readCur[i]], i);
			//PQL0.push(mndata0[i]);
			victoryTree0.setVal(i, recvData0[readCur[i]]);
		}

		//MNData<unsigned char> temp, minVal;
		//int32 minPos = 0;
		victoryTree0.CreateVictoryTree();

		RecvData<unsigned char> minVal;
		int64 minNodeId = 0;
		while (!victoryTree0.isFinish()) {
			i = victoryTree0.getMaxIndx();
			//cout << "当前最大值是：" << i << endl;
			//recvData0[readCur[i]].printA();
			offset = i * commSize;
			if (recvData0[readCur[i]].prePos == commSize) {
				sTime = clock();
				MPI_Send((char *)(sendBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				commColumn += commSize * sizeof(uint_type);
				commTime += clock() - sTime;
				victoryTree0.setPosFinsih(i);
			}
			else {
				if (minVal != recvData0[readCur[i]]) {
					gPos--;
					minVal.setMember(recvData0[readCur[i]]);
					minNodeId = i;
				}
				else if (minNodeId == i) {
					markData[i].isSelfCom = 1;
				}
				else {
					markData[i].isCycle = 1;
				}

				j = recvData0[readCur[i]].prePos;
				//cout << "hhh  " <<j<< endl;
				if (j != -1)recvData0[offset + j].suffixGIndex = gPos;
				sendBuf[readCur[i]++] = gPos;
				//PQL0.pop();
				//cout << ",,," << endl;
				if (readCur[i] == (offset + commSize)) {
					//发送数据
					sTime = clock();
					MPI_Send((char *)(sendBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					commColumn += commSize * sizeof(uint_type);

#ifdef SASDEBUG
					if (level > LEVELL) {
						cout << "level is " << level << " computeSAS send data to  " << i << " ：" << endl;
						for (j = 0; j < commSize; j++) {
							cout << sendBuf[i * commSize + j] << ",";
						}
						cout << endl;
					}

#endif
					readCur[i] = offset;
					//接收数据
					MPI_Recv((char*)(recvData0 + offset), commSize * sizeof(RecvData<unsigned char>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
					commColumn += commSize * sizeof(RecvData<unsigned char>);
					commTime += clock() - sTime;

#ifdef SASDEBUG
					if (level > LEVELL) {
						cout << "level is " << level << " computeSAS get data from " << i << " ：" << endl;
						for (int64 t = 0; t < numOfChild; t++) {
							for (j = 0; j < commSize; j++) {
								recvData0[t * commSize + j].printA();
							}
						}
					}

#endif

				}
				victoryTree0.setAndFresh(i, recvData0[readCur[i]]);
				//cout << "i is " << i << endl;
				//mndata0[i].setMember(recvData0[readCur[i]], i);
				//mndata0[i].printA();
				//MNData<unsigned char> temp(recvData0[readCur[i]], i);
				//PQL0.push(mndata0[i]);
			}
		}

	}
	void computeGSAsTree1(int32 level, uint_type gPos) {
		uint_type i, j, offset = 0;
		//cout << "computeGSAsTree0 is start" << endl;
		victoryTree1.initTree(numOfChild);
		for (i = 0; i < numOfChild; i++) {
			readCur[i] = i * commSize;
			sTime = clock();
			MPI_Recv((char*)(recvData1 + readCur[i]), commSize * sizeof(RecvData<uint_type>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			commColumn += commSize * sizeof(RecvData<uint_type>);
			commTime += clock() - sTime;
			/*SendData<unsigned char> aa[commSize];
			MPI_Recv((char*)aa, commSize * sizeof(SendData<unsigned char>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			cout << i << " 接收到的数据是：" << endl;
			for (j = 0; j < commSize; j++) {
			aa[j].printA();
			}*/
		}

#ifdef SASDEBUG
		if (level > LEVELL) {

			for (i = 0; i < numOfChild; i++) {
				cout << "level is " << level << " computeSAS get data from " << i << " ：" << endl;
				for (j = 0; j < commSize; j++) {
					recvData0[i * commSize + j].printA();
				}
			}
		}

#endif

		for (i = 0; i < numOfChild; i++) {
			//mndata0[i].setMember(recvData0[readCur[i]], i);
			//PQL0.push(mndata0[i]);
			victoryTree1.setVal(i, recvData1[readCur[i]]);
		}

		//MNData<unsigned char> temp, minVal;
		//int32 minPos = 0;
		victoryTree1.CreateVictoryTree();

		RecvData<uint_type> minVal;
		int64 minNodeId = 0;
		while (!victoryTree1.isFinish()) {
			i = victoryTree1.getMaxIndx();
			//cout << "当前最大值是：" << i << endl;
			//recvData0[readCur[i]].printA();
			offset = i * commSize;
			if (recvData1[readCur[i]].prePos == commSize) {
				sTime = clock();
				MPI_Send((char *)(sendBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				commColumn += commSize * sizeof(uint_type);
				commTime += clock() - sTime;
				victoryTree1.setPosFinsih(i);
			}
			else {
				if (minVal != recvData1[readCur[i]]) {
					gPos--;
					minVal.setMember(recvData1[readCur[i]]);
					minNodeId = i;
				}
				else if (minNodeId == i) {
					markData[i].isSelfCom = 1;
				}
				else {
					markData[i].isCycle = 1;
				}

				j = recvData1[readCur[i]].prePos;
				//cout << "hhh  " <<j<< endl;
				if (j != -1)recvData1[offset + j].suffixGIndex = gPos;
				sendBuf[readCur[i]++] = gPos;
				//PQL0.pop();
				//cout << ",,," << endl;
				if (readCur[i] == (offset + commSize)) {
					//发送数据
					sTime = clock();
					MPI_Send((char *)(sendBuf + offset), commSize * sizeof(uint_type), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					commColumn += commSize * sizeof(uint_type);

#ifdef SASDEBUG
					if (level > LEVELL) {
						cout << "level is " << level << " computeSAS send data to  " << i << " ：" << endl;
						for (j = 0; j < commSize; j++) {
							cout << sendBuf[i * commSize + j] << ",";
						}
						cout << endl;
					}

#endif
					readCur[i] = offset;
					//接收数据
					MPI_Recv((char*)(recvData1 + offset), commSize * sizeof(RecvData<uint_type>), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
					commColumn += commSize * sizeof(RecvData<uint_type>);
					commTime += clock() - sTime;

#ifdef SASDEBUG
					if (level > LEVELL) {
						cout << "level is " << level << " computeSAS get data from " << i << " ：" << endl;
						for (int64 t = 0; t < numOfChild; t++) {
							for (j = 0; j < commSize; j++) {
								recvData0[t * commSize + j].printA();
							}
						}
					}

#endif

				}
				victoryTree1.setAndFresh(i, recvData1[readCur[i]]);
				//cout << "i is " << i << endl;
				//mndata0[i].setMember(recvData0[readCur[i]], i);
				//mndata0[i].printA();
				//MNData<unsigned char> temp(recvData0[readCur[i]], i);
				//PQL0.push(mndata0[i]);
			}
		}

	}
	~MainNode()
	{
		delete[] sizeOfChild;

		delete[] recvBuf;
		delete[] sendBuf;
		delete[] recvData0;
		delete[] recvData1;
		delete[] readCur;
	}
};
#endif