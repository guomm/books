#ifndef  __MAINNODE
#define  __MAINNODE
#include <iostream>
#include <unordered_map>
#include "mycommon.h"
//#include "mystack.h"

#define tget(j,i) ( (recvBuf[j+i/8]&mask[(i)%8]) ? 1 : 0 )
//#define isGget(i) ( (recvBuf[(i)/8]&mask[(i)%8]) ? 1 : 0 )

#define SADEBUG
//#define SALDEBUG
//#define SASDEBUG
//#define RENAMEDEBUG

#define LEVELL 7
//#define PQQ
using namespace std;

class MainNode {
private:
	int32 numOfChild;
	MPI_Status status;
	int64 *sizeOfChild;
	
	unsigned char *recvBuf;
	int64 *sendBuf;
	unordered_map<int64, int64> *map;
	RecvData *recvData;
	priority_queue<RecvData, vector<RecvData>, greater<RecvData>>PQL;
	priority_queue<RecvData> PQS;
	priority_queue<ReNameData, vector<ReNameData>, greater<ReNameData>>rnPQ;
	ReNameData *rnData;
	int32 *readCur;
public:
	MainNode(int32 num_child_node){
		numOfChild = num_child_node;
		int64 recvLen = 2 * commSize * sizeof(int64) + commSize + commSize * sizeof(int64);
		int32 sumLen = commSize * numOfChild;
		recvBuf = new unsigned char[recvLen * numOfChild];
		sendBuf = new int64[sumLen];
		map = new unordered_map<int64, int64>[numOfChild];
		recvData = new RecvData[sumLen];
		rnData = new ReNameData[sumLen];
		readCur = new int32[sumLen];
		sizeOfChild = new int64[numOfChild];
		for (int32 i = 0; i < numOfChild; i++) {
			sizeOfChild[i] = 0;
			unordered_map<int64, int64> mapTemp;
			map[i] = mapTemp;
		}
		for (int32 i = 0; i < sumLen; i++) {
			readCur[i] = 0;
			sendBuf[i] = 0;
			RecvData temp;
			recvData[i] = temp;
			ReNameData rnTemp;
			rnData[i] = rnTemp;
		}
	}

	void compute() {
		computeGSA(0);
	}

	void reName(int32 level, int64 &rank, bool &isCycle) {
#ifdef RENAMEDEBUG
		cout << "主节点开始重命名: numOfChild is " << numOfChild << endl;
#endif
		int64  offset = 0;
		//rank = 0;

		while (!rnPQ.empty())rnPQ.pop();

		int64 *curBuf;
		int64 i, j = 0;
		//MPI_Status status;
		for (i = 0; i < numOfChild; i++) {
			//	sumChildSize[i] = sumSize;
			//	sumSize += sizeOfChild[i];
			
			offset = i * commSize;
			readCur[i] = offset;
			//sendCur[i] = offset;
			MPI_Recv((char *)(recvBuf + offset * sizeof(int64)), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			curBuf = (int64 *)(recvBuf + offset * sizeof(int64));
			for (j = 0; j < commSize; j++) {
				rnData[offset + j].setMember(curBuf[j], i);
				if (curBuf[j] == EMPTY) break;
			}

#ifdef RENAMEDEBUG
			cout << "重命名：从 " << i << "接收到的数据：" << endl;
			for (j = 0; j < commSize; j++) {
				rnData[offset + j].printA();
			}
#endif
		}

		for (i = 0; i < numOfChild; i++) {
			rnPQ.push(rnData[readCur[i]]);
		}
		ReNameData temp, minVal;
		while (!rnPQ.empty()) {
			temp = rnPQ.top();
			i = temp.nodeId;
			offset = i * commSize;
			if (temp.index == EMPTY) {
	
				//发送数据
				MPI_Send(((char *)(sendBuf + offset)), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				rnPQ.pop();
#ifdef RENAMEDEBUG
				std::cout << "重命名:最后向 " << i << "发送的数据：" << endl;
				for (j = 0; j < commSize; j++) {
					cout << sendBuf[i * commSize + j] << ",";
				}
				cout << endl;
#endif
			}
			else {

				if (minVal != temp) {
					rank++;
					minVal.setMember(temp);
				}
				else {
					isCycle = 1;
				}

				sendBuf[readCur[i]++] = rank;
				rnPQ.pop();

				if (readCur[i] == (offset + commSize)) {
					//发送数据
					MPI_Send((char *)(sendBuf + offset), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);

#ifdef RENAMEDEBUG
					std::cout << "重命名:向 " << i << "发送的数据：" << endl;
					for (j = 0; j < commSize; j++) {
						cout << sendBuf[offset + j] << ",";
					}
					cout << endl;
#endif

					MPI_Recv((char *)(recvBuf + offset * sizeof(int64)), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
					readCur[i] = offset;
					curBuf = (int64 *)(recvBuf + offset * sizeof(int64));
					for (j = 0; j < commSize; j++) {
						rnData[offset + j].setMember(curBuf[j], i);
						if (curBuf[j] == EMPTY) break;
					}
#ifdef RENAMEDEBUG
					std::cout << "重命名：从 " << i << "接收到的数据：" << endl;
					for (j = 0; j < commSize; j++) {
						rnData[offset + j].printA();
					}
#endif
				}
				rnPQ.push(rnData[readCur[i]]);

			}
		}
	}

	void computeGSA(int32 level) {

		while (!rnPQ.empty())rnPQ.pop();
		int64 n = 0;
		//接受每个子节点的大小
		int32 i, j = 0;
		for (i = 0; i < numOfChild; i++) {
			MPI_Recv(&(sizeOfChild[i]), 1, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD, &status);
			n += sizeOfChild[i];
		}
#ifdef	SADEBUG
		for (i = 0; i < numOfChild; i++) {
			std::cout << "the size of " << i << "is " << sizeOfChild[i] << endl;
		}

		std::cout << "n is " << n << endl;
#endif

		//开始计算SA
		int64 gPos = numOfChild;

		computeGSAl(level, gPos);

		gPos = n + 1;

		computeGSAs(level, gPos);

		bool isCycle = 0;
		int64 rank = 0;
		//rename
		reName(level, rank, isCycle);

		if (isCycle) {
			for (i = 0; i < numOfChild;i++)MPI_Send(&rank, 1, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD);
			computeGSA(level + 1);
		}
		else {
			rank = 0;
			for (i = 0; i < numOfChild; i++)MPI_Send(&rank, 1, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD);
		}

		gPos = numOfChild;

		computeGSAl(level, gPos);

		gPos = n + 1;

		computeGSAs(level, gPos);
	}

	void computeGSAl(int32 level, int64 gPos) {
		int64 recvLen= 2 * commSize * sizeof(int64) + commSize / 8 + ((level == 0) ? commSize : (commSize * sizeof(int64)));
		int64 i, j,offset = 0;

		while (!PQL.empty())PQL.pop();
		for (i = 0; i < numOfChild; i++) {
			readCur[i] = i * commSize;
			MPI_Recv((char*)recvBuf + i * recvLen, recvLen, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
		}

		for (i = 0; i < numOfChild; i++) {
			for (j = 0; j < commSize; j++) {
				
					if (level == 0)
						recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
							((int64 *)(recvBuf + i * recvLen))[commSize + j],
							tget(i * recvLen + 2 * commSize * sizeof(int64), j),
							(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8)[j],
							i);
					else {
						recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
							((int64 *)(recvBuf + i * recvLen))[commSize + j],
							tget(i * recvLen + 2 * commSize * sizeof(int64), j),
							((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 ))[j],
							i);
					}
					map[i][recvData[i * commSize + j].SA] = i * commSize + j;
					if (((int64 *)(recvBuf + i * recvLen))[j] == EMPTY) break;
			}

#ifdef SALDEBUG
			if (level > LEVELL) {
				std::cout << "level is " << level << " computeSAL getData from " << i << "：" << endl;
				for (j = 0; j < commSize; j++) {
					recvData[i * commSize + j].printA();
				}
			}
			
#endif
		}
		
		for (i = 0; i < numOfChild; i++) {
			PQL.push(recvData[readCur[i]]);
		}
		//
#ifdef PQQ
		cout << "比较是否相等(recvData[readCur[0]] < recvData[readCur[1] + 1]):" << (recvData[readCur[0] + 1] < recvData[readCur[1]]) << " PQ size is " << PQ.size() << endl;
		cout << "目前PQ中的值是：" << endl;
		vector<RecvData>tmp;
		while (!PQ.empty())
		{
			RecvData a = PQ.top(); PQ.pop();
			//对a进行操作
			a.printA();
			tmp.push_back(a);
		}
		for (i = 0; i<tmp.size(); i++) PQ.push(tmp[i]);
#endif // PQQ

		RecvData temp,minVal;		
		while (!PQL.empty()) {
			temp = PQL.top();
			i = temp.nodeId;
			offset = i * commSize;
			if (temp.SA == EMPTY) {
				//发送数据
				MPI_Send((char *)(sendBuf + offset), commSize * sizeof(int64) , MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				PQL.pop();

				//清除对应map中的数据
				map[i].clear();

#ifdef SALDEBUG
				if (level > LEVELL) {
					cout << "level is " << level << " computeSAL last send data to " << i << " ：" << endl;
					for (j = 0; j < commSize; j++) {
						cout << sendBuf[i * commSize + j] << ",";
						//recvData[i * commSize + j].printA();
					}
					cout << endl;
				}
				
#endif
			}
			else {
				if (minVal != temp) {
					gPos++;
					minVal.setMember(temp);
				}
				j = temp.SA - 1;
				

				sendBuf[readCur[i]++] = gPos;
				if (j >= 0 && map[i].find(j) != map[i].end() && recvData[map[i][j]].type == L_TYPE) {
					recvData[map[i][j]].suffixGIndex = gPos;
				}

				PQL.pop();
				
				if (readCur[i] == (offset + commSize)) {
					//发送数据
					MPI_Send((char *)(sendBuf + offset), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					//清除对应map中的数据
					map[i].clear();
					//time = 0;
#ifdef SALDEBUG
					if (level > LEVELL) {
						cout << "level is " << level << " computeSAL send data to  " << i << " ：" << endl;
						for (j = 0; j < commSize; j++) {
							cout << sendBuf[i * commSize + j] << ",";
						}
						cout << endl;
					}
					
#endif
					//接收数据
					//if (!markFinish[i]) {
						MPI_Recv((char*)recvBuf + i * recvLen, recvLen, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
						//转化数据
						for (j = 0; j < commSize; j++) {
							if (level == 0)
								recvData[offset + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64), j),
								//	isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8)[j],
									i);
							else {
								recvData[offset + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64), j),
									((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8))[j],
									i);
							}
							map[i][recvData[offset + j].SA] = offset + j;
							if (((int64 *)(recvBuf + i * recvLen))[j] == EMPTY) break;
						}
						readCur[i] = offset;
#ifdef SALDEBUG
						if (level > LEVELL) {
							cout << "level is " << level << " computeSAL get data from " << i << " ：" << endl;
							for (i = 0; i < numOfChild; i++) {
								for (j = 0; j < commSize; j++) {
									recvData[i * commSize + j].printA();
								}
							}
						}
						
#endif
				//	}
					
				}
				PQL.push(recvData[readCur[i]]);
			}
			
#ifdef PQQ
			cout << "目前PQ中的值是：" << endl;
			vector<RecvData>tmp;
			while (!PQ.empty())
			{
				RecvData a = PQ.top(); PQ.pop();
				//对a进行操作
				a.printA();
				tmp.push_back(a);
			}
			for (i = 0; i<tmp.size(); i++) PQ.push(tmp[i]);
#endif // PQQ
		}
		
	}
	
	void computeGSAs(int32 level, int64 &gPos) {
		int64 recvLen = 2 * commSize * sizeof(int64) + commSize / 8 + ((level == 0) ? commSize : (commSize * sizeof(int64)));
		int64 i, j, offset = 0;

#ifdef SASDEBUG
		cout << "gPos is " << gPos << endl;
#endif
		while (!PQS.empty())PQS.pop();
		for (i = 0; i < numOfChild; i++) {
			readCur[i] = i * commSize;
			MPI_Recv((char*)recvBuf + i * recvLen, recvLen, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
		}
		for (i = 0; i < numOfChild; i++) {
			for (j = 0; j < commSize; j++) {
				if (level == 0)
					recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
						((int64 *)(recvBuf + i * recvLen))[commSize + j],
						tget(i * recvLen + 2 * commSize * sizeof(int64), j),
						(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8)[j],
						i);
				else {
					recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
						((int64 *)(recvBuf + i * recvLen))[commSize + j],
						tget(i * recvLen + 2 * commSize * sizeof(int64), j),
						((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8))[j],
						i);
				}
				map[i][recvData[i * commSize + j].SA] = i * commSize + j;
				if (((int64 *)(recvBuf + i * recvLen))[j] == MAX) break;
			}

#ifdef SASDEBUG
			if (level > LEVELL) {
				cout << "level is " << level << " computeSAS get data from " << i << "：" << endl;
				for (j = 0; j < commSize; j++) {
					recvData[i * commSize + j].printA();
				}
			}
			
#endif
		}

		for (i = 0; i < numOfChild; i++) {
			PQS.push(recvData[readCur[i]]);
		}

		//int32 time = 0, minNodeId = 0;
		RecvData temp, minVal;

		while (!PQS.empty()) {
			temp = PQS.top();
			i = temp.nodeId;
			offset = i * commSize;
			if (temp.SA == MAX) {
				//发送数据
				MPI_Send((char*)(sendBuf + offset), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				PQS.pop();

				//清除对应map中的数据
				map[i].clear();
#ifdef SASDEBUG
				if (level > LEVELL) {
					cout << "level is " << level << " computeSAS last send data to " << i << " ：" << endl;
					for (j = 0; j < commSize; j++) {
						cout << sendBuf[i * commSize + j] << ",";
					}
					cout << endl;
			}
#endif

			}
			else {

				if (minVal != temp) {
					gPos--;
					minVal.setMember(temp);
				}
				j = temp.SA - 1;

				sendBuf[readCur[i]++] = gPos;

				if (temp.SA > 0 && map[i].find(j) != map[i].end() && recvData[map[i][j]].type == S_TYPE) {
					recvData[map[i][j]].suffixGIndex = gPos;
				}
				PQS.pop();
				
				if (readCur[i] == (offset + commSize)) {
					//发送数据
					MPI_Send((char*)(sendBuf + offset), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					//清除对应map中的数据
					map[i].clear();

#ifdef SASDEBUG
					if (level > LEVELL) {
						cout << "level is " << level << " computeSAS send data to " << i << " ：" << endl;
						for (j = 0; j < commSize; j++) {
							cout << sendBuf[i * commSize + j] << ",";
						}
						cout << endl;
					}
					
#endif
					//接收数据
				//	if (!markFinish[i]) {
						MPI_Recv((char*)recvBuf + i * recvLen, recvLen, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
						//转化数据
						for (j = 0; j < commSize; j++) {
							if (level == 0)
								recvData[offset + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64), j),
									//isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									(recvBuf + i * recvLen + 2 * commSize * sizeof(int64)  + commSize / 8)[j],
									i);
							else {
								recvData[offset + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64), j),
									//isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8))[j],
									i);
							}
							map[i][recvData[offset + j].SA] = offset + j;
							if (((int64 *)(recvBuf + i * recvLen))[j] == MAX) break;
						}
						readCur[i] = offset;
						//sendCur[i] = offset;
						//time = 0;
						//PQ.push(recvData[readCur[i]++]);
#ifdef SASDEBUG
						if (level > LEVELL) {
							cout << "level is " << level << " computeSAS get data from " << i << "：" << endl;
							for (j = 0; j < commSize; j++) {
								recvData[i * commSize + j].printA();
							}
						}
						
#endif
		//			}

				}
				PQS.push(recvData[readCur[i]]);
			}

		}

	}

	~MainNode()
	{
		delete[] sizeOfChild;
		delete[] recvBuf;
		delete[] sendBuf;
		delete[] map;
		delete[] recvData;
		delete[] readCur;
	}
};
#endif