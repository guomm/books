#ifndef  __MAINNODE
#define  __MAINNODE
#include <iostream>
#include <unordered_map>
#include "mycommon.h"
//#include "mystack.h"

#define tget(j,i) ( (recvBuf[j+i/8]&mask[(i)%8]) ? 1 : 0 )
//#define isGget(i) ( (recvBuf[(i)/8]&mask[(i)%8]) ? 1 : 0 )

//#define SADEBUG
//#define SALDEBUG
//#define SASDEBUG
//#define RENAMEDEBUG

using namespace std;

class MainNode {
private:
	int32 numOfChild;
	MPI_Status status;
	int64 *sizeOfChild;
	
	bool *markFinish;
	unsigned char *recvBuf;
	int64 *sendBuf;
	unordered_map<int64, int64> *map;
	RecvData *recvData;
	priority_queue<RecvData> PQ;
	priority_queue<ReNameData> rnPQ;
	ReNameData *rnData;
	int32 *readCur;
	int32 *sendCur;
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
		sendCur = new int32[sumLen];
		markFinish = new bool[numOfChild];
		sizeOfChild = new int64[numOfChild];
		for (int32 i = 0; i < numOfChild; i++) {
			markFinish[i] = false;
			sizeOfChild[i] = 0;
			unordered_map<int64, int64> mapTemp;
			map[i] = mapTemp;
		}
		for (int32 i = 0; i < sumLen; i++) {
			sendCur[i] = 0;
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

	void reName(int32 level) {
#ifdef RENAMEDEBUG
		cout << "主节点开始重命名: numOfChild is " << numOfChild << endl;
#endif
		int64 rank = 0, offset = 0;
		rank = 0;

		int64 *curBuf;
		int64 i, j = 0;
		//MPI_Status status;
		for (i = 0; i < numOfChild; i++) {
			//	sumChildSize[i] = sumSize;
			//	sumSize += sizeOfChild[i];
			
			offset = i * commSize;
			readCur[i] = offset;
			sendCur[i] = offset;
			MPI_Recv((char *)recvBuf + offset * sizeof(int64), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
			curBuf = (int64 *)(recvBuf + offset);
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
			rnPQ.push(rnData[readCur[i]++]);
			markFinish[i] = 0;
		}
		ReNameData temp, minVal;
		int32 time = 0;
		while (!rnPQ.empty()) {
			temp = rnPQ.top();
			i = temp.nodeId;
			if (temp.index == EMPTY) {
				markFinish[i] = 1;
				//发送数据
				MPI_Send(((char *)(sendBuf + i * commSize)), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				rnPQ.pop();
#ifdef RENAMEDEBUG
				std::cout << "重命名:最后向 " << i << "发送的数据：" << endl;
				for (j = 0; j < commSize; j++) {
					cout << sendBuf[i * commSize + j] << endl;
				}
#endif
			}
			else {

				if (minVal != temp) {
					rank++;
					minVal.setMember(temp);
				}

				sendBuf[i * commSize + sendCur[i]] = rank;
				sendCur[i]++;
				time++;
				rnPQ.pop();
				//rnPQ.push(rnData[i * commSize + readCur[i]]);
				//readCur[i]++;

				if (time == commSize) {
					//发送数据
					offset = i * commSize;
					MPI_Send((char *)(sendBuf + offset), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);

#ifdef RENAMEDEBUG
					std::cout << "重命名:向 " << i << "发送的数据：" << endl;
					for (j = 0; j < commSize; j++) {
						cout << sendBuf[offset + j] << endl;
					}
#endif
					MPI_Recv((char *)(recvBuf + offset * sizeof(int64)), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
					readCur[i] = i * commSize;
					sendCur[i] = i * commSize;

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
				rnPQ.push(rnData[readCur[i]++]);

			}
		}
	}

	void computeGSA(int32 level) {

		//cout << "....mainNode..." << endl;
		while (!rnPQ.empty())rnPQ.pop();
		int64 n = 0;
		//接受每个子节点的大小
		int32 i, j = 0;
		//int64 sumSize = 0;
		for (i = 0; i < numOfChild; i++) {
			MPI_Recv(&(sizeOfChild[i]), 1, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD, &status);
			//cout << "the size of " << i << "is " << sizeOfChild[i] << endl;
			n += sizeOfChild[i];
		}
#ifdef	SADEBUG
		for (i = 0; i < numOfChild; i++) {
			//MPI_Recv(&(sizeOfChild[i]), 1, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD, &status);
			std::cout << "the size of " << i << "is " << sizeOfChild[i] << endl;
		}

		std::cout << "n is " << n << endl;
#endif

		//开始计算SA
		int64 gPos = numOfChild;

		computeGSAl(level, gPos);

		gPos = n + 1;

		computeGSAs(level, gPos);

		//rename
		reName(level);
		//MPI_Recv((char *)recvBuf, commSize * sizeof(int64), MPI_CHAR, 1, 1, MPI_COMM_WORLD, &status);
		//cout << "hhhhhhgggg" << endl;
		//判断是否递归
		bool isCycle = 0;
		if (gPos != numOfChild + 1) {
			//需要递归,接收每个节点传送过来的lms字符的全局排名
			isCycle = 1;
			for (i = 0; i < numOfChild;i++)MPI_Send(&isCycle, 1, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
			computeGSA(level + 1);
		}
		else {
			for (i = 0; i < numOfChild; i++)MPI_Send(&isCycle, 1, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
		}

		gPos = numOfChild;

		computeGSAl(level, gPos);

		gPos = n + 1;

		computeGSAs(level, gPos);
	}

	void computeGSAl(int32 level, int64 gPos) {
		int64 recvLen= 2 * commSize * sizeof(int64) + commSize / 8 + ((level == 0) ? commSize : (commSize * sizeof(int64)));
		int64 i, j, numFinised = 0;
		//int32 chrOffset = 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8;
		//int64 gPos = numOfChild + 1;
		int64 sumSize = 0;
		//int64 *sumChildSize = new int64[numOfChild];
		while (!PQ.empty())PQ.pop();
		for (i = 0; i < numOfChild; i++) {
		//	sumChildSize[i] = sumSize;
		//	sumSize += sizeOfChild[i];
			markFinish[i] = 0;
			readCur[i] = i * commSize;
			sendCur[i] = i * commSize;
			MPI_Recv((char*)recvBuf + i * recvLen, recvLen, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
		}
		//cout << "接收到初始化数据长度:"<< recvLen<<endl;
		for (i = 0; i < numOfChild; i++) {
			for (j = 0; j < commSize; j++) {
				
					if (level == 0)
						recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
							((int64 *)(recvBuf + i * recvLen))[commSize + j],
							tget(i * recvLen + 2 * commSize * sizeof(int64), j),
							//isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
							(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8)[j],
							i);
					else {
						recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
							((int64 *)(recvBuf + i * recvLen))[commSize + j],
							tget(i * recvLen + 2 * commSize * sizeof(int64), j),
						//	isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
							((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 ))[j],
							i);
					}
					map[i][recvData[j].SA] = i * commSize + j;
					if (((int64 *)(recvBuf + i * recvLen))[j] == EMPTY) break;
			}

#ifdef SALDEBUG
			if (level > 0)cout <<"level is "<<level<<" computeSAL 从 "<<i<< "接收到的数据：" << endl;
			for (i = 0; i < numOfChild; i++) {
				for (j = 0; j < commSize; j++) {
					if(level > 0)recvData[i * commSize + j].printA();
				}
			}
#endif
		}
		
		for (i = 0; i < numOfChild; i++) {
			PQ.push(recvData[readCur[i]++]);
		}

		int32 time = 0;
		RecvData temp,minVal;
		
		while (!PQ.empty()) {
			//cout << "111" << endl;
			temp = PQ.top();
			i = temp.nodeId;
			if (temp.SA == EMPTY) {
				//cout << "向 " << i << " 发送的数据：" << endl;
				markFinish[i] = 1;
				//发送数据
				MPI_Send((char *)(sendBuf + i * commSize), commSize * sizeof(int64) , MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				PQ.pop();

				//清除对应map中的数据
				map[i].clear();

#ifdef SALDEBUG
				if (level > 0)cout << "level is " << level << " computeSAL 最后向 " << i << " 发送的数据：" << endl;
				for (j = 0; j < commSize; j++) {
					if (level > 0)cout<<sendBuf[i * commSize + j]<<",";
				}
				cout << endl;
#endif
			}
			else {
				//cout << "222" <<endl;
				if (minVal != temp) {
					gPos++;
					minVal.setMember(temp);
				}
				j = temp.SA - 1;
				

				sendBuf[sendCur[i]++] = gPos;
				//sendCur[i]++;
				time++;
				
				if (temp.SA > 0 && map[i].find(j) != map[i].end()) {
					recvData[map[i][j]].suffixGIndex = gPos;
					
				}
				
				//cout <<"gPos is" << gPos<< " minVal is:" << endl;
				//minVal.printA();
				//cout << "gPos is" << gPos << endl;
				PQ.pop();
				
				//readCur[i]++;
				if (time == commSize) {
					//发送数据
					MPI_Send((char *)(sendBuf + i * commSize), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					//清除对应map中的数据
					map[i].clear();
					time = 0;
#ifdef SALDEBUG
					if (level > 0)cout << "level is " << level << " computeSAL 向 "<<i << " 发送的数据：" << endl;
					for (j = 0; j < commSize; j++) {
						if (level > 0)cout << sendBuf[i * commSize + j] << ",";
					}
					cout << endl;
#endif
					//接收数据
					if (!markFinish[i]) {
						MPI_Recv((char*)recvBuf + i * recvLen, recvLen, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
						//转化数据
						for (j = 0; j < commSize; j++) {
							if (level == 0)
								recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64), j),
								//	isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8)[j],
									i);
							else {
								recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64), j),
								//	isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8))[j],
									i);
							}
							map[i][recvData[j].SA] = i * commSize + j;
						}
						readCur[i] = i * commSize;
						sendCur[i] = i * commSize;
						
						PQ.push(recvData[readCur[i]++]);
#ifdef SALDEBUG
						if (level > 0)cout << "level is " << level << " computeSAL 从 " << i << "接收到的数据：" << endl;
						for (i = 0; i < numOfChild; i++) {
							for (j = 0; j < commSize; j++) {
								if (level > 0)recvData[i * commSize + j].printA();
							}
						}
#endif
					}
					
				}
				else {
					PQ.push(recvData[readCur[i]++]);
				}
			}
			
		}
		
	}
	
	void computeGSAs(int32 level, int64 &gPos) {
		int64 recvLen = 2 * commSize * sizeof(int64) + commSize / 8 + ((level == 0) ? commSize : (commSize * sizeof(int64)));
		int64 i, j, sumSize = 0;
		int32 chrOffset = 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8;
		//int64 gPos = numOfChild + 1;
		//int32 sumSize = 0;
		//int64 *sumChildSize = new int64[numOfChild];
#ifdef SASDEBUG
		cout << "gPos is " << gPos << endl;
#endif
		while (!PQ.empty())PQ.pop();
		for (i = 0; i < numOfChild; i++) {
			//	sumChildSize[i] = sumSize;
			//	sumSize += sizeOfChild[i];
			markFinish[i] = 0;
			readCur[i] = i * commSize;
			sendCur[i] = i * commSize;
			MPI_Recv((char*)recvBuf + i * recvLen, recvLen, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
		}
		for (i = 0; i < numOfChild; i++) {
			for (j = 0; j < commSize; j++) {
				if (((int64 *)(recvBuf + i * recvLen))[j] != EMPTY) {
					if (level == 0)
						recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
							((int64 *)(recvBuf + i * recvLen))[commSize + j],
							tget(i * recvLen + 2 * commSize * sizeof(int64), j),
							//isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
							(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8)[j],
							i);
					else {
						recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
							((int64 *)(recvBuf + i * recvLen))[commSize + j],
							tget(i * recvLen + 2 * commSize * sizeof(int64), j),
							//isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
							((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8))[j],
							i);
					}
					map[i][recvData[j].SA] = i * commSize + j;
				}
			}
		}
#ifdef SASDEBUG
		if (level > 0)cout << "level is " << level << " computeSAS 从 " << i << "接收到的数据：" << endl;
		for (i = 0; i < numOfChild; i++) {
			for (j = 0; j < commSize; j++) {
				if (level > 0)recvData[i * commSize + j].printA();
			}
		}
#endif
		for (i = 0; i < numOfChild; i++) {
			PQ.push(recvData[readCur[i]++]);
		}

		int32 time = 0, minNodeId = 0;
		RecvData temp, minVal;

		while (!PQ.empty()) {
			temp = PQ.top();
			i = temp.nodeId;
			if (temp.SA == EMPTY) {
				markFinish[i] == 1;
				//发送数据
				MPI_Send((char*)(sendBuf + i * commSize), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				PQ.pop();

				//清除对应map中的数据
				map[i].clear();
#ifdef SASDEBUG
				if (level > 0)cout << "level is " << level << " computeSAS 最后向 " << i << " 发送的数据：" << endl;
				for (j = 0; j < commSize; j++) {
					if (level > 0)cout << sendBuf[i * commSize + j] << ",";
				}
				cout << endl;
#endif

			}
			else {

				if (minVal != temp) {
					gPos--;
					minVal.setMember(temp);
				}
				j = temp.SA - 1;

				sendBuf[sendCur[i]++] = gPos;
				//sendCur[i]++;
				time++;
				if (temp.SA > 0 && map[i].find(j) != map[i].end()) {
					recvData[map[i][j]].suffixGIndex = gPos;
				}
				PQ.pop();
				
			//	readCur[i]++;
				if (time == commSize) {
					//发送数据
					MPI_Send((char*)(sendBuf + i * commSize), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					//清除对应map中的数据
					map[i].clear();

#ifdef SASDEBUG
					if (level > 0)cout << "level is " << level << " computeSAS 向 " << i << " 发送的数据：" << endl;
					for (j = 0; j < commSize; j++) {
						if (level > 0)cout << sendBuf[i * commSize + j] << ",";
					}
					cout << endl;
#endif
					//接收数据
					if (!markFinish[i]) {
						MPI_Recv((char*)recvBuf + i * recvLen, recvLen, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
						//转化数据
						for (j = 0; j < commSize; j++) {
							if (level == 0)
								recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64), j),
									//isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									(recvBuf + i * recvLen + 2 * commSize * sizeof(int64)  + commSize / 8)[j],
									i);
							else {
								recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64), j),
									//isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8))[j],
									i);
							}
							map[i][recvData[j].SA] = i * commSize + j;
						}
						readCur[i] = i * commSize;
						sendCur[i] = i * commSize;
						time = 0;
						PQ.push(recvData[readCur[i]++]);
#ifdef SASDEBUG
						if (level > 0)cout << "level is " << level << " computeSAS 从 " << i << "接收到的数据：" << endl;
						for (i = 0; i < numOfChild; i++) {
							for (j = 0; j < commSize; j++) {
								if (level > 0)recvData[i * commSize + j].printA();
							}
						}
#endif
					}

				}
				else {
					PQ.push(recvData[readCur[i]++]);
				}
			}

		}

	}

	~MainNode()
	{
		delete[] sizeOfChild;
		delete[] markFinish;
		delete[] recvBuf;
		delete[] sendBuf;
		delete[] map;
		delete[] recvData;
		delete[] readCur;
		delete[] sendCur;
	}
};
#endif