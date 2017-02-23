#ifndef  __MAINNODE
#define  __MAINNODE
#include <iostream>
#include <hash_map>
#include "mycommon.h"
#include "mystack.h"

#define tget(i) ( (recvBuf[(i)/8]&mask[(i)%8]) ? 1 : 0 )
#define isGget(i) ( (recvBuf[(i)/8]&mask[(i)%8]) ? 1 : 0 )

class MainNode {
private:
	int32 numOfChild;
	MPI_Status status;
	int64 *sizeOfChild;
	int64 n;
	bool *markFinish;
	unsigned char * recvBuf;
	int64 *sendBuf;
	hash_map<int64, RecvData> *map;
	RecvData *recvData;
	priority_queue<RecvData> PQ;
	int32 *readCur;
	int32 *sendCur;
public:
	MainNode(int32 num_child_node) :numOfChild(num_child_node){
		int64 recvLen = 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8 + commSize * sizeof(int64);
		int32 sumLen = commSize * numOfChild;
		recvBuf = new unsigned char[recvLen * numOfChild];
		sendBuf = new int64[sumLen];
		map = new hash_map<int64, RecvData>[numOfChild];
		recvData = new RecvData[sumLen];
		readCur = new int32[sumLen];
		sendCur = new int32[sumLen];
		for (int32 i = 0; i < numOfChild; i++) {
			markFinish[i] = false;
		}
		for (int32 i = 0; i < sumLen; i++) {
			sendCur[i] = 0;
			readCur[i] = 0;
			sendBuf[i] = 0;
			RecvData temp;
			recvData[i] = temp;
		}
	}

	void computeGSA(unsigned char *recvBuf, int64 *sendBuf, hash_map<int64, RecvData> *map, RecvData *recvData, priority_queue<RecvData> &PQ, priority_queue<ReNameData> &rnPQ,ReNameData *rnData, int32 * readCur, int32 * sendCur, bool *markFinish,  int32 level) {

		
		while (!rnPQ.empty())rnPQ.pop();

		//接受每个子节点的大小
		int32 i, j = 0;
		//int64 sumSize = 0;
		for (i = 0; i < numOfChild; i++) {
			MPI_Recv(&(sizeOfChild[i]), 1, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD, &status);
			cout << "the size of " << i << "is " << sizeOfChild[i] << endl;
			n += sizeOfChild[i];
		}
		//开始计算SA
		int64 gPos = numOfChild + 1;

		computeGSAl(recvBuf, sendBuf, map, recvData, PQ, readCur, sendCur, markFinish, level, gPos);

		gPos = n;

		computeGSAs(recvBuf, sendBuf, map, recvData, PQ, readCur, sendCur, markFinish, level, gPos);

		//rename
		int64 rank = 0, offset = 0;
		for (i = 0; i < numOfChild; i++) {
			//	sumChildSize[i] = sumSize;
			//	sumSize += sizeOfChild[i];
			offset = i * commSize;
			readCur[i] = offset;
			sendCur[i] = offset;
			MPI_Recv((int64 *)recvBuf + offset, commSize, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD, &status);
			int64 *curBuf = (int64 *)recvBuf + offset;
			for (j = 0; j < commSize; j++) {
				ReNameData temp1(curBuf[j], i);
				rnData[offset + j] = temp1;
				if (curBuf[j] == EMPTY) break;
			}
		}

		for (i = 0; i < numOfChild; i++) {
			rnPQ.push(rnData[readCur[i]++]);
			markFinish[i] = 0;
		}
		ReNameData temp, minVal;
		while (!rnPQ.empty()) {
			temp = rnPQ.top();
			i = temp.nodeId;
			if (temp.index == EMPTY) {
				markFinish[i] = 1;
				//发送数据
				MPI_Send(sendBuf + i * commSize, commSize, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD);
				rnPQ.pop();
			}
			else {
				sendBuf[i * commSize + sendCur[i]] = rank++;
				sendCur[i]++;
				rnPQ.pop();
				rnPQ.push(rnData[i * commSize + readCur[i]]);
				readCur[i]++;

				if (readCur[i] == commSize) {
					//发送数据
					offset = i * commSize;
					MPI_Send(sendBuf + offset, commSize, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD);
					MPI_Recv((int64 *)recvBuf + offset, commSize, MPI_LONG, i + 1, i + 1, MPI_COMM_WORLD, &status);
					readCur[i] = 0;
					sendCur[i] = 0;

					int64 *curBuf = (int64 *)recvBuf + offset;
					for (j = 0; j < commSize; j++) {
						ReNameData temp1(curBuf[j], i);
						rnData[offset + j] = temp1;
						if (curBuf[j] == EMPTY) break;
					}
				}
			}
		}
		//判断是否递归
		bool isCycle = 0;
		if (gPos != numOfChild + 1) {
			//需要递归,接收每个节点传送过来的lms字符的全局排名
			isCycle = 1;
			for (i = 0; i < numOfChild;i++)MPI_Send(&isCycle, 1, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
			computeGSA(recvBuf,sendBuf,map,recvData,PQ,rnPQ,rnData,readCur,sendCur,markFinish,level + 1);
		}
		else {
			for (i = 0; i < numOfChild; i++)MPI_Send(&isCycle, 1, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
		}
	}

	void computeGSAl(unsigned char *recvBuf, int64 *sendBuf, hash_map<int64, RecvData> *map, RecvData *recvData, priority_queue<RecvData> &PQ, int32 * readCur, int32 * sendCur, bool *markFinish, int32 level, int64 &gPos) {
		int64 recvLen= 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8 + (level == 0) ? commSize : (commSize * sizeof(int64));
		int32 i, j, numFinised = 0;
		int32 chrOffset = 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8;
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
		for (i = 0; i < numOfChild; i++) {
			for (j = 0; j < commSize; j++) {
				
					if (level == 0)
						recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
							((int64 *)(recvBuf + i * recvLen))[commSize + j],
							tget(i * recvLen + 2 * commSize * sizeof(int64) + j),
							isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
							(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8)[j],
							i);
					else {
						recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
							((int64 *)(recvBuf + i * recvLen))[commSize + j],
							tget(i * recvLen + 2 * commSize * sizeof(int64) + j),
							isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
							((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8))[j],
							i);
					}
					map[i][recvData[j].SA] = recvData[i * commSize + j];
					if (((int64 *)(recvBuf + i * recvLen))[j] == EMPTY) break;
			}
		}
		
		for (i = 0; i < numOfChild; i++) {
			PQ.push(recvData[readCur[i]++]);
		}

		int32 minNodeId = 0;
		RecvData temp,minVal;
		
		while (!PQ.empty) {
			temp = PQ.top();
			i = temp.nodeId;
			if (temp.SA == EMPTY) {
				markFinish[i] = 1;
				//发送数据
				MPI_Send((char*)(sendBuf + i * commSize), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				PQ.pop();

				//清除对应map中的数据
				map[i].clear();
			}
			else {
				
				j = temp.SA - 1;
				sendBuf[i * commSize + sendCur[i]] = gPos;
				sendCur[i]++;
				if (minVal != temp)gPos++;
				minVal = temp;
				if (temp.SA > 0 && map[i].find(j) != map[i].end()) {
					map[i][j].suffixGIndex = gPos;
					map[i][j].isGlobal = 1;
				}
				PQ.pop();
				PQ.push(recvData[readCur[i]]);
				readCur[i]++;
				if (readCur[i] == commSize) {
					//发送数据
					MPI_Send((char*)(sendBuf + i * commSize), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					//清除对应map中的数据
					map[i].clear();
					//接收数据
					if (!markFinish[i]) {
						MPI_Recv((char*)recvBuf + i * recvLen, recvLen, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
						//转化数据
						for (j = 0; j < commSize; j++) {
							if (level == 0)
								recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64) + j),
									isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8)[j],
									i);
							else {
								recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64) + j),
									isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8))[j],
									i);
							}
							map[i][recvData[j].SA] = recvData[i * commSize + j];
						}
						readCur[i] = 0;
						sendCur[i] = 0;
					}
					
				}
			}
			
		}
		
	}
	
	void computeGSAs(unsigned char *recvBuf, int64 *sendBuf, hash_map<int64, RecvData> *map, RecvData *recvData, priority_queue<RecvData> &PQ, int32 * readCur, int32 * sendCur, bool *markFinish, int32 level, int64 &gPos) {
		int64 recvLen = 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8 + (level == 0) ? commSize : (commSize * sizeof(int64));
		int32 i, j, numFinised = 0;
		int32 chrOffset = 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8;
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
		for (i = 0; i < numOfChild; i++) {
			for (j = 0; j < commSize; j++) {
				if (((int64 *)(recvBuf + i * recvLen))[j] != EMPTY) {
					if (level == 0)
						recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
							((int64 *)(recvBuf + i * recvLen))[commSize + j],
							tget(i * recvLen + 2 * commSize * sizeof(int64) + j),
							isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
							(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8)[j],
							i);
					else {
						recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
							((int64 *)(recvBuf + i * recvLen))[commSize + j],
							tget(i * recvLen + 2 * commSize * sizeof(int64) + j),
							isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
							((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8))[j],
							i);
					}
					map[i][recvData[j].SA] = recvData[i * commSize + j];
				}
			}
		}

		for (i = 0; i < numOfChild; i++) {
			PQ.push(recvData[readCur[i]++]);
		}

		int32 minNodeId = 0;
		RecvData temp, minVal;

		while (!PQ.empty) {
			temp = PQ.top();
			i = temp.nodeId;
			if (temp.SA == EMPTY) {
				markFinish[i] == 1;
				//发送数据
				MPI_Send((char*)(sendBuf + i * commSize), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
				PQ.pop();

				//清除对应map中的数据
				map[i].clear();
			}
			else {

				j = temp.SA - 1;
				sendBuf[i * commSize + sendCur[i]] = gPos;
				sendCur[i]++;
				if (minVal != temp)gPos--;
				minVal = temp;
				if (temp.SA > 0 && map[i].find(j) != map[i].end()) {
					map[i][j].suffixGIndex = gPos;
					map[i][j].isGlobal = 1;
				}
				PQ.pop();
				PQ.push(recvData[readCur[i]]);
				readCur[i]++;
				if (readCur[i] == commSize) {
					//发送数据
					MPI_Send((char*)(sendBuf + i * commSize), commSize * sizeof(int64), MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD);
					//清除对应map中的数据
					map[i].clear();
					//接收数据
					if (!markFinish[i]) {
						MPI_Recv((char*)recvBuf + i * recvLen, recvLen, MPI_CHAR, i + 1, i + 1, MPI_COMM_WORLD, &status);
						//转化数据
						for (j = 0; j < commSize; j++) {
							if (level == 0)
								recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64) + j),
									isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8)[j],
									i);
							else {
								recvData[i * commSize + j].setMember(((int64 *)(recvBuf + i * recvLen))[j],
									((int64 *)(recvBuf + i * recvLen))[commSize + j],
									tget(i * recvLen + 2 * commSize * sizeof(int64) + j),
									isGget(i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + j),
									((int64 *)(recvBuf + i * recvLen + 2 * commSize * sizeof(int64) + commSize / 8 + commSize / 8))[j],
									i);
							}
							map[i][recvData[j].SA] = recvData[i * commSize + j];
						}
						readCur[i] = 0;
						sendCur[i] = 0;
					}

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