#ifndef  __CHILDNODE
#define  __CHILDNODE

#include <iostream>
#include <string>
#include "mycommon.h"
#include "sais.h"
#include <vector>
#include <unordered_map>
#include <functional> 
#include <time.h>

using namespace std;
//using SendData<T> = RecvData<T>;

#define tget(i) ( (t[(i)/8]&mask[(i)%8]) ? 1 : 0 )
#define tset(i, b) t[(i)/8]=(b) ? (mask[(i)%8]|t[(i)/8]) : ((~mask[(i)%8])&t[(i)/8])
#define chr(i) (cs==sizeof(uint_type)?((uint_type*)s)[i]:((unsigned char *)s)[i])
#define isLMS(i) (i>0 && tget(i) && !tget(i-1))

#define bufset(i,j, b) bufSend[i+j/8]=(b) ? (mask[(j)%8]|bufSend[i+j/8]) : ((~mask[(j)%8])&bufSend[i+j/8])

//#define SALDEBUG

//#define SADEBUG
//#define SALDEBUGC
//#define SASDEBUG
//#define RENAMEDEBUGC
//#define LL -1

class ChildNode {
private:
	int32 nodeId;
	MPI_Status status;
	uint_type n;
	uint_type *bufRecv = new uint_type[commSize];
	uint_type *bufSend = new uint_type[commSize];
	uint_type *bufSendSA = new uint_type[commSize];
	unordered_map<uint_type, uint_type> map;
	SendData<unsigned char> *bufSend0 = new SendData<unsigned char>[commSize];
	SendData<uint_type> *bufSend1 = new SendData<uint_type>[commSize];
public:
	clock_t commTime, sTime;
	ChildNode(int32 nodeIdC) {
		nodeId = nodeIdC;

		for (int32 i = 0; i < commSize; i++) {
			SendData<unsigned char> tt;
			bufSend0[i] = tt;
			SendData<uint_type> hh;
			bufSend1[i] = hh;
			bufRecv[i] = 0;
			bufSend[i] = 0;
			bufSendSA[i] = 0;
		}
	}


	void getBuckets(unsigned char *s, uint_type *bkt, uint_type n, uint_type K, int32 cs, bool end) {
		uint_type i, sum = 0;
		for (i = 0; i <= K; i++) bkt[i] = 0; // clear all buckets
		for (i = 0; i<n; i++) bkt[chr(i)]++; // compute the size of each bucket
		for (i = 0; i <= K; i++) { sum += bkt[i]; bkt[i] = end ? sum - 1 : sum - bkt[i]; }
	}


	void computeSA(unsigned char *s, uint_type *SA, uint_type *suffixGIndex, uint_type n, int32 cs, uint_type K, uint_type level) {
		//发送大小
		int64 size = n;
		MPI_Send(&size, 1, MPI_LONG, MAINNODEID, nodeId, MPI_COMM_WORLD);
		if (n == 1) {
			//只有一个字符的情况
			induceSAl1();
			induceSAs1();
			reName1();
			int64 isLoop = 0;
			MPI_Recv(&(isLoop), 1, MPI_LONG, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
			if (isLoop) {
				uint_type *SA1 = SA, *s1 = SA + n - 1;
				s1[0] = nodeId;
				computeSA((unsigned char *)s1, SA1, suffixGIndex, 1, sizeof(int64), 1, level + 1);
			}
			induceSAl1();
			induceSAs1();
			SA[0] = 0;
			suffixGIndex[0] = nodeId;
			return;
		}
		//初始化数据
		uint_type *bkt = new uint_type[K + 1]; // bucket counters

		unsigned char *t = new unsigned char[n / 8 + 1];

		uint_type i, j = 0;
		tset(n - 2, 0); tset(n - 1, 1); // the sentinel must be in s1, important!!!
		for (i = n - 3; i >= 0; i--) tset(i, (chr(i) < chr(i + 1) || (chr(i) == chr(i + 1) && tget(i + 1) == 1)) ? 1 : 0);
		//uint_type maxVal = 0;
		//for (i = n - 1; i >= 0; i--)if (chr(i) > maxVal)maxVal = chr(i);
		//cout << "level is " << level << " nodeId is " << nodeId << " K is " <<K<<"cs is "<<cs<<" maxVal is "<<maxVal<< endl;
		getBuckets(s, bkt, n, K, cs, true); // find ends of buckets



		for (i = 0; i < n; i++) {
			SA[i] = EMPTY; suffixGIndex[i] = EMPTY;
		}
		/*cout << "level is " << level << " nodeId is " << nodeId << " 开始......." << endl;
		for (i = 1; i < n - 1; i++) {
		cout << bkt[i] << ",";
		}
		cout << endl;*/
		//cout << "level is " << level << " nodeId is " << nodeId << " 初始化lms..." << endl;
		for (i = 1; i < n - 1; i++)
			if (isLMS(i)) {
				SA[bkt[chr(i)]--] = i;
				suffixGIndex[i + 1] = 0;
			}
		SA[0] = n - 1; // set the single sentinel LMS-substring
		suffixGIndex[n - 1] = nodeId;

		//cout << "level is "<<level <<" nodeId is "<<nodeId<<" 开始induceL..." << endl;
		induceSAl(t, SA, s, bkt, n, K, cs, level, suffixGIndex, nodeId, status);

		//cout << "level is " << level << " nodeId is " << nodeId << " 开始induceS..." << endl;
		induceSAs(t, SA, s, bkt, n, K, cs, suffixGIndex, nodeId, status, level);

#ifdef RENAMEDEBUGC
		cout << "level is " << level << " 子节点 " << nodeId << " 开始重命名:" << std::endl;
#endif

		MarkData markData;
		MPI_Recv((char *)(&markData), 1 * sizeof(MarkData), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		//markData.print();
		uint_type n1 = 0;
		if (markData.isCycle == 0 && markData.isSelfCom == 1) {
			cout << "level is "<<level<< " nodeId is "<<nodeId<<" need sais，not recycle" << endl;
			//全局不需要递归，自己SAIS求解本地SA
			//int64 n1 = 0;
			for (i = 0; i < n; i++)
				if (isLMS(SA[i]))
					SA[n1++] = SA[i];

			// Init the name array buffer
			for (i = n1; i < n; i++) SA[i] = EMPTY;
			// find the lexicographic names of all substrings
			uint_type name = 0, prev = -1;
			for (i = 0; i < n1; i++) {
				int64 pos = SA[i]; bool diff = false;
				for (int64 d = 0; d < n; d++)
					if (prev == -1 || pos + d == n - 1 || prev + d == n - 1 ||
						chr(pos + d) != chr(prev + d) ||
						tget(pos + d) != tget(prev + d))
					{
						diff = true; break;
					}
					else
						if (d > 0 && (isLMS(pos + d) || isLMS(prev + d)))
							break;

				if (diff)
				{
					name++; prev = pos;
				}
				pos = pos / 2; //(pos%2==0)?pos/2:(pos-1)/2;
				SA[n1 + pos] = name - 1;
			}
			for (i = n - 1, j = n - 1; i >= n1; i--)
				if (SA[i] != EMPTY) SA[j--] = SA[i];
			uint_type *s1 = SA + n - n1, *SA1 = SA;
			SA_IS((unsigned char *)s1, SA1, n1, name - 1, sizeof(uint_type), 0);
			for (j = 0, i = 1; i < n; i++)
				if (isLMS(i)) s1[j++] = i; // get p1
			for (i = 0; i < n1; i++) SA1[i] = s1[SA1[i]]; // get index in s1
			for (i = 0, j = 0; i < n1; i++)
				s1[j++] = suffixGIndex[SA1[i]];
			for (i = 0, j = 0; i < n1; i++)
				suffixGIndex[j++] = s1[i];
			uint_type sendPos = 0, offset = 0;
			suffixGIndex[n1] = 0;
			j = n1 + 1;
			uint_type p = 0;
			while (sendPos < j) {
				sTime = clock();
				MPI_Send((char *)(suffixGIndex + sendPos), commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
				MPI_Recv((char *)(suffixGIndex + sendPos), commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
				commTime += clock() - sTime;
				sendPos += commSize;
			}

			
		}
		else if (markData.isCycle == 1) {
			cout << "level is " << level << " nodeId is " << nodeId << " need recycle" << endl;
			//全局需要递归
			n1 = 0;
			for (i = 0; i < n; i++)
				if (isLMS(SA[i])) SA[n1++] = SA[i];
			uint_type *s1 = SA + n - n1;
			for (i = 0, j = 0; i < n1; i++)
				s1[j++] = suffixGIndex[SA[i]];
			for (i = 0, j = 0; i < n1; i++)
				suffixGIndex[j++] = s1[i];
			suffixGIndex[n1] = 0;
			uint_type sendPos = 0, offset = 0;
			j = n1 + 1;
			
			while (sendPos < j) {
				sTime = clock();
				MPI_Send((char *)(suffixGIndex + sendPos), commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
				MPI_Recv((char *)(suffixGIndex + sendPos), commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
				commTime += clock() - sTime;
				sendPos += commSize;
			}

			for (i = n - n1, j = 0; i < n; i++) {
				SA[i] = suffixGIndex[j];
				//SA[suffixGIndex[i]] = suffixGIndex[j++];
				suffixGIndex[i] = SA[j++];
			}
			//j =  n - n1;
			//cout << ",,,," << endl;
			uint_type tempD = 0;
			/*i = 0; j = n1 - 1;
			while (i < j) {
				tempD = SA[i];
				SA[i] = SA[j];
				SA[j] = tempD;
				i++; j--;
			}*/

			i = n-n1; j = n - 1;
			while (i < j) {
				tempD = SA[i];
				SA[i] = SA[j];
				SA[j] = tempD;

				tempD = suffixGIndex[i];
				suffixGIndex[i] = suffixGIndex[j];
				suffixGIndex[j] = tempD;
				i++; j--;
			}
			//cout << "...." << endl;
			s1 = SA + n - n1;
			computeSA((unsigned char *)s1, SA, suffixGIndex, n1, sizeof(uint_type), s1[0], level + 1);
			//cout << "level is " << level << " nodeId is " << nodeId << " 递归结束" << endl;
			offset = n - n1;
			for (i = 0; i < n1; i++) {
				SA[i] = suffixGIndex[offset+SA[i]];
			}	
			i = 0; j = n1 - 1;
			while (i < j) {
			tempD = suffixGIndex[i];
			suffixGIndex[i] = suffixGIndex[j];
			suffixGIndex[j] = tempD;
			i++; j--;
			}
		}
		else {
			cout << "level is " << level << " nodeId is " << nodeId << "don't need sais，also recycle" << endl;
			//本地不需要SAIS求解，全局不需要递归
			/*cout << "level is " << level << " 子节点 " << nodeId << " n1 is :" << n1 << "重命名:S：" << endl;
			for (i = 0; i < n; i++) {
				cout << chr(i) << ",";
			}
			cout << endl;
			cout << "level is " << level << " 子节点 " << nodeId << " n1 is :" << n1 << "重命名:SA：" << endl;
			for (i = 0; i < n; i++) {
				cout << SA[i] << ",";
			}
			cout << endl;*/
			n1 = 0;
			for (i = 0; i < n; i++)
				if (isLMS(SA[i])) SA[n1++] = SA[i];
			uint_type *s1 = SA + n - n1;
			for (i = 0, j = 0; i < n1; i++)
				s1[j++] = suffixGIndex[SA[i]];
			for (i = 0, j = 0; i < n1; i++)
				suffixGIndex[j++] = s1[i];
			suffixGIndex[n1] = 0;
			uint_type sendPos = 0, offset = 0;
			j = n1 + 1;
			uint_type p = 0;
			/*cout << "level is " << level << " 子节点 " << nodeId << " n1 is :" << n1 << "重命名:发送的的数据是：" << endl;
			for (i = 0; i < n1; i++) {
			cout << suffixGIndex[i] << ",";
			}
			cout << endl;*/
			while (sendPos < j) {
				sTime = clock();
				MPI_Send((char *)(suffixGIndex + sendPos), commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
				MPI_Recv((char *)(suffixGIndex + sendPos), commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
				commTime += clock() - sTime;
				sendPos += commSize;
			}
			/*cout << "level is " << level << " 子节点 " << nodeId << " n1 is :" << n1 << "重命名:接收的的数据是：" << endl;
			for (i = 0; i < n1; i++) {
				cout << suffixGIndex[i] << ",";
			}
			cout << endl;*/
		}
		/*cout << "level is " << level << " 子节点 " << nodeId << " n1 is :" << n1 << "重命名:SA：" << endl;
		for (i = 0; i < n1; i++) {
		cout << SA[i] << ",";
		}
		cout << endl;
		cout << "level is " << level << " 子节点 " << nodeId << " n1 is :" << n1 << "重命名:GRank：" << endl;
		for (i = 0; i < n1; i++) {
			cout << suffixGIndex[i] << ",";
		}
		cout << endl;*/
		for (i = n - n1, j = 0; i < n; i++) {
			SA[i] = suffixGIndex[j++];
			//现在suffixgindex前n1个字符是服务器重命名后的值，最后n1个字符是lms字符的sa
		}
		for (i = 0; i < n; i++)suffixGIndex[i] = EMPTY;


		for (i = n - n1, j = 0; i < n; i++) {
			suffixGIndex[SA[j++]] = SA[i];
			//SA[suffixGIndex[i]] = suffixGIndex[j++];
		}

		for (i = n1; i < n; i++) SA[i] = EMPTY;

		getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
		for (i = n1 - 1; i >= 0; i--) {
			j = SA[i]; SA[i] = EMPTY;
			if (level == 0 && i == 0)
				SA[0] = n - 1;
			else
				SA[bkt[chr(j)]--] = j;
		}

		suffixGIndex[n - 1] = nodeId;
		//cout << "start to indcueL " << endl;
		//cout << "level is " << level << " 子节点 " << nodeId << " 开始重命名:" << std::endl;

		//offset = n - n1;
		//for (i = 0; i < n1; i++)SA[offset++] = suffixGIndex[SA[i]];
		//SA[n - n1] = nodeId;

		////利用suffixGIndex发送全局排名，重命名全局排名
		//for (i = n - n1, j = 0; i < n; i++)suffixGIndex[j++] = SA[i];

#ifdef RENAMEDEBUGC
		if (level > LL) {
			cout << "重命名:发送的的数据是：" << endl;
			for (i = 0; i < n1; i++) {
				cout << suffixGIndex[i] << endl;
			}
		}
#endif
		/*cout << "level is " << level << " 子节点 " << nodeId << " n1 is :" << n1 << "重命名:发送的的数据是：" << endl;
		for (i = 0; i < n1; i++) {
		cout << suffixGIndex[i] << ",";
		}
		cout << endl;*/
	
#ifdef RENAMEDEBUGC
		cout << "level is " << level << " nodeId is " << nodeId << " lms字符已放好，开始排L n is " << n << " n1 is " << n1 << endl;
#endif // RENAMEDEBUGC
		//
		//cout << "level is " << level << " nodeId is " << nodeId << " lms字符已放好，开始排L n is " << n << " n1 is " << n1 << endl;
		induceSAl(t, SA, s, bkt, n, K, cs, level, suffixGIndex, nodeId, status);

		/*j = suffixGIndex[n - 1];
		for (i = n - 2; i >= 0; i--)suffixGIndex[i + 1] = suffixGIndex[i];
		suffixGIndex[0] = j;
		suffixGIndex[n - 1] = nodeId;*/
		//cout << "start to indcueS " << endl;
#ifdef RENAMEDEBUGC
		cout << "level is " << level << " nodeId is " << nodeId << " lms字符已放好，开始排S" << endl;
#endif
		//cout << "level is " << level << " nodeId is " << nodeId << " lms字符已放好，开始排S" << endl;
		induceSAs(t, SA, s, bkt, n, K, cs, suffixGIndex, nodeId, status, level);
#ifdef RENAMEDEBUGC
		cout << "level is " << level << " nodeId is " << nodeId << " 排序结束" << endl;
#endif
	/*	j = suffixGIndex[n - 1];
		for (i = n - 2; i >= 0; i--) {
			suffixGIndex[i + 1] = suffixGIndex[i];
		}
		suffixGIndex[0] = j;

		suffixGIndex[n - 1] = nodeId;*/
		//cout << "nodeId is " << nodeId << " level " << level << " is over." << endl;
		/*cout <<"nodeId is "<<nodeId<<" level "<<level << " 排完序的s，SA，GRank是：" << endl;
		for (i = 0; i < n; i++)
		cout << chr(i) << ",," << SA[i] <<",,"<< suffixGIndex[i]<<",,"<<nodeId<< endl;*/
		delete[] bkt;
		delete[] t;
	}

	// compute SAl
	void induceSAl(unsigned char *t, uint_type *SA, unsigned char *s, uint_type *bkt,
		uint_type n, uint_type K, uint_type cs, uint_type level, uint_type *suffixGIndex, int32 nodeId, MPI_Status &status) {
		uint_type i, j, sendPos = 0, prePos = 0;
		map.clear();
		//uint_type sendLength = 2 * commSize * sizeof(uint_type) + commSize / 8 + ((level != 0) ? (commSize * sizeof(uint_type)) : commSize);
		getBuckets(s, bkt, n, K, cs, false); // find heads of buckets
		if (level == 0) bkt[0]++;
		//cout << "SA:"  <<n<< endl;
		//for (i = 0; i < n; i++)cout << SA[i] << ",";
		//cout << endl;
		//sentinal 不发送
		j = SA[0] - 1;
		SA[bkt[chr(j)]++] = j;
		for (i = 1; i<n; i++)
			if (SA[i] != EMPTY) {
				j = SA[i] - 1;
				if (j >= 0 && !tget(j)) SA[bkt[chr(j)]++] = j;
				//((uint_type *)bufSend)[sendPos] = SA[i];

				//((uint_type *)bufSend)[commSize + sendPos] = suffixGIndex[SA[i]];

				//bufset(commSize * sizeof(uint_type) * 2, sendPos, tget(SA[i]));//type
				//if (cs == sizeof(uint_type)) {
				//	((uint_type *)(bufSend + commSize * sizeof(uint_type) * 2 + commSize / 8))[sendPos] = ((uint_type *)s)[SA[i]];
				//}
				//else {
				//	bufSend[commSize * sizeof(uint_type) * 2 + commSize / 8 + sendPos] = s[SA[i]];
				//}
				if (map.find(SA[i] + 1) != map.end()) {
					if (level == 0)bufSend0[map[SA[i] + 1]].prePos = sendPos;
					else bufSend1[map[SA[i] + 1]].prePos = sendPos;
				}
					
				//else prePos = -1;

				bufSendSA[sendPos] = SA[i];
				if (level == 0) {
					bufSend0[sendPos].prePos = -1;
					bufSend0[sendPos].suffixGIndex = suffixGIndex[SA[i]+1];
					bufSend0[sendPos].type = tget(SA[i]);
					bufSend0[sendPos].data = s[SA[i]];
				}
				else {
					bufSend1[sendPos].prePos = -1;
					bufSend1[sendPos].suffixGIndex = suffixGIndex[SA[i]+1];
					bufSend1[sendPos].type = tget(SA[i]);
					bufSend1[sendPos].data = chr(SA[i]);
				}
				map[SA[i]] = sendPos;
				sendPos++;
				//如果缓冲区满了向服务器发送数据
				if (sendPos == commSize) {
					sTime = clock();
					if (level == 0) {
						//cout << "node " << nodeId << "发送的值是：" << endl;
						//for (j = 0; j < commSize; j++)bufSend0[j].printA();
						MPI_Send((char*)bufSend0, commSize * sizeof(SendData<unsigned char>), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);

					}
					else MPI_Send((char*)bufSend1, commSize * sizeof(SendData<uint_type>), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
					map.clear();
					MPI_Recv((char*)bufRecv, commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
					commTime += clock() - sTime;
					//更新suffixGIndex
					for (j = 0; j < commSize; j++) {
						suffixGIndex[bufSendSA[j]] = bufRecv[j];
						//if (bufSendSA[j] == 0)suffixGIndex[n - 1] = bufRecv[j];
						//else suffixGIndex[bufSendSA[j] - 1] = bufRecv[j];
					}
					//if (level == 0) {
					//	for (j = 0; j < commSize; j++) {
					//		suffixGIndex[bufSendSA[j]]= bufRecv[j];
					//		//if (bufSendSA[j] == 0)suffixGIndex[n - 1] = bufRecv[j];
					//		//else suffixGIndex[bufSendSA[j] - 1] = bufRecv[j];
					//	}
					//}
					//else {
					//	for (j = 0; j < commSize; j++) {
					//		suffixGIndex[bufSendSA[j]] = bufRecv[j];
					//		//if (bufSendSA[j] == 0)suffixGIndex[n - 1] = bufRecv[j];
					//		//else suffixGIndex[bufSendSA[j] - 1] = bufRecv[j];
					//	}
					//}

					sendPos = 0;
				}
			}

		//向服务器发送剩余的数据(加一个结束位)
		if (level == 0) {
			bufSend0[sendPos].prePos = commSize;
			bufSend0[sendPos].data = 0;
		}
		else {
			bufSend1[sendPos].prePos = commSize;
			bufSend1[sendPos].data = 0;
		}
		//((uint_type *)bufSend)[sendPos] = EMPTY;
		sendPos++;
		sTime = clock();
		if (level == 0) {
			//cout << "node " << nodeId << "发送的值是：" << endl;
			//for (j = 0; j < commSize; j++)bufSend0[j].printA();
			MPI_Send((char*)bufSend0, commSize * sizeof(SendData<unsigned char>), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);

		}
		else MPI_Send((char*)bufSend1, commSize * sizeof(SendData<uint_type>), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		map.clear();

		MPI_Recv((char*)bufRecv, commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		commTime += clock() - sTime;
		sendPos--;
		for (j = 0; j < sendPos; j++) {
			suffixGIndex[bufSendSA[j]] = bufRecv[j];
		}
		/*if (level == 0) {
			for (j = 0; j < sendPos; j++) {
				if (bufSendSA[j] == 0)suffixGIndex[n - 1] = bufRecv[j];
				else suffixGIndex[bufSendSA[j] - 1] = bufRecv[j];
			}
		}
		else {
			for (j = 0; j < sendPos; j++) {
				if (bufSendSA[j] == 0)suffixGIndex[n - 1] = bufRecv[j];
				else suffixGIndex[bufSendSA[j] - 1] = bufRecv[j];
			}
		}*/

	}

	void induceSAl1() {
		//int64 sendLength = 2 * commSize * sizeof(uint_type) + commSize / 8 + commSize * sizeof(uint_type);
		//int64 sendPos = 0;
		bufSend1[0].prePos = commSize;
		sTime = clock();
		MPI_Send((char*)bufSend1, commSize * sizeof(SendData<uint_type>), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		MPI_Recv((char*)bufRecv, commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		commTime += clock() - sTime;
	}

	// compute SAs
	void induceSAs(unsigned char *t, uint_type *SA, unsigned char *s, uint_type *bkt,
		uint_type n, uint_type K, uint_type cs, uint_type *suffixGIndex, int32 nodeId, MPI_Status &status, int32 level) {
		uint_type i, j, sendPos = 0, prePos = 0;
		map.clear();
		//uint_type sendLength = 2 * commSize * sizeof(uint_type) + commSize / 8 + ((cs == sizeof(uint_type)) ? (commSize * sizeof(uint_type)) : commSize);
		getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
		for (i = n - 1; i > 0; i--)
			if (SA[i] != EMPTY) {
				j = SA[i] - 1;
				if (j >= 0 && tget(j)) SA[bkt[chr(j)]--] = j;

				if (map.find(SA[i] + 1) != map.end()) {
					if (level == 0)bufSend0[map[SA[i] + 1]].prePos = sendPos;
					else bufSend1[map[SA[i] + 1]].prePos = sendPos;
				}
					//bufSend0[map[SA[i] + 1]].prePos = sendPos;

				bufSendSA[sendPos] = SA[i];
				if (level == 0) {
					bufSend0[sendPos].prePos = -1;
					bufSend0[sendPos].suffixGIndex = suffixGIndex[SA[i]+1];
					bufSend0[sendPos].type = tget(SA[i]);
					bufSend0[sendPos].data = s[SA[i]];
				}
				else {
					bufSend1[sendPos].prePos = -1;
					bufSend1[sendPos].suffixGIndex = suffixGIndex[SA[i]+1];
					bufSend1[sendPos].type = tget(SA[i]);
					bufSend1[sendPos].data = chr(SA[i]);
				}
				map[SA[i]] = sendPos;
				sendPos++;
				//如果缓冲区满了向服务器发送数据
				if (sendPos == commSize) {
					sTime = clock();
					if (level == 0)MPI_Send((char*)bufSend0, commSize * sizeof(SendData<unsigned char>), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
					else MPI_Send((char*)bufSend1, commSize * sizeof(SendData<uint_type>), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
					map.clear();

					MPI_Recv((char*)bufRecv, commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
					commTime += clock() - sTime;
					//recvPos += commSize;
					//更新suffixGIndex
					for (j = 0; j < commSize; j++) {
						suffixGIndex[bufSendSA[j]] = bufRecv[j];
					}
					/*if (level == 0) {
						for (j = 0; j < commSize; j++) {
							if (bufSendSA[j] == 0)suffixGIndex[n - 1] = bufRecv[j];
							else suffixGIndex[bufSendSA[j] - 1] = bufRecv[j];
						}
					}
					else {
						for (j = 0; j < commSize; j++) {
							if (bufSendSA[j] == 0)suffixGIndex[n - 1] = bufRecv[j];
							else suffixGIndex[bufSendSA[j] - 1] = bufRecv[j];
						}
					}*/
					sendPos = 0;
				}
			}

		//向服务器发送剩余的数据(加一个结束位)
		if (level == 0) {
			bufSend0[sendPos].prePos = commSize;
			bufSend0[sendPos].data = n;
		}
		else {
			bufSend1[sendPos].prePos = commSize;
			bufSend1[sendPos].data = n;
		}
		//((uint_type *)bufSend)[sendPos] = MAX;
		sendPos++;
		//向服务器发送剩余的数据
		sTime = clock();
		if (level == 0)MPI_Send((char*)bufSend0, commSize * sizeof(SendData<unsigned char>), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		else MPI_Send((char*)bufSend1, commSize * sizeof(SendData<uint_type>), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);

		MPI_Recv((char *)bufRecv, commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		commTime += clock() - sTime;
		sendPos--;
		for (j = 0; j < sendPos; j++) {
			suffixGIndex[bufSendSA[j]] = bufRecv[j];
		}
		/*if (level == 0) {
			for (j = 0; j < sendPos; j++) {
				if (bufSendSA[j] == 0)suffixGIndex[n - 1] = bufRecv[j];
				else suffixGIndex[bufSendSA[j] - 1] = bufRecv[j];
			}
		}
		else {
			for (j = 0; j < sendPos; j++) {
				if (bufSendSA[j] == 0)suffixGIndex[n - 1] = bufRecv[j];
				else suffixGIndex[bufSendSA[j] - 1] = bufRecv[j];
			}
		}*/

	}

	void induceSAs1() {
		//int64 sendLength = 2 * commSize * sizeof(uint_type) + commSize / 8 + commSize * sizeof(uint_type);
		//int64 sendPos = 0;
		bufSend1[0].prePos = -1;
		bufSend1[0].suffixGIndex = nodeId;
		bufSend1[0].type = 1;
		bufSend1[0].data = nodeId;
		//((uint_type *)bufSend)[commSize + sendPos] = nodeId;
		//bufset(commSize * sizeof(uint_type) * 2, sendPos, 1);//type
		//((uint_type *)(bufSend + commSize * sizeof(uint_type) * 2 + commSize / 8))[sendPos] = nodeId;
		//sendPos++;
		//((uint_type *)bufSend)[sendPos] = MAX;
		bufSend1[1].prePos = commSize;
		sTime = clock();
		MPI_Send((char*)bufSend1, commSize * sizeof(SendData<uint_type>), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		MPI_Recv((char*)bufRecv, commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		commTime += clock() - sTime;
	}

	void reName1() {
		bufSend[0] = nodeId;
		bufSend[1] = EMPTY;
		sTime = clock();
		MPI_Send((char *)(bufSend), commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		MPI_Recv((char *)(bufRecv), commSize * sizeof(uint_type), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		commTime += clock() - sTime;
	}
	void compute(string fileName) {
		char buf[100];
		sprintf_s(buf, "%d", nodeId);
		fileName = fileName + buf + ".txt";
		//cout << "fileName " << fileName<<endl;
		n = MyIO<char>::getFileLenCh(fileName);
		n++;
		unsigned char *data = new unsigned char[n];
		MyIO<unsigned char>::read(data, fileName, n - 1, std::ios_base::in | std::ios_base::binary, 0);
		data[n - 1] = 0;//add sentinal
		uint_type *SA = new uint_type[n];
		uint_type *suffixGIndex = new uint_type[n];
		clock_t start, finish;
		start = clock();
		computeSA(data, SA, suffixGIndex, n, sizeof(unsigned char), 256, 0);
		finish = clock();
		cout << "node " << nodeId << " completed, total time is" << (double)(finish - start) / CLOCKS_PER_SEC << ", communicate time is " << (double)(commTime) / CLOCKS_PER_SEC << endl;
		string writefn = "D://test//result";
		char buf1[100];
		sprintf_s(buf1, "%d", nodeId);
		writefn = writefn + buf1 + ".txt";
		MyIO<uint_type>::write(suffixGIndex, writefn, n, std::ios_base::out | std::ios_base::binary, 0);

		delete[] suffixGIndex;
		delete[] SA;
		delete[] data;

	}


	~ChildNode()
	{

		delete[] bufRecv;
		delete[] bufSend;
	}
};
#endif