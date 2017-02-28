#ifndef  __CHILDNODE
#define  __CHILDNODE

#include <iostream>
#include <string>
#include "mycommon.h"
#include <vector>
#include <queue>
#include <functional> 


using namespace std;

#define tget(i) ( (t[(i)/8]&mask[(i)%8]) ? 1 : 0 )
#define tset(i, b) t[(i)/8]=(b) ? (mask[(i)%8]|t[(i)/8]) : ((~mask[(i)%8])&t[(i)/8])
#define chr(i) (cs==sizeof(int64)?((int64*)s)[i]:((unsigned char *)s)[i])
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
	int64 n;
	unsigned char *bufRecv = new unsigned char[commSize * sizeof(int64)];
	unsigned char *bufSend = new unsigned char[commSize * sizeof(int64) * 4];
	
public:
	ChildNode(int32 nodeIdC){
		nodeId = nodeIdC;
	}


	void getBuckets(unsigned char *s, int64 *bkt, int64 n, int64 K, int32 cs, bool end) {
		int64 i, sum = 0;
		for (i = 0; i <= K; i++) bkt[i] = 0; // clear all buckets
		for (i = 0; i<n; i++) bkt[chr(i)]++; // compute the size of each bucket
		for (i = 0; i <= K; i++) { sum += bkt[i]; bkt[i] = end ? sum - 1 : sum - bkt[i]; }
	}


	void computeSA(unsigned char *s, int64 *SA, int64 *suffixGIndex, int64 n, int32 cs, int64 K, int64 level) {
		//发送大小
		MPI_Send(&n, 1, MPI_LONG, MAINNODEID, nodeId, MPI_COMM_WORLD);

		//初始化数据
		int64 *bkt = new int64[K + 1]; // bucket counters
		
		unsigned char *t = new unsigned char[n / 8 + 1];
		
		int64 i,j = 0;
		tset(n - 2, 0); tset(n - 1, 1); // the sentinel must be in s1, important!!!
		for (i = n - 3; i >= 0; i--) tset(i, (chr(i)<chr(i + 1) || (chr(i) == chr(i + 1) && tget(i + 1) == 1)) ? 1 : 0);
		int64 maxVal = 0;
		for (i = n - 1; i >= 0; i--)if (chr(i) > maxVal)maxVal = chr(i);
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
		for (i = 1; i < n-1; i++)
			if (isLMS(i)) {
				SA[bkt[chr(i)]--] = i;
				suffixGIndex[i] = 0;
			}
		SA[0] = n - 1; // set the single sentinel LMS-substring
		suffixGIndex[n - 2] = nodeId;
		
		//cout << "level is "<<level <<" nodeId is "<<nodeId<<" 开始induceL..." << endl;
		induceSAl(t,SA,s,bkt,n,K,cs,level, suffixGIndex, nodeId, status);
		//将后缀排名向后移一位，得到本字符的全局排名
		j = suffixGIndex[n-1];
		for (i = n - 2; i >= 0; i--)suffixGIndex[i+1] = suffixGIndex[i];
		suffixGIndex[0] = j;
		suffixGIndex[n - 1] = nodeId;
		//cout << "level is " << level << " nodeId is " << nodeId << " 开始induceS..." << endl;
		induceSAs(t, SA, s, bkt, n, K, cs, suffixGIndex, nodeId, status);

#ifdef RENAMEDEBUGC
		cout << "level is " << level << " 子节点 "<<nodeId<<" 开始重命名:" << std::endl;
#endif
		//因为每次迭代最少减少一半，所以用SA前半部分存储本地SA，后半部分存储全局排名
		int64 n1 = 0, sendPos = 0, offset = 0;
		for (i = 0; i<n; i++)
			if (isLMS(SA[i])) {
				SA[n1++] = SA[i];
			}
		//cout << "level is " << level << " 子节点 " << nodeId << " 开始重命名:" << std::endl;

		offset = n - n1;
		for (i = 0; i < n1;i++)SA[offset++] = suffixGIndex[SA[i] - 1];
		SA[n - n1] = nodeId;

		//利用suffixGIndex发送全局排名，重命名全局排名
		for (i = n - n1, j = 0; i < n; i++)suffixGIndex[j++] = SA[i];
		suffixGIndex[n1] = EMPTY;
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
		sendPos = 0;
		j = n1 + 1;
		int64 p = 0;
		while (sendPos < j) {
			//cout << "开始发送数据：nodeId is "<<nodeId <<" n1 is "<< n1 << endl;
			/*cout << "发送的数据：nodeId is " << nodeId << endl;
			for (i = sendPos; i < sendPos + commSize; i++) {
				cout << suffixGIndex[i] << ",";
			}
			cout << endl;*/
			if ((j - sendPos) < commSize) {
				for (i = sendPos,p=0; i < j; i++)((int64 *)bufSend)[p++] = suffixGIndex[i];
				MPI_Send((char *)(bufSend), commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
				MPI_Recv((char *)(bufRecv), commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
				for (i = sendPos,p=0; i < j; i++)suffixGIndex[i] = ((int64 *)bufRecv)[p++];
				sendPos = j;
			}
			else {
				MPI_Send((char *)(suffixGIndex + sendPos), commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
				MPI_Recv((char *)(suffixGIndex + sendPos), commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
				sendPos += commSize;
			}
			
			
		}
#ifdef RENAMEDEBUGC
		if (level > LL) {
			cout << "接收到的重命名后的数据是：" << endl;
			for (i = 0; i < n1; i++) {
				cout << suffixGIndex[i] << endl;
			}
		}

#endif
		//将lms的SA存入suffixGIndex的末尾
		j = 0;
		for (i = n - n1; i < n; i++) {
			suffixGIndex[i] = SA[j++];
			//现在suffixGIndex前n1个字符是服务器重命名后的值，最后n1个字符是lms字符的SA
		}
		for (i = 0; i < n; i++)SA[i] = EMPTY;

		j = 0;
		for (i = n - n1; i < n; i++) {
			SA[suffixGIndex[i]]= suffixGIndex[j++];
		}

		for (i = n - 1, j = n - 1; i >= 0; i--)
			if (SA[i] != EMPTY) SA[j--] = SA[i];//此时SA最后n1个字符就是lms字符串的全局排名的顺序

		

		int64 *SA1 = SA, *s1 = SA + n - n1;

		int64 isLoop = 0;
		MPI_Recv(&(isLoop), 1, MPI_LONG, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		if (isLoop) {
			cout << "level is " << level << " nodeId is  " << nodeId << " start recycle " << endl;
			/*cout << "level is " << level << "nodeId is  "<<nodeId<<" 开始递归 s1是：" << endl;
			for (i =0; i < n1; i++) {
				cout << s1[i] << ",";
			}
			cout << endl;*/
			computeSA((unsigned char *)s1,SA1, suffixGIndex,n1,sizeof(int64), isLoop,level+1);
		}
		else {
			cout << "level is " << level << " don't need recycle" << endl;

			//不需要递归时才求本地排名

			//将全局排名先存入suffixGIndex前n1个位置
			
			for (i = n - n1, j = 0; i < n; i++) {
				suffixGIndex[j++] = SA[i];
			}
			//此时suffixGIndex前n1是lms字符的全局排名，后n1字符是lms字符的SA
			//利用SA数组求本地排名
			for (i = 0; i < n; i++)SA[i] = EMPTY;
			
			for (i = n - n1, j = 0; i < n; i++) {
				SA[suffixGIndex[i]] = j++;
			}

			for (i = 0, j = 0; i < n; i++)
				if (SA[i] != EMPTY) SA[j++] = SA[i];//此时SA前n1个字符就是本地排名的顺序

		
		}
		//cout << "level is " << level << " nodeId is "<<nodeId<<" start to put lms character.." << endl;
		//suffixGIndex是全局排名 SA1是本地排名 s1是全局重命名后的字符串
		
		//调整全局排名的位置到SA的最后端
		for (i = n - n1, j = 0; i < n; i++)SA[i] = suffixGIndex[j++];

		/*cout << "level is " << level <<"nodeId is "<<nodeId  << " 全局重命名后lms 字符的全局排名和本地排名是:" << endl;
		for (i = 0; i < n1; i++) {
			cout <<suffixGIndex[i] << "," << SA1[i]<<","<< nodeId << endl;
		}*/

		 // put all left-most S characters int64o their buckets
		getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
		
			 
		//利用SuffixGIndex来保存本地的lms字符的SA
		
		int64 *SAT = new int64[n1];
		s1 = SAT;
		for (j = 0,i = 1; i < n; i++) {
			if (isLMS(i)) {
				s1[j++] = i; // get p1
			}
		}
		/*cout << "level is " << level << " lsm字符的SA是:" << endl;
		for (i =0; i < n1; i++) {
			cout << s1[i] << endl;
		}*/
		for (i = 0; i < n; i++)suffixGIndex[i] = EMPTY;
		
		for (i = n - n1, j = 0; i<n; i++) suffixGIndex[s1[j++]] = SA[i]; // lms字符的suffixGIndex已放好

		
		for (i = n-n1, j = 0; i<n; i++) SA[i] = SA1[j++]; // 
		SA1 = SA + n - n1;

		
		for (i = 0; i < n1; i++)SA[i] = s1[SA1[i]];

		for (i = 0; i < n1; i++) {
			SAT[i]=SA[i];
		}
#ifdef RENAMEDEBUGC
		if (level > LL) {
			cout <<"nodeID is "<<nodeId<< "SAT is :" << endl;
			for (i = 0; i < n1; i++) {
			cout<<SAT[i]<<",";
			}
			cout << endl;
		}
#endif // RENAMEDEBUGC

		
		for (i = 0; i<n; i++) SA[i] = EMPTY; // init SA[n1..n-1]
		//j = 0;
		for (i = n1 - 1; i >= 0; i--) {
			j = SAT[i]; 
			SA[bkt[chr(j)]--] = j;
		}
		SA[0] = n - 1;
		delete[] SAT;

		
		suffixGIndex[n - 2] = nodeId;
		
#ifdef RENAMEDEBUGC
		cout << "level is " << level << " nodeId is " << nodeId << " lms字符已放好，开始排L n is "<<n<<" n1 is "<< n1 << endl;
#endif // RENAMEDEBUGC
		//
		//cout << "level is " << level << " nodeId is " << nodeId << " lms字符已放好，开始排L n is " << n << " n1 is " << n1 << endl;
		induceSAl(t, SA, s, bkt, n, K, cs, level, suffixGIndex, nodeId, status);

		j = suffixGIndex[n - 1];
		for (i = n - 2; i >= 0; i--)suffixGIndex[i + 1] = suffixGIndex[i];
		suffixGIndex[0] = j;
		suffixGIndex[n - 1] = nodeId;

#ifdef RENAMEDEBUGC
		cout << "level is " << level << " nodeId is " << nodeId << " lms字符已放好，开始排S" << endl;
#endif
		//cout << "level is " << level << " nodeId is " << nodeId << " lms字符已放好，开始排S" << endl;
		induceSAs(t, SA, s, bkt, n, K, cs, suffixGIndex, nodeId, status);
#ifdef RENAMEDEBUGC
		cout << "level is " << level << " nodeId is " << nodeId << " 排序结束" << endl;
#endif
		j = suffixGIndex[n - 1];
		for (i = n-2; i >=0; i--) {
			suffixGIndex[i+1] = suffixGIndex[i];
		}
		suffixGIndex[0] = j;

		suffixGIndex[n - 1] = nodeId; 
		//cout << "nodeId is " << nodeId << " level " << level << " is over." << endl;
		/*cout <<"nodeId is "<<nodeId<<" level "<<level << " 排完序的s，SA，GRank是：" << endl;
		for (i = 0; i < n; i++)
			cout << chr(i) << ",," << SA[i] <<",,"<< suffixGIndex[i]<<",,"<<nodeId<< endl;*/
		delete[] bkt;
		delete[] t;
	}

	// compute SAl
	void induceSAl(unsigned char *t, int64 *SA, unsigned char *s, int64 *bkt,
		int64 n, int64 K, int64 cs, int64 level, int64 *suffixGIndex,int32 nodeId, MPI_Status &status) {
		int64 i, j, sendPos = 0;
		int64 sendLength = 2 * commSize * sizeof(int64) + commSize / 8  + ((level != 0) ? (commSize * sizeof(int64)) : commSize);
		getBuckets(s, bkt, n, K, cs, false); // find heads of buckets
		if (level == 0) bkt[0]++;
		//sentinal 不发送
		j = SA[0] - 1;
		SA[bkt[chr(j)]++] = j;
		for (i = 1; i<n; i++)
			if (SA[i] != EMPTY) {
				j = SA[i] - 1;
				if (j >= 0 && !tget(j)) SA[bkt[chr(j)]++] = j;
				((int64 *)bufSend)[sendPos] = SA[i];
				
				((int64 *)bufSend)[commSize + sendPos] = suffixGIndex[SA[i]];

				bufset(commSize * sizeof(int64) * 2, sendPos, tget(SA[i]));//type
				if (cs == sizeof(int64)) {
					((int64 *)(bufSend + commSize * sizeof(int64) * 2 + commSize / 8 ))[sendPos] = ((int64 *)s)[SA[i]];
				}
				else {
					bufSend[commSize * sizeof(int64) * 2 + commSize / 8  + sendPos] = s[SA[i]];
				}
				sendPos++;
				//如果缓冲区满了向服务器发送数据
				if (sendPos == commSize) {
					MPI_Send((char*)bufSend, sendLength, MPI_CHAR, MAINNODEID, nodeId , MPI_COMM_WORLD);
					MPI_Recv((char*)bufRecv, commSize * sizeof(int64) , MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
					//更新suffixGIndex
					for (j = 0; j < commSize; j++) {
						if(((int64 *)bufSend)[j] == 0)suffixGIndex[n - 1] = ((int64 *)bufRecv)[j];
						else {
							suffixGIndex[((int64 *)bufSend)[j] - 1] = ((int64 *)bufRecv)[j];
						}
					}
					sendPos = 0;
				}
			}
		
		//向服务器发送剩余的数据(加一个结束位)
		((int64 *)bufSend)[sendPos] = EMPTY;
		sendPos++;
		MPI_Send((char*)bufSend, sendLength, MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		MPI_Recv((char*)bufRecv, commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		sendPos--;
		for (j = 0; j < sendPos; j++) {
			if (((int64 *)bufSend)[j] == 0)suffixGIndex[n - 1] = ((int64 *)bufRecv)[j];
			else {
				suffixGIndex[((int64 *)bufSend)[j] - 1] = ((int64 *)bufRecv)[j];
			}
		}
		
	}

	// compute SAs
	void induceSAs(unsigned char *t, int64 *SA, unsigned char *s, int64 *bkt,
		int64 n, int64 K, int64 cs, int64 *suffixGIndex, int32 nodeId, MPI_Status &status) {
		int64 i, j, sendPos = 0, recvPos = 0;
		int64 sendLength = 2 * commSize * sizeof(int64) + commSize / 8  + ((cs == sizeof(int64)) ? (commSize * sizeof(int64)) : commSize);
		getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
		for (i = n - 1; i > 0; i--)
			if (SA[i] != EMPTY) {
				j = SA[i] - 1;
				if (j >= 0 && tget(j)) SA[bkt[chr(j)]--] = j;

				((int64 *)bufSend)[sendPos] = SA[i];
				((int64 *)bufSend)[commSize + sendPos] = suffixGIndex[SA[i]];
				bufset(commSize * sizeof(int64) * 2, sendPos, tget(SA[i]));//type
				if (cs == sizeof(int64)) {
					((int64 *)(bufSend + commSize * sizeof(int64) * 2 + commSize / 8))[sendPos] = ((int64 *)s)[SA[i]];
				}
				else {
					bufSend[commSize * sizeof(int64) * 2 + commSize / 8 + sendPos] = s[SA[i]];
				}
				sendPos++;
				//如果缓冲区满了向服务器发送数据
				if (sendPos == commSize) {
					MPI_Send((char*)bufSend, sendLength, MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
					MPI_Recv((char*)bufRecv, commSize * sizeof(int64) , MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
					//recvPos += commSize;
					//更新suffixGIndex
					for (j = 0; j < commSize; j++) {
						if (((int64 *)bufSend)[j] == 0)suffixGIndex[n - 1] = ((int64 *)bufRecv)[j];
						else {
							suffixGIndex[((int64 *)bufSend)[j] - 1] = ((int64 *)bufRecv)[j];
						}
					}
					sendPos = 0;
				}
			}

		//向服务器发送剩余的数据(加一个结束位)
		((int64 *)bufSend)[sendPos] = MAX;
		sendPos++;
		//向服务器发送剩余的数据
		MPI_Send((char *)bufSend, sendLength, MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		MPI_Recv((char *)bufRecv, commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		sendPos--;
		for (j = 0; j < sendPos; j++) {
			if (((int64 *)bufSend)[j] == 0)suffixGIndex[n - 1] = ((int64 *)bufRecv)[j];
			else {
				suffixGIndex[((int64 *)bufSend)[j] - 1] = ((int64 *)bufRecv)[j];
			}
		}
		
	}

	void compute(string fileName) {
		char buf[100];
		sprintf_s(buf, "%d", nodeId);
		fileName = fileName + buf + ".txt";
		n = MyIO<char>::getFileLenCh(fileName);
		n++;
		unsigned char *data = new unsigned char[n];
		MyIO<unsigned char>::read(data, fileName, n - 1, std::ios_base::in | std::ios_base::binary, 0);
		data[n - 1] = 0;//add sentinal
		int64 *SA = new int64[n];
		int64 *suffixGIndex = new int64[n];
	
		computeSA(data,SA, suffixGIndex,n,sizeof(unsigned char),256,0);

		string writefn = "D://test//result";
		char buf1[100];
		sprintf_s(buf1, "%d", nodeId);
		writefn = writefn + buf1 + ".txt";
		MyIO<int64>::write(suffixGIndex, writefn, n, std::ios_base::out | std::ios_base::binary, 0);

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