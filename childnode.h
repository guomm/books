#ifndef  __CHILDNODE
#define  __CHILDNODE

#include <iostream>
#include <string>
#include "mycommon.h"
//#include "sais.h"
#include <vector>
#include <queue>
#include <functional> 


using namespace std;
//#include "mycommon.h"


//#define setchr(i,chr) (cs==sizeof(int64)?((int64*)bufSend)[i]:((unsigned char *)bufSend)[i])

#define tget(i) ( (t[(i)/8]&mask[(i)%8]) ? 1 : 0 )
#define tset(i, b) t[(i)/8]=(b) ? (mask[(i)%8]|t[(i)/8]) : ((~mask[(i)%8])&t[(i)/8])

//#define isGget(i) ( (isG[(i)/8]&mask[(i)%8]) ? 1 : 0 )
//#define isGset(i, b) isG[(i)/8]=(b) ? (mask[(i)%8]|isG[(i)/8]) : ((~mask[(i)%8])&isG[(i)/8])

#define chr(i) (cs==sizeof(int64)?((int64*)s)[i]:((unsigned char *)s)[i])
#define isLMS(i) (i>0 && tget(i) && !tget(i-1))

#define bufset(i,j, b) bufSend[i+j/8]=(b) ? (mask[(j)%8]|bufSend[i+j/8]) : ((~mask[(j)%8])&bufSend[i+j/8])

//#define SALDEBUG

//#define SADEBUG
//#define SALDEBUG
//#define SASDEBUG
#define RENAMEDEBUG


class ChildNode {
private:
	int32 nodeId;
	MPI_Status status;
	int64 n;
	unsigned char *bufRecv = (unsigned char *)malloc(commSize * sizeof(int64));//GRank;;
	unsigned char *bufSend = (unsigned char *)malloc(1024*8);//SA+suffixGIndex+Type+isG+s;

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


	void computeSA(unsigned char *s, int64 *SA, int64 *GRank, int64 *suffixGIndex, int64 n, int32 cs, int64 K, int64 level) {
		//���ʹ�С
		MPI_Send(&n, 1, MPI_LONG, MAINNODEID, nodeId, MPI_COMM_WORLD);

		//��ʼ������
		int64 *bkt = (int64 *)malloc(sizeof(int64)*(K + 1)); // bucket counters
		//unsigned char *data = (unsigned char *)malloc(2 * n * sizeof(int64) + n / 8 + 1 + n / 8 + 1 + (cs == sizeof(int64)) ? (n * sizeof(int64)) : n);//SA+suffixGIndex+Type+isG+s
		
		int32 bufPos=0;
		unsigned char *t = (unsigned char *)malloc(n / 8 + 1);
		unsigned char *isG = (unsigned char *)malloc(n / 8 + 1);
		//int64 *SA = (int64 *)data;
		
		int64 i = 0;
		tset(n - 2, 0); tset(n - 1, 1); // the sentinel must be in s1, important!!!
		//isGset(n - 2, 0); isGset(n - 1, 0);
		for (i = n - 3; i >= 0; i--) tset(i, (chr(i)<chr(i + 1) || (chr(i) == chr(i + 1) && tget(i + 1) == 1)) ? 1 : 0);
		getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
		for (i = 0; i < n; i++) {
			SA[i] = EMPTY; GRank[i] = EMPTY; suffixGIndex[i] = EMPTY;
		}
		for (i = 1; i < n-1; i++)
			if (isLMS(i)) {
				SA[bkt[chr(i)]--] = i;
				suffixGIndex[i] = 0;
			//	isGset(i, 1);
			}
		SA[0] = n - 1; // set the single sentinel LMS-substring
		GRank[0] = nodeId;
		suffixGIndex[SA[0] - 1] = nodeId;
		//isGset(SA[0] - 1, 1);
		/*for (i = 0; i < n; i++) {
			cout << SA[i] <<"," << suffixGIndex[SA[i]]<< endl;
		}*/
		induceSAl(t,SA,s,bkt,n,K,cs,level, isG, GRank, suffixGIndex, nodeId, status);

		induceSAs(t, SA, s, bkt, n, K, cs, isG, GRank, suffixGIndex, nodeId, status);

		free(bkt);
#ifdef RENAMEDEBUG
		cout << "�ӽڵ� "<<nodeId<<" ��ʼ������:" << std::endl;
#endif
		//��Ϊÿ�ε������ټ���һ�룬������SAǰ�벿�ִ洢����SA����벿�ִ洢ȫ������
		int64 n1 = 0, sendPos = 0, offset = n / 2, desOffset = 0;
		for (i = 0; i<n; i++)
			if (isLMS(SA[i])) {
				SA[n1++] = SA[i];
				SA[offset++] = suffixGIndex[SA[i] - 1];
			}
		//		GRank[n1++] = GRank[i];
		SA[offset] = EMPTY;
		//����suffixGIndex����ȫ��������������ȫ������
		offset = n / 2;
		desOffset = offset + n1;
		int64 j = 0;
		for (i = offset; i < desOffset; i++)suffixGIndex[j++] = SA[i];
		suffixGIndex[n1] = EMPTY;
#ifdef RENAMEDEBUG
		cout << "������:���͵ĵ������ǣ�" << endl;
		for (i = 0; i < n1; i++) {
			cout << suffixGIndex[i] << endl;
		}

#endif
		sendPos = 0;
		//MPI_Send((char *)(suffixGIndex + sendPos), commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		//MPI_Recv((char *)(suffixGIndex + sendPos), commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		//cout<<
		while (sendPos <= n1) {
			//cout << "��ʼ�������ݣ�nodeId is "<<nodeId <<" n1 is "<< n1 << endl;
			MPI_Send((char *)(suffixGIndex + sendPos), commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		//	cout << "------" << endl;
			MPI_Recv((char *)(suffixGIndex + sendPos), commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
			//cout << "���յ����ݣ�" << endl;
			sendPos += commSize;
		}
		
		//��lms��SA����suffixGIndex��ĩβ
		j = 0;
		for (i = n - n1; i < n; i++) {
			suffixGIndex[i] = SA[j++];
			//����suffixGIndexǰn1���ַ��Ƿ��������������ֵ�����n1���ַ���lms�ַ���SA
		}
		for (i = 0; i < n; i++)SA[i] = EMPTY;

		j = 0;
		for (i = n - n1; i < n; i++) {
			SA[suffixGIndex[i]]= suffixGIndex[j++];
		}

		for (i = n - 1, j = n - 1; i >= 0; i--)
			if (SA[i] != EMPTY) SA[j--] = SA[i];//��ʱSA���n1���ַ�����lms�ַ�����ȫ��������˳��

		

		int64 *SA1 = SA, *s1 = SA + n - n1;
		////��ȫ��������������ݷ���SA�����
		//j = 0;
		//for (i = n - n1; i < n; i++) {
		//	SA[i] = suffixGIndex[j++];//ÿ�ε������ټ���һ��
		//}

		////��lms�ַ���SA����suffixGIndex�����(Ϊ���󱾵ص�SA����)
		//j = 0;
		//for (i = n - n1; i < n; i++) {
		//	suffixGIndex[i] = SA[j++];//ÿ�ε������ټ���һ��
		//}

#ifdef RENAMEDEBUG
		cout << "���յ�����������������ǣ�" << endl;
		for (i = n - n1; i < n; i++) {
			cout << SA[i] << endl;
		}
		

#endif
//		int64 *SA1 = SA, *s1 = SA + n - n1;


		bool isLoop = 0;
		MPI_Recv(&(isLoop), 1, MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		if (isLoop) {
			cout << "��ʼ�ݹ�" << endl;
			computeSA((unsigned char *)s1,SA1, GRank, suffixGIndex,n1,sizeof(int64), n1,level+1);
		}
		else {
			cout << "����Ҫ�ݹ�" << endl;

			//����Ҫ�ݹ�ʱ���󱾵�����
			//int64 *SATemp = new int64[n1];

			//��ȫ�������ȴ���suffixGIndexǰn1��λ��
			j = 0;
			for (i = n - n1; i < n; i++) {
				suffixGIndex[j++] = SA[i];
			}
			//��ʱsuffixGIndexǰn1��lms�ַ���ȫ����������n1�ַ���lms�ַ���SA
			//����SA�����󱾵�����
			for (i = 0; i < n; i++)SA[i] = EMPTY;
			j = 0;
			for (i = n - n1; i < n; i++) {
				SA[suffixGIndex[i]] = j++;
			}

			for (i = 0, j = 0; i < n; i++)
				if (SA[i] != EMPTY) SA[j++] = SA[i];//��ʱSAǰn1���ַ����Ǳ���������˳��

			//����ȫ��������λ�õ�SA������
			for (i = n - n1, j = 0; i < n; i++)SA[i] = suffixGIndex[j++];

			//��suffixGIndexĩβ�����lms�ַ���SA����SA����
			//j = 0;
			//for (i = n - n1; i < n; i++) {
			//	SA[j++] = suffixGIndex[i];
			//}
			////��s1����suffixGIndex�б���,��s1��SA1�ֱ�洢
			//for (i = 0; i < n1; i++) {
			//	suffixGIndex[i] = s1[i];
			//}
			////s1 = suffixGIndex;

			////��SA1������������
			//for (i = n1; i<n; i++) SA[i] = EMPTY;
			//int64 name = 0, prev = -1;
			//for (i = 0; i<n1; i++) {
			//	int64 pos = SA[i]; bool diff = false;
			//	for (int64 d = 0; d<n; d++)
			//		if (prev == -1 || pos + d == n - 1 || prev + d == n - 1 ||
			//			chr(pos + d) != chr(prev + d) ||
			//			tget(pos + d) != tget(prev + d))
			//		{
			//			diff = true; break;
			//		}
			//		else
			//			if (d>0 && (isLMS(pos + d) || isLMS(prev + d)))
			//				break;

			//	if (diff)
			//	{
			//		name++; prev = pos;
			//	}
			//	pos = pos / 2; //(pos%2==0)?pos/2:(pos-1)/2;
			//	SA[n1 + pos] = name - 1;
			//}
			//for (i = n - 1, j = n - 1; i >= n1; i--)
			//	if (SA[i] != EMPTY) SA[j--] = SA[i];
			//SA1= SA + n - n1;

			////����SA1��λ�ã�SA1(�ֲ�����)��SA�����ǰ�벿��,s1��ȫ����������SA����ĺ�벿��
			//for (i = 0; i < n1; i++) {
			//	SA[i] = SA1[i];
			//}
			//SA1 = SA;
			//j = 0;
			//for (i = n- n1; i < n; i++) {
			//	SA[i] = suffixGIndex[j++];
			//}
			//s1 = SA + n - n1;
		}

		//s1��ȫ������ SA1�Ǳ�������
		
		cout << " lms �ַ���ȫ�������ͱ���������:" << endl;
		for (i = 0; i < n1; i++) {
			cout << s1[i] << "," << SA1[i] << endl;
		}

		bkt = (int64 *)malloc(sizeof(int64)*(K + 1)); // bucket counters

														  // put all left-most S characters int64o their buckets
		getBuckets(s, bkt, n, K, cs, true); // find ends of buckets
		j = 0;
			 
		//����SuffixGIndex�����汾�ص�lms�ַ���SA
		s1 = suffixGIndex;

		for (i = 1; i < n; i++) {
			if (isLMS(i)) {
				s1[j++] = i; // get p1
			}
		}
		j = 0;
		for (i = 0; i<n1; i++) SA1[i] = s1[SA1[i]]; // get index in s1
		//for (i = 0; i<n1; i++) SA1[i] = s1[SA1[i]]; // get index in s1
		for (i = 0; i < n; i++)suffixGIndex[i] = EMPTY;
		//��ȫ��������ֵ��suffixGIndex
		for (i = n - n1; i<n; i++) suffixGIndex[SA1[j++]] = SA[i]; // init SA[n1..n-1]
		for (i = n1; i<n; i++) SA[i] = EMPTY; // init SA[n1..n-1]
		j = 0;
		for (i = n1 - 1; i >= 0; i--) {
			j = SA[i]; SA[i] = EMPTY;
			if (level == 0 && i == 0) {
				SA[0] = n - 1;
			}
			else {
				SA[bkt[chr(j)]--] = j;
			}
		}
		induceSAl(t, SA, s, bkt, n, K, cs, level, isG, GRank, suffixGIndex, nodeId, status);

		induceSAs(t, SA, s, bkt, n, K, cs, isG, GRank, suffixGIndex, nodeId, status);

		free(bkt);
		free(t);
		free(isG);
		free(bufSend);
	}

	// compute SAl
	void induceSAl(unsigned char *t, int64 *SA, unsigned char *s, int64 *bkt,
		int64 n, int64 K, int64 cs, int64 level, unsigned char *isG, int64 *GRank , int64 *suffixGIndex,int32 nodeId, MPI_Status &status) {
		int64 i, j, sendPos = 0;
		int64 sendLength = 2 * commSize * sizeof(int64) + commSize / 8  + ((cs == sizeof(int64)) ? (commSize * sizeof(int64)) : commSize);
		getBuckets(s, bkt, n, K, cs, false); // find heads of buckets
		if (level == 0) bkt[0]++;

		//sentinal ������
		j = SA[0] - 1;
		SA[bkt[chr(j)]++] = j;

		for (i = 1; i<n; i++)
			if (SA[i] != EMPTY) {
				j = SA[i] - 1;
				if (j >= 0 && !tget(j)) SA[bkt[chr(j)]++] = j;
				((int64 *)bufSend)[sendPos] = SA[i];
				
				((int64 *)bufSend)[commSize + sendPos] = suffixGIndex[SA[i]];

				//cout << SA[i] << "," << ((int64 *)bufSend)[sendPos] << endl;

				bufset(commSize * sizeof(int64) * 2, sendPos, tget(SA[i]));//type

				//cout << SA[i] << "," << ((int64 *)bufSend)[sendPos] << endl;
				//bufset(commSize * sizeof(int64) * 2 + commSize / 8 + sendPos, isGget(i));//type
				if (cs == sizeof(int64)) {
					((int64 *)bufSend)[commSize * 2 + commSize / 8  / 8 + sendPos] = ((int64 *)s)[SA[i]];
				}
				else {
					bufSend[commSize * sizeof(int64) * 2 + commSize / 8  + sendPos] = s[SA[i]];
				}
				//cout <<SA[i]<<","<< ((int64 *)bufSend)[sendPos] << "," << ((int64 *)bufSend)[commSize + sendPos] << "," << tget(SA[i]) << ","<< bufSend[commSize * sizeof(int64) * 2 + commSize / 8 + sendPos] << endl;
				sendPos++;
				//����������������������������
				if (sendPos == commSize) {
					//cout << "���ͳ�ʼ������:" << sendLength << " sizeof64: "<<sizeof(int64)<< endl;
					MPI_Send((char*)bufSend, sendLength, MPI_CHAR, MAINNODEID, nodeId , MPI_COMM_WORLD);
					MPI_Recv((char*)GRank, commSize * sizeof(int64) , MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
					//recvPos += commSize;
					//����suffixGIndex
					for (j = 0; j < commSize; j++) {
						if(((int64 *)bufSend)[j] == 0)suffixGIndex[n - 1] = GRank[j];
						else {
							suffixGIndex[((int64 *)bufSend)[j] - 1] = GRank[j];
						}
							
						//isGset(((int64 *)bufSend)[j] - 1, 1);
						//cout << GRank[j] << endl;
					}
					sendPos = 0;
				}
			}

		//�����������ʣ�������(��һ������λ)
		((int64 *)bufSend)[sendPos] = EMPTY;
		sendPos++;
		MPI_Send((char*)bufSend, sendLength, MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		MPI_Recv((char*)GRank, commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		sendPos--;
		for (j = 0; j < sendPos; j++) {
			suffixGIndex[((int64 *)bufSend)[j] - 1] = GRank[j];
		}
		
	}

	// compute SAs
	void induceSAs(unsigned char *t, int64 *SA, unsigned char *s, int64 *bkt,
		int64 n, int64 K, int64 cs, unsigned char *isG, int64 *GRank, int64 *suffixGIndex, int32 nodeId, MPI_Status &status) {
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
				//bufset(commSize * sizeof(int64) * 2 + commSize / 8 + sendPos, isGget(i));//type
				if (cs == sizeof(int64)) {
					((int64 *)bufSend)[commSize * 2 + commSize / 8  / 8 + sendPos] = ((int64 *)s)[SA[i]];
				}
				else {
					bufSend[commSize * sizeof(int64) * 2 + commSize / 8 + sendPos] = s[SA[i]];
				}
				sendPos++;
				//����������������������������
				if (sendPos == commSize) {
					MPI_Send((char*)bufSend, sendLength, MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
					MPI_Recv((char*)GRank, commSize * sizeof(int64) , MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
					//recvPos += commSize;
					//����suffixGIndex
					for (j = 0; j < commSize; j++) {
						suffixGIndex[((int64 *)bufSend)[j] - 1] = GRank[j];
						//isGset(((int64 *)bufSend)[j] - 1, 1);
					}
					sendPos = 0;
				}
			}

		//�����������ʣ�������(��һ������λ)
		((int64 *)bufSend)[sendPos] = EMPTY;
		sendPos++;
		//�����������ʣ�������
		MPI_Send((char *)bufSend, sendLength, MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD);
		MPI_Recv((char *)GRank, commSize * sizeof(int64), MPI_CHAR, MAINNODEID, nodeId, MPI_COMM_WORLD, &status);
		sendPos--;
		for (j = 0; j < sendPos; j++) {
			suffixGIndex[((int64 *)bufSend)[j] - 1] = GRank[j];
		}
		
	}

	void compute(string fileName) {
		char buf[100];
		sprintf_s(buf, "%d", nodeId);
		fileName = fileName + buf + ".txt";
		//cout << "....fileName..." << endl;
		//this->nodeId = nodeId;
		n = MyIO<char>::getFileLenCh(fileName);
		n++;
		unsigned char *data = new unsigned char[n];
		MyIO<unsigned char>::read(data, fileName, n - 1, std::ios_base::in | std::ios_base::binary, 0);
		data[n - 1] = 0;//add sentinal
		int64 *SA = new int64[n];
		int64 *GRank = new int64[n];
		int64 *suffixGIndex = new int64[n];
		
		//cout << "n is " << n << endl;
		//printf("%s", data);

		computeSA(data,SA, GRank, suffixGIndex,n,sizeof(unsigned char),256,0);

		delete[] SA;
		delete[] GRank;
		delete[] suffixGIndex;
		//computeSA<unsigned char>(data,n,128,0);
	}


	~ChildNode()
	{
		free(bufRecv);
		free(bufSend);
	}
};
#endif