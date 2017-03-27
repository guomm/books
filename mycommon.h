#ifndef  __MYCOMMON
#define  __MYCOMMON
//#include<io.h>
#include <fstream> 
#include <iostream> 
#include<vector>
#include <string>
#include <algorithm>
#include "mpi.h"  

#ifdef	_MSC_VER
using int8 = __int8;
using uint8 = unsigned __int8;
using int16 = __int16;
using uint16 = unsigned __int16;
using int32 = __int32;
using uint32 = unsigned __int32;
using int64 = __int64;
using uint64 = unsigned __int64;
#else
using int8 = char;
using uint8 = unsigned char;
using int16 = short;
using uint16 = unsigned short;
using int32 = int;
using uint32 = unsigned int;
using int64 = long long int;
using uint64 = unsigned long long int;
#endif

using uint_type = int32;

struct MyFile {
	int64 size;
	std::string name;
};

template<typename T>
class MyIO {
public:
	static void write(T * _ta, std::string _fn, int64 _n, std::ios_base::openmode _mode, int64 _offset);
	static void read(T * _ta, std::string _fn, int64 _n, std::ios_base::openmode _mode, int64 _offset);
	static int64 getFileLenCh(std::string _fn);
	static void getFilesUnderDir(std::string _fn, std::vector<MyFile>& files);
	static bool compareFile(MyFile &file1, MyFile &file2);

	static void getFilesInOrder(std::string _fn, std::vector<MyFile>& files);
};

template<typename T>
void MyIO<T>::write(T * _ta, std::string _fn, int64 _n, std::ios_base::openmode _mode, int64 _offset) {
	std::ofstream fout(_fn, _mode);
	fout.seekp(_offset, std::ios_base::beg);
	fout.write((char *)_ta, _n * sizeof(T) / sizeof(char));
}

template<typename T>
void MyIO<T>::read(T * _ta, std::string _fn, int64 _n, std::ios_base::openmode _mode, int64 _offset) {
	std::ifstream fin(_fn, _mode);
	fin.seekg(_offset, std::ios_base::beg);
	fin.read((char *)_ta, _n * sizeof(T) / sizeof(char));
}

template<typename T>
int64 MyIO<T>::getFileLenCh(std::string _fn) {
	std::ifstream fin(_fn);
	fin.seekg(0, std::ios_base::end);
	return fin.tellg();
}

//template<typename T>
//void MyIO<T>::getFilesUnderDir(std::string path, std::vector<MyFile>& files) {
//	//文件句柄
//	long   hFile = 0;
//	//文件信息
//	struct _finddata_t fileinfo;
//	std::string p;
//	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
//	{
//		do
//		{
//			//如果是目录,迭代之
//			//如果不是,加入列表
//			if ((fileinfo.attrib &  _A_SUBDIR))
//			{
//				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
//					MyIO<T>::getFilesUnderDir(p.assign(path).append("\\").append(fileinfo.name), files);
//			}
//			else
//			{
//				//struct _finddata_t fileTemp = fileinfo;
//				struct MyFile fileTemp = { fileinfo.size, p.assign(path).append("\\").append(fileinfo.name) };
//				files.push_back(fileTemp);
//			}
//		} while (_findnext(hFile, &fileinfo) == 0);
//		_findclose(hFile);
//	}
//}
//
//
//template<typename T>
//bool  MyIO<T>::compareFile(MyFile &file1, MyFile &file2) {
//	return file1.size < file2.size;
//}
//
//template<typename T>
//void MyIO<T>::getFilesInOrder(std::string _fn, std::vector<MyFile>& files) {
//	MyIO<T>::getFilesUnderDir(_fn, files);
//	std::sort(files.begin(), files.end(), MyIO<T>::compareFile);
//}
//

constexpr int32 commSize = 4;
constexpr bool L_TYPE = 0;
constexpr bool S_TYPE = 1;
constexpr int32 MAINNODEID = 0;
unsigned char mask[] = { 0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01 };
const int64 EMPTY = 0xffffffffffffffff;
const uint_type MAXU = (std::numeric_limits<uint_type>::max)();
const uint_type MAXC = (std::numeric_limits<unsigned char>::max)();
const uint_type MINU = (std::numeric_limits<uint_type>::min)();
const uint_type MINC = (std::numeric_limits<unsigned char>::min)();

#if _MSC_VER
#pragma pack(push, 1)
#endif
template <typename T>
struct RecvData {
	uint_type suffixGIndex;
	int32 prePos;
	bool  type;
	T data;


	RecvData(int32 prePosC, uint_type suffixGIndexT, bool typeT, T dataT) {
		type = typeT;
		data = dataT;
		//isGlobal = isGlobalT;
		suffixGIndex = suffixGIndexT;
		prePos = prePosC;
	}

	RecvData() {
		suffixGIndex = 0;
		//isGlobal = 0;
		type = 0;
		data = 0;
		prePos = 0;
	}

	void setMember(int64 prePosC, int64 suffixGIndexT, bool typeT, T dataT) {
		type = typeT;
		data = dataT;
		//isGlobal = isGlobalT;
		suffixGIndex = suffixGIndexT;
		prePos = prePosC;
	}

	void setMember(const RecvData<T> & other) {
		type = other.type;
		data = other.data;
		//isGlobal = isGlobalT;
		suffixGIndex = other.suffixGIndex;
		prePos = other.prePos;
	}

	bool operator<(const RecvData<T> &other) const
	{
		//if(!isGlobal || !other.isGlobal)std::cout << "不是全局变量不能比较<" << std::endl;
		//if (this == nullptr) return true;
		if (data < other.data)return true;
		if (data > other.data)return false;
		if (type < other.type)return true;
		if (type > other.type)return false;
		if (suffixGIndex < other.suffixGIndex)return true;
		return false;
	}

	bool operator>(const RecvData<T> &other) const
	{
		//if (this == nullptr) return true;
		//if (!isGlobal || !other.isGlobal)std::cout << "不是全局变量不能比较>" << std::endl;
		if (data > other.data)return true;
		if (data < other.data)return false;
		if (type > other.type)return true;
		if (type < other.type)return false;
		if (suffixGIndex > other.suffixGIndex)return true;
		return false;
	}

	bool operator==(const RecvData<T> &other) const
	{
		//if (this == nullptr) return true;
		//if (!isGlobal || !other.isGlobal)std::cout << "不是全局变量不能比较==" << std::endl;
		if ((data == other.data) && (type == other.type) && (suffixGIndex == other.suffixGIndex))return true;
		return false;
	}

	bool operator!=(const RecvData<T> &other) const
	{
		//if (this == nullptr) return true;
		//if (!isGlobal || !other.isGlobal)std::cout << "不是全局变量不能比较==" << std::endl;
		if ((data != other.data) || (type != other.type) || (suffixGIndex != other.suffixGIndex))return true;
		return false;
	}

	void printA() {
		std::cout << "<" << (int32)data << "," << prePos << "," << type << "," << suffixGIndex << ">" << std::endl;
	}
}



#if _MSC_VER
;
#pragma pack(pop)
#pragma pack(push, 1)
#else
__attribute__((packed));
#endif


template <typename T>
struct SendData {
	uint_type suffixGIndex;
	int32 prePos;
	bool  type;
	T data;


	SendData(int64 prePosC, int64 suffixGIndexT, bool typeT, T dataT) {
		type = typeT;
		data = dataT;
		//isGlobal = isGlobalT;
		suffixGIndex = suffixGIndexT;
		prePos = prePosC;
	}

	SendData() {
		suffixGIndex = 0;
		type = 0;
		data = 0;
		prePos = 0;
	}

	void setMember(int64 prePosC, int64 suffixGIndexT, bool typeT, T dataT) {
		type = typeT;
		data = dataT;
		suffixGIndex = suffixGIndexT;
		prePos = prePosC;
	}

	void setMember(const SendData<T> & other) {
		type = other.type;
		data = other.data;
		suffixGIndex = other.suffixGIndex;
		prePos = other.prePos;
	}


	void printA() {
		std::cout << "<" << data << "," << prePos << "," << type << "," << suffixGIndex << ">" << std::endl;
	}
}

#if _MSC_VER
;
#pragma pack(pop)
#else
__attribute__((packed));
#endif


template <typename T>
class MPIComm {
public:
	static void send(T* data, int32 size, int32 dest, int32 tag);
	static void recv(T* data, int32 size, int32 dest, int32 tag, MPI_Status &status);
};

template <typename T>
void MPIComm<T>::send(T* data, int32 size, int32 dest, int32 tag) {
	MPI_Send((char*)data, size*sizeof(T) / sizeof(char), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
}

template <typename T>
void MPIComm<T>::recv(T* data, int32 size, int32 dest, int32 tag, MPI_Status &status) {
	MPI_Recv((char*)data, size * sizeof(T) / sizeof(char), MPI_CHAR, dest, tag, MPI_COMM_WORLD, &status);
}

struct MarkData {
public:
	int64 isCycle;
	int64 isSelfCom;

	MarkData() {
		isCycle = 0;
		isSelfCom = 0;
	}

	void print() {
		std::cout <<"<"<< isCycle << "," << isSelfCom <<">"<< std::endl;
	}
};

template <typename T>
class LoserTree {
public:
	T *myData = nullptr;
	int32 *ls = nullptr;
	int32 k;
	int32 finishNum;
	T MaxData;
	T MinData;
	LoserTree() {
		//MaxData = MaxData1;
		//MinData = MinData1;
		//MaxData.setMember(0, MAXU, 1, (std::numeric_limits<T>::max)());
		//MinData.setMember(0, MINU, 0, (std::numeric_limits<T>::min)());
	}

	void setMAXMIN(T MinData1, T MaxData1) {
		MaxData = MaxData1;
		MinData = MinData1;
	}
	void initTree(int32 numOfChild) {
		k = numOfChild;
		finishNum = 0;
		if (myData == nullptr) {
			myData = new T[numOfChild + 1];
			/*for (int32 i = numOfChild ; i >= 0 ; i--) {
			T temp;
			myData[i] = temp;
			}*/
		}
		myData[numOfChild] = MinData;
		if (ls == nullptr) ls = new int32[numOfChild];

		for (int32 i = 0; i < numOfChild; i++) {
			ls[i] = numOfChild;
		}

	}

	void setVal(int32 i, const T &val) {
		myData[i] = val;
	}

	void setAndFresh(int32 i, const T &val) {
		setVal(i, val);
		adjust(i);
	}

	void setPosFinsih(int32 i) {
		finishNum++;
		//std::cout << "finishNum is " << finishNum<<std::endl;
		//MaxData.printA();
		setVal(i, MaxData);
		adjust(i);
		//myData[i].printA();
		//cout << "目前容器中的值是：" << endl;
		//for (int32 j = 0; j < k;j++)myData[j].printA();
	}

	void adjust(int32 s)
	{
		int32 t = (s + k) / 2;
		int32 temp = 0;
		while (t>0)
		{
			if (myData[s] > myData[ls[t]])
			{
				temp = s;
				s = ls[t];
				ls[t] = temp;
			}
			t = t / 2;
		}
		ls[0] = s;
	}

	void CreateLoserTree()
	{
		int32 i;
		for (i = k - 1; i >= 0; i--)adjust(i);
	}

	bool isFinish() {
		return finishNum == k;
	}
	int32 getMinIndx() {
		return ls[0];
	}
};

template <typename T>
class VictoryTree {
public:
	RecvData<T> *myData = nullptr;
	int32 *ls = nullptr;
	int32 k;
	int32 finishNum;
	RecvData<T> MinData;
	RecvData<T> MaxData;
	VictoryTree() {
		MinData.setMember(0, MINU, 0, (std::numeric_limits<T>::min)());
		MaxData.setMember(0, MAXU, 1, (std::numeric_limits<T>::max)());
	}

	void initTree(int32 numOfChild) {
		k = numOfChild;
		finishNum = 0;
		if (myData == nullptr) {
			myData = new RecvData<T>[numOfChild + 1];
			for (int32 i = numOfChild; i >= 0; i--) {
				RecvData<T> temp;
				myData[i] = temp;
			}
		}
		myData[numOfChild].setMember(MaxData);
		if (ls == nullptr) ls = new int32[numOfChild];

		for (int32 i = 0; i < numOfChild; i++) {
			ls[i] = numOfChild;
		}

	}

	void setVal(int32 i, const RecvData<T> &val) {
		myData[i].setMember(val);
	}

	void setAndFresh(int32 i, const RecvData<T> &val) {
		setVal(i, val);
		adjust(i);
	}

	void setPosFinsih(int32 i) {
		finishNum++;
		setVal(i, MinData);
		adjust(i);
	}

	void adjust(int32 s)
	{
		int32 t = (s + k) / 2;
		int32 temp = 0;
		while (t>0)
		{
			if (myData[s] < myData[ls[t]])
			{
				temp = s;
				s = ls[t];
				ls[t] = temp;
			}
			t = t / 2;
		}
		ls[0] = s;
	}

	void CreateVictoryTree()
	{
		int32 i;
		for (i = k - 1; i >= 0; i--)adjust(i);
		//cout << "目前容器中的值是：" << endl;
		//for (int32 j = 0; j < k;j++)myData[j].printA();
		//cout << "目前容器中的最大值下标值是：" << getMaxIndx()<< endl;
	}

	bool isFinish() {
		return finishNum == k;
	}
	int32 getMaxIndx() {
		return ls[0];
	}
};
#endif