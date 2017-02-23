#ifndef  __MYCOMMON
#define  __MYCOMMON
#include<io.h>
#include<vector>
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

constexpr int32 my_pointer_size = sizeof(void *);

template<int32 ptrSize>
struct choose_int_types {};

template<>
struct choose_int_types<4> { //32-bit
	using int_type = int32;
	using unsigned_type = uint32;
};

template<>
struct choose_int_types<8> { //64-bit
	using int_type = int64;
	using unsigned_type = uint64;
};

using int_type = choose_int_types<my_pointer_size>::int_type;
using uint_type = choose_int_types<my_pointer_size>::unsigned_type;
using ptype = uint32;

template <typename HighType>
class uint_pair
{
public:
	//! lower part type, always 32-bit
	typedef uint32 low_type;
	//! higher part type, currently either 8-bit or 16-bit
	typedef HighType high_type;

private:
	low_type low;

	high_type high;

	static uint_type low_max() {
		return std::numeric_limits<low_type>::max();
	}

	static const size_t low_bits = 8 * sizeof(low_type);

	static uint_type high_max() {
		return std::numeric_limits<high_type>::max();
	}

	static const size_t high_bits = 8 * sizeof(high_type);

public:
	static const size_t digits = low_bits + high_bits;

	static const size_t bytes = sizeof(low_type) + sizeof(high_type);

	inline uint_pair() {}

	inline uint_pair(const low_type& l, const high_type& h) : low(l), high(h) {}

	inline uint_pair(const uint_pair& a) : low(a.low), high(a.high) {}

	inline uint_pair(const uint32& a) : low(a), high(0) {}

	//! assumption: a >= 0
	inline uint_pair(const int32& a) : low(a), high(0) {}


	//!assumption: no overflow
	inline uint_pair(const uint64& a) : low((low_type)(a & low_max())), high((high_type)((a >> low_bits) & high_max())) {}

	//!assumption: a >= 0 and no overflow
	inline uint_pair(const int64& a) : low((low_type)(a & low_max())), high((high_type)((a >> low_bits) & high_max())) {}

	inline uint64 u64() const {
		return ((uint64)high) << low_bits | (uint64)low;
	}

	inline operator uint64 () const {
		return u64();
	}

	inline uint_pair& operator ++ () {
		if (low == low_max())
			++high, low = 0;
		else
			++low;
		return *this;
	}

	inline uint_pair& operator -- () {
		if (low == 0)
			--high, low = (low_type)low_max();
		else
			--low;
		return *this;
	}

	inline uint_pair& operator += (const uint_pair& b) {
		uint64 add = (uint64)low + b.low;
		low = (low_type)(add & low_max());
		high = (high_type)(high + b.high + ((add >> low_bits) & high_max()));
		return *this;
	}

	inline bool operator == (const uint_pair& b) const {
		return (low == b.low) && (high == b.high);
	}

	inline bool operator != (const uint_pair& b) const {
		return (low != b.low) || (high != b.high);
	}

	inline bool operator < (const uint_pair& b) const {
		return (high < b.high) || (high == b.high && low < b.low);
	}

	inline bool operator <= (const uint_pair& b) const {
		return (high < b.high) || (high == b.high && low <= b.low);
	}

	inline bool operator > (const uint_pair& b) const {
		return (high > b.high) || (high == b.high && low > b.low);
	}

	// !assumption: b >= 0
	inline bool operator > (const int32 & b) const {
		return (high > 0 || high == 0 && low > b);
	}
	inline bool operator >= (const uint_pair& b) const
	{
		return (high > b.high) || (high == b.high && low >= b.low);
	}

	friend std::ostream& operator << (std::ostream& os, const uint_pair& a)
	{
		return os << a.u64();
	}

	static uint_pair min()
	{
		return uint_pair(std::numeric_limits<low_type>::min(),
			std::numeric_limits<high_type>::min());
	}

	static uint_pair max()
	{
		return uint_pair(std::numeric_limits<low_type>::max(),
			std::numeric_limits<high_type>::max());
	}
};

namespace std {
	//! template class providing some numeric_limits fields for uint_pair types.
	template <typename HighType>
	class numeric_limits<uint_pair<HighType> > {
	public:
		//! yes we have information about uint_pair
		static const bool is_specialized = true;

		//! return an uint_pair instance containing the smallest value possible
		static uint_pair<HighType> min() {
			return uint_pair<HighType>::min();
		}

		//! return an uint_pair instance containing the largest value possible
		static uint_pair<HighType> max() {
			return uint_pair<HighType>::max();
		}

		//! return an uint_pair instance containing the smallest value possible
		static uint_pair<HighType> lowest() {
			return min();
		}

		//! unit_pair types are unsigned
		static const bool is_signed = false;

		//! uint_pair types are integers
		static const bool is_integer = true;

		//! unit_pair types contain exact integers
		static const bool is_exact = true;

		//! unit_pair radix is binary
		static const int radix = 2;

		//! number of binary digits (bits) in uint_pair
		static const int digits = uint_pair<HighType>::digits;

		//! epsilon is zero
		static const uint_pair<HighType> epsilon() {
			return uint_pair<HighType>(0, 0);
		}

		//! rounding error is zero
		static const uint_pair<HighType> round_error() {
			return uint_pair<HighType>(0, 0);
		}

		//! no exponent
		static const int min_exponent = 0;

		//! no exponent
		static const int min_exponent10 = 0;

		//! no exponent
		static const int max_exponent = 0;

		//! no exponent
		static const int max_exponent10 = 0;

		//! no infinity
		static const bool has_infinity = false;
	};
}

struct MyFile {
	int64 size;
	std::string name;
};

template<typename T>
class MyIO {
public:
	static void write(T * _ta, std::string _fn, uint_type _n, std::ios_base::openmode _mode, uint_type _offset);
	static void read(T * _ta, std::string _fn, uint_type _n, std::ios_base::openmode _mode, uint_type _offset);
	static uint_type getFileLenCh(std::string _fn);
	static void getFilesUnderDir(std::string _fn,std::vector<MyFile>& files);
	static bool compareFile(MyFile &file1, MyFile &file2);

	static void getFilesInOrder(std::string _fn, std::vector<MyFile>& files);
};

template<typename T>
void MyIO<T>::write(T * _ta, std::string _fn, uint_type _n, std::ios_base::openmode _mode, uint_type _offset) {
	std::ofstream fout(_fn, _mode);
	fout.seekp(_offset, std::ios_base::beg);
	fout.write((char *)_ta, _n * sizeof(T) / sizeof(char));
}

template<typename T>
void MyIO<T>::read(T * _ta, std::string _fn, uint_type _n, std::ios_base::openmode _mode, uint_type _offset) {
	std::ifstream fin(_fn, _mode);
	fin.seekg(_offset, std::ios_base::beg);
	fin.read((char *)_ta, _n * sizeof(T) / sizeof(char));
}

template<typename T>
uint_type MyIO<T>::getFileLenCh(std::string _fn) {
	std::ifstream fin(_fn);
	fin.seekg(0, std::ios_base::end);
	return fin.tellg();
}

template<typename T>
void MyIO<T>::getFilesUnderDir(std::string path, std::vector<MyFile>& files) {
	//文件句柄
	long   hFile = 0;
	//文件信息
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*").c_str(), &fileinfo)) != -1)
	{
		do
		{
			//如果是目录,迭代之
			//如果不是,加入列表
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					MyIO<T>::getFilesUnderDir(p.assign(path).append("\\").append(fileinfo.name), files);
			}
			else
			{
				//struct _finddata_t fileTemp = fileinfo;
				struct MyFile fileTemp = { fileinfo.size, p.assign(path).append("\\").append(fileinfo.name) };
				files.push_back(fileTemp);
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}


template<typename T>
bool  MyIO<T>::compareFile(MyFile &file1, MyFile &file2) {
	return file1.size < file2.size;
}

template<typename T>
void MyIO<T>::getFilesInOrder(std::string _fn, std::vector<MyFile>& files) {
	MyIO<T>::getFilesUnderDir(_fn, files);
	std::sort(files.begin(), files.end(), MyIO<T>::compareFile);
}



using uint40 = uint_pair<uint8>;

using uint48 = uint_pair<uint16>;

constexpr uint64 TBUFSIZE = 1 * 1024 * 1024 * 1024; //total buffer size in bytes
constexpr uint64 SBUFSIZE = 1 * 1024 * 1024;  //small buffer size
constexpr uint64 MBUFSIZE = 64 * 1024 * 1024;  //middle buffer size
constexpr uint64 LBUFSIZE = 512 * 1024 * 1024;  //large buffer size

constexpr int32 K1 = 1024; //main node buffer size
constexpr int32 K2 = 1024; //child node buffer size
//constexpr int64 EMPTY = 0xffffffffffffffff;
constexpr int64 CHILDCAPACITY = 1024  * 20;//200M

constexpr int32 commSize=5;
constexpr bool L_TYPE = 0;
constexpr bool S_TYPE = 1;
constexpr int32 MAINNODEID = 0;
#if _MSC_VER
#pragma pack(push, 1)
#endif

//template<typename T>
//struct CommunicateData {
//	int64 suffixGIndex, pos;
//	bool isGlobal, type;
//	T data;
//
//	CommunicateData(bool typeT, T dataT, bool isGlobalT, int64 suffixGIndexT, int64 posC) {
//		pos = posC;
//		type = typeT;
//		data = dataT;
//		isGlobal = isGlobalT;
//		suffixGIndex = suffixGIndexT;
//	}
//
//	/*CommunicateData(ChildMetaData<T> dataC) {
//	type = dataC.type;
//	data = dataC.data;
//	isGlobal = dataC.isGlobal;
//	suffixGIndex = dataC.suffixGIndex;
//	}*/
//
//	/*void setMember(ChildMataData &val) {
//	lastCharIndex = val.lastCharGIndex;
//	data = val.data;
//	globalIndex = val.globalIndex;
//	}*/
//
//	CommunicateData() {
//		suffixGIndex = 0;
//		isGlobal = 0;
//		type = 0;
//		data = 0;
//		pos = 0;
//	}
//
//	bool operator < (CommunicateData<T> &other) const
//	{
//		if (data < other.data)return true;
//		if (data > other.data)return false;
//		if (type < other.type)return true;
//		if (type > other.type)return false;
//		if (suffixGIndex < other.suffixGIndex)return true;
//		return false;
//	}
//
//	bool operator > (CommunicateData<T> &other) const
//	{
//		if (data > other.data)return true;
//		if (data < other.data)return false;
//		if (type > other.type)return true;
//		if (type < other.type)return false;
//		if (suffixGIndex > other.suffixGIndex)return true;
//		return false;
//	}
//
//	bool operator == (CommunicateData<T> &other) const
//	{
//		if (data == other.data && type == other.type && suffixGIndex == other.suffixGIndex)return true;
//		return false;
//	}
//
//	void printA() {
//		std::cout << "<" << isGlobal << "," << suffixGIndex << "," << data << ">" << std::endl;
//	}
//}

struct ReNameData {
	int64 index;
	int32 nodeId;

	ReNameData(int64 indexC, int64 nodeIdC) {
		index = indexC;
		nodeId = nodeIdC;
	}

	/*CommunicateData(ChildMetaData<T> dataC) {
	type = dataC.type;
	data = dataC.data;
	isGlobal = dataC.isGlobal;
	suffixGIndex = dataC.suffixGIndex;
	}*/

	/*void setMember(ChildMataData &val) {
	lastCharIndex = val.lastCharGIndex;
	data = val.data;
	globalIndex = val.globalIndex;
	}*/

	ReNameData() {
		index = 0;
		nodeId = 0;
	}

	bool operator < (ReNameData &other) const
	{
		return index < other.index;
	}

	bool operator > (ReNameData &other) const
	{
		return index > other.index;
	}

	bool operator == (ReNameData &other) const
	{
		return index == other.index;
	}

	void printA() {
		std::cout << "<" << index << "," << nodeId << std::endl;
	}
}

#if _MSC_VER
;
#pragma pack(pop)
#pragma pack(push, 1)
#else
__attribute__((packed));
#endif


struct RecvData {
	int64 suffixGIndex, SA;
	bool isGlobal, type;
	int64 data;
	int32 nodeId;

	RecvData(int64 SAC, int64 suffixGIndexT, bool typeT, bool isGlobalT, int64 dataT, int32 nodeIdC) {
		type = typeT;
		data = dataT;
		isGlobal = isGlobalT;
		suffixGIndex = suffixGIndexT;
		SA = SAC;
		nodeId = nodeIdC;
	}

	RecvData() {
		suffixGIndex = 0;
		isGlobal = 0;
		type = 0;
		data = 0;
		SA = 0;
		nodeId = 0;
	}

	void setMember(int64 SAC, int64 suffixGIndexT, bool typeT, bool isGlobalT, int64 dataT, int32 nodeIdC) {
		type = typeT;
		data = dataT;
		isGlobal = isGlobalT;
		suffixGIndex = suffixGIndexT;
		SA = SAC;
		nodeId = nodeIdC;
	}

	bool operator < (RecvData &other) const
	{
		if(!isGlobal || !other.isGlobal)std::cout << "不是全局变量不能比较<" << std::endl;
		if (data < other.data)return true;
		if (data > other.data)return false;
		if (type < other.type)return true;
		if (type > other.type)return false;
		if (suffixGIndex < other.suffixGIndex)return true;
		return false;
	}

	bool operator > (RecvData &other) const
	{
		if (!isGlobal || !other.isGlobal)std::cout << "不是全局变量不能比较>" << std::endl;
		if (data > other.data)return true;
		if (data < other.data)return false;
		if (type > other.type)return true;
		if (type < other.type)return false;
		if (suffixGIndex > other.suffixGIndex)return true;
		return false;
	}

	bool operator == (RecvData &other) const
	{
		if (!isGlobal || !other.isGlobal)std::cout << "不是全局变量不能比较==" << std::endl;
		if (data == other.data && type == other.type && suffixGIndex == other.suffixGIndex)return true;
		return false;
	}

	bool operator != (RecvData &other) const
	{
		if (!isGlobal || !other.isGlobal)std::cout << "不是全局变量不能比较==" << std::endl;
		if (data != other.data || type != other.type || suffixGIndex != other.suffixGIndex)return true;
		return false;
	}

	void printA() {
		std::cout << "<" << isGlobal << "," << suffixGIndex << "," << data << ">" << std::endl;
	}
}
//#if _MSC_VER
//;
//#pragma pack(pop)
//#pragma pack(push, 1)
//#else
//__attribute__((packed));
//#endif
//
//
//struct MainNodeMetaData {
//	int64 globalIndex, lastCharIndex,fromNodeId,currentIndex;
//	unsigned char data;
//
//	void setMember(ChildMataData &val,int32 nodeId) {
//		lastCharIndex = val.lastCharGIndex;
//		data = val.data;
//		globalIndex = val.globalIndex;
//		fromNodeId = nodeId;
//		currentIndex = val.currentIndex;
//	}
//
//	MainNodeMetaData() {
//		fromNodeId = 0;
//		globalIndex = 0;
//		lastCharIndex = 0;
//		data = 0;
//		currentIndex = 0;
//	}
//
//	MainNodeMetaData(int64 lastCI, unsigned char dataI) {
//		lastCharIndex = lastCI;
//		data = dataI;
//	}
//
//	bool operator < (MainNodeMetaData &other) const
//	{
//		if (data < other.data)return true;
//		if (data > other.data)return false;
//		if (lastCharIndex < other.lastCharIndex)return true;
//		return false;
//	}
//
//	bool operator > (MainNodeMetaData &other) const
//	{
//		if (data > other.data)return true;
//		if (data < other.data)return false;
//		if (lastCharIndex > other.lastCharIndex)return true;
//		return false;
//	}
//
//
//	void printA() {
//		std::cout << "<" << fromNodeId <<","<< globalIndex << "," << lastCharIndex << "," << data << ">" << std::endl;
//	}
//}

#if _MSC_VER
;
#pragma pack(pop)
#else
__attribute__((packed));
#endif

//class MPIDataType {
//public:
//	static void  constructMPIDataType(MPI_Datatype &myDataType) {
//		CommunicateData b;
//		//下面构建一个自定义变量，利用MPI函数构建，若用C语言的结构体
//		//则在Send和Recv是不能识别数据类型
//
//		MPI_Datatype old_types[2];
//		MPI_Aint indices[2];
//		//指定每个块中的变量个数,这里只有2个块，其中包含4个MPI_INT,6个MPI_FLOAT
//		int32 blocklens[2];
//		blocklens[0] = 2;
//		blocklens[1] = 1;
//		//指定原来旧的数据类型
//		old_types[0] = MPI_LONG;
//		old_types[1] = MPI_UNSIGNED_CHAR;
//		//指定每个块中变量的偏移量，需要依靠一个结构体的实例，这里是b。
//		MPI_Address(&b, &indices[0]);
//		MPI_Address(&b.data, &indices[1]);
//		indices[1] -= indices[0];
//		//std::cout << indices[1] << std::endl;
//		indices[0] = 0;
//		//看来indices[0]也可以一开始就应当赋值为0
//		//创建新数据于myDataType之中
//		MPI_Type_struct(2, blocklens, indices, old_types, &myDataType);
//		//注册新数据
//		MPI_Type_commit(&myDataType);
//	}
//};

template <typename T>
class MPIComm {
public:
	static void send(T* data,int32 size,int32 dest,int32 tag);
	static void recv(T* data, int32 size, int32 dest, int32 tag, MPI_Status &status);
};

template <typename T>
void MPIComm<T>::send(T* data, int32 size, int32 dest, int32 tag) {
	MPI_Send((char*)data, size*sizeof(T) / sizeof(char), MPI_CHAR, dest, tag, MPI_COMM_WORLD);
}

template <typename T>
void MPIComm<T>::recv(T* data, int32 size, int32 dest, int32 tag, MPI_Status &status){
	MPI_Recv((char*)data, size * sizeof(T) / sizeof(char), MPI_CHAR, dest, tag, MPI_COMM_WORLD, &status);
}

#endif