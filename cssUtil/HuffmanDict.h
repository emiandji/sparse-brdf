#pragma once

#include "defs.h"


namespace boost
{
	namespace serialization
	{
		// --------------------------------------------------------------------

		template < class Archive, typename Block, typename Allocator >
		inline void save(Archive & ar, boost::dynamic_bitset< Block, Allocator > const & t, const unsigned int /* version */)
		{
			// Serialize bitset size
			std::size_t size = t.size();
			ar << size;

			// Convert bitset into a vector
			std::vector< Block > v(t.num_blocks());
			to_block_range(t, v.begin());

			// Serialize vector
			ar & v;
		}

		// --------------------------------------------------------------------

		template < class Archive, typename Block, typename Allocator >
		inline void load(Archive & ar, boost::dynamic_bitset< Block, Allocator > & t, const unsigned int /* version */)
		{
			std::size_t size;
			ar & size;

			t.resize(size);

			// Load vector
			std::vector< Block > v;
			ar & v;

			// Convert vector into a bitset
			from_block_range(v.begin(), v.end(), t);
		}

		// --------------------------------------------------------------------

		template <class Archive, typename Block, typename Allocator>
		inline void serialize(Archive & ar, boost::dynamic_bitset<Block, Allocator> & t, const unsigned int version)
		{
			boost::serialization::split_free(ar, t, version);
		}

		// --------------------------------------------------------------------

	}	//	namespace serialization
}	//	namespace boost



using namespace std;


class ClNode
{
public:
	const float m_fProb;
	virtual ~ClNode() {}
protected:
	ClNode(float _fProb) :m_fProb(_fProb) {}

};

struct SNodeCompare
{
	bool operator()(const ClNode* _lNode, const ClNode* _rNode) const { return _lNode->m_fProb > _rNode->m_fProb; }
};

class CInternalNode : public ClNode
{
public:
	ClNode *const m_cLeft;
	ClNode *const m_cRight;

	CInternalNode(ClNode *_c0, ClNode *_c1) : ClNode(_c0->m_fProb + _c1->m_fProb), m_cLeft(_c0), m_cRight(_c1) {}
	~CInternalNode()
	{
		delete m_cLeft;
		delete m_cRight;
	}


};


template <typename T>
class CLeafNode :public ClNode
{
public:
	const T m_fC;
	CLeafNode(float _prob, T _fC) : ClNode(_prob), m_fC(_fC) {}

};


void injectLoopDyn(boost::dynamic_bitset<unsigned long>& bs1, const boost::dynamic_bitset<unsigned long>& bs2)
{

	size_t b1size = bs1.size();
	size_t b2size = bs2.size();
	bs1.resize(b1size + b2size);
	for (size_t i = 0; i < b2size; i++)
		bs1[i + b1size] = bs2[i];
};



template <typename T>
class CHuffmanDict
{
public:

	std::map<T, boost::dynamic_bitset<unsigned long> > m_mHuffmanDictMap;

	CHuffmanDict();
	CHuffmanDict(std::map<T, boost::dynamic_bitset<unsigned long> > & _mHuffmanDictMap);
	~CHuffmanDict();

	void CreateTree(std::vector<T>& _Sig, std::vector<float>& _fProb);

	void GenerateDict(const ClNode* _node, const  boost::dynamic_bitset<unsigned long>& _prefix);
	void Decode(const  boost::dynamic_bitset<unsigned long>& _encodedSig, size_t iOrigSignalSize, vector<T>& _decodedSig);
	void Encode(const std::vector<T> &_testSignal, boost::dynamic_bitset<unsigned long>& _Huffencoded);
	inline ClNode* GetRoot() { return m_cRoot; };

	bool ReadDict(const std::string _filename);
	bool WriteDict(const std::string _filename);

	bool WriteCodedSignal(boost::dynamic_bitset<unsigned long> &_encodedSig, const std::string _filename, size_t iOrigSignalSize);
	bool ReadCodedSignal(boost::dynamic_bitset<unsigned long> &_encodedSig, const std::string _filename, size_t& iOrigSignalSize);

private:

	ClNode* m_cRoot;
};

#include "HuffmanDict.inl"
