/** $lic$
 * Copyright (C) 2012-2015 by Massachusetts Institute of Technology
 * Copyright (C) 2010-2013 by The Board of Trustees of Stanford University
 *
 * This file is part of zsim.
 *
 * zsim is free software; you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, version 2.
 *
 * If you use this software in your research, we request that you reference
 * the zsim paper ("ZSim: Fast and Accurate Microarchitectural Simulation of
 * Thousand-Core Systems", Sanchez and Kozyrakis, ISCA-40, June 2013) as the
 * source of the simulator in any publications that use this software, and that
 * you send us a citation of your work.
 *
 * zsim is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CACHE_ARRAYS_H_
#define CACHE_ARRAYS_H_

#include "cache_compressor.h"
#include "memory_hierarchy.h"
#include "stats.h"
#include "hash.h"
#include "zsim.h"
#include <iostream>

typedef int64_t compressionFactor;

// Mask for the compression factor. 
// Currently accepts 4 CF values: 0, 1, 2, 3.
// Each value represents, respectively, 64, 32, 16, 8 bits in the compressed data
const uint64_t CF_MASK = 0x3;

/*
 * How to access an address and its contents:
 * std::cout << "Address  - " << std::hex << (lineAddr << lineBits) << std::dec << std::endl;
 * std::cout << "Contents - " << std::hex << *(Address*)(lineAddr << lineBits) << std::dec << std::endl;
 */

/* General interface of a cache array. The array is a fixed-size associative container that
 * translates addresses to line IDs. A line ID represents the position of the tag. The other
 * cache components store tag data in non-associative arrays indexed by line ID.
 */
class CacheArray : public GlobAlloc {
    public:
        /* Returns tag's ID if present, -1 otherwise. If updateReplacement is set, call the replacement policy's update() on the line accessed*/
        virtual int32_t lookup(const Address lineAddr, const MemReq* req, bool updateReplacement) = 0;

        /* Runs replacement scheme, returns tag ID of new pos, address of line to write back, and tag ID of lines to evict*/
        virtual uint32_t preinsert(const Address lineAddr, const MemReq* req, std::vector<Address>* wbLineAddr, std::vector<int32_t>* evictId) = 0;

        /* Actually do the replacement, writing the new address in lineId.
         * NOTE: This method is guaranteed to be called after preinsert, although
         * there may be some intervening calls to lookup. The implementation is
         * allowed to keep internal state in preinsert() and use it in postinsert()
         */
        virtual void postinsert(const Address lineAddr, const MemReq* req, uint32_t lineId) = 0;

        virtual void initStats(AggregateStat* parent) {}

        /* Checks if the lineAddr matches the array entry at id  */
        virtual bool matchArray(const Address lineAddr, uint32_t id) { return false; };
};

class ReplPolicy;
class HashFamily;
class CPack;

/* Set-associative cache array */
class SetAssocArray : public CacheArray {
    protected:
        Address* array;
        ReplPolicy* rp;
        HashFamily* hf;
        uint32_t numLines;
        uint32_t numSets;
        uint32_t assoc;
        uint32_t setMask;

    public:
        SetAssocArray(uint32_t _numLines, uint32_t _assoc, ReplPolicy* _rp, HashFamily* _hf);

        int32_t lookup(const Address lineAddr, const MemReq* req, bool updateReplacement);
        uint32_t preinsert(const Address lineAddr, const MemReq* req, std::vector<Address>* wbLineAddr, std::vector<int32_t>* evictId);
        void postinsert(const Address lineAddr, const MemReq* req, uint32_t candidate);

        bool matchArray(const Address lineAddr, uint32_t id);
};

/* Skewed Compressed Cache array 
 * As we do not keep a data array, data insertion and decompression is not done.
 */
class SCCArray : public CacheArray {
    protected:
        Address* array; // This holds the superblock tags, not the usual line tags
        ReplPolicy* rp;
        HashFamily* hf;
        uint32_t numLines;
        uint32_t numSets;
        uint32_t assoc;
        uint32_t setMask;
        uint32_t wayMask;
        uint32_t superBlockSize;
        uint32_t superBlockBits;

        // Keeps the super-block validation tags. Up to 64 valid bits per super-block.
        uint64_t* validTag;

        // Cache compressor
        CacheCompressor* compressor;

        // On a lookup the accessing block's compressibility is not known,
        // so this function calculates the compression factor based on the
        // address of the block and the way being accessed:
        // {CF1,CF0} = {A_(lineBits+superBlockBits+1),A_(lineBits+superBlockBits)} ^ {W1,W0}
        compressionFactor calcCF(const Address superBlockAddr, uint32_t way){
            return (superBlockAddr ^ way) & CF_MASK;
        }

        // Calculate way group using the following formula:
        // {W1,W0} = {A_(lineBits+superBlockBits+1),A_(lineBits+superBlockBits)} ^ {CF1,CF0}
        uint32_t calcWay(const Address superBlockAddr, compressionFactor cf){
            return (superBlockAddr ^ cf) & wayMask;
        }

        // Calculate set index based on the following equation:
        // Index = / h({A47-A11,A8,A7,A6}) , CF = 0
        //         | h({A47-A11,A8,A7})    , CF = 1
        //         | h({A47-A11,A8})       , CF = 2
        //         L h({A47-A11})          , CF = 3
        // The hash function used is different for each way
        uint32_t calcSetIndex(const Address lineAddr, compressionFactor cf, uint32_t way){
            return hf->hash(way, ((lineAddr & (0x00000FFFFFFFF800 >> lineBits)) >> (superBlockSize + cf - lineBits)) | 
                                 ((lineAddr >> cf) & ((superBlockSize-1) >> cf))) & setMask;
        }

        // Calculate valid bit offset
        uint64_t calcBitOffset(const uint64_t lineAddr, int cf){
            return (superBlockSize-1) - ((lineAddr & ((superBlockSize-1) >> (CF_MASK-cf))) << (CF_MASK - cf));
        }

        // Shift bit offset in order to make a mask for the tag array
        uint64_t calcTagMask(const uint64_t lineAddr, int cf){
            return 0x1 << calcBitOffset(lineAddr, cf);
        }

        // Calculate lineAddr given a super-block address, and the data entry's set, way and byte offset.
        Address calcAddress(const Address superBlockAddr, uint32_t setIndex, uint32_t way, uint32_t offset);
        // Compress cache line to calculate compression factor
        compressionFactor getLineCF(const Address lineAddr);

        // Map an array id to a candidate
        inline uint32_t arrayIdToCand(uint32_t id, uint32_t offset);

        // Map a candidate to its corresponding array entry
        inline uint32_t candToArrayId(uint32_t candidate);

        inline Address toSuperBlockAddr(const Address lineAddr);

    public:
        SCCArray(uint32_t _numLines, uint32_t _assoc, uint32_t _superBlockSize, ReplPolicy* _rp, HashFamily* _hf, CacheCompressor* _compressor);

        virtual int32_t lookup(const Address lineAddr, const MemReq* req, bool updateReplacement);
        uint32_t preinsert(const Address lineAddr, const MemReq* req, std::vector<Address>* wbLineAddr, std::vector<int32_t>* evictId);
        void postinsert(const Address lineAddr, const MemReq* req, uint32_t candidate);

        bool matchArray(const Address lineAddr, uint32_t id);
};

/* Dynamically Skewed Compressed Cache array 
 * As we do not keep a data array, data insertion and decompression is not done.
 */
// TODO Only works for LRU and random replacement policies
class DSCCArray : public SCCArray {
    protected:
        Address* array; // This holds the superblock tags, not the usual line tags
        ReplPolicy* rp;
        HashFamily* hf;
        uint32_t numLines;
        uint32_t numSets;
        uint32_t assoc;
        uint32_t setMask;
        uint32_t wayMask;
        uint32_t superBlockSize;
        uint32_t superBlockBits;

        // Keeps the super-block validation tags. Up to 64 valid bits per super-block.
        uint64_t* validTag;

        // Cache compressor
        CacheCompressor* compressor;

    public:
        DSCCArray(uint32_t _numLines, uint32_t _assoc, uint32_t _superBlockSize, ReplPolicy* _rp, HashFamily* _hf, CacheCompressor* _compressor);

        int32_t lookup(const Address lineAddr, const MemReq* req, bool updateReplacement);
        uint32_t preinsert(const Address lineAddr, const MemReq* req, std::vector<Address>* wbLineAddr, std::vector<int32_t>* evictId);
        void postinsert(const Address lineAddr, const MemReq* req, uint32_t candidate);

        bool matchArray(const Address lineAddr, uint32_t id);
};

/* The cache array that started this simulator :) */
class ZArray : public CacheArray {
    private:
        Address* array; //maps line id to address
        uint32_t* lookupArray; //maps physical position to lineId
        ReplPolicy* rp;
        HashFamily* hf;
        uint32_t numLines;
        uint32_t numSets;
        uint32_t ways;
        uint32_t cands;
        uint32_t setMask;

        //preinsert() stores the swaps that must be done here, postinsert() does the swaps
        uint32_t* swapArray; //contains physical positions
        uint32_t swapArrayLen; //set in preinsert()

        uint32_t lastCandIdx;

        Counter statSwaps;

    public:
        ZArray(uint32_t _numLines, uint32_t _ways, uint32_t _candidates, ReplPolicy* _rp, HashFamily* _hf);

        int32_t lookup(const Address lineAddr, const MemReq* req, bool updateReplacement);
        uint32_t preinsert(const Address lineAddr, const MemReq* req, std::vector<Address>* wbLineAddr, std::vector<int32_t>* evictId);
        void postinsert(const Address lineAddr, const MemReq* req, uint32_t candidate);

        //zcache-specific, since timing code needs to know the number of swaps, and these depend on idx
        //Should be called after preinsert(). Allows intervening lookups
        uint32_t getLastCandIdx() const {return lastCandIdx;}

        void initStats(AggregateStat* parentStat);

        bool matchArray(const Address lineAddr, uint32_t id);
};

// Simple wrapper classes and iterators for candidates in each case; simplifies replacement policy interface without sacrificing performance
// NOTE: All must implement the same interface and be POD (we pass them by value)
struct SetAssocCands {
    struct iterator {
        uint32_t x;
        explicit inline iterator(uint32_t _x) : x(_x) {}
        inline void inc() {x++;} //overloading prefix/postfix too messy
        inline uint32_t operator*() const { return x; }
        inline bool operator==(const iterator& it) const { return it.x == x; }
        inline bool operator!=(const iterator& it) const { return it.x != x; }
    };

    uint32_t b, e;
    inline SetAssocCands(uint32_t _b, uint32_t _e) : b(_b), e(_e) {}
    inline iterator begin() const {return iterator(b);}
    inline iterator end() const {return iterator(e);}
    inline uint32_t numCands() const { return e-b; }
};


struct ZWalkInfo {
    uint32_t pos;
    uint32_t lineId;
    int32_t parentIdx;

    inline void set(uint32_t p, uint32_t i, int32_t x) {pos = p; lineId = i; parentIdx = x;}
};

struct ZCands {
    struct iterator {
        ZWalkInfo* x;
        explicit inline iterator(ZWalkInfo* _x) : x(_x) {}
        inline void inc() {x++;} //overloading prefix/postfix too messy
        inline uint32_t operator*() const { return x->lineId; }
        inline bool operator==(const iterator& it) const { return it.x == x; }
        inline bool operator!=(const iterator& it) const { return it.x != x; }
    };

    ZWalkInfo* b;
    ZWalkInfo* e;
    inline ZCands(ZWalkInfo* _b, ZWalkInfo* _e) : b(_b), e(_e) {}
    inline iterator begin() const {return iterator(b);}
    inline iterator end() const {return iterator(e);}
    inline uint32_t numCands() const { return e-b; }
};

struct DSCCWalkInfo {
    compressionFactor cf;
    uint32_t way;
    uint32_t setIndex;
    uint32_t offset;
    uint32_t lineId;

    inline void set(compressionFactor c, uint32_t w, uint32_t s, uint32_t o, uint32_t id){cf = c; way = w; setIndex = s; offset = o; lineId = id;}
};

struct DSCCCands {
    struct iterator {
        DSCCWalkInfo* x;
        explicit inline iterator(DSCCWalkInfo* _x) : x(_x) {}
        inline void inc() {x++;} //overloading prefix/postfix too messy
        inline uint32_t operator*() const { return x->lineId; }
        inline bool operator==(const iterator& it) const { return it.x == x; }
        inline bool operator!=(const iterator& it) const { return it.x != x; }
    };

    DSCCWalkInfo* b;
    DSCCWalkInfo* e;
    inline DSCCCands(DSCCWalkInfo* _b, DSCCWalkInfo* _e) : b(_b), e(_e) {}
    inline iterator begin() const {return iterator(b);}
    inline iterator end() const {return iterator(e);}
    inline uint32_t numCands() const { return e-b; }
};

#endif  // CACHE_ARRAYS_H_
