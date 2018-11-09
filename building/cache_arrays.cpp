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

#include "cache_arrays.h"
#include "repl_policies.h"
#include "zsim.h"

uint64_t cacheHits;
size_t otherCF;

/* Set-associative array implementation */

SetAssocArray::SetAssocArray(uint32_t _numLines, uint32_t _assoc, ReplPolicy* _rp, HashFamily* _hf) : rp(_rp), hf(_hf), numLines(_numLines), assoc(_assoc)  {
    array = gm_calloc<Address>(numLines);
    numSets = numLines/assoc;
    setMask = numSets - 1;
    assert_msg(isPow2(numSets), "must have a power of 2 # sets, but you specified %d", numSets);
}

bool SetAssocArray::matchArray(const Address lineAddr, uint32_t id){
    return lineAddr == array[id];
}

int32_t SetAssocArray::lookup(const Address lineAddr, const MemReq* req, bool updateReplacement) {
    uint32_t set = hf->hash(0, lineAddr) & setMask;
    uint32_t first = set*assoc;
    for (uint32_t id = first; id < first + assoc; id++) {
        if (array[id] == lineAddr) {
            if (updateReplacement) rp->update(id, req);
            return id;
        }
    }
    return -1;
}

uint32_t SetAssocArray::preinsert(const Address lineAddr, const MemReq* req, std::vector<Address>* wbLineAddr, std::vector<int32_t>* evictId) { //TODO: Give out valid bit of wb cand?
    uint32_t set = hf->hash(0, lineAddr) & setMask;
    uint32_t first = set*assoc;

    uint32_t candidate = rp->rankCands(req, SetAssocCands(first, first+assoc));

    wbLineAddr->push_back(array[candidate]);
    evictId->push_back(candidate);

    return candidate;
}

void SetAssocArray::postinsert(const Address lineAddr, const MemReq* req, uint32_t candidate) {
    rp->replaced(candidate);
    array[candidate] = lineAddr;
    rp->update(candidate, req);
}

/* SCC implementation */

SCCArray::SCCArray(uint32_t _numLines, uint32_t _assoc, uint32_t _superBlockSize, ReplPolicy* _rp, HashFamily* _hf, CacheCompressor* _compressor) : rp(_rp), hf(_hf), numLines(_numLines), assoc(_assoc), superBlockSize(_superBlockSize), compressor(_compressor)  {
    assert_msg(isPow2(superBlockSize), "must have a power of 2 # super-block size, but you specified %d", superBlockSize);

    // Number of lines is superblocksize times bigger due to a workaround in init.cpp, as ZSim doesn't handle compressed caches
    numLines /= superBlockSize;
    superBlockBits = __builtin_ctz(superBlockSize);
    lineBits = ilog2(zinfo->lineSize);

    array    = gm_calloc<Address>(numLines);
    validTag = gm_calloc<uint64_t>(numLines);
    numSets  = numLines/assoc;
    setMask  = numSets - 1;
    wayMask  = assoc - 1;

    assert_msg(isPow2(numSets), "must have a power of 2 # sets, but you specified %d", numSets);
    assert_msg(isPow2(assoc), "must have a power of 2 # ways, but you specified %d", assoc);
}

inline Address SCCArray::toSuperBlockAddr(const Address lineAddr){
    return lineAddr >> superBlockBits;
}

inline uint32_t SCCArray::candToArrayId(uint32_t candidate){
    return candidate/superBlockSize;
}

inline uint32_t SCCArray::arrayIdToCand(uint32_t id, uint32_t offset){
    return id*superBlockSize + offset;
}

bool SCCArray::matchArray(const Address lineAddr, uint32_t id){
    return toSuperBlockAddr(lineAddr) == array[candToArrayId(id)];
}

Address SCCArray::calcAddress(const Address superBlockAddr, uint32_t setIndex, uint32_t way, uint32_t offset) {
    // Calc CF
    compressionFactor cf = calcCF(superBlockAddr, way);

    // Only iterate through possible addresses
    for (Address lineAddr = (superBlockAddr << superBlockBits) + (superBlockSize - 1) - (offset/(superBlockSize/(1<<cf))); lineAddr >= (superBlockAddr << superBlockBits); lineAddr -= 1<<cf)
        if ((offset == calcBitOffset(lineAddr, cf)) && (setIndex == calcSetIndex(lineAddr, cf, way))) return lineAddr;

    // If address hasn't been found, something is wrong
    panic("Address could not be calculated!");

    return 0;
}

compressionFactor SCCArray::getLineCF(const Address lineAddr){
    // Minimum compressed size is 8 (and thus the hard coded 3s)
    uint64_t compWords[2<<(lineBits-3)];
    uint64_t cacheLine[1<<(lineBits-3)];

    // Get data to be compressed
    for (int i = 0; i < 1<<(lineBits-3); i++) cacheLine[i] = *(Address*)((lineAddr<<lineBits)+(i<<3));

    // Apply compression
    int compressedSize = compressor->compress(cacheLine, compWords)/(1<<(lineBits-3));

    // Check if compression is not good (compressed size is bigger than uncompressed)
    if (compressedSize > (1<<lineBits)) compressedSize = 1<<lineBits;
    // Check if needs to round up compression due to the minimum compressed size
    else if (compressedSize < 8) compressedSize = 8;

    // Calculate compression factor. The compressed size is rounded up to the closest power of 2
    return CF_MASK - ((!isPow2(compressedSize) ? 32-__builtin_clz(compressedSize) : __builtin_ctz(compressedSize)) - 3);
}

int32_t SCCArray::lookup(const Address lineAddr, const MemReq* req, bool updateReplacement) {
    // A super-block contains superBlockSize 64-bit blocks
    Address superBlockAddr = toSuperBlockAddr(lineAddr);

    // Search for the address in all ways
    for (uint32_t way = 0; way < assoc; way++){
        // Calculate compression factor
        compressionFactor cf = calcCF(superBlockAddr, way);

        // Calculate set index and respective array id
        uint32_t id = calcSetIndex(lineAddr, cf, way)*assoc + way;

        // The mask is the block's position within the super-block
        uint64_t tagMask = calcTagMask(lineAddr, cf);

        // Check if block is valid and tags match
        if (((validTag[id] & tagMask) != 0) && (array[id] == superBlockAddr)){
            cacheHits++;
            std::cout << "Hits: " << std::dec << cacheHits << std::endl;
            if (updateReplacement) rp->update(id, req);

            return arrayIdToCand(id, calcBitOffset(lineAddr, cf));
        }
    }

    return -1;
}

uint32_t SCCArray::preinsert(const Address lineAddr, const MemReq* req, std::vector<Address>* wbLineAddr, std::vector<int32_t>* evictId) { //TODO: Give out valid bit of wb cand?
    uint32_t evictStart, evictEnd;

    // A super-block contains superBlockSize 64-bit blocks
    Address superBlockAddr = toSuperBlockAddr(lineAddr);

    // Calculate compression factor
    compressionFactor cf = getLineCF(lineAddr);

    // Calculate way, set index and offset
    uint32_t way = calcWay(superBlockAddr, cf);
    uint32_t set = calcSetIndex(lineAddr, cf, way);
    uint32_t offset = calcBitOffset(lineAddr, cf);

    // Find replacement candidate
    uint32_t id = set*assoc + way;
    uint32_t candidate = rp->rankCands(req, SetAssocCands(id, id+1));

    // Check if candidate belongs to the same super block. If it doesn't, evict ALL blocks
    if (array[candidate] == superBlockAddr){
        evictStart = offset + 1 - superBlockSize/(0x1 << cf);
        evictEnd = offset;
    } else {
        evictStart = 0;
        evictEnd = superBlockSize-1;
    }

    // Get all addresses to be evicted
    for (uint32_t i = evictStart; i <= evictEnd; i++){
        if (validTag[candidate] & (0x1 << i)) {
            Address address = calcAddress(array[candidate], set, way, i);
            wbLineAddr->push_back(address);
            evictId->push_back(arrayIdToCand(candidate, i));
        }
    }

    // Return the location where the block will be placed
    return arrayIdToCand(candidate, offset);
}

void SCCArray::postinsert(const Address lineAddr, const MemReq* req, uint32_t candidate) {
    // A super-block contains superBlockSize 64-bit blocks
    Address superBlockAddr = toSuperBlockAddr(lineAddr);

    // Correct candidate as the simulator thinks the cache is superBlockSize bigger thank it really is
    candidate = candToArrayId(candidate);

    // Calculate compression factor
    compressionFactor cf = calcCF(superBlockAddr, candidate%assoc);

    // Calculate tag offset and mask
    uint64_t offset = calcBitOffset(lineAddr, cf);
    uint64_t tagMask = calcTagMask(lineAddr, cf);

    rp->replaced(candidate);

    // Update validation tag. Remove evicted blocks bits, but keep others
    if (array[candidate] == superBlockAddr){
        // Create removal mask
        // First create a mask where all bits are set if CF_0, half are set if CF_1, and so on
        uint64_t removalMask = ((0x1 << superBlockSize) - 1) >> ((superBlockSize-1) - ((superBlockSize-1) >> cf));
        // Then shift mask to its proper position
        removalMask <<= offset - ((superBlockSize-1) >> cf);
        // Last, reverse it
        removalMask = ~removalMask;
        // Inneficient method: set tagMask bit and keep shifting right until all required bits are set
        //for (int i = 1, uint64_t removalMask = tagMask; i < (superBlockSize>>cf); i++) removalMask = (removalMask >> 1) | tagMask;

        validTag[candidate] = (validTag[candidate] & removalMask) | tagMask;
    } else validTag[candidate] = tagMask;

    array[candidate] = superBlockAddr;

    rp->update(candidate, req);
}

DSCCArray::DSCCArray(uint32_t _numLines, uint32_t _assoc, uint32_t _superBlockSize, ReplPolicy* _rp, HashFamily* _hf, CacheCompressor* _compressor) : SCCArray(_numLines, _assoc, _superBlockSize, _rp, _hf, _compressor), rp(_rp), hf(_hf), numLines(_numLines), assoc(_assoc), superBlockSize(_superBlockSize), compressor(_compressor) {
    assert_msg(isPow2(superBlockSize), "must have a power of 2 # super-block size, but you specified %d", superBlockSize);

    // Number of lines is superblocksize times bigger due to a workaround in init.cpp, as ZSim doesn't handle compressed caches
    numLines /= superBlockSize;
    superBlockBits = __builtin_ctz(superBlockSize);

    array    = gm_calloc<Address>(numLines);
    validTag = gm_calloc<uint64_t>(numLines);
    numSets  = numLines/assoc;
    setMask  = numSets - 1;
    wayMask  = assoc - 1;

    rp->setCacheArray(this);

    assert_msg(isPow2(numSets), "must have a power of 2 # sets, but you specified %d", numSets);
    assert_msg(isPow2(assoc), "must have a power of 2 # ways, but you specified %d", assoc);
}

bool DSCCArray::matchArray(const Address lineAddr, uint32_t id){
    return toSuperBlockAddr(lineAddr) == array[candToArrayId(id)];
}

int32_t DSCCArray::lookup(const Address lineAddr, const MemReq* req, bool updateReplacement) {
    // A super-block contains superBlockSize 64-bit blocks
    Address superBlockAddr = toSuperBlockAddr(lineAddr);

    // Search for the address in all ways
    for (uint32_t way = 0; way < assoc; way++){
        // Calculate compression factor
        compressionFactor cf = calcCF(superBlockAddr, way);

        // Calculate set index and respective array id
        uint32_t id = calcSetIndex(lineAddr, cf, way)*assoc + way;

        // The mask is the block's position within the super-block
        uint64_t tagMask = calcTagMask(lineAddr, cf);

        // Check if block is valid and tags match
        if (((validTag[id] & tagMask) != 0) && (array[id] == superBlockAddr)){
            cacheHits++;
            std::cout << "Hits: " << std::dec << cacheHits << std::endl;
            if (updateReplacement) rp->update(id, req);

            return arrayIdToCand(id, calcBitOffset(lineAddr, cf));
        }
    }

    return -1;
}

uint32_t DSCCArray::preinsert(const Address lineAddr, const MemReq* req, std::vector<Address>* wbLineAddr, std::vector<int32_t>* evictId) { //TODO: Give out valid bit of wb cand?
    uint32_t numCandidates = 0;
    uint32_t evictStart, evictEnd;
    compressionFactor cf;
    uint32_t way = 0, set = 0, offset = 0;

    // A super-block contains superBlockSize 64-bit blocks
    Address superBlockAddr = toSuperBlockAddr(lineAddr);

    // Calculate all possible compression factors and placement locations
    // The compressed size is rounded up to the closest power of 2. Minimum compressed size is 8 (the hard coded 3)
    DSCCWalkInfo candidates[CF_MASK+1];
    bool gotCF = false;
    for (cf = getLineCF(lineAddr); cf >= 0; cf--){
        // Calculate way and set index
        way = calcWay(superBlockAddr, cf);
        set = calcSetIndex(lineAddr, cf, way);

        // Push valid values to vectors
        // If a block can be compressed, don't leave it uncompressed
        if (!gotCF or cf != 0) candidates[numCandidates++].set(cf, way, set, calcBitOffset(lineAddr, cf), set*assoc + way);

        gotCF = true;
    }

    // Find replacement candidate
    uint32_t bestCandidate = rp->rankCands(req, DSCCCands(&candidates[0], &candidates[numCandidates]));

    // Find out the chosen candidate cf, way, set and offset
    for (uint32_t i = 0; i < numCandidates; i++){
        if (candidates[i].lineId == bestCandidate){
            cf = candidates[i].cf;
            way = candidates[i].way;
            set = candidates[i].setIndex;
            offset = candidates[i].offset;

            break;
        }
    }

    // Check if candidate belongs to the same super block. If it doesn't, evict ALL blocks
    if (array[bestCandidate] == superBlockAddr){
        evictStart = offset + 1 - superBlockSize/(0x1 << cf);
        evictEnd = offset;
    } else {
        evictStart = 0;
        evictEnd = superBlockSize-1;
    }

    // Get all addresses to be evicted
    for (uint32_t i = evictStart; i <= evictEnd; i++){
        if (validTag[bestCandidate] & (0x1 << i)) {
            Address address = calcAddress(array[bestCandidate], set, way, i);
            wbLineAddr->push_back(address);
            evictId->push_back(arrayIdToCand(bestCandidate, i));
        }
    }

    // Update counter of number of CF different from max CF
    if (cf != getLineCF(lineAddr)) {
        otherCF++;
        std::cout << otherCF << std::endl;
    }

    // Return the location where the block will be placed
    return arrayIdToCand(bestCandidate, offset);
}

void DSCCArray::postinsert(const Address lineAddr, const MemReq* req, uint32_t candidate) {
    // A super-block contains superBlockSize 64-bit blocks
    Address superBlockAddr = toSuperBlockAddr(lineAddr);

    // Correct candidate as the simulator thinks the cache is superBlockSize bigger thank it really is
    candidate = candToArrayId(candidate);

    // Calculate compression factor
    compressionFactor cf = calcCF(superBlockAddr, candidate%assoc);

    // Calculate tag offset and mask
    uint64_t offset = calcBitOffset(lineAddr, cf);
    uint64_t tagMask = calcTagMask(lineAddr, cf);

    rp->replaced(candidate);

    // Update validation tag. Remove evicted blocks bits, but keep others
    if (array[candidate] == superBlockAddr){
        // Create removal mask
        // First create a mask where all bits are set if CF_0, half are set if CF_1, and so on
        uint64_t removalMask = ((0x1 << superBlockSize) - 1) >> ((superBlockSize-1) - ((superBlockSize-1) >> cf));
        // Then shift mask to its proper position
        removalMask <<= offset - ((superBlockSize-1) >> cf);
        // Last, reverse it
        removalMask = ~removalMask;
        // Inneficient method: set tagMask bit and keep shifting right until all required bits are set
        //for (int i = 1, uint64_t removalMask = tagMask; i < (superBlockSize>>cf); i++) removalMask = (removalMask >> 1) | tagMask;

        validTag[candidate] = (validTag[candidate] & removalMask) | tagMask;
    } else validTag[candidate] = tagMask;

    array[candidate] = superBlockAddr;

    rp->update(candidate, req);
}

/* ZCache implementation */
ZArray::ZArray(uint32_t _numLines, uint32_t _ways, uint32_t _candidates, ReplPolicy* _rp, HashFamily* _hf) //(int _size, int _lineSize, int _assoc, int _zassoc, ReplacementPolicy<T>* _rp, int _hashType)
    : rp(_rp), hf(_hf), numLines(_numLines), ways(_ways), cands(_candidates)
{
    assert_msg(ways > 1, "zcaches need >=2 ways to work");
    assert_msg(cands >= ways, "candidates < ways does not make sense in a zcache");
    assert_msg(numLines % ways == 0, "number of lines is not a multiple of ways");

    //Populate secondary parameters
    numSets = numLines/ways;
    assert_msg(isPow2(numSets), "must have a power of 2 # sets, but you specified %d", numSets);
    setMask = numSets - 1;

    lookupArray = gm_calloc<uint32_t>(numLines);
    array = gm_calloc<Address>(numLines);
    for (uint32_t i = 0; i < numLines; i++) {
        lookupArray[i] = i;  // start with a linear mapping; with swaps, it'll get progressively scrambled
    }
    swapArray = gm_calloc<uint32_t>(cands/ways + 2);  // conservative upper bound (tight within 2 ways)
}

bool ZArray::matchArray(const Address lineAddr, uint32_t id){
    return lineAddr == array[id];
}

void ZArray::initStats(AggregateStat* parentStat) {
    AggregateStat* objStats = new AggregateStat();
    objStats->init("array", "ZArray stats");
    statSwaps.init("swaps", "Block swaps in replacement process");
    objStats->append(&statSwaps);
    parentStat->append(objStats);
}

int32_t ZArray::lookup(const Address lineAddr, const MemReq* req, bool updateReplacement) {
    /* Be defensive: If the line is 0, panic instead of asserting. Now this can
     * only happen on a segfault in the main program, but when we move to full
     * system, phy page 0 might be used, and this will hit us in a very subtle
     * way if we don't check.
     */
    if (unlikely(!lineAddr)) panic("ZArray::lookup called with lineAddr==0 -- your app just segfaulted");

    for (uint32_t w = 0; w < ways; w++) {
        uint32_t lineId = lookupArray[w*numSets + (hf->hash(w, lineAddr) & setMask)];
        if (array[lineId] == lineAddr) {
            if (updateReplacement) {
                rp->update(lineId, req);
            }
            return lineId;
        }
    }
    return -1;
}

uint32_t ZArray::preinsert(const Address lineAddr, const MemReq* req, std::vector<Address>* wbLineAddr, std::vector<int32_t>* evictId) {
    ZWalkInfo candidates[cands + ways]; //extra ways entries to avoid checking on every expansion

    bool all_valid = true;
    uint32_t fringeStart = 0;
    uint32_t numCandidates = ways; //seeds

    //info("Replacement for incoming 0x%lx", lineAddr);

    //Seeds
    for (uint32_t w = 0; w < ways; w++) {
        uint32_t pos = w*numSets + (hf->hash(w, lineAddr) & setMask);
        uint32_t lineId = lookupArray[pos];
        candidates[w].set(pos, lineId, -1);
        all_valid &= (array[lineId] != 0);
        //info("Seed Candidate %d addr 0x%lx pos %d lineId %d", w, array[lineId], pos, lineId);
    }

    //Expand fringe in BFS fashion
    while (numCandidates < cands && all_valid) {
        uint32_t fringeId = candidates[fringeStart].lineId;
        Address fringeAddr = array[fringeId];
        assert(fringeAddr);
        for (uint32_t w = 0; w < ways; w++) {
            uint32_t hval = hf->hash(w, fringeAddr) & setMask;
            uint32_t pos = w*numSets + hval;
            uint32_t lineId = lookupArray[pos];

            // Logically, you want to do this...
#if 0
            if (lineId != fringeId) {
                //info("Candidate %d way %d addr 0x%lx pos %d lineId %d parent %d", numCandidates, w, array[lineId], pos, lineId, fringeStart);
                candidates[numCandidates++].set(pos, lineId, (int32_t)fringeStart);
                all_valid &= (array[lineId] != 0);
            }
#endif
            // But this compiles as a branch and ILP sucks (this data-dependent branch is long-latency and mispredicted often)
            // Logically though, this is just checking for whether we're revisiting ourselves, so we can eliminate the branch as follows:
            candidates[numCandidates].set(pos, lineId, (int32_t)fringeStart);
            all_valid &= (array[lineId] != 0);  // no problem, if lineId == fringeId the line's already valid, so no harm done
            numCandidates += (lineId != fringeId); // if lineId == fringeId, the cand we just wrote will be overwritten
        }
        fringeStart++;
    }

    //Get best candidate (NOTE: This could be folded in the code above, but it's messy since we can expand more than zassoc elements)
    assert(!all_valid || numCandidates >= cands);
    numCandidates = (numCandidates > cands)? cands : numCandidates;

    //info("Using %d candidates, all_valid=%d", numCandidates, all_valid);

    uint32_t bestCandidate = rp->rankCands(req, ZCands(&candidates[0], &candidates[numCandidates]));
    assert(bestCandidate < numLines);

    //Fill in swap array

    //Get the *minimum* index of cands that matches lineId. We need the minimum in case there are loops (rare, but possible)
    uint32_t minIdx = -1;
    for (uint32_t ii = 0; ii < numCandidates; ii++) {
        if (bestCandidate == candidates[ii].lineId) {
            minIdx = ii;
            break;
        }
    }
    assert(minIdx >= 0);
    //info("Best candidate is %d lineId %d", minIdx, bestCandidate);

    lastCandIdx = minIdx; //used by timing simulation code to schedule array accesses

    int32_t idx = minIdx;
    uint32_t swapIdx = 0;
    while (idx >= 0) {
        swapArray[swapIdx++] = candidates[idx].pos;
        idx = candidates[idx].parentIdx;
    }
    swapArrayLen = swapIdx;
    assert(swapArrayLen > 0);

    //Write address of line we're replacing
    wbLineAddr->push_back(array[bestCandidate]);
    evictId->push_back(bestCandidate);

    return bestCandidate;
}

void ZArray::postinsert(const Address lineAddr, const MemReq* req, uint32_t candidate) {
    //We do the swaps in lookupArray, the array stays the same
    assert(lookupArray[swapArray[0]] == candidate);
    for (uint32_t i = 0; i < swapArrayLen-1; i++) {
        //info("Moving position %d (lineId %d) <- %d (lineId %d)", swapArray[i], lookupArray[swapArray[i]], swapArray[i+1], lookupArray[swapArray[i+1]]);
        lookupArray[swapArray[i]] = lookupArray[swapArray[i+1]];
    }
    lookupArray[swapArray[swapArrayLen-1]] = candidate; //note that in preinsert() we walk the array backwards when populating swapArray, so the last elem is where the new line goes
    //info("Inserting lineId %d in position %d", candidate, swapArray[swapArrayLen-1]);

    rp->replaced(candidate);
    array[candidate] = lineAddr;
    rp->update(candidate, req);

    statSwaps.inc(swapArrayLen-1);
}

