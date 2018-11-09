#ifndef CPACK_HPP
#define CPACK_HPP

#include "galloc.h"
#include <cstdint>
#include <vector>

#define DICTIONARY_WIDTH 4
#define DICTIONARY_SIZE (1<<4)

class CacheCompressor : public GlobAlloc {
    public:
        virtual ~CacheCompressor() {};

        // After compression is done, the compressed data is returned, along with its size in bits
        virtual std::size_t compress(uint64_t* cacheLine, uint64_t* compWords) = 0;
        // The decompression step returns the uncompressed data
        virtual void decompress(uint64_t* compWords, uint64_t* cacheLine) = 0;
};

class CPack : public CacheCompressor {
    private:
        const uint32_t code[6]   = {0x0, 0x1, 0x2, 0xC, 0xD, 0xE};
        const uint32_t length[6] = {2, 34, 6, 24, 12, 16};
        enum Pattern {ZZZZ, XXXX, MMMM, MMXX, ZZZX, MMMX};

        // The dictionary
        std::vector<uint32_t> dictionary;
        // Dictionary size
        std::size_t dictionarySize;
        // Number of entries in the dictionary
        std::size_t numEntries;

        // Returns size of compressed word
        std::size_t compressWord(uint32_t data, uint64_t& compWord);
        // Returns decompressed word
        uint32_t decompressWord(uint64_t compData);

        void resetDictionary();

    public:
        CPack(std::size_t wordsPerCacheLine);

        std::size_t compress(uint64_t* cacheLine, uint64_t* compWords);
        void decompress(uint64_t* compWords, uint64_t* cacheLine);

        void operator=(const CPack &c);
};

#endif
