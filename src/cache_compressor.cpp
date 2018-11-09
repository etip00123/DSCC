#include "cache_compressor.h"
#include <iostream>

CPack::CPack(std::size_t wordsPerCacheLine){
    this->dictionarySize = wordsPerCacheLine;
    dictionary.resize(dictionarySize);

    resetDictionary();
}

void CPack::operator=(const CPack &c){ 
    dictionarySize = c.dictionarySize;
    numEntries = c.numEntries;

    dictionary = c.dictionary;
}

void CPack::resetDictionary(){
    numEntries = 0;

    for (std::size_t i = 0; i < dictionarySize; i++) dictionary[i] = 0;
}

std::size_t CPack::compressWord(uint32_t data, uint64_t& compWord){
    std::size_t i;
    std::size_t matchPos  = 0;    // Position of best dictionary match  
    std::size_t numBMatch = 0;    // Number of bytes matched so far

    // Pattern matching (zzzz and zzzx)
    if ((data & 0xFFFFFF00) == 0x00000000){
        // Matched zzzx
        if (data & 0x000000FF){
            compWord = (code[ZZZX]<<8) | (data & 0x000000FF);
            return length[ZZZX];
        // Matched zzzz
        } else {
            compWord = code[ZZZZ];
            return length[ZZZZ];
        }
    }

    // Dictionary matching
    for (i = 0; i < numEntries; i++){
        // Check if matched mmmm
        if (dictionary[i] == data){
            matchPos = i;
            numBMatch = 4;
            break;
        }

        if (numBMatch < 3){
            // Check if matched mmmx
            if ((dictionary[i] & 0xFFFFFF00) == (data & 0xFFFFFF00)){
                matchPos = i;
                numBMatch = 3;
            // Check if matched mmxx
            } else if ((numBMatch < 2) && ((dictionary[i] & 0xFFFF0000) == (data & 0xFFFF0000))){
                matchPos = i;
                numBMatch = 2;
            }
        }
    }

    // Push into dictionary
    if (numEntries < dictionarySize) dictionary[numEntries++] = data;

    // Generate output
    switch (numBMatch){
        case 4:
            compWord = (code[MMMM] << DICTIONARY_WIDTH) | matchPos;
            return length[MMMM];
        case 3:
            compWord = (((code[MMMX] << DICTIONARY_WIDTH) | matchPos) << 8) | (data & 0x000000FF);
            return length[MMMX];
        case 2:
            compWord = (((code[MMXX] << DICTIONARY_WIDTH) | matchPos) << 16) | (data & 0x0000FFFF);
            return length[MMXX];
        default:
            compWord = ((uint64_t)code[XXXX] << 32) | data;
            return length[XXXX];
    }
}

std::size_t CPack::compress(uint64_t* cacheLine, uint64_t* compWords){
//std::cout << "Debug 0" << std::endl;
    std::size_t compSize = 0;

//std::cout << "Debug 1" << std::endl;
    // Reset dictionary
    resetDictionary();

//std::cout << "Debug 2" << std::endl;
    // Compress every word sequentially
    for (std::size_t i = 0; i < dictionarySize/2; i++) {
        uint32_t firstWord = (cacheLine[i]&0xFFFFFFFF00000000) >> 32;
//std::cout << "Debug 3" << std::endl;
        uint32_t secondWord = cacheLine[i]&0x00000000FFFFFFFF;
//std::cout << "Debug 4" << std::endl;
        compSize += compressWord(firstWord, compWords[2*i]);
//std::cout << "Debug 5" << std::endl;
        compSize += compressWord(secondWord, compWords[2*i+1]);
    }
//std::cout << "Debug 6" << std::endl;

    // Return compressed block size
    return compSize;
}

uint32_t CPack::decompressWord(uint64_t compData){
    // Check if matches xxxx
    if (compData & 0x0000000100000000) {
        // Rebuild dictionary
        dictionary[numEntries++] = compData;
        return (uint32_t) compData;
    // Check if matches mmxx
    } else if (compData & 0xF00000) return (dictionary[(compData & 0xF0000) >> 16] & 0xFFFF0000) | (compData & 0x0000FFFF); 
    // Check if matches mmmx
    else if (compData & 0xF000) return (dictionary[(compData & 0xF00) >> 8] & 0xFFFFFF00) | (compData & 0x000000FF); 
    // Check if matches zzzx
    else if (compData & 0xF00) return compData & 0xFF;
    // Check if matches mmmm
    else if (compData & 0xF0) return dictionary[compData & 0xF];
    // Matches zzzz
    else return 0;
}

void CPack::decompress(uint64_t* compWords, uint64_t* cacheLine){
    // Reset dictionary
    resetDictionary();

    // Decompress every word sequentially
    for (std::size_t i = 0; i < dictionarySize/2; i++) {
        cacheLine[i] = (uint64_t)decompressWord(compWords[2*i]) << 32;
        cacheLine[i] |= (uint64_t)decompressWord(compWords[2*i+1]);
    }
}
/*
int main(){
  size_t lineSize = 256;
  size_t wordsPerLine = 256/32;
  CPack cpack(wordsPerLine);
  uint64_t cacheLine[] = {0x12345678AAAAAAAA, 0x123400003527894E, 0x000000ABBBBB2022, 0x123456AABB1011D4};
  uint64_t compWords[wordsPerLine], decompDWords[wordsPerLine/2];

  // Test compression
  size_t compDataSize = cpack.compress(cacheLine, compWords);
  std::cout << "Original size: " << lineSize << std::endl << "Compressed size: " << compDataSize << std::endl << std::endl;
  for (int i = 0; i < wordsPerLine/2; i++)
    std::cout << "Compressed " << std::hex << cacheLine[i] << " to " << compWords[2*i] << " and " << compWords[2*i+1] << std::endl;

  // Test decompression
  cpack.decompress(compWords, decompDWords);
  for (int i = 0; i < wordsPerLine/2; i++)
    std::cout << "Decompressed " << std::hex << compWords[2*i] << " and " << compWords[2*i+1] << " to " << decompDWords[i] << std::endl;

  return 0;
}*/
