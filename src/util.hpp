#ifndef UTIL_H
#define UTIL_H

#include <string>
#include <algorithm>
#include <functional>
#include <vector>
#include <sstream>

// UTILITY FUNCTIONS -------------------------------------------------
// trim taken from stack overflow
// http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
static inline std::string &ltrim(std::string &s) {
    s.erase(s.begin(),
            std::find_if(s.begin(), s.end(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::not1(std::ptr_fun<int, int>(std::isspace))).base(),
            s.end());
    return s;
}

// trim from both ends
static inline std::string &trim(std::string &s) {
    return ltrim(rtrim(s));
}


static inline std::vector<std::string> &split(const std::string &s, char delim,
                                              std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

static inline bool ends_with(std::string const &value,
                             std::string const &ending) {
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

static inline char bits_to_char(unsigned value){
    switch(value) {
    case 0:
        return 'A';
    case 1:
        return 'C';
    case 2:
        return 'G';
    case 3:
        return 'T';
    default:
        return 'N';
    }
}

static inline int char_to_bits(char c) {
    int cvalue = -1;
    switch (c) {
        case 'A':
        case 'a':
            cvalue = 0;
            break;
        case 'C':
        case 'c':
            cvalue = 1;
            break;
        case 'G':
        case 'g':
            cvalue = 2;
            break;
        case 'T':
        case 't':
            cvalue = 3;
            break;
    }
    return cvalue;
}

static inline std::string reverse_complement(std::string read){
    // add reverse complementary
    std::reverse(read.begin(), read.end());
    for(auto rb = read.begin();rb != read.end(); rb++){
        switch(*rb){
        case 'A':
            *rb = 'T';
            break;
        case 'C':
            *rb = 'G';
            break;
        case 'G':
            *rb = 'C';
            break;
        case 'T':
            *rb = 'A';
            break;
        default:
            *rb = 'N';
            break;
        }
    }
    return read;
}


template<typename RecordCompartor, typename RecordType, typename DataType,
         typename IndexType=unsigned>
void sort_by(std::vector<RecordType>& records, std::vector<DataType>& srtby){
    std::vector<IndexType> srt_idx, cur_idx, inv_cur_idx;

    // Initialize sorted index
    srt_idx.resize(records.size());
    cur_idx.resize(records.size());
    inv_cur_idx.resize(records.size());
    for(IndexType i = 0u; i < records.size(); i++)
        srt_idx[i] = cur_idx[i] = inv_cur_idx[i] = i;
    //sort the index
    RecordCompartor dcomp(srtby);
    std::sort(srt_idx.begin(), srt_idx.end(), dcomp);

    for(IndexType i = 0u; i < srt_idx.size(); i++){
        IndexType& x = srt_idx[i];
        IndexType& y = inv_cur_idx[x];
        IndexType& p = cur_idx[i];
        IndexType& q = cur_idx[y];

        std::swap(records[i], records[y]); // swap data
        std::swap(cur_idx[i], cur_idx[y]); // Swap cur_idx x and y
        std::swap(inv_cur_idx[p], inv_cur_idx[q]); // Swap for inv_cur_idx
    }
}

#endif /* UTIL_H */
