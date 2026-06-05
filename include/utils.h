#pragma once

#include <algorithm>
#include <bitset>
#include <blake3.h>
#include <cmath>
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Common/block.h>
#include <cryptoTools/Crypto/AES.h>
#include <cstddef>
#include <cstdint>
#include <emmintrin.h>
#include <libOTe/TwoChooseOne/ConfigureCode.h>
#include <stdexcept>
#include <tuple>
#include <vector>

using u32 = oc::u32;
using u64 = oc::u64;
using u8 = oc::u8;
using block = oc::block;

std::tuple<std::vector<std::vector<u64>>, std::vector<std::vector<u64>>> genInputs(u64 n, u64 d);

inline u64 low(oc::block &blk)
{
    u64 low64 = _mm_extract_epi64(blk, 0);

    return low64;
}

inline u64 high(oc::block &blk)
{
    u64 high64 = _mm_extract_epi64(blk, 1);

    return high64;
}

inline bool lsb(oc::block &blk)
{
    u8 lsb = _mm_extract_epi8(blk, 0) & 1;

    return lsb;
}

// NOTE: This function doesn't return bitlen. It returns log_alpha
inline int bitlen(int x)
{ // x & (-x)
    if (x < 1)
        return 0;
    for (int i = 0; i < 32; i++) {
        int curr = 1 << i;
        if (curr >= x)
            return i;
    }
    return 0;
}

inline uint8_t bool_to_uint8(const uint8_t *data, size_t len)
{
    if (len != 0)
        len = (len > 8 ? 8 : len);
    else
        len = 8;
    uint8_t res = 0;
    for (size_t i = 0; i < len; ++i) {
        if (data[i])
            res |= (1ULL << i);
    }
    return res;
}

inline void uint8_to_bool(uint8_t *data, uint8_t input, int length)
{
    for (int i = 0; i < length; ++i) {
        data[i] = (input & 1) == 1;
        input >>= 1;
    }
}

// Decompose the interval [start, end] using an improved method in appendix
inline std::vector<block> getIntervalPrefix(u64 start, u64 end, int shift = 0)
{
    if (start > end) {
        throw std::runtime_error("decompose improve: end should >= start");
    }

    if (start == end) {
        return { block(0, start) };
    }

    u64 interval_len = end - start + 1; // interval length
    u64 bit_width = static_cast<u64>(std::log2(interval_len)) + 1;
    u64 aligned_start = start;

    // step 1
    // find aligned_start >= start, s.t. 2^(bit_width-1) | aligned_start
    while (aligned_start <= end && (aligned_start & ((1 << (bit_width - 1)) - 1))) {
        aligned_start++;
    }

    if (aligned_start > end) {
        throw std::runtime_error("decompose improve: can't find aligned_start");
    }

    // step 2
    // right_len = end - aligned_start + 1; left_len = aligned_start - start
    u64 right_len = end - aligned_start + 1; // right side length
    u64 left_len = aligned_start - start;    // left side length

    // convert to binary representation, bitset access is from low bit to high bit
    std::bitset<64> right_bits(right_len); // right_len + left_len = 2\delta + 1
    std::bitset<64> left_bits(left_len);

    std::vector<block> results;       // result set
    u64 right_pos = aligned_start;    // move right
    u64 left_pos = aligned_start - 1; // move left

    std::vector<block> temp_results;

    // traverse from high bit to low bit
    for (u64 i = bit_width; i >= 1; i--) {
        if (right_bits[i - 1]) {
            // if right_len's i-th bit is 1
            if (i - 1 > (bit_width - 1 - shift)) {
                temp_results.push_back(block(i - 1, right_pos >> (i - 1)));
            } else {
                results.push_back(block(i - 1, right_pos >> (i - 1)));
            }
            right_pos += (1 << (i - 1));
        }
        if (left_bits[i - 1]) {
            // if left_len's i-th bit is 1
            if (i - 1 > (bit_width - 1 - shift)) {
                temp_results.push_back(block(i - 1, left_pos >> (i - 1)));
            } else {
                results.push_back(block(i - 1, left_pos >> (i - 1)));
            }
            left_pos -= (1 << (i - 1));
        }
    }

    if (shift != 0) {
        for (int i = 0; i < temp_results.size(); i++) {
            int len = high(temp_results[i]);
            u64 base = low(temp_results[i]);
            for (int j = 0; j < (1 << (len - (bit_width - 1 - shift))); j++) {
                results.push_back(block(bit_width - 1 - shift, (base << (len - (bit_width - 1 - shift))) + j));
            }
        }
    }

    return results;
}

inline u64 firstLessThan(u64 x, const std::vector<u64> &U)
{
    auto it = std::lower_bound(U.begin(), U.end(), x); // find the first element not less than x

    if (it == U.begin()) {
        // No element is less than x, return 0 or other appropriate value
        return 0;
    }

    --it; // Move to the largest element less than x
    return *it;
}

inline std::vector<block> getIntervalPrefixSet(u64 start, u64 end, std::vector<u64> U)
{
    std::vector<block> finalPrefixes;
    auto prefixes = getIntervalPrefix(start, end);
    for (auto &p : prefixes) {
        u64 len = high(p);
        u64 base = low(p);
        if (std::find(U.begin(), U.end(), len) != U.end()) {
            finalPrefixes.push_back(p);
        } else {
            u64 newLen = firstLessThan(len, U);

            for (int j = 0; j < (1 << (len - newLen)); j++) {
                finalPrefixes.push_back(block(newLen, (base << (len - newLen)) + j));
            }
        }
    }

    return finalPrefixes;
}

inline std::vector<block> getPrefix(u64 x, int maxLen)
{
    std::vector<block> res;

    for (int len = 0; len < maxLen; len++) {
        res.push_back(block(len, x >> (len)));
    }

    return res;
}

inline std::vector<block> getPrefixSet(u64 x, std::vector<u64> U)
{
    std::vector<block> res;

    for (auto &len : U) {
        res.push_back(block(len, x >> (len)));
    }

    return res;
}

inline u64 upBound(block prefix)
{
    u64 len = high(prefix);
    u64 base = low(prefix);
    u64 upper = (1 << len) - 1 + (base << len);

    return upper;
}

inline u64 lowBound(block prefix)
{
    u64 len = high(prefix);
    u64 base = low(prefix);
    u64 lower = base << len;

    return lower;
}

void inline Hash(std::vector<block> &input)
{
    auto n8 = input.size() / 8 * 8;

    // notOneBlock
    // block mask = OneBlock ^ AllOneBlock;

    auto r = input.data();

    for (u64 i = 0; i < n8; i += 8) {
        // r[0] = r[0] & mask;
        // r[1] = r[1] & mask;
        // r[2] = r[2] & mask;
        // r[3] = r[3] & mask;
        // r[4] = r[4] & mask;
        // r[5] = r[5] & mask;
        // r[6] = r[6] & mask;
        // r[7] = r[7] & mask;

        oc::mAesFixedKey.hashBlocks<8>(r, r);
        r += 8;
    }
    for (u64 i = n8; i < input.size(); i++) {
        // input[i] = input[i] & mask;
        input[i] = oc::mAesFixedKey.hashBlock(input[i]);
    }
}

inline uint64_t combination(uint64_t n, uint64_t k)
{
    if (k > n)
        return 0;
    if (k == 0 || k == n)
        return 1;

    // C(n, k) = C(n, n-k)
    if (k > n - k)
        k = n - k;

    uint64_t result = 1;
    for (uint64_t i = 1; i <= k; ++i) {
        result = result * (n - i + 1) / i;
    }

    return result;
}

inline std::vector<u64> cell(std::vector<u64> data, u64 sidelen)
{
    std::vector<u64> res(data.size(), 0);
    for (size_t i = 0; i < data.size(); ++i) {
        res[i] = data[i] / sidelen;
    }
    return res;
}

inline block blake3_hash(const std::vector<u64> &cell, u64 dim, block val)
{
    blake3_hasher hasher;
    block hash_out;
    blake3_hasher_init(&hasher);

    auto merge_val = block(high(val) + (dim << 32), low(val));
    blake3_hasher_update(&hasher, merge_val.data(), sizeof(block));
    blake3_hasher_update(&hasher, cell.data(), cell.size() * sizeof(u64));

    blake3_hasher_finalize(&hasher, hash_out.data(), 16);

    return hash_out;
}

inline block blake3_hash(const std::vector<u64> &cell, u64 dim, u64 val)
{
    blake3_hasher hasher;
    block hash_out;
    blake3_hasher_init(&hasher);

    auto merge_val = block(dim, val);
    blake3_hasher_update(&hasher, merge_val.data(), sizeof(block));
    blake3_hasher_update(&hasher, cell.data(), cell.size() * sizeof(u64));

    blake3_hasher_finalize(&hasher, hash_out.data(), 16);

    return hash_out;
}

inline block blake3_hash(const std::vector<u64> &cell, u64 hash_id)
{
    blake3_hasher hasher;
    block hash_out;
    blake3_hasher_init(&hasher);

    auto merge_val = block(0, hash_id);
    blake3_hasher_update(&hasher, merge_val.data(), sizeof(block));
    blake3_hasher_update(&hasher, cell.data(), cell.size() * sizeof(u64));

    blake3_hasher_finalize(&hasher, hash_out.data(), 16);

    return hash_out;
}

inline block blake3_hash(const std::vector<u64> &cell)
{
    blake3_hasher hasher;
    block hash_out;
    blake3_hasher_init(&hasher);

    blake3_hasher_update(&hasher, cell.data(), cell.size() * sizeof(u64));

    blake3_hasher_finalize(&hasher, hash_out.data(), 16);

    return hash_out;
}

inline std::vector<std::vector<u64>> neigh(std::vector<u64> point, u64 delta)
{
    int d = point.size();
    u64 cell_num = 1 << d;
    std::vector<u64> k0(d);
    std::vector<u64> k1(d);
    std::vector<std::vector<u64>> neighbors;

    for (int i = 0; i < d; ++i) {
        k0[i] = (point[i] - delta) / (2 * delta);
        k1[i] = (point[i] + delta) / (2 * delta);
        if (k0[i] == k1[i]) {
            throw std::runtime_error("neigh error");
        }
    }

    for (size_t i = 0; i < cell_num; i++) {
        std::vector<u64> neighbor(d);
        for (int j = 0; j < d; j++) {
            if ((i >> j) & 1) {
                neighbor[j] = k1[j];
            } else {
                neighbor[j] = k0[j];
            }
        }
        neighbors.push_back(neighbor);
    }
    return neighbors;
}

const osuCrypto::MultType type = osuCrypto::MultType::Tungsten;

struct NoHash {
    inline size_t operator()(const block &v) const
    {
        return v.get<size_t>(0);
    }
};
