#include "cmp.h"
#include <coproto/Common/macoro.h>
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Common/Defines.h>
#include <libOTe/Tools/Coproto.h>
#include <libOTe/Triple/SilentOtTriple/SilentOtTriple.h>
#include <vector>
#include <volePSI/Defines.h>

// void CmpSender::compare(std::vector<u64> &data, std::vector<u8> &out, Socket &chl)
// {
//     coproto::sync_wait(triple->genBaseOts(*prng, chl));

//     int numDigits = mBitlen / M;

//     std::vector<u8> leaf_messages_cmp(mNum * numDigits * (1 << M));

//     std::vector<u8> leaf_res_cmp(mNum * numDigits);

//     std::vector<u8> leaf_messages_eq(mNum * numDigits * (1 << M));

//     std::vector<u8> leaf_res_eq(mNum * numDigits);

//     for (u64 i = 0; i < mNum; ++i) {
//         u64 val = data[i];
//         for (int d = 0; d < numDigits; ++d) {
//             u64 digit = (val >> (d * M)) & ((1 << M) - 1);
//             for (u64 k = 0; k < (1 << M); ++k) {
//                 leaf_messages_cmp[i * numDigits * (1 << M) + d * (1 << M) + k] = leaf_res_cmp[i * numDigits + d] ^ ((k < digit) ? 1 : 0);
//                 leaf_messages_eq[i * numDigits * (1 << M) + d * (1 << M) + k] = leaf_res_eq[i * numDigits + d] ^ ((k == digit) ? 1 : 0);
//             }
//         }
//     }

//     mOt->send(leaf_messages_cmp, chl);

//     mOt->send(leaf_messages_eq, chl);

//     auto numAND = (numDigits - 1) * mNum;
//     std::vector<block> ai(numAND);
//     std::vector<block> bi(numAND);
//     std::vector<block> ci(numAND);

//     coproto::sync_wait(triple->expand(ai, bi, ci, *prng, chl));

//     u8 *ei = new u8[numAND];
//     u8 *fi = new u8[numAND];
//     u8 *e = new u8[mNum];
//     u8 *f = new u8[mNum];

// }

void NcoOTSender::send(u8 **messages, Socket &chl)
{
    // coproto::sync_wait(mOt->genSilentBaseOts(*prng, chl));

    std::vector<std::array<block, 2>> otMessages(mNum * 4);

    coproto::sync_wait(mOt->send(otMessages, *prng, chl));

    BitVector correctMessages0(mNum * 16);
    BitVector correctMessages1(mNum * 16);

    for (u64 i = 0; i < mNum; ++i) {
        for (int j = 0; j < 16; ++j) {
            correctMessages0[i * 16 + j] = messages[i][j] & 1;
            correctMessages1[i * 16 + j] = (messages[i][j] >> 1) & 1;
            for (int k = 0; k < 4; ++k) {
                correctMessages0[i * 16 + j] ^= lsb(otMessages[i * 4 + k][(j >> k) & 1]);
                correctMessages1[i * 16 + j] ^= lsb(otMessages[i * 4 + k][(j >> k) & 1]);
            }
        }
    }

    coproto::sync_wait(chl.send(correctMessages0));
    coproto::sync_wait(chl.send(correctMessages1));
}

void NcoOTRecver::recv(u8 *outs, u8 *choices, Socket &chl)
{
    // coproto::sync_wait(mOt->genSilentBaseOts(*prng, chl));

    std::vector<block> otMessages(mNum * 4);

    BitVector choiceBits(mNum * 4);

    for (u64 i = 0; i < mNum; ++i) {
        for (int j = 0; j < 4; ++j) {
            choiceBits[i * 4 + j] = (choices[i] >> j) & 1;
        }
    }

    coproto::sync_wait(mOt->receive(choiceBits, otMessages, *prng, chl));

    BitVector correctMessages0(mNum * 16);
    BitVector correctMessages1(mNum * 16);

    coproto::sync_wait(chl.recv(correctMessages0));
    coproto::sync_wait(chl.recv(correctMessages1));

    for (u64 i = 0; i < mNum; ++i) {
        u8 out0 = correctMessages0[i * 16 + choices[i]];
        u8 out1 = correctMessages1[i * 16 + choices[i]];
        for (int j = 0; j < 4; ++j) {
            out0 ^= lsb(otMessages[i * 4 + j]);
            out1 ^= lsb(otMessages[i * 4 + j]);
        }
        outs[i] = (out1 << 1) | out0;
    }
}

MillionaireProtocolSender::MillionaireProtocolSender(int num_cmps, int bitlength, int radix_base)
{
    this->num_cmps = num_cmps;

    configure(bitlength, radix_base);

    otpack = new NcoOTSender(this->num_digits * roundUpTo(this->num_cmps, 8));
    num_triples_round = roundUpTo(this->num_triples * roundUpTo(this->num_cmps, 8), 128);
    triple_gen = new SilentOtTriple();
    triple_gen->init(1, num_triples_round);
}

void MillionaireProtocolSender::configure(int bitlength, int radix_base)
{
    assert(radix_base <= 8);
    assert(bitlength <= 64);
    this->l = bitlength;
    this->beta = radix_base;

    this->num_digits = ceil((double)l / beta);
    this->r = l % beta;
    this->log_alpha = bitlen(num_digits) - 1;
    this->log_num_digits = log_alpha + 1;
    this->num_triples_corr = 2 * num_digits - 2 - 2 * log_num_digits;
    this->num_triples_std = log_num_digits;
    this->num_triples = num_triples_std + num_triples_corr;
    if (beta == 8)
        this->mask_beta = -1;
    else
        this->mask_beta = (1 << beta) - 1;
    this->mask_r = (1 << r) - 1;
    this->beta_pow = 1 << beta;
    this->prng = new PRNG(oc::sysRandomSeed());
}

void MillionaireProtocolSender::drelu(uint8_t *res, uint64_t *data, osuCrypto::Socket &chl)
{
    int input_size = num_cmps;

    u8 *carry = new u8[num_cmps];
    u8 *local_msb = new u8[num_cmps];
    u64 *values = new u64[num_cmps];
    u64 mask = (1ull << (l - 1)) - 1;

    for (u64 i = 0; i < num_cmps; i++) {
        values[i] = data[i] & mask;
        local_msb[i] = (data[i] >> (l - 1)) & 1;
    }

    compare(carry, values, chl);

    for (u64 i = 0; i < input_size; i++) {
        res[i] = carry[i] ^ local_msb[i] ^ 1;
        res[i] &= 1;
    }

    delete[] carry;
    delete[] local_msb;
    delete[] values;
}

void MillionaireProtocolSender::compare(uint8_t *res, uint64_t *data, osuCrypto::Socket &chl, bool greater_than, bool equality, int radix_base)
{
    coproto::sync_wait(triple_gen->genBaseOts(*prng, chl));

    int origin_num_cmps = num_cmps;
    // num_cmps should be a multiple of 8
    num_cmps = ceil(num_cmps / 8.0) * 8;

    uint64_t *data_ext;
    if (origin_num_cmps == num_cmps)
        data_ext = data;
    else {
        data_ext = new uint64_t[num_cmps];
        memcpy(data_ext, data, origin_num_cmps * sizeof(uint64_t));
        memset(data_ext + origin_num_cmps, 0, (num_cmps - origin_num_cmps) * sizeof(uint64_t));
    }

    uint8_t *digits;       // num_digits * num_cmps
    uint8_t *leaf_res_cmp; // num_digits * num_cmps
    uint8_t *leaf_res_eq;  // num_digits * num_cmps

    digits = new uint8_t[num_digits * num_cmps];
    leaf_res_cmp = new uint8_t[num_digits * num_cmps];
    leaf_res_eq = new uint8_t[num_digits * num_cmps];

    // Extract radix-digits from data
    for (int i = 0; i < num_digits; i++) // Stored from LSB to MSB
        for (int j = 0; j < num_cmps; j++)
            digits[i * num_cmps + j] = (uint8_t)(data_ext[j] >> i * beta) & mask_beta;

    {
        uint8_t **leaf_ot_messages; // (num_digits * num_cmps) X beta_pow (=2^beta)
        leaf_ot_messages = new uint8_t *[num_digits * num_cmps];
        for (int i = 0; i < num_digits * num_cmps; i++)
            leaf_ot_messages[i] = new uint8_t[beta_pow];

        // Set Leaf OT messages
        prng->get((uint8_t *)leaf_res_cmp, num_digits * num_cmps);
        prng->get((uint8_t *)leaf_res_eq, num_digits * num_cmps);

        for (int i = 0; i < num_digits * num_cmps; i++) {
            leaf_res_cmp[i] &= 1;
            leaf_res_eq[i] &= 1;
        }

        for (int i = 0; i < num_digits; i++) {
            for (int j = 0; j < num_cmps; j++) {
                if (i == 0) {
                    set_leaf_ot_messages(
                        leaf_ot_messages[i * num_cmps + j], digits[i * num_cmps + j], beta_pow, leaf_res_cmp[i * num_cmps + j], 0, greater_than, false);
                } else {
                    set_leaf_ot_messages(
                        leaf_ot_messages[i * num_cmps + j],
                        digits[i * num_cmps + j],
                        beta_pow,
                        leaf_res_cmp[i * num_cmps + j],
                        leaf_res_eq[i * num_cmps + j],
                        greater_than);
                }
            }
        }

        // Perform Leaf OTs

        // otpack->kkot_beta->send(leaf_ot_messages, num_cmps*(num_digits), 2);
        otpack->send(leaf_ot_messages, chl);

        // Cleanup
        for (int i = 0; i < num_digits * num_cmps; i++)
            delete[] leaf_ot_messages[i];
        delete[] leaf_ot_messages;
    }

    traverse_and_compute_ANDs(*triple_gen, num_cmps, leaf_res_eq, leaf_res_cmp, chl);

    for (int i = 0; i < origin_num_cmps; i++)
        res[i] = leaf_res_cmp[i];

    // Cleanup
    if (origin_num_cmps != num_cmps)
        delete[] data_ext;
    delete[] digits;
    delete[] leaf_res_cmp;
    delete[] leaf_res_eq;
}

void MillionaireProtocolSender::set_leaf_ot_messages(uint8_t *ot_messages, uint8_t digit, int N, uint8_t mask_cmp, uint8_t mask_eq, bool greater_than, bool eq)
{
    for (int i = 0; i < N; i++) {
        if (greater_than) {
            ot_messages[i] = ((digit > i) ^ mask_cmp);
        } else {
            ot_messages[i] = ((digit < i) ^ mask_cmp);
        }
        if (eq) {
            ot_messages[i] = (ot_messages[i] << 1) | ((digit == i) ^ mask_eq);
        }
    }
}

void MillionaireProtocolSender::traverse_and_compute_ANDs(
    SilentOtTriple &triple_gen, int num_cmps, uint8_t *leaf_res_eq, uint8_t *leaf_res_cmp, osuCrypto::Socket &chl)
{
    std::vector<block> A(num_triples_round / 128);
    std::vector<block> B(num_triples_round / 128);
    std::vector<block> C(num_triples_round / 128);

    // Generate required Bit-Triples
#if USE_CHEETAH
    coproto::sync_wait(triple_gen.expand(A, B, C, *prng, chl));

    u8 *ai = reinterpret_cast<uint8_t *>(A.data());
    u8 *bi = reinterpret_cast<uint8_t *>(B.data());
    u8 *ci = reinterpret_cast<uint8_t *>(C.data());

#elif defined(WAN_EXEC)
    // std::cout<<"Running on WAN_EXEC; Skipping correlated triples"<<std::endl;
    triple_gen->generate(party, &triples_std, _16KKOT_to_4OT);
#else
    triple_gen->generate(party, &triples_corr, _8KKOT);
    triple_gen->generate(party, &triples_std, _16KKOT_to_4OT);
#endif
    // std::cout << "Bit Triples Generated" << std::endl;

    // Combine leaf OT results in a bottom-up fashion
    int counter_std = 0, old_counter_std = 0;
    int counter_corr = 0, old_counter_corr = 0;
    int counter_combined = 0, old_counter_combined = 0;
    uint8_t *ei = new uint8_t[(num_triples * num_cmps) / 8];
    uint8_t *fi = new uint8_t[(num_triples * num_cmps) / 8];
    uint8_t *e = new uint8_t[(num_triples * num_cmps) / 8];
    uint8_t *f = new uint8_t[(num_triples * num_cmps) / 8];

    for (int i = 1; i < num_digits; i *= 2) {
        for (int j = 0; j < num_digits and j + i < num_digits; j += 2 * i) {
            if (j == 0) {
                AND_step_1(
                    ei + (counter_std * num_cmps) / 8,
                    fi + (counter_std * num_cmps) / 8,
                    leaf_res_cmp + j * num_cmps,
                    leaf_res_eq + (j + i) * num_cmps,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_std++;
                counter_combined++;

            } else {
                AND_step_1(
                    ei + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    fi + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    leaf_res_cmp + j * num_cmps,
                    leaf_res_eq + (j + i) * num_cmps,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_combined++;
                AND_step_1(
                    ei + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    fi + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    leaf_res_eq + j * num_cmps,
                    leaf_res_eq + (j + i) * num_cmps,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_combined++;
                counter_corr++;
            }
        }
        int offset_std = (old_counter_std * num_cmps) / 8;
        int size_std = ((counter_std - old_counter_std) * num_cmps) / 8;
        int offset_corr = ((num_triples_std + 2 * old_counter_corr) * num_cmps) / 8;
        int size_corr = (2 * (counter_corr - old_counter_corr) * num_cmps) / 8;
        if (size_corr != 0) {
            coproto::sync_wait(chl.send(oc::span<u8>(ei + offset_std, size_std)));
            coproto::sync_wait(chl.send(oc::span<u8>(ei + offset_corr, size_corr)));
            coproto::sync_wait(chl.send(oc::span<u8>(fi + offset_std, size_std)));
            coproto::sync_wait(chl.send(oc::span<u8>(fi + offset_corr, size_corr)));

            coproto::sync_wait(chl.recv(oc::span<u8>(e + offset_std, size_std)));
            coproto::sync_wait(chl.recv(oc::span<u8>(e + offset_corr, size_corr)));
            coproto::sync_wait(chl.recv(oc::span<u8>(f + offset_std, size_std)));
            coproto::sync_wait(chl.recv(oc::span<u8>(f + offset_corr, size_corr)));
        } else {
            coproto::sync_wait(chl.send(oc::span<u8>(ei + offset_std, size_std)));
            coproto::sync_wait(chl.send(oc::span<u8>(fi + offset_std, size_std)));

            coproto::sync_wait(chl.recv(oc::span<u8>(e + offset_std, size_std)));
            coproto::sync_wait(chl.recv(oc::span<u8>(f + offset_std, size_std)));
        }
        for (int i = 0; i < size_std; i++) {
            e[i + offset_std] ^= ei[i + offset_std];
            f[i + offset_std] ^= fi[i + offset_std];
        }
        for (int i = 0; i < size_corr; i++) {
            e[i + offset_corr] ^= ei[i + offset_corr];
            f[i + offset_corr] ^= fi[i + offset_corr];
        }

        counter_std = old_counter_std;
        counter_corr = old_counter_corr;

        counter_combined = old_counter_combined;

        for (int j = 0; j < num_digits and j + i < num_digits; j += 2 * i) {
            if (j == 0) {
                AND_step_2(
                    leaf_res_cmp + j * num_cmps,
                    e + (counter_std * num_cmps) / 8,
                    f + (counter_std * num_cmps) / 8,
                    ei + (counter_std * num_cmps) / 8,
                    fi + (counter_std * num_cmps) / 8,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    (ci) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_combined++;

                for (int k = 0; k < num_cmps; k++)
                    leaf_res_cmp[j * num_cmps + k] ^= leaf_res_cmp[(j + i) * num_cmps + k];
                counter_std++;
            } else {
                AND_step_2(
                    leaf_res_cmp + j * num_cmps,
                    e + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    f + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    ei + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    fi + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    (ci) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_combined++;
                AND_step_2(
                    leaf_res_eq + j * num_cmps,
                    e + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    f + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    ei + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    fi + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    (ci) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_combined++;

                for (int k = 0; k < num_cmps; k++)
                    leaf_res_cmp[j * num_cmps + k] ^= leaf_res_cmp[(j + i) * num_cmps + k];
                counter_corr++;
            }
        }
        old_counter_std = counter_std;
        old_counter_corr = counter_corr;

        old_counter_combined = counter_combined;
    }

#if defined(WAN_EXEC) || USE_CHEETAH
    assert(counter_combined == num_triples);
#else
    assert(counter_std == num_triples_std);
    assert(2 * counter_corr == num_triples_corr);
#endif

    // cleanup
    delete[] ei;
    delete[] fi;
    delete[] e;
    delete[] f;
}

void MillionaireProtocolSender::AND_step_1(
    uint8_t *ei, // evaluates batch of 8 ANDs
    uint8_t *fi,
    uint8_t *xi,
    uint8_t *yi,
    uint8_t *ai,
    uint8_t *bi,
    int num_ANDs)
{
    assert(num_ANDs % 8 == 0);
    for (int i = 0; i < num_ANDs; i += 8) {
        ei[i / 8] = ai[i / 8];
        fi[i / 8] = bi[i / 8];
        ei[i / 8] ^= bool_to_uint8(xi + i, 8);
        fi[i / 8] ^= bool_to_uint8(yi + i, 8);
    }
}
void MillionaireProtocolSender::AND_step_2(
    uint8_t *zi, // evaluates batch of 8 ANDs
    uint8_t *e,
    uint8_t *f,
    uint8_t *ei,
    uint8_t *fi,
    uint8_t *ai,
    uint8_t *bi,
    uint8_t *ci,
    int num_ANDs)
{
    assert(num_ANDs % 8 == 0);
    for (int i = 0; i < num_ANDs; i += 8) {
        uint8_t temp_z;

        temp_z = e[i / 8] & f[i / 8];

        temp_z ^= f[i / 8] & ai[i / 8];
        temp_z ^= e[i / 8] & bi[i / 8];
        temp_z ^= ci[i / 8];
        uint8_to_bool(zi + i, temp_z, 8);
    }
}

MillionaireProtocolSender::~MillionaireProtocolSender()
{
    delete triple_gen;
    delete otpack;
    delete prng;
}

MillionaireProtocolRecver::MillionaireProtocolRecver(int num_cmps, int bitlength, int radix_base)
{
    this->num_cmps = num_cmps;

    configure(bitlength, radix_base);

    otpack = new NcoOTRecver(this->num_digits * roundUpTo(this->num_cmps, 8));
    num_triples_round = roundUpTo(this->num_triples * roundUpTo(this->num_cmps, 8), 128);
    triple_gen = new SilentOtTriple();
    triple_gen->init(0, num_triples_round);
}

void MillionaireProtocolRecver::configure(int bitlength, int radix_base)
{
    assert(radix_base <= 8);
    assert(bitlength <= 64);
    this->l = bitlength;
    this->beta = radix_base;

    this->num_digits = ceil((double)l / beta);
    this->r = l % beta;
    this->log_alpha = bitlen(num_digits) - 1;
    this->log_num_digits = log_alpha + 1;
    this->num_triples_corr = 2 * num_digits - 2 - 2 * log_num_digits;
    this->num_triples_std = log_num_digits;
    this->num_triples = num_triples_std + num_triples_corr;
    if (beta == 8)
        this->mask_beta = -1;
    else
        this->mask_beta = (1 << beta) - 1;
    this->mask_r = (1 << r) - 1;
    this->beta_pow = 1 << beta;
    this->prng = new PRNG(oc::sysRandomSeed());
}

void MillionaireProtocolRecver::drelu(uint8_t *res, uint64_t *data, osuCrypto::Socket &chl)
{
    int input_size = num_cmps;

    u8 *carry = new u8[num_cmps];
    u8 *local_msb = new u8[num_cmps];
    u64 *values = new u64[num_cmps];
    u64 mask = (1ull << (l - 1)) - 1;

    for (u64 i = 0; i < num_cmps; i++) {
        values[i] = data[i] & mask;
        values[i] = mask - values[i];
        local_msb[i] = (data[i] >> (l - 1)) & 1;
    }

    compare(carry, values, chl);

    for (u64 i = 0; i < input_size; i++) {
        res[i] = carry[i] ^ local_msb[i] ^ 0;
        res[i] &= 1;
    }

    delete[] carry;
    delete[] local_msb;
    delete[] values;
}

void MillionaireProtocolRecver::compare(uint8_t *res, uint64_t *data, osuCrypto::Socket &chl, bool greater_than, bool equality, int radix_base)
{
    coproto::sync_wait(triple_gen->genBaseOts(*prng, chl));

    int origin_num_cmps = num_cmps;
    // num_cmps should be a multiple of 8
    num_cmps = ceil(num_cmps / 8.0) * 8;

    uint64_t *data_ext;
    if (origin_num_cmps == num_cmps)
        data_ext = data;
    else {
        data_ext = new uint64_t[num_cmps];
        memcpy(data_ext, data, origin_num_cmps * sizeof(uint64_t));
        memset(data_ext + origin_num_cmps, 0, (num_cmps - origin_num_cmps) * sizeof(uint64_t));
    }

    uint8_t *digits;       // num_digits * num_cmps
    uint8_t *leaf_res_cmp; // num_digits * num_cmps
    uint8_t *leaf_res_eq;  // num_digits * num_cmps

    digits = new uint8_t[num_digits * num_cmps];
    leaf_res_cmp = new uint8_t[num_digits * num_cmps];
    leaf_res_eq = new uint8_t[num_digits * num_cmps];

    // Extract radix-digits from data
    for (int i = 0; i < num_digits; i++) // Stored from LSB to MSB
        for (int j = 0; j < num_cmps; j++)
            digits[i * num_cmps + j] = (uint8_t)(data_ext[j] >> i * beta) & mask_beta;

    { // party = sci::BOB
      // Perform Leaf OTs

        otpack->recv(leaf_res_cmp, digits, chl);

        // Extract equality result from leaf_res_cmp
        for (int i = num_cmps; i < num_digits * num_cmps; i++) {
            leaf_res_eq[i] = leaf_res_cmp[i] & 1;
            leaf_res_cmp[i] >>= 1;
        }
    }

    traverse_and_compute_ANDs(*triple_gen, num_cmps, leaf_res_eq, leaf_res_cmp, chl);

    for (int i = 0; i < origin_num_cmps; i++)
        res[i] = leaf_res_cmp[i];

    // Cleanup
    if (origin_num_cmps != num_cmps)
        delete[] data_ext;
    delete[] digits;
    delete[] leaf_res_cmp;
    delete[] leaf_res_eq;
}

void MillionaireProtocolRecver::set_leaf_ot_messages(uint8_t *ot_messages, uint8_t digit, int N, uint8_t mask_cmp, uint8_t mask_eq, bool greater_than, bool eq)
{
    for (int i = 0; i < N; i++) {
        if (greater_than) {
            ot_messages[i] = ((digit > i) ^ mask_cmp);
        } else {
            ot_messages[i] = ((digit < i) ^ mask_cmp);
        }
        if (eq) {
            ot_messages[i] = (ot_messages[i] << 1) | ((digit == i) ^ mask_eq);
        }
    }
}

void MillionaireProtocolRecver::traverse_and_compute_ANDs(
    SilentOtTriple &triple_gen, int num_cmps, uint8_t *leaf_res_eq, uint8_t *leaf_res_cmp, osuCrypto::Socket &chl)
{
#if defined(WAN_EXEC) || USE_CHEETAH
    std::vector<block> A(num_triples_round / 128);
    std::vector<block> B(num_triples_round / 128);
    std::vector<block> C(num_triples_round / 128);
#else
    Triple triples_corr(num_triples_corr * num_cmps, true, num_cmps);
    Triple triples_std(num_triples_std * num_cmps, true);
#endif
    // Generate required Bit-Triples
#if USE_CHEETAH

    coproto::sync_wait(triple_gen.expand(A, B, C, *prng, chl));

    u8 *ai = reinterpret_cast<uint8_t *>(A.data());
    u8 *bi = reinterpret_cast<uint8_t *>(B.data());
    u8 *ci = reinterpret_cast<uint8_t *>(C.data());

#elif defined(WAN_EXEC)
    // std::cout<<"Running on WAN_EXEC; Skipping correlated triples"<<std::endl;
    triple_gen->generate(party, &triples_std, _16KKOT_to_4OT);
#else
    triple_gen->generate(party, &triples_corr, _8KKOT);
    triple_gen->generate(party, &triples_std, _16KKOT_to_4OT);
#endif
    // std::cout << "Bit Triples Generated" << std::endl;

    // Combine leaf OT results in a bottom-up fashion
    int counter_std = 0, old_counter_std = 0;
    int counter_corr = 0, old_counter_corr = 0;
    int counter_combined = 0, old_counter_combined = 0;
    uint8_t *ei = new uint8_t[(num_triples * num_cmps) / 8];
    uint8_t *fi = new uint8_t[(num_triples * num_cmps) / 8];
    uint8_t *e = new uint8_t[(num_triples * num_cmps) / 8];
    uint8_t *f = new uint8_t[(num_triples * num_cmps) / 8];

    for (int i = 1; i < num_digits; i *= 2) {
        for (int j = 0; j < num_digits and j + i < num_digits; j += 2 * i) {
            if (j == 0) {
                AND_step_1(
                    ei + (counter_std * num_cmps) / 8,
                    fi + (counter_std * num_cmps) / 8,
                    leaf_res_cmp + j * num_cmps,
                    leaf_res_eq + (j + i) * num_cmps,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_std++;
                counter_combined++;
            } else {
#if defined(WAN_EXEC) || USE_CHEETAH
                AND_step_1(
                    ei + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    fi + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    leaf_res_cmp + j * num_cmps,
                    leaf_res_eq + (j + i) * num_cmps,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_combined++;
                AND_step_1(
                    ei + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    fi + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    leaf_res_eq + j * num_cmps,
                    leaf_res_eq + (j + i) * num_cmps,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_combined++;
                counter_corr++;
#else
                AND_step_1(
                    ei + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    fi + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    leaf_res_cmp + j * num_cmps,
                    leaf_res_eq + (j + i) * num_cmps,
                    (triples_corr.ai) + (2 * counter_corr * num_cmps) / 8,
                    (triples_corr.bi) + (2 * counter_corr * num_cmps) / 8,
                    num_cmps);
                AND_step_1(
                    ei + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    fi + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    leaf_res_eq + j * num_cmps,
                    leaf_res_eq + (j + i) * num_cmps,
                    (triples_corr.ai) + ((2 * counter_corr + 1) * num_cmps) / 8,
                    (triples_corr.bi) + ((2 * counter_corr + 1) * num_cmps) / 8,
                    num_cmps);
                counter_corr++;
#endif
            }
        }
        int offset_std = (old_counter_std * num_cmps) / 8;
        int size_std = ((counter_std - old_counter_std) * num_cmps) / 8;
        int offset_corr = ((num_triples_std + 2 * old_counter_corr) * num_cmps) / 8;
        int size_corr = (2 * (counter_corr - old_counter_corr) * num_cmps) / 8;

        // party = sci::BOB
        if (size_corr != 0) {
            coproto::sync_wait(chl.recv(oc::span<u8>(e + offset_std, size_std)));
            coproto::sync_wait(chl.recv(oc::span<u8>(e + offset_corr, size_corr)));
            coproto::sync_wait(chl.recv(oc::span<u8>(f + offset_std, size_std)));
            coproto::sync_wait(chl.recv(oc::span<u8>(f + offset_corr, size_corr)));

            coproto::sync_wait(chl.send(oc::span<u8>(ei + offset_std, size_std)));
            coproto::sync_wait(chl.send(oc::span<u8>(ei + offset_corr, size_corr)));
            coproto::sync_wait(chl.send(oc::span<u8>(fi + offset_std, size_std)));
            coproto::sync_wait(chl.send(oc::span<u8>(fi + offset_corr, size_corr)));
        } else {
            coproto::sync_wait(chl.recv(oc::span<u8>(e + offset_std, size_std)));
            coproto::sync_wait(chl.recv(oc::span<u8>(f + offset_std, size_std)));

            coproto::sync_wait(chl.send(oc::span<u8>(ei + offset_std, size_std)));
            coproto::sync_wait(chl.send(oc::span<u8>(fi + offset_std, size_std)));
        }
        for (int i = 0; i < size_std; i++) {
            e[i + offset_std] ^= ei[i + offset_std];
            f[i + offset_std] ^= fi[i + offset_std];
        }
        for (int i = 0; i < size_corr; i++) {
            e[i + offset_corr] ^= ei[i + offset_corr];
            f[i + offset_corr] ^= fi[i + offset_corr];
        }

        counter_std = old_counter_std;
        counter_corr = old_counter_corr;
#if defined(WAN_EXEC) || USE_CHEETAH
        counter_combined = old_counter_combined;
#endif
        for (int j = 0; j < num_digits and j + i < num_digits; j += 2 * i) {
            if (j == 0) {
#if defined(WAN_EXEC) || USE_CHEETAH
                AND_step_2(
                    leaf_res_cmp + j * num_cmps,
                    e + (counter_std * num_cmps) / 8,
                    f + (counter_std * num_cmps) / 8,
                    ei + (counter_std * num_cmps) / 8,
                    fi + (counter_std * num_cmps) / 8,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    (ci) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_combined++;
#else
                AND_step_2(
                    leaf_res_cmp + j * num_cmps,
                    e + (counter_std * num_cmps) / 8,
                    f + (counter_std * num_cmps) / 8,
                    ei + (counter_std * num_cmps) / 8,
                    fi + (counter_std * num_cmps) / 8,
                    (triples_std.ai) + (counter_std * num_cmps) / 8,
                    (triples_std.bi) + (counter_std * num_cmps) / 8,
                    (triples_std.ci) + (counter_std * num_cmps) / 8,
                    num_cmps);
#endif
                for (int k = 0; k < num_cmps; k++)
                    leaf_res_cmp[j * num_cmps + k] ^= leaf_res_cmp[(j + i) * num_cmps + k];
                counter_std++;
            } else {
#if defined(WAN_EXEC) || USE_CHEETAH
                AND_step_2(
                    leaf_res_cmp + j * num_cmps,
                    e + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    f + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    ei + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    fi + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    (ci) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_combined++;
                AND_step_2(
                    leaf_res_eq + j * num_cmps,
                    e + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    f + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    ei + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    fi + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    (ai) + (counter_combined * num_cmps) / 8,
                    (bi) + (counter_combined * num_cmps) / 8,
                    (ci) + (counter_combined * num_cmps) / 8,
                    num_cmps);
                counter_combined++;
#else
                AND_step_2(
                    leaf_res_cmp + j * num_cmps,
                    e + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    f + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    ei + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    fi + ((num_triples_std + 2 * counter_corr) * num_cmps) / 8,
                    (triples_corr.ai) + (2 * counter_corr * num_cmps) / 8,
                    (triples_corr.bi) + (2 * counter_corr * num_cmps) / 8,
                    (triples_corr.ci) + (2 * counter_corr * num_cmps) / 8,
                    num_cmps);
                AND_step_2(
                    leaf_res_eq + j * num_cmps,
                    e + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    f + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    ei + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    fi + ((num_triples_std + (2 * counter_corr + 1)) * num_cmps) / 8,
                    (triples_corr.ai) + ((2 * counter_corr + 1) * num_cmps) / 8,
                    (triples_corr.bi) + ((2 * counter_corr + 1) * num_cmps) / 8,
                    (triples_corr.ci) + ((2 * counter_corr + 1) * num_cmps) / 8,
                    num_cmps);
#endif
                for (int k = 0; k < num_cmps; k++)
                    leaf_res_cmp[j * num_cmps + k] ^= leaf_res_cmp[(j + i) * num_cmps + k];
                counter_corr++;
            }
        }
        old_counter_std = counter_std;
        old_counter_corr = counter_corr;
#if defined(WAN_EXEC) || USE_CHEETAH
        old_counter_combined = counter_combined;
#endif
    }

#if defined(WAN_EXEC) || USE_CHEETAH
    assert(counter_combined == num_triples);
#else
    assert(counter_std == num_triples_std);
    assert(2 * counter_corr == num_triples_corr);
#endif

    // cleanup
    delete[] ei;
    delete[] fi;
    delete[] e;
    delete[] f;
}

void MillionaireProtocolRecver::AND_step_1(
    uint8_t *ei, // evaluates batch of 8 ANDs
    uint8_t *fi,
    uint8_t *xi,
    uint8_t *yi,
    uint8_t *ai,
    uint8_t *bi,
    int num_ANDs)
{
    assert(num_ANDs % 8 == 0);
    for (int i = 0; i < num_ANDs; i += 8) {
        ei[i / 8] = ai[i / 8];
        fi[i / 8] = bi[i / 8];
        ei[i / 8] ^= bool_to_uint8(xi + i, 8);
        fi[i / 8] ^= bool_to_uint8(yi + i, 8);
    }
}
void MillionaireProtocolRecver::AND_step_2(
    uint8_t *zi, // evaluates batch of 8 ANDs
    uint8_t *e,
    uint8_t *f,
    uint8_t *ei,
    uint8_t *fi,
    uint8_t *ai,
    uint8_t *bi,
    uint8_t *ci,
    int num_ANDs)
{
    assert(num_ANDs % 8 == 0);
    for (int i = 0; i < num_ANDs; i += 8) {
        uint8_t temp_z;
        temp_z = 0;
        temp_z ^= f[i / 8] & ai[i / 8];
        temp_z ^= e[i / 8] & bi[i / 8];
        temp_z ^= ci[i / 8];
        uint8_to_bool(zi + i, temp_z, 8);
    }
}

MillionaireProtocolRecver::~MillionaireProtocolRecver()
{
    delete triple_gen;
    delete otpack;
    delete prng;
}