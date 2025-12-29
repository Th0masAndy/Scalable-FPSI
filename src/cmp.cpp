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
    coproto::sync_wait(mOt->genSilentBaseOts(*prng, chl));

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
    coproto::sync_wait(mOt->genSilentBaseOts(*prng, chl));

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
