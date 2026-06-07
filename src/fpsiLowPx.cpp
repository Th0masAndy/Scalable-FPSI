#include <algorithm>
#include <coproto/Socket/AsioSocket.h>
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Common/CuckooIndex.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtSender.h>
#include <thread>
#include "and.h"
#include "cmp.h"
#include "eq.h"
#include "mul.h"
#include "mux.h"
#include "opprf.h"
#include "params.h"
#include "permute.h"
#include "sparsehash/dense_hash_map"
#include "utils.h"

using namespace volePSI;

namespace {
std::vector<block> tagsToBlocks(const std::vector<u64> &tags)
{
    std::vector<block> blocks(tags.size());
    for (u64 i = 0; i < tags.size(); ++i) {
        blocks[i] = block(tags[i], tags[i]);
    }
    return blocks;
}

std::vector<u8> bitVectorToBytes(const BitVector &bits)
{
    std::vector<u8> out(bits.size());
    for (u64 i = 0; i < bits.size(); ++i) {
        out[i] = static_cast<u8>(bits[i]);
    }
    return out;
}
} // namespace

void normL1(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl);

void normL2(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl);

void normL0(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl);

void preProcessPrefix(std::vector<std::vector<u64>> &inputs, std::vector<block> &listKey, std::vector<block> &listVal, std::vector<block> &r_R, int delta)
{
    PRNG prng(oc::sysRandomSeed());
    int d = inputs[0].size();
    for (size_t i = 0; i < inputs.size(); i++) {
        auto neighbors = neigh(inputs[i], delta);
        for (int j = 0; j < d; j++) {
            auto prefixes = getIntervalPrefixSet(inputs[i][j] - delta, inputs[i][j] + delta, prefixLenMapNaive.at(2 * delta));
            for (auto prefix : prefixes) {
                for (auto neighbor : neighbors) {
                    listKey.push_back(blake3_hash(neighbor, j, prefix));
                    listVal.push_back(r_R[i * d + j]);
                }
            }
        }
    }
}

void sampleData(std::vector<std::vector<u64>> &sendSet, std::vector<std::vector<u64>> &recvSet, int delta, u64 n, size_t d, PRNG &prng)
{
    for (u64 i = 0; i < n; i++) {
        std::vector<u64> tmp;
        for (u64 j = 0; j < d; j++) {
            tmp.push_back(prng.get<u64>() * 4 * delta + prng.get<u64>() % (2 * delta + 1) - delta); // 2 delta apart
        }
        sendSet.push_back(tmp);
    }

    for (u64 i = 0; i < n; i++) {
        std::vector<u64> tmp;
        for (u64 j = 0; j < d; j++) {
            tmp.push_back(prng.get<u64>() * 8 * delta + prng.get<u64>() % (4 * delta + 1) - 2 * delta); // 4 delta apart
        }
        recvSet.push_back(tmp);
    }
}

std::vector<u64> shift(const std::vector<u64> &input, int delta, bool isCell = false)
{
    std::vector<u64> cellId = cell(input, delta * 2);
    if (isCell) {
        cellId = input;
    }
    std::vector<u64> corner(cellId.size());
    for (u64 i = 0; i < cellId.size(); i++) {
        corner[i] = cellId[i] - 1;
        corner[i] = corner[i] * (delta * 2);
    }
    return corner;
}

void fpsiLowLpPx(const oc::CLP &cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    size_t d = cmd.getOr("d", 2);
    int delta = cmd.getOr("delta", 2);
    int verbose = cmd.getOr("v", 0);

    int lp = cmd.getOr("p", 2);

    int interSize = log2ceil(n);

    int u = 64;

    if (lp == 0) {
        u = log2ceil(6 * delta) + 2;
    } else if (lp == 1) {
        u = log2ceil(6 * delta * d) + 1;
    } else if (lp == 2) {
        u = log2ceil(u64(d) * u64(3 * delta) * u64(3 * delta));
    } else {
        throw std::runtime_error("unsupported lp norm");
    }

    int bytesLen = std::max<int>(2, divCeil(u, 8));

    if (verbose) {
        std::cout << "bytesLen: " << bytesLen << std::endl;
    }

    int numTry = cmd.getOr("try", 1);

    PRNG prng(sysRandomSeed());
    std::vector<std::vector<u64>> recvSet;
    std::vector<block> recvListKey;
    std::vector<block> recvListVal;
    std::vector<block> rand_R(n * (1 << d), ZeroBlock);
    std::vector<block> s_R(n * d * (1 << d));
    prng.get(s_R.data(), s_R.size());

    std::vector<std::vector<u64>> sendSet;
    std::vector<block> sendListKey;
    std::vector<block> sendListVal;
    std::vector<block> rand_S(n, ZeroBlock);
    std::vector<block> r_S(n * d);
    prng.get(r_S.data(), r_S.size());

    sampleData(sendSet, recvSet, delta, n, d, prng);

    for (u64 i = 0; i < interSize; i++) {
        u64 averageDiff = 0;
        if (lp == 0) {
            averageDiff = delta;
        } else if (lp == 1) {
            averageDiff = floor(delta / d);
        } else if (lp == 2) {
            averageDiff = floor(delta / std::sqrt(d));
        }
        for (size_t j = 0; j < d; j++) {
            recvSet[i][j] = sendSet[i][j] + prng.get<u64>() % (2 * averageDiff + 1) - averageDiff;
        }
    }

    auto socket = coproto::AsioSocket::makePair();

    oc::Timer time;
    auto s = time.setTimePoint("begin");

    for (size_t t = 0; t < numTry; t++) {
        std::vector<u64> x;
        std::vector<u64> y;
        std::vector<u64> tag_s;
        std::vector<u64> tag_r;
        std::vector<u64> cuckooMap(n);
        std::vector<u64> reverseCuckooMap;

        std::thread sendThr([&]() {
            std::vector<block> keys(n);
            for (u64 i = 0; i < n; ++i) {
                keys[i] = blake3_hash(cell(sendSet[i], 2 * delta));
            }

            auto cuckoo = oc::CuckooIndex<>{};
            cuckoo.init(keys.size(), 40, 0, 3);
            cuckoo.insert(keys, ZeroBlock);

            std::vector<block> Tx(cuckoo.mNumBins);
            auto r = oc::Matrix<u8>(Tx.size(), 8 + d * bytesLen);
            u64 numBins = cuckoo.mBins.size();
            std::vector<u64> mMapping(numBins);
            reverseCuckooMap.resize(numBins);

            for (u64 i = 0; i < numBins; ++i) {
                auto &bin = cuckoo.mBins[i];
                if (bin.isEmpty() == false) {
                    auto j = bin.hashIdx();
                    auto b = bin.idx();

                    Tx[i] = block(0, j) ^ (keys[b] << 2);
                    mMapping[i] = b;
                    cuckooMap[b] = i;
                    reverseCuckooMap[i] = b;
                } else {
                    Tx[i] = block(i, 0);
                    mMapping[i] = keys.size(); // invalid
                    reverseCuckooMap[i] = n;   // invalid
                }
            }

            PRNG prng(sysRandomSeed());

            RsOpprfReceiver recver;
            recver.setTimer(time);
            coproto::sync_wait(recver.receive(n * (1 << d) * 3, Tx, r, prng, 1, socket[1]));

            x.resize(numBins * d);
            for (u64 i = 0; i < numBins; i++) {
                auto b = mMapping[i];
                auto ptr = r.data(i);
                if (b != keys.size()) {
                    auto corner = shift(sendSet[b], delta);
                    for (u64 k = 0; k < d; ++k) {
                        u64 val = 0;
                        memcpy(&val, ptr + 8 + k * bytesLen, bytesLen);
                        x[i * d + k] = sendSet[b][k] - corner[k] - val;
                    }
                }
                tag_s.push_back(reinterpret_cast<uint64_t *>(ptr)[0]);
            }
        });

        std::thread recvThr([&]() {
            std::vector<block> keys;
            std::vector<std::vector<u64>> all_cells;
            std::vector<u64> mMapping(n * (1 << d), 0);
            u64 count = 0;
            for (u64 i = 0; i < n; ++i) {
                auto neighbors = neigh(recvSet[i], delta);
                for (size_t j = 0; j < neighbors.size(); j++) {
                    all_cells.push_back(neighbors[j]);
                    mMapping[count++] = i;
                }
            }

            for (auto cell : all_cells) {
                keys.push_back(blake3_hash(cell));
            }

            while (keys.size() != (n * (1 << d))) {
                keys.push_back(prng.get<block>());
            }

            auto sIdx = volePSI::SimpleIndex{};
            auto params = oc::CuckooIndex<>::selectParams(n, 40, 0, 3);
            auto numBins = params.numBins();
            sIdx.init(numBins, keys.size(), 40, 3);
            sIdx.insertItems(keys, ZeroBlock);

            auto r = std::vector<u8>(numBins * d * bytesLen);
            auto v = std::vector<u64>(numBins);
            prng.get(r.data(), r.size());
            prng.get(v.data(), v.size());

            tag_r = v;
            y.resize(numBins * d);

            std::vector<block> Ty(keys.size() * 3);
            auto ctr = 0;
            auto Tv = oc::Matrix<u8>(Ty.size(), 8 + d * bytesLen);

            for (u64 i = 0; i < numBins; ++i) {
                auto bin = sIdx.mBins[i];
                auto size = sIdx.mBinSizes[i];

                for (u64 p = 0; p < size; ++p) {
                    auto j = bin[p].hashIdx();
                    auto b = bin[p].idx();
                    Ty[ctr] = block(0, j) ^ (keys[b] << 2);
                    auto ptr = Tv.data(ctr);
                    auto maskPtr = r.data() + i * d * bytesLen;
                    memcpy(ptr, &v[i], sizeof(u64));

                    auto corner = shift(all_cells[b], delta, 1);
                    for (u64 k = 0; k < d; ++k) {
                        u64 mask_i = 0;
                        memcpy(&mask_i, maskPtr + k * bytesLen, bytesLen);
                        y[i * d + k] = mask_i;
                        u64 val = recvSet[mMapping[b]][k] - corner[k] + mask_i;
                        memcpy(ptr + 8 + k * bytesLen, &val, bytesLen);
                    }
                    ctr++;
                }
            }
            PRNG prng(sysRandomSeed());

            RsOpprfSender sender;
            sender.setTimer(time);
            coproto::sync_wait(sender.send(numBins, Ty, Tv, prng, 1, socket[0]));
        });

        sendThr.join();
        recvThr.join();

        auto t_matching = time.setTimePoint("matching done");

        // std::cout << "matching time: " << std::chrono::duration_cast<std::chrono::microseconds>(t_matching - s).count() / double(1000 * 1000) << " seconds"
        //           << std::endl;

        // std::cout << "comm after matching: " << (socket[0].bytesReceived() + socket[0].bytesSent()) * 1.0 / 1024 / 1024 << " MB" << std::endl;

        std::vector<u8> resBits0(x.size() / d, 0);
        std::vector<u8> resBits1(y.size() / d, 0);

        // todo: replace mux with all and

        if (lp == 1) {
            normL1(oc::span<u64>(x.data(), x.size()), oc::span<u64>(y.data(), y.size()), resBits0, resBits1, d, delta, bytesLen, socket);
        } else if (lp == 2) {
            normL2(oc::span<u64>(x.data(), x.size()), oc::span<u64>(y.data(), y.size()), resBits0, resBits1, d, delta, bytesLen, socket);
        } else if (lp == 0) {
            normL0(oc::span<u64>(x.data(), x.size()), oc::span<u64>(y.data(), y.size()), resBits0, resBits1, d, delta, bytesLen, socket);
        } else {
            throw std::runtime_error("lp type not supported");
        }

        time.setTimePoint("norm done");

        std::vector<u8> andRes0;
        std::vector<u8> andRes1;

        std::thread andSendThr([&]() {
            auto tagBlocks = tagsToBlocks(tag_s);
            BitVector eqResBit;
            ssPEQT(1, tagBlocks, eqResBit, socket[1], 1);
            auto eqRes = bitVectorToBytes(eqResBit);

            AndSender sender(resBits0.size(), &socket[1]);
            sender.andBits(eqRes, resBits1, andRes1);
        });

        std::thread andRecvThr([&]() {
            auto tagBlocks = tagsToBlocks(tag_r);
            BitVector eqResBit;
            ssPEQT(0, tagBlocks, eqResBit, socket[0], 1);
            auto eqRes = bitVectorToBytes(eqResBit);

            AndRecver recver(resBits0.size(), &socket[0]);
            recver.andBits(eqRes, resBits0, andRes0);
        });

        andSendThr.join();
        andRecvThr.join();

        time.setTimePoint("AND done");

        BitVector finalBits(andRes0.size());
        BitVector output0;
        BitVector output1;

        std::vector<u32> permute;

        std::thread permuteSendThr([&]() {
            BitVector bits1(andRes1.size());
            for (u64 i = 0; i < andRes1.size(); i++) {
                bits1[i] = andRes1[i] & 1;
            }

            Opermute1(0, bits1, output1, permute, socket[1], 1);

            coproto::sync_wait(socket[1].send(output1));
        });

        std::thread permuteRecvThr([&]() {
            BitVector bits0(andRes0.size());
            for (u64 i = 0; i < andRes0.size(); i++) {
                bits0[i] = andRes0[i] & 1;
            }
            std::vector<u32> pi;
            Opermute1(1, bits0, output0, pi, socket[0], 1);

            coproto::sync_wait(socket[0].recv(finalBits));
            finalBits = finalBits ^ output0;
        });

        permuteSendThr.join();
        permuteRecvThr.join();

        time.setTimePoint("permute done");

        if (verbose) {
            u64 counter = 0;
            for (u64 i = 0; i < output0.size(); i++) {
                if (finalBits[i]) {
                    std::cout << "final match at index " << i << " oringal cuckoo index " << permute[i] << " oringal index " << reverseCuckooMap[permute[i]]
                              << std::endl;
                    counter++;
                }
            }
            std::cout << "total final match num: " << counter << std::endl;
        }

        std::thread sendOT([&]() {
            auto num = finalBits.size();
            SilentOtExtSender sender;
            PRNG prng(oc::sysRandomSeed());
            sender.configure(num);

            coproto::sync_wait(sender.genSilentBaseOts(prng, socket[1]));

            std::vector<std::array<block, 2>> messages(num);

            coproto::sync_wait(sender.send(messages, prng, socket[1]));

            std::vector<u64> correctMessages(num * d);
            for (u64 i = 0; i < num; i++) {
                u64 real_idx = reverseCuckooMap[permute[i]];
                if (real_idx < n) {
                    prng.SetSeed(messages[i][1]);
                    for (u64 j = 0; j < d; j++) {
                        correctMessages[i * d + j] = prng.get<u64>() ^ sendSet[real_idx][j];
                    }
                } else {
                    prng.SetSeed(messages[i][0]);
                    for (u64 j = 0; j < d; j++) {
                        correctMessages[i * d + j] = prng.get<u64>();
                    }
                }
            }
            coproto::sync_wait(socket[1].send(correctMessages));
        });

        std::vector<std::vector<u64>> result;
        std::thread recvOT([&]() {
            oc::BitVector choiceVec = finalBits;
            auto num = choiceVec.size();
            SilentOtExtReceiver receiver;
            PRNG prng(oc::sysRandomSeed());
            receiver.configure(num);

            coproto::sync_wait(receiver.genSilentBaseOts(prng, socket[0]));

            std::vector<block> recvMessages(num);

            coproto::sync_wait(receiver.receive(choiceVec, recvMessages, prng, socket[0]));

            std::vector<u64> recvCorrectMessages(num * d);
            coproto::sync_wait(socket[0].recv(recvCorrectMessages));

            for (u64 i = 0; i < num; i++) {
                if (choiceVec[i] & 1) {
                    prng.SetSeed(recvMessages[i]);
                    std::vector<u64> tmp(d);
                    for (u64 j = 0; j < d; j++) {
                        tmp[j] = prng.get<u64>() ^ recvCorrectMessages[i * d + j];
                    }
                    result.push_back(tmp);
                }
            }
        });

        sendOT.join();
        recvOT.join();

        time.setTimePoint("OT done");
    }

    auto e = time.setTimePoint("all done");

    if (verbose) {
        std::cout << time << std::endl;
    }

    std::cout << (socket[0].bytesReceived() + socket[0].bytesSent()) * 1.0 / double(numTry) / 1024 / 1024 << " MB" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(e - s).count() * 1.0 / double(numTry) / double(1000 * 1000) << " seconds" << std::endl;
}

void normL0(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl)
{
    int bitsLen = bytesLen * 8;
    u64 mask = (1ull << bitsLen) - 1;

    auto n = x.size() / d;

    std::thread cmpSendThr([&]() {
        MillionaireProtocolSender sender(x.size(), bitsLen);

        std::vector<u8> cmpShare(x.size());
        sender.drelu(cmpShare.data(), x.data(), chl[1]);

        MuxSender mux(x.size(), &chl[1]);

        std::vector<u64> res(x.size());
        std::vector<u64> x_vec(x.begin(), x.end());

        mux.muxA(cmpShare, x_vec, res, bitsLen);

        std::vector<u64> abs(x.size());
        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = 2 * res[i] - x[i];
        }

        std::vector<u64> dis_max(n, 0);

        MillionaireProtocolSender sender2(n, bitsLen);
        MuxSender mux2(n, &chl[1]);

        for (u64 i = 0; i < d; i++) {
            std::vector<u8> compare_res(n);
            std::vector<u64> curr(n, 0);
            std::vector<u64> res(n, 0);
            for (u64 j = 0; j < n; ++j) {
                curr[j] = abs[j * d + i] - dis_max[j];
            }
            sender2.drelu(compare_res.data(), curr.data(), chl[1]);
            mux2.muxA(compare_res, curr, res, bitsLen);
            for (u64 j = 0; j < n; ++j) {
                dis_max[j] += res[j];
            }
        }

        for (u64 i = 0; i < dis_max.size(); ++i) {
            dis_max[i] = delta - dis_max[i];
        }

        sender2.drelu(resBits1.data(), dis_max.data(), chl[1]);
    });

    std::thread cmpRecvThr([&]() {
        MillionaireProtocolRecver recver(y.size(), bitsLen);

        std::vector<u8> cmpShare(y.size());

        recver.drelu(cmpShare.data(), y.data(), chl[0]);

        MuxRecver mux(y.size(), &chl[0]);
        std::vector<u64> res(y.size());
        std::vector<u64> y_vec(y.begin(), y.end());

        mux.muxA(cmpShare, y_vec, res, bitsLen);

        std::vector<u64> abs(y.size());

        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = 2 * res[i] - y[i];
        }

        std::vector<u64> dis_max(n, 0);

        MillionaireProtocolRecver recver2(n, bitsLen);
        MuxRecver mux2(n, &chl[0]);

        for (u64 i = 0; i < d; i++) {
            std::vector<u8> compare_res(n);
            std::vector<u64> curr(n, 0);
            std::vector<u64> res(n, 0);
            for (u64 j = 0; j < n; ++j) {
                curr[j] = abs[j * d + i] - dis_max[j];
            }
            recver2.drelu(compare_res.data(), curr.data(), chl[0]);
            mux2.muxA(compare_res, curr, res, bitsLen);
            for (u64 j = 0; j < n; ++j) {
                dis_max[j] += res[j];
            }
        }

        for (u64 i = 0; i < dis_max.size(); ++i) {
            dis_max[i] = (-dis_max[i]) & mask;
        }

        recver2.drelu(resBits0.data(), dis_max.data(), chl[0]);
    });

    cmpSendThr.join();
    cmpRecvThr.join();
}

void normL1(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl)
{
    int bitsLen = bytesLen * 8;
    u64 mask = (1ull << bitsLen) - 1;

    auto n = x.size() / d;

    std::thread cmpSendThr([&]() {
        MillionaireProtocolSender sender(x.size(), bitsLen);

        std::vector<u8> cmpShare(x.size());
        sender.drelu(cmpShare.data(), x.data(), chl[1]);

        MuxSender mux(x.size(), &chl[1]);

        std::vector<u64> res(x.size());
        std::vector<u64> x_vec(x.begin(), x.end());

        mux.muxA(cmpShare, x_vec, res, bitsLen);

        std::vector<u64> abs(x.size());
        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = 2 * res[i] - x[i];
        }

        std::vector<u64> dis(n, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += abs[i * d + j];
            }
            dis[i] = delta - dis[i];
        }

        MillionaireProtocolSender sender2(dis.size(), bitsLen);

        sender2.drelu(resBits1.data(), dis.data(), chl[1]);
    });

    std::thread cmpRecvThr([&]() {
        MillionaireProtocolRecver recver(y.size(), bitsLen);

        std::vector<u8> cmpShare(y.size());
        recver.drelu(cmpShare.data(), y.data(), chl[0]);

        MuxRecver mux(y.size(), &chl[0]);
        std::vector<u64> res(y.size());
        std::vector<u64> y_vec(y.begin(), y.end());

        mux.muxA(cmpShare, y_vec, res, bitsLen);

        std::vector<u64> abs(y.size());

        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = 2 * res[i] - y[i];
        }

        std::vector<u64> dis(y.size() / d, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += abs[i * d + j];
            }
            dis[i] = (-dis[i]) & mask;
        }

        MillionaireProtocolRecver recver2(dis.size(), bitsLen);

        recver2.drelu(resBits0.data(), dis.data(), chl[0]);
    });

    cmpSendThr.join();
    cmpRecvThr.join();
}

void normL2(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl)
{
    u64 delta_p = delta * delta;
    int bitsLen = bytesLen * 8;

    u64 mask = (1ull << bitsLen) - 1;

    auto n = x.size() / d;

    std::thread cmpSendThr([&]() {
        MulSender sender(x.size(), &chl[1], bitsLen);

        std::vector<u64> x_vec(x.begin(), x.end());
        std::vector<u64> dots(x.size());
        sender.mul(x_vec, dots);

        std::vector<u64> absSquare(x.size());
        for (u64 i = 0; i < absSquare.size(); ++i) {
            u64 square = (__uint128_t(x[i]) * __uint128_t(x[i])) & mask;
            absSquare[i] = square + 2 * dots[i];
            absSquare[i] = absSquare[i] & mask;
        }

        std::vector<u64> dis(n, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += absSquare[i * d + j];
            }
            dis[i] = (delta_p - dis[i]) & mask;
        }

        MillionaireProtocolSender sender2(dis.size(), bitsLen);

        sender2.drelu(resBits1.data(), dis.data(), chl[1]);
    });

    std::thread cmpRecvThr([&]() {
        MulRecver recver(y.size(), &chl[0], bitsLen);

        std::vector<u64> y_vec(y.begin(), y.end());
        std::vector<u64> dots(y.size());
        recver.mul(y_vec, dots);

        std::vector<u64> absSquare(y.size());

        for (u64 i = 0; i < absSquare.size(); ++i) {
            u64 square = (__uint128_t(y[i]) * __uint128_t(y[i])) & mask;
            absSquare[i] = square + 2 * dots[i];
            absSquare[i] = absSquare[i] & mask;
        }

        std::vector<u64> dis(y.size() / d, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += absSquare[i * d + j];
            }
            dis[i] = (-dis[i]) & mask;
        }

        MillionaireProtocolRecver recver2(dis.size(), bitsLen);

        recver2.drelu(resBits0.data(), dis.data(), chl[0]);
    });

    cmpSendThr.join();
    cmpRecvThr.join();
}
