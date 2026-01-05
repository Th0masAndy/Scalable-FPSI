#include <cmath>
#include <coproto/Socket/AsioSocket.h>
#include <coproto/Socket/LocalAsyncSock.h>
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Common/CuckooIndex.h>
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Common/block.h>
#include <cstring>
#include <macoro/sync_wait.h>
#include <thread>
#include <vector>
#include <volePSI/Defines.h>
#include <volePSI/GMW/SilentTripleGen.h>
#include "cmp.h"
#include "eq.h"
#include "mul.h"
#include "mux.h"
#include "opprf.h"
#include "params.h"
#include "permute.h"
#include "sparsehash/dense_hash_map"
#include "utils.h"

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
            tmp.push_back(prng.get<u64>() * 2 * delta); // 2 delta apart
        }
        sendSet.push_back(tmp);
    }

    for (u64 i = 0; i < n; i++) {
        std::vector<u64> tmp;
        for (u64 j = 0; j < d; j++) {
            tmp.push_back(prng.get<u64>() * 4 * delta); // 4 delta apart
        }
        recvSet.push_back(tmp);
    }
}

void fpsiLowPx(const oc::CLP &cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    size_t d = cmd.getOr("d", 2);
    int delta = cmd.getOr("delta", 2);
    int verbose = cmd.getOr("v", 0);

    int numTry = cmd.getOr("try", 1);

    int prefixNum = static_cast<int>(std::ceil(std::log2(delta * 2 + 1)));
    int prefixLen = static_cast<int>(std::floor(std::log2(delta * 2 + 1))) + 1;
    int interSize = cmd.getOr("nn", 4);

    PRNG prng(sysRandomSeed());
    std::vector<std::vector<u64>> recvSet;
    std::vector<block> recvListKey;
    std::vector<block> recvListVal;
    std::vector<block> rand_R(n * (1 << d), ZeroBlock);
    std::vector<block> g(n * d * (1 << d) * prefixNum);
    std::vector<block> s_R(n * d * (1 << d));
    prng.get(s_R.data(), s_R.size());

    std::vector<std::vector<u64>> sendSet;
    std::vector<block> sendListKey;
    std::vector<block> sendListVal;
    std::vector<block> rand_S(n, ZeroBlock);
    std::vector<block> r_S(n * d);
    prng.get(r_S.data(), r_S.size());

    sampleData(sendSet, recvSet, delta, n, d, prng);

    oc::Timer time;
    time.setTimePoint("begin");

    for (size_t i = 0; i < recvSet.size(); i++) {
        auto neighbors = neigh(recvSet[i], delta);
        for (int j = 0; j < d; j++) {
            auto prefixes = getPrefixSet(recvSet[i][j], prefixLenMapNaive.at(2 * delta));
            for (auto neighbor : neighbors) {
                for (auto prefix : prefixes) {
                    recvListKey.push_back(blake3_hash(neighbor, j, prefix));
                }
            }
        }
    }

    for (size_t i = 0; i < sendSet.size(); i++) {
        for (int j = 0; j < d; j++) {
            auto prefixes = getIntervalPrefixSet(sendSet[i][j] - delta, sendSet[i][j] + delta, prefixLenMapNaive.at(2 * delta));
            for (auto prefix : prefixes) {
                sendListKey.push_back(blake3_hash(cell(sendSet[i], 2 * delta), j, prefix));
                sendListVal.push_back(r_S[i * d + j]);
            }
        }
    }

    while (sendListKey.size() < n * d * prefixLen) {
        sendListKey.push_back(prng.get<block>());
        sendListVal.push_back(prng.get<block>());
    }

    auto s = time.setTimePoint("preprocess done");

    auto socket = coproto::AsioSocket::makePair();

    std::thread recvThr([&]() {
        OpprfRevcer recver(n * d * (1 << d) * prefixNum, n * d * prefixLen);
        recver.setTimer(time);
        recver.recv(recvListKey, g, socket[0]);
    });

    std::thread sendThr([&]() {
        OpprfSender sender(n * d * (1 << d) * prefixNum, n * d * prefixLen);
        sender.setTimer(time);
        sender.send(sendListKey, sendListVal, socket[1]);
    });

    recvThr.join();
    sendThr.join();

    time.setTimePoint("first OPPRF done");
    std::cout << (socket[0].bytesReceived() + socket[0].bytesSent()) * 1.0 / 1024 / 1024 << " MB" << std::endl;

    std::thread recvThrReverse([&]() {
        std::vector<block> inputKeys(g.size());
        std::vector<block> outputVals(g.size());
        for (u64 i = 0; i < g.size(); i++) {
            inputKeys[i] = g[i];
            outputVals[i] = s_R[i / prefixNum];
        }

        OpprfSender sender(n * d, n * d * (1 << d) * prefixNum);
        sender.setTimer(time);
        sender.send(inputKeys, outputVals, socket[0]);

        for (int i = 0; i < rand_R.size(); i++) {
            for (int j = 0; j < d; j++) {
                rand_R[i] += s_R[i / (1 << d) * d * (1 << d) + i % (1 << d) + j * (1 << d)];
            }
        }
    });

    std::thread sendThrReverse([&]() {
        std::vector<block> out(r_S.size());

        OpprfRevcer recver(n * d, n * d * (1 << d) * prefixNum);
        recver.setTimer(time);
        recver.recv(r_S, out, socket[1]);

        for (int i = 0; i < rand_S.size(); i++) {
            for (int j = 0; j < d; j++) {
                rand_S[i] += out[i * d + j];
            }
        }
    });

    recvThrReverse.join();
    sendThrReverse.join();

    time.setTimePoint("second OPPRF done");
    std::cout << (socket[0].bytesReceived() + socket[0].bytesSent()) * 1.0 / 1024 / 1024 << " MB" << std::endl;

    std::thread sendMaskThr([&]() {
        PRNG prng;
        auto len = d * sizeof(u64);
        auto numBlocks = (len + sizeof(block) - 1) / sizeof(block);

        std::vector<u8> buffer(len * n);
        for (int i = 0; i < n; i++) {
            prng.SetSeed(rand_S[i], numBlocks);
            std::vector<u8> tmp(len);
            prng.get(tmp.data(), len);
            for (int j = 0; j < len; j++) {
                buffer[i * len + j] = tmp[j] ^ ((u8 *)&sendSet[i][0])[j];
            }
        }

        Hash(rand_S);

        coproto::sync_wait(socket[1].send(buffer));
        coproto::sync_wait(socket[1].send(rand_S));
    });

    std::thread recvMaskThr([&]() {
        PRNG prng;
        auto len = d * sizeof(u64);
        auto numBlocks = (len + sizeof(block) - 1) / sizeof(block);
        std::vector<block> rand_S(n, ZeroBlock);

        std::vector<u8> buffer(len * n);
        coproto::sync_wait(socket[0].recv(buffer));
        coproto::sync_wait(socket[0].recv(rand_S));

        auto map = google::dense_hash_map<block, u64, NoHash>{};
        map.resize(rand_R.size());
        map.set_empty_key(oc::ZeroBlock);

        Hash(rand_R);

        for (auto i = 0; i < rand_R.size(); i++) {
            map.insert({ rand_R[i], i });
        }

        for (auto i = 0; i < rand_S.size(); i++) {
            if (map.find(rand_S[i]) != map.end()) {
                std::vector<u8> tmp(len);
                prng.SetSeed(rand_S[i], numBlocks);
                prng.get(tmp.data(), len);
                for (int j = 0; j < len; j++) {
                    tmp[j] ^= buffer[i * len + j];
                }
            }
        }
    });

    sendMaskThr.join();
    recvMaskThr.join();

    auto e = time.setTimePoint("wLPSI done");

    // std::cout << time << std::endl;

    std::cout << (socket[0].bytesReceived() + socket[0].bytesSent()) * 1.0 / 1024 / 1024 << " MB" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(e - s).count() / double(1000 * 1000) << " seconds" << std::endl;
}

std::vector<u64> shift(const std::vector<u64> &input, int delta)
{
    auto cellId = cell(input, delta * 2);
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

    int u = std::ceil(std::log2(6 * delta + 1));
    int bytesLen = divCeil(u, 8);

    int numTry = cmd.getOr("try", 1);

    int prefixNum = static_cast<int>(std::ceil(std::log2(delta * 2 + 1)));
    int prefixLen = static_cast<int>(std::floor(std::log2(delta * 2 + 1))) + 1;
    int interSize = cmd.getOr("nn", 4);

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

    auto socket = coproto::AsioSocket::makePair();

    oc::Timer time;
    auto s = time.setTimePoint("begin");

    std::vector<u64> x;
    std::vector<u64> y;
    std::vector<u64> tag_s;
    std::vector<u64> tag_r;

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

        for (u64 i = 0; i < numBins; ++i) {
            auto &bin = cuckoo.mBins[i];
            if (bin.isEmpty() == false) {
                auto j = bin.hashIdx();
                auto b = bin.idx();

                Tx[i] = block(0, j) ^ (keys[b] << 2);
                mMapping[i] = b;
            } else {
                Tx[i] = block(i, 0);
                mMapping[i] = keys.size(); // invalid
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
                    u64 val;
                    memcpy(&val, ptr + 8 + k * bytesLen, bytesLen);
                    x[i * d + k] = val - (sendSet[b][k] - corner[k]);
                }
            }
            tag_s.push_back(reinterpret_cast<uint64_t *>(ptr)[0]);
        }
    });

    std::thread recvThr([&]() {
        std::vector<block> keys;
        for (u64 i = 0; i < n; ++i) {
            auto neighbors = neigh(recvSet[i], delta);
            for (size_t j = 0; j < neighbors.size(); j++) {
                keys.push_back(blake3_hash(neighbors[j]));
            }
        }
        while (keys.size() != (n * (1 << d))) {
            keys.push_back(prng.get<block>());
        }

        auto sIdx = SimpleIndex{};
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
                std::memcpy(ptr, &v[i], sizeof(u64));

                u64 idx = b / (1 << d);
                auto corner = shift(recvSet[idx], delta);
                for (u64 k = 0; k < d; ++k) {
                    u64 mask_i = 0;
                    memcpy(&mask_i, maskPtr + 8 + k * bytesLen, bytesLen);
                    y[i * d + k] = mask_i;
                    u64 val = recvSet[idx][k] - corner[k] - mask_i;
                    std::memcpy(ptr + 8 + k * bytesLen, &val, bytesLen);
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

    time.setTimePoint("matching done");

    std::vector<u8> resBits0(x.size() / d, 0);
    std::vector<u8> resBits1(y.size() / d, 0);

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
        PRNG prng(sysRandomSeed());
        SilentOtTriple tripleGen;
        u64 numTriples = roundUpTo(resBits0.size(), 128);
        tripleGen.init(1, numTriples);
        coproto::sync_wait(tripleGen.genBaseOts(prng, socket[1]));

        std::vector<block> A(numTriples / 128);
        std::vector<block> B(numTriples / 128);
        std::vector<block> C(numTriples / 128);

        coproto::sync_wait(tripleGen.expand(A, B, C, prng, socket[1]));

        u8 *ai = reinterpret_cast<uint8_t *>(A.data());
        u8 *bi = reinterpret_cast<uint8_t *>(B.data());
        u8 *ci = reinterpret_cast<uint8_t *>(C.data());
        u8 *ei = new u8[numTriples];
        u8 *fi = new u8[numTriples];
        u8 *e = new u8[numTriples];
        u8 *f = new u8[numTriples];

        std::vector<block> tag_s_block(tag_s.size());
        for (u64 i = 0; i < tag_s.size(); i++) {
            tag_s_block[i] = block(0, tag_s[i]);
        }
        BitVector eqResBit;
        ssPEQT(1, tag_s_block, eqResBit, socket[1], 1);

        std::vector<u8> eqRes(eqResBit.size());
        for (u64 i = 0; i < eqResBit.size(); i++) {
            eqRes[i] = (u8)eqResBit[i];
        }

        MillionaireProtocolSender::AND_step_1(ei, fi, eqRes.data(), resBits1.data(), ai, bi, numTriples);

        coproto::sync_wait(socket[1].send(oc::span<u8>(ei, numTriples)));
        coproto::sync_wait(socket[1].send(oc::span<u8>(fi, numTriples)));

        coproto::sync_wait(socket[1].recv(oc::span<u8>(e, numTriples)));
        coproto::sync_wait(socket[1].recv(oc::span<u8>(f, numTriples)));

        for (u64 i = 0; i < numTriples; i++) {
            e[i] ^= ei[i];
            f[i] ^= fi[i];
        }

        andRes1.resize(numTriples);
        MillionaireProtocolSender::AND_step_2(andRes1.data(), e, f, ei, fi, ai, bi, ci, numTriples);

        delete[] ei;
        delete[] fi;
        delete[] e;
        delete[] f;
    });

    std::thread andRecvThr([&]() {
        PRNG prng(sysRandomSeed());
        SilentOtTriple tripleGen;
        u64 numTriples = roundUpTo(resBits0.size(), 128);
        tripleGen.init(0, numTriples);
        coproto::sync_wait(tripleGen.genBaseOts(prng, socket[0]));

        std::vector<block> A(numTriples / 128);
        std::vector<block> B(numTriples / 128);
        std::vector<block> C(numTriples / 128);

        coproto::sync_wait(tripleGen.expand(A, B, C, prng, socket[0]));

        u8 *ai = reinterpret_cast<uint8_t *>(A.data());
        u8 *bi = reinterpret_cast<uint8_t *>(B.data());
        u8 *ci = reinterpret_cast<uint8_t *>(C.data());
        u8 *ei = new u8[numTriples];
        u8 *fi = new u8[numTriples];
        u8 *e = new u8[numTriples];
        u8 *f = new u8[numTriples];

        std::vector<block> tag_r_block(tag_r.size());
        for (u64 i = 0; i < tag_r.size(); i++) {
            tag_r_block[i] = block(0, tag_r[i]);
        }
        BitVector eqResBit;
        ssPEQT(0, tag_r_block, eqResBit, socket[0], 1);

        std::vector<u8> eqRes(eqResBit.size());
        for (u64 i = 0; i < eqResBit.size(); i++) {
            eqRes[i] = (u8)eqResBit[i];
        }

        MillionaireProtocolRecver::AND_step_1(ei, fi, eqRes.data(), resBits0.data(), ai, bi, numTriples);

        coproto::sync_wait(socket[0].send(oc::span<u8>(ei, numTriples)));
        coproto::sync_wait(socket[0].send(oc::span<u8>(fi, numTriples)));

        coproto::sync_wait(socket[0].recv(oc::span<u8>(e, numTriples)));
        coproto::sync_wait(socket[0].recv(oc::span<u8>(f, numTriples)));

        for (u64 i = 0; i < numTriples; i++) {
            e[i] ^= ei[i];
            f[i] ^= fi[i];
        }

        andRes0.resize(numTriples);
        MillionaireProtocolRecver::AND_step_2(andRes0.data(), e, f, ei, fi, ai, bi, ci, numTriples);

        delete[] ei;
        delete[] fi;
        delete[] e;
        delete[] f;
    });

    andSendThr.join();
    andRecvThr.join();

    time.setTimePoint("AND done");

    std::thread permuteSendThr([&]() {
        BitVector bits1(andRes1.data(), andRes1.size());
        BitVector output;
        std::vector<u32> pi;
        Opermute1(0, bits1, output, pi, socket[1], 1);
    });

    std::thread permuteRecvThr([&]() {
        BitVector bits0(andRes0.data(), andRes0.size());
        BitVector output;
        std::vector<u32> pi;
        Opermute1(1, bits0, output, pi, socket[0], 1);
    });

    permuteSendThr.join();
    permuteRecvThr.join();

    auto e = time.setTimePoint("permute done");

    // std::cout << time << std::endl;

    std::cout << (socket[0].bytesReceived() + socket[0].bytesSent()) * 1.0 / 1024 / 1024 << " MB" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(e - s).count() / double(1000 * 1000) << " seconds" << std::endl;
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
    int extBitsLen = bitsLen + static_cast<int>(std::ceil(std::log2(d)));

    auto n = x.size() / d;

    std::thread cmpSendThr([&]() {
        MillionaireProtocolSender sender(x.size(), bitsLen);

        std::vector<u8> cmpShare(x.size());
        sender.compare(cmpShare.data(), x.data(), chl[1]);

        MuxSender mux(x.size(), &chl[1]);

        std::vector<u64> res(x.size());
        std::vector<u64> x_vec(x.begin(), x.end());

        mux.muxA(cmpShare, x_vec, res, bitsLen);

        std::vector<u64> abs(x.size());
        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = x[i] - 2 * res[i];
        }

        std::vector<u64> dis_max(n, 0);

        MillionaireProtocolSender sender2(n, bitsLen);
        MuxSender mux2(n, &chl[1]);

        for (u64 i = 0; i < d; i++) {
            std::vector<u8> compare_res(n);
            std::vector<u64> curr(n, 0);
            for (u64 j = 0; j < n; ++j) {
                curr[j] = abs[j * d + i] - dis_max[j];
            }
            sender2.compare(compare_res.data(), curr.data(), chl[1]);
            mux2.muxA(compare_res, curr, curr, bitsLen);
            for (u64 j = 0; j < n; ++j) {
                dis_max[j] += curr[j];
            }
        }

        sender2.compare(resBits1.data(), dis_max.data(), chl[1]);
    });

    std::thread cmpRecvThr([&]() {
        MillionaireProtocolRecver recver(y.size(), bitsLen);

        std::vector<u8> cmpShare(y.size());
        recver.compare(cmpShare.data(), y.data(), chl[0]);

        MuxRecver mux(y.size(), &chl[0]);
        std::vector<u64> res(y.size());
        std::vector<u64> y_vec(y.begin(), y.end());

        for (u64 i = 0; i < y_vec.size(); ++i) {
            y_vec[i] = -y_vec[i];
        }

        mux.muxA(cmpShare, y_vec, res, bitsLen);

        std::vector<u64> abs(y.size());

        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = y[i] - 2 * res[i];
        }

        std::vector<u64> dis_max(n, 0);

        MillionaireProtocolRecver recver2(n, bitsLen);
        MuxRecver mux2(n, &chl[0]);

        for (u64 i = 0; i < d; i++) {
            std::vector<u8> compare_res(n);
            std::vector<u64> curr(n, 0);
            for (u64 j = 0; j < n; ++j) {
                curr[j] = abs[j * d + i] - dis_max[j];
            }
            recver2.compare(compare_res.data(), curr.data(), chl[0]);
            mux2.muxA(compare_res, curr, curr, bitsLen);
            for (u64 j = 0; j < n; ++j) {
                dis_max[j] += curr[j];
            }
        }

        recver2.compare(resBits0.data(), dis_max.data(), chl[0]);
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
    u64 delta_p = delta;
    int bitsLen = bytesLen * 8;
    int extBitsLen = bitsLen + static_cast<int>(std::ceil(std::log2(d)));
    int prefixLenIfmat = static_cast<int>(std::ceil(std::log2(delta_p * 2 + 1)));

    auto n = x.size() / d;

    std::thread cmpSendThr([&]() {
        MillionaireProtocolSender sender(x.size(), bitsLen);

        std::vector<u8> cmpShare(x.size());
        sender.compare(cmpShare.data(), x.data(), chl[1]);

        MuxSender mux(x.size(), &chl[1]);

        std::vector<u64> res(x.size());
        std::vector<u64> x_vec(x.begin(), x.end());

        mux.muxA(cmpShare, x_vec, res, bitsLen);

        std::vector<u64> abs(x.size());
        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = x[i] - 2 * res[i];
        }

        std::vector<u64> dis(n, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += abs[i * d + j];
            }
        }

        MillionaireProtocolSender sender2(dis.size(), extBitsLen);

        sender2.compare(resBits1.data(), dis.data(), chl[1]);
    });

    std::thread cmpRecvThr([&]() {
        MillionaireProtocolRecver recver(y.size(), bitsLen);

        std::vector<u8> cmpShare(y.size());
        recver.compare(cmpShare.data(), y.data(), chl[0]);

        MuxRecver mux(y.size(), &chl[0]);
        std::vector<u64> res(y.size());
        std::vector<u64> y_vec(y.begin(), y.end());

        for (u64 i = 0; i < y_vec.size(); ++i) {
            y_vec[i] = -y_vec[i];
        }

        mux.muxA(cmpShare, y_vec, res, bitsLen);

        std::vector<u64> abs(y.size());

        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = y[i] - 2 * res[i];
        }

        std::vector<u64> dis(y.size() / d, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += abs[i * d + j];
            }
        }

        MillionaireProtocolRecver recver2(dis.size(), extBitsLen);

        recver2.compare(resBits0.data(), dis.data(), chl[0]);
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
    int extBitsLen = 2 * bitsLen + static_cast<int>(std::ceil(std::log2(d)));
    int prefixLenIfmat = static_cast<int>(std::ceil(std::log2(delta_p * 2 + 1)));

    auto n = x.size() / d;

    std::thread cmpSendThr([&]() {
        MulSender sender(x.size(), &chl[1]);

        std::vector<u64> x_vec(x.begin(), x.end());
        std::vector<u64> dots(x.size());
        sender.mul(x_vec, dots);

        std::vector<u64> absSquare(x.size());
        for (u64 i = 0; i < absSquare.size(); ++i) {
            absSquare[i] = x[i] * x[i] + 2 * dots[i];
        }

        std::vector<u64> dis(n, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += absSquare[i * d + j];
            }
        }

        MillionaireProtocolSender sender2(dis.size(), extBitsLen);

        sender2.compare(resBits1.data(), dis.data(), chl[1]);
    });

    std::thread cmpRecvThr([&]() {
        MulRecver recver(y.size(), &chl[0]);

        std::vector<u64> y_vec(y.begin(), y.end());
        std::vector<u64> dots(y.size());
        recver.mul(y_vec, dots);

        std::vector<u64> absSquare(y.size());

        for (u64 i = 0; i < absSquare.size(); ++i) {
            absSquare[i] = y[i] * y[i] + 2 * dots[i];
        }

        std::vector<u64> dis(y.size() / d, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += absSquare[i * d + j];
            }
        }

        MillionaireProtocolRecver recver2(dis.size(), extBitsLen);

        recver2.compare(resBits0.data(), dis.data(), chl[0]);
    });

    cmpSendThr.join();
    cmpRecvThr.join();
}