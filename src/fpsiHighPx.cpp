#include <coproto/Common/macoro.h>
#include <coproto/Socket/AsioSocket.h>
#include <coproto/Socket/LocalAsyncSock.h>
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Common/CuckooIndex.h>
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Common/Matrix.h>
#include <cryptoTools/Common/block.h>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <macoro/sync_wait.h>
#include <thread>
#include <vector>
#include <volePSI/Defines.h>
#include <volePSI/RsOpprf.h>
#include "cmp.h"
#include "eq.h"
#include "mul.h"
#include "mux.h"
#include "opprf.h"
#include "params.h"
#include "sparsehash/dense_hash_map"
#include "utils.h"

void normL1(std::vector<u64> x, std::vector<u64> y, std::vector<u8> &choiceBit, u64 d, int delta, std::array<coproto::AsioSocket, 2> &chl);
void normL2(oc::span<u64> x, oc::span<u64> y, std::vector<u8> &choiceBit, u64 d, int delta, std::array<coproto::AsioSocket, 2> &chl);

void sampleDataGlobalDisj(std::vector<std::vector<u64>> &sendSet, std::vector<std::vector<u64>> &recvSet, int delta, u64 n, size_t d, PRNG &prng)
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

void fpsiHighPx(const oc::CLP &cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    size_t d = cmd.getOr("d", 2);
    int delta = cmd.getOr("delta", 2);
    int verbose = cmd.getOr("v", 0);

    int numTry = cmd.getOr("try", 1);

    int prefixNum = static_cast<int>(std::ceil(std::log2(delta * 2 + 1)));
    int prefixLen = static_cast<int>(std::floor(std::log2(delta * 2 + 1))) + 1;

    int interSize = log2ceil(n);

    PRNG prng(sysRandomSeed());
    std::vector<std::vector<u64>> recvSet;
    std::vector<block> recvListKey;
    std::vector<block> recvListVal;
    std::vector<block> rand_R(n, ZeroBlock);
    std::vector<block> g(n * d * prefixNum);
    std::vector<block> s_R(n * d);
    prng.get(s_R.data(), s_R.size());

    std::vector<std::vector<u64>> sendSet;
    std::vector<block> sendListKey;
    std::vector<block> sendListVal;
    std::vector<block> rand_S(n, ZeroBlock);
    std::vector<block> r_S(n * d);
    prng.get(r_S.data(), r_S.size());

    sampleDataGlobalDisj(sendSet, recvSet, delta, n, d, prng);

    for (u64 i = 0; i < interSize; i++) {
        for (int j = 0; j < d; j++) {
            recvSet[i][j] = sendSet[i][j] + prng.get<u64>() % (2 * delta + 1) - delta;
        }
    }

    oc::Timer time;
    time.setTimePoint("begin");

    for (size_t i = 0; i < recvSet.size(); i++) {
        for (int j = 0; j < d; j++) {
            auto prefixes = getPrefixSet(recvSet[i][j], prefixLenMapNaive.at(2 * delta));
            for (auto prefix : prefixes) {
                recvListKey.push_back(prefix ^ block(j << 4, 0));
            }
        }
    }

    for (size_t i = 0; i < sendSet.size(); i++) {
        for (int j = 0; j < d; j++) {
            auto prefixes = getIntervalPrefixSet(sendSet[i][j] - delta, sendSet[i][j] + delta, prefixLenMapNaive.at(2 * delta));
            for (auto prefix : prefixes) {
                sendListKey.push_back(prefix ^ block(j << 4, 0));
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

    for (int t = 0; t < numTry; t++) {
        std::thread recvThr([&]() {
            OpprfRevcer recver(n * d * prefixNum, n * d * prefixLen);
            recver.setTimer(time);
            recver.recv(recvListKey, g, socket[0]);
        });

        std::thread sendThr([&]() {
            OpprfSender sender(n * d * prefixNum, n * d * prefixLen);
            sender.setTimer(time);
            sender.send(sendListKey, sendListVal, socket[1]);
        });

        recvThr.join();
        sendThr.join();

        time.setTimePoint("first OPPRF done");

        std::thread recvThrReverse([&]() {
            std::vector<block> inputKeys(g.size());
            std::vector<block> outputVals(g.size());
            for (u64 i = 0; i < g.size(); i++) {
                inputKeys[i] = g[i];
                outputVals[i] = s_R[i / prefixNum];
            }

            OpprfSender sender(n * d, n * d * prefixNum);
            sender.setTimer(time);
            sender.send(inputKeys, outputVals, socket[0]);

            for (int i = 0; i < rand_R.size(); i++) {
                for (int j = 0; j < d; j++) {
                    rand_R[i] += s_R[i * d + j];
                }
            }
        });

        std::thread sendThrReverse([&]() {
            std::vector<block> out(r_S.size());

            OpprfRevcer recver(n * d, n * d * prefixNum);
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

        if (verbose) {
            for (u64 i = 0; i < rand_R.size(); i++) {
                if (rand_R[i] == rand_S[i]) {
                    std::cout << "match at index " << i << std::endl;
                }
            }
        }

        time.setTimePoint("second OPPRF done");

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

        time.setTimePoint("wLPSI done");
    }

    auto e = time.setTimePoint("all done");

    if (verbose) {
        std::cout << time << std::endl;
    }

    std::cout << (socket[0].bytesReceived() + socket[0].bytesSent()) * 1.0 / numTry / 1024 / 1024 << " MB" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(e - s).count() * 1.0 / numTry / double(1000 * 1000) << " seconds" << std::endl;
}

void fpsiHighLpPx(const oc::CLP &cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    size_t d = cmd.getOr("d", 2);
    int delta = cmd.getOr("delta", 2);
    int verbose = cmd.getOr("v", 0);

    int lp = cmd.getOr("p", 2);

    int numTry = cmd.getOr("try", 1);

    int prefixNum = static_cast<int>(std::ceil(std::log2(delta * 2 + 1)));
    int prefixLen = static_cast<int>(std::floor(std::log2(delta * 2 + 1))) + 1;

    int interSize = log2ceil(n);

    PRNG prng(sysRandomSeed());
    std::vector<std::vector<u64>> recvSet;
    std::vector<block> recvListKey;
    std::vector<block> recvListVal;
    std::vector<block> rand_R(n, ZeroBlock);
    std::vector<block> g(n * d * prefixNum);
    std::vector<block> s_R(n * d);
    prng.get(s_R.data(), s_R.size());

    std::vector<std::vector<u64>> sendSet;
    std::vector<block> sendListKey;
    std::vector<block> sendListVal;
    std::vector<block> rand_S(n, ZeroBlock);
    std::vector<block> r_S(n * d);
    prng.get(r_S.data(), r_S.size());

    sampleDataGlobalDisj(sendSet, recvSet, delta, n, d, prng);

    for (u64 i = 0; i < interSize; i++) {
        u64 averageDiff = 0;
        if (lp == 1) {
            averageDiff = floor(delta * 1.0 / d);
        } else if (lp == 2) {
            averageDiff = floor(delta / std::sqrt(d));
        }
        for (size_t j = 0; j < d; j++) {
            recvSet[i][j] = sendSet[i][j] + prng.get<u64>() % (2 * averageDiff + 1) - averageDiff;
        }
    }

    oc::Timer time;
    time.setTimePoint("begin");

    for (size_t i = 0; i < recvSet.size(); i++) {
        for (int j = 0; j < d; j++) {
            auto prefixes = getPrefixSet(recvSet[i][j], prefixLenMapNaive.at(2 * delta));
            for (auto prefix : prefixes) {
                recvListKey.push_back(prefix ^ block(j << 4, 0));
            }
        }
    }

    for (size_t i = 0; i < sendSet.size(); i++) {
        for (int j = 0; j < d; j++) {
            auto prefixes = getIntervalPrefixSet(sendSet[i][j] - delta, sendSet[i][j] + delta, prefixLenMapNaive.at(2 * delta));
            for (auto prefix : prefixes) {
                sendListKey.push_back(prefix ^ block(j << 4, 0));
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

    {
        std::thread recvThr([&]() {
            OpprfRevcer recver(n * d * prefixNum, n * d * prefixLen);
            recver.setTimer(time);
            recver.recv(recvListKey, g, socket[0]);
        });

        std::thread sendThr([&]() {
            OpprfSender sender(n * d * prefixNum, n * d * prefixLen);
            sender.setTimer(time);
            sender.send(sendListKey, sendListVal, socket[1]);
        });

        recvThr.join();
        sendThr.join();

        time.setTimePoint("first OPPRF done");

        std::thread recvThrReverse([&]() {
            std::vector<block> inputKeys(g.size());
            std::vector<block> outputVals(g.size());
            for (u64 i = 0; i < g.size(); i++) {
                inputKeys[i] = g[i];
                outputVals[i] = s_R[i / prefixNum];
            }

            OpprfSender sender(n * d, n * d * prefixNum);
            sender.setTimer(time);
            sender.send(inputKeys, outputVals, socket[0]);

            for (int i = 0; i < rand_R.size(); i++) {
                for (int j = 0; j < d; j++) {
                    rand_R[i] += s_R[i * d + j];
                }
            }
        });

        std::thread sendThrReverse([&]() {
            std::vector<block> out(r_S.size());

            OpprfRevcer recver(n * d, n * d * prefixNum);
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

        for (u64 i = 0; i < rand_R.size(); i++) {
            if (rand_R[i] == rand_S[i]) {
                std::cout << "match at index " << i << std::endl;
            }
        }

        std::vector<u64> x;
        std::vector<u64> y;
        std::vector<u32> mMapping;
        std::vector<u32> cuckooMap(n);

        std::thread matchSendThr([&]() {
            auto sIdx = SimpleIndex{};
            auto params = oc::CuckooIndex<>::selectParams(rand_R.size(), 40, 0, 3);
            auto numBins = params.numBins();
            sIdx.init(numBins, rand_S.size(), 40, 3);
            sIdx.insertItems(rand_S, ZeroBlock);

            auto r = std::vector<u64>(numBins * d);
            prng.get(r.data(), r.size());
            x = r;

            std::vector<block> Ty(rand_S.size() * 3);
            auto ctr = 0;
            auto Tv = oc::Matrix<u8>(Ty.size(), d * sizeof(u64));

            for (u64 i = 0; i < numBins; ++i) {
                auto bin = sIdx.mBins[i];
                auto size = sIdx.mBinSizes[i];

                auto mask = r.data() + i * d;
                for (u64 p = 0; p < size; ++p) {
                    auto j = bin[p].hashIdx();
                    auto b = bin[p].idx();
                    Ty[ctr] = block(0, j) ^ (rand_S[b] << 2);
                    auto ptr = Tv.data(ctr);
                    for (u64 k = 0; k < d; ++k) {
                        auto val = sendSet[b][k] - mask[k];
                        std::memcpy(ptr + k * sizeof(u64), &val, sizeof(u64));
                    }
                    ctr++;
                }
            }
            PRNG prng(sysRandomSeed());

            RsOpprfSender sender;
            sender.setTimer(time);
            coproto::sync_wait(sender.send(numBins, Ty, Tv, prng, 1, socket[1]));
        });

        std::thread matchRecvThr([&]() {
            auto cuckoo = oc::CuckooIndex<>{};
            cuckoo.init(rand_R.size(), 40, 0, 3);
            cuckoo.insert(rand_R, ZeroBlock);

            std::vector<block> Tx(cuckoo.mNumBins);
            auto r = oc::Matrix<u8>(Tx.size(), d * sizeof(u64));
            u64 numBins = cuckoo.mBins.size();
            mMapping.resize(numBins);

            for (u64 i = 0; i < numBins; ++i) {
                auto &bin = cuckoo.mBins[i];
                if (bin.isEmpty() == false) {
                    auto j = bin.hashIdx();
                    auto b = bin.idx();

                    Tx[i] = block(0, j) ^ (rand_R[b] << 2);
                    mMapping[i] = b;
                    cuckooMap[b] = i;
                } else {
                    Tx[i] = block(i, 0);
                    mMapping[i] = rand_R.size(); // invalid
                }
            }

            PRNG prng(sysRandomSeed());

            RsOpprfReceiver recver;
            recver.setTimer(time);
            coproto::sync_wait(recver.receive(rand_R.size() * 3, Tx, r, prng, 1, socket[0]));

            auto r_ptr = reinterpret_cast<uint64_t *>(r.data());
            y = std::vector<u64>(numBins * d);
            for (auto i = 0; i < y.size(); i++) {
                y[i] = r_ptr[i];
            }

            for (u64 i = 0; i < numBins; ++i) {
                auto b = mMapping[i];
                if (b != rand_R.size()) {
                    for (u64 k = 0; k < d; ++k) {
                        y[i * d + k] -= recvSet[b][k];
                    }
                } // else empty bin
            }
        });

        matchSendThr.join();
        matchRecvThr.join();

        for (u64 i = 0; i < interSize; i++) {
            auto b = cuckooMap[i];
            for (u64 j = 0; j < d; j++) {
                std::cout << int64_t(x[b * d + j] + y[b * d + j]) << ", ";
            }
            std::cout << "at index " << b << " ";
            std::cout << std::endl;
        }

        time.setTimePoint("matching done");

        std::vector<u8> choiceBit(x.size() / d, 0);

        if (lp == 1) {
            normL1(x, y, choiceBit, d, delta, socket);
        } else if (lp == 2) {
            normL2(x, y, choiceBit, d, delta, socket);
        } else {
            throw std::runtime_error("lp type not supported");
        }

        time.setTimePoint("norm done");

        for (u64 i = 0; i < choiceBit.size(); i++) {
            if ((choiceBit[i] & 1)) {
                std::cout << "fianal match at index " << i << std::endl;
            }
        }

        // todo: OT

        time.setTimePoint("OT done");
    }

    auto e = time.setTimePoint("all done");

    if (verbose) {
        std::cout << time << std::endl;
    }

    std::cout << (socket[0].bytesReceived() + socket[0].bytesSent()) * 1.0 / 1024 / 1024 << " MB" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(e - s).count() / double(1000 * 1000) << " seconds" << std::endl;
}

void normL1(std::vector<u64> x, std::vector<u64> y, std::vector<u8> &choiceBit, u64 d, int delta, std::array<coproto::AsioSocket, 2> &chl)
{
    int bitsLen = 56;
    u64 mask = (1ull << bitsLen) - 1;
    int extBitlen = roundUpTo(bitsLen + log2ceil(d), 8);

    auto n = x.size() / d;

    for (auto i = 0; i < x.size(); ++i) {
        x[i] = x[i] & mask;
        y[i] = y[i] & mask;
    }

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
            abs[i] = abs[i] & mask;
        }

        std::vector<u64> dis(n, 0);
        std::vector<u64> ext_dis(n, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += abs[i * d + j];
            }
            dis[i] = dis[i] & mask;
            ext_dis[i] = dis[i];
            dis[i] = dis[i] - (1ULL << bitsLen);
        }

        MillionaireProtocolSender sender2(dis.size(), extBitlen);
        MuxSender mux2(dis.size(), &chl[1]);

        std::vector<u8> carry(dis.size());
        std::vector<u64> mod(dis.size(), 0);
        std::vector<u64> res_mux(dis.size());
        sender2.drelu(carry.data(), dis.data(), chl[1]);

        mux2.muxA(carry, mod, res_mux);

        for (u64 i = 0; i < dis.size(); ++i) {
            ext_dis[i] = ext_dis[i] - res_mux[i];
        }

        std::vector<u8> resBits1(dis.size());

        for (u64 i = 0; i < dis.size(); ++i) {
            ext_dis[i] = (u64)delta - ext_dis[i];
        }

        sender2.drelu(resBits1.data(), ext_dis.data(), chl[1]);

        coproto::sync_wait(chl[1].send(resBits1));
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
            abs[i] = abs[i] & mask;
        }

        std::vector<u64> abs_recv(x.size());

        std::vector<u64> dis(n, 0);
        std::vector<u64> ext_dis(n, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += abs[i * d + j];
            }
            dis[i] = dis[i] & mask;
            ext_dis[i] = dis[i];
        }

        MillionaireProtocolRecver recver2(dis.size(), extBitlen);
        MuxRecver mux2(dis.size(), &chl[0]);

        std::vector<u8> carry(dis.size());
        std::vector<u64> mod(dis.size(), (1ULL << bitsLen));
        std::vector<u64> res_mux(dis.size());

        recver2.drelu(carry.data(), dis.data(), chl[0]);
        mux2.muxA(carry, mod, res_mux);

        std::vector<u8> carry_recv(dis.size());

        for (u64 i = 0; i < dis.size(); ++i) {
            ext_dis[i] = ext_dis[i] - res_mux[i];
        }

        std::vector<u8> resBits0(dis.size());

        for (u64 i = 0; i < dis.size(); ++i) {
            ext_dis[i] = (u64)0 - ext_dis[i];
        }

        recver2.drelu(resBits0.data(), ext_dis.data(), chl[0]);

        std::vector<u8> resBits1(dis.size());
        coproto::sync_wait(chl[0].recv(resBits1));
        for (u64 i = 0; i < dis.size(); ++i) {
            choiceBit[i] = (resBits1[i] ^ resBits0[i]) & 1;
        }
    });

    cmpSendThr.join();
    cmpRecvThr.join();

    // std::vector<u64> x_copy(x.size());
    // std::vector<u64> y_copy(y.size());

    // for (auto i = 0; i < x.size(); ++i) {
    //     if (std::abs(int64_t(x[i] + y[i])) < delta) {
    //         std::cout << "match at index " << i << std::endl;
    //     }
    // }

    // u64 delta_p = delta;
    // int prefixLenIfmat = static_cast<int>(std::ceil(std::log2(delta_p * 2 + 1)));
    // PRNG prng(sysRandomSeed());
    // auto n = x.size() / d;

    // std::vector<u64> dis_s(x.size() / d, 0);
    // std::vector<u64> dis_r(y.size() / d, 0);

    // std::thread cmpSendThr([&]() {
    //     MillionaireProtocolSender sender(x.size(), 32);

    //     std::vector<u8> cmpShare(x.size());
    //     sender.drelu(cmpShare.data(), x_copy.data(), chl[1]);

    //     MuxSender mux(x.size(), &chl[1]);

    //     std::vector<u64> res(x.size());
    //     std::vector<u64> x_vec(x.begin(), x.end());

    //     mux.muxA(cmpShare, x_vec, res);

    //     std::vector<u64> abs(x.size());
    //     for (u64 i = 0; i < abs.size(); ++i) {
    //         abs[i] = 2 * res[i] - x[i];
    //     }

    //     for (u64 i = 0; i < dis_s.size(); ++i) {
    //         for (u64 j = 0; j < d; ++j) {
    //             dis_s[i] += abs[i * d + j];
    //         }
    //     }

    //     std::vector<block> prefixS;

    //     for (u64 i = 0; i < n; i++) {
    //         auto pre = getIntervalPrefix(0ULL - dis_s[i], delta_p - dis_s[i]);
    //         for (auto &p : pre) {
    //             p = p ^ block(i << 32, 0);
    //             prefixS.push_back(p);
    //         }
    //     }
    //     while (prefixS.size() != (n * prefixLenIfmat)) {
    //         prefixS.push_back(prng.get<block>());
    //     }

    //     PEqTSender eqSend(n * prefixLenIfmat, 1, false, &chl[1]);

    //     eqSend.eq(prefixS);
    // });

    // std::thread cmpRecvThr([&]() {
    //     MillionaireProtocolRecver recver(y.size(), 32);

    //     std::vector<u8> cmpShare(y.size());
    //     recver.drelu(cmpShare.data(), y_copy.data(), chl[0]);

    //     MuxRecver mux(y.size(), &chl[0]);
    //     std::vector<u64> res(y.size());
    //     std::vector<u64> y_vec(y.begin(), y.end());

    //     mux.muxA(cmpShare, y_vec, res);

    //     std::vector<u64> abs(y.size());

    //     for (u64 i = 0; i < abs.size(); ++i) {
    //         abs[i] = 2 * res[i] - y[i];
    //     }

    //     for (u64 i = 0; i < dis_r.size(); ++i) {
    //         for (u64 j = 0; j < d; ++j) {
    //             dis_r[i] += abs[i * d + j];
    //         }
    //     }

    //     std::vector<block> prefixR;

    //     for (u64 i = 0; i < n; i++) {
    //         auto pre = getPrefix(dis_r[i], prefixLenIfmat);
    //         for (auto &p : pre) {
    //             p = p ^ block(i << 32, 0);
    //             prefixR.push_back(p);
    //         }
    //     }

    //     PEqTRecver eqRecv(n * prefixLenIfmat, 1, false, &chl[0]);

    //     std::vector<u64> intersection;

    //     eqRecv.eq(prefixR, intersection);

    //     for (auto &v : intersection) {
    //         choiceBit[v / prefixLenIfmat] = 1;
    //     }
    // });

    // cmpSendThr.join();
    // cmpRecvThr.join();

    // for (auto i = 0; i < dis_r.size(); ++i) {
    //     if (std::abs(int64_t(dis_r[i] + dis_s[i])) < delta) {
    //         std::cout << "dis match at index " << i << std::endl;
    //     }
    // }
}

void normL2(oc::span<u64> x, oc::span<u64> y, std::vector<u8> &choiceBit, u64 d, int delta, std::array<coproto::AsioSocket, 2> &chl)
{
    u64 delta_p = delta * delta;
    int prefixLenIfmat = static_cast<int>(std::ceil(std::log2(delta_p * 2 + 1)));
    PRNG prng(sysRandomSeed());
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

        std::vector<block> prefixS;

        for (u64 i = 0; i < n; i++) {
            auto pre = getIntervalPrefix(0ULL - dis[i], delta_p - dis[i]);
            for (auto &p : pre) {
                p = p ^ block(i << 32, 0);
                prefixS.push_back(p);
            }
        }
        while (prefixS.size() != (n * prefixLenIfmat)) {
            prefixS.push_back(prng.get<block>());
        }

        PEqTSender eqSend(n * prefixLenIfmat, 1, false, &chl[1]);

        eqSend.eq(prefixS);
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

        std::vector<block> prefixR;

        for (u64 i = 0; i < n; i++) {
            auto pre = getPrefix(dis[i], prefixLenIfmat);
            for (auto &p : pre) {
                p = p ^ block(i << 32, 0);
                prefixR.push_back(p);
            }
        }

        PEqTRecver eqRecv(n * prefixLenIfmat, 1, false, &chl[0]);

        std::vector<u64> intersection;

        eqRecv.eq(prefixR, intersection);

        for (auto &v : intersection) {
            choiceBit[v / prefixLenIfmat] = 1;
        }
    });

    cmpSendThr.join();
    cmpRecvThr.join();
}