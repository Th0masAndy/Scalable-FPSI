#include <coproto/Common/macoro.h>
#include <coproto/Socket/AsioSocket.h>
#include <coproto/Socket/LocalAsyncSock.h>
#include <cryptoTools/Common/CLP.h>
#include <cryptoTools/Common/CuckooIndex.h>
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Common/Matrix.h>
#include <cryptoTools/Common/block.h>
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

void normL1(oc::span<u64> x, oc::span<u64> y, std::vector<u8> &choiceBit, u64 d, int delta, std::array<coproto::AsioSocket, 2> &chl);
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
    int interSize = cmd.getOr("nn", 4);

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

    std::cout << time << std::endl;

    std::cout << (socket[0].bytesReceived() + socket[0].bytesSent()) / 1024 / 1024 << " MB" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(e - s).count() / double(1000 * 1000) << " seconds" << std::endl;
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
    int interSize = cmd.getOr("nn", 4);

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

    oc::span<u64> x;
    oc::span<u64> y;

    std::thread matchSendThr([&]() {
        auto sIdx = SimpleIndex{};
        auto params = oc::CuckooIndex<>::selectParams(rand_R.size(), 40, 0, 3);
        auto numBins = params.numBins();
        sIdx.init(numBins, rand_S.size(), 40, 3);
        sIdx.insertItems(rand_S, ZeroBlock);

        auto r = std::vector<u64>(numBins * d);
        prng.get(r.data(), r.size());
        x = oc::span<u64>(r.data(), r.size());

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
        std::vector<u64> mMapping(numBins);

        for (u64 i = 0; i < numBins; ++i) {
            auto &bin = cuckoo.mBins[i];
            if (bin.isEmpty() == false) {
                auto j = bin.hashIdx();
                auto b = bin.idx();

                Tx[i] = block(0, j) ^ (rand_R[b] << 2);
                mMapping[i] = b;
            } else {
                Tx[i] = block(i, 0);
                mMapping[i] = rand_R.size(); // invalid
            }
        }

        PRNG prng(sysRandomSeed());

        RsOpprfReceiver recver;
        recver.setTimer(time);
        coproto::sync_wait(recver.receive(rand_R.size() * 3, Tx, r, prng, 1, socket[0]));

        y = oc::span<u64>(reinterpret_cast<uint64_t *>(r.data()), r.size() / (sizeof(u64)));

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

    auto e = time.setTimePoint("OT done");
    std::cout << time << std::endl;

    std::cout << (socket[0].bytesReceived() + socket[0].bytesSent()) / 1024 / 1024 << " MB" << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::microseconds>(e - s).count() / double(1000 * 1000) << " seconds" << std::endl;
}

void normL1(oc::span<u64> x, oc::span<u64> y, std::vector<u8> &choiceBit, u64 d, int delta, std::array<coproto::AsioSocket, 2> &chl)
{
    u64 delta_p = delta;
    int prefixLenIfmat = static_cast<int>(std::ceil(std::log2(delta_p * 2 + 1)));
    PRNG prng(sysRandomSeed());
    auto n = x.size() / d;

    std::thread cmpSendThr([&]() {
        MillionaireProtocolSender sender(x.size(), 32);

        std::vector<u8> cmpShare(x.size());
        sender.compare(cmpShare.data(), x.data(), chl[1]);

        MuxSender mux(x.size(), &chl[1]);

        std::vector<u64> res(x.size());
        std::vector<u64> x_vec(x.begin(), x.end());

        mux.muxA(cmpShare, x_vec, res);

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
        MillionaireProtocolRecver recver(y.size(), 32);

        std::vector<u8> cmpShare(y.size());
        recver.compare(cmpShare.data(), y.data(), chl[0]);

        MuxRecver mux(y.size(), &chl[0]);
        std::vector<u64> res(y.size());
        std::vector<u64> y_vec(y.begin(), y.end());

        for (u64 i = 0; i < y_vec.size(); ++i) {
            y_vec[i] = -y_vec[i];
        }

        mux.muxA(cmpShare, y_vec, res);

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