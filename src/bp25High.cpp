#include <coproto/Socket/AsioSocket.h>
#include <cryptoTools/Common/CuckooIndex.h>
#include "opprf.h"
#include "sparsehash/dense_hash_map"
#include "utils.h"

using namespace oc;

void sampleDataBp25High(std::vector<std::vector<u64>> &sendSet, std::vector<std::vector<u64>> &recvSet, int delta, u64 n, size_t d, PRNG &prng)
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

void bp25HighLp(const oc::CLP &cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    size_t d = cmd.getOr("d", 2);
    int delta = cmd.getOr("delta", 2);
    int verbose = cmd.getOr("v", 0);

    int lp = cmd.getOr("p", 2);

    u64 delta_p = std::pow(delta, lp);

    int prefixLenIfmat = static_cast<int>(std::ceil(std::log2(delta_p + 1)));

    int numTry = cmd.getOr("try", 1);

    int interSize = cmd.getOr("nn", 4);

    PRNG prng(sysRandomSeed());
    std::vector<std::vector<u64>> recvSet;
    std::vector<block> recvListKey;
    std::vector<block> recvListVal;
    std::vector<u64> rand_R(n, 0);
    std::vector<block> rand_R_prefix;
    std::vector<block> g(n * d);

    std::vector<std::vector<u64>> sendSet;
    std::vector<block> sendListKey;
    std::vector<block> sendListVal;
    std::vector<u64> rand_S(n, 0);
    std::vector<block> rand_S_prefix;
    std::vector<u64> r_S(n * d);
    prng.get(r_S.data(), r_S.size());
    for (int i = 0; i < n; i++) {
        u64 sum = 0;
        for (int j = 0; j < d; j++) {
            sum += r_S[i * d + j];
        }
        rand_S[i] = sum;
    }

    sampleDataBp25High(sendSet, recvSet, delta, n, d, prng);

    oc::Timer time;
    time.setTimePoint("begin");

    for (size_t i = 0; i < recvSet.size(); i++) {
        for (int j = 0; j < d; j++) {
            recvListKey.push_back(block(j, recvSet[i][j]));
        }
    }

    for (size_t i = 0; i < sendSet.size(); i++) {
        for (int j = 0; j < d; j++) {
            for (int k = -delta; k <= delta; k++) {
                sendListKey.push_back(block(j, sendSet[i][j] + k));
                if (lp == 1) {
                    sendListVal.push_back(block(0, r_S[i * d + j] + std::abs(k)));
                } else if (lp == 2) {
                    sendListVal.push_back(block(0, r_S[i * d + j] + k * k));
                } else {
                    throw std::runtime_error("unsupported lp");
                }
            }
        }
    }

    while (sendListKey.size() < n * d * (2 * delta + 1)) {
        sendListKey.push_back(prng.get<block>());
        sendListVal.push_back(prng.get<block>());
    }

    auto s = time.setTimePoint("preprocess done");

    auto socket = coproto::AsioSocket::makePair();

    std::thread recvThr([&]() {
        OpprfRevcer recver(n * d, n * d * (2 * delta + 1));
        recver.setTimer(time);
        recver.recv(recvListKey, g, socket[0]);
    });

    std::thread sendThr([&]() {
        OpprfSender sender(n * d, n * d * (2 * delta + 1));
        sender.setTimer(time);
        sender.send(sendListKey, sendListVal, socket[1]);
    });

    recvThr.join();
    sendThr.join();

    time.setTimePoint("first OPPRF done");

    std::thread sendMaskThr([&]() {
        PRNG prng(sysRandomSeed());
        auto len = d * sizeof(u64);
        auto numBlocks = (len + sizeof(block) - 1) / sizeof(block);

        std::vector<u8> buffer(len * prefixLenIfmat * n);
        for (int i = 0; i < n; i++) {
            auto prefixes = getIntervalPrefix(rand_S[i], rand_S[i] + delta_p);
            for (auto prefix : prefixes) {
                rand_S_prefix.push_back(prefix);
            }
            for (int i = 0; i < prefixLenIfmat - prefixes.size(); i++) {
                rand_S_prefix.push_back(prng.get<block>());
            }

            for (auto p = 0; p < prefixLenIfmat; p++) {
                prng.SetSeed(rand_S_prefix[p], numBlocks);
                std::vector<u8> tmp(len);
                prng.get(tmp.data(), len);
                for (int j = 0; j < len; j++) {
                    buffer[(i * prefixLenIfmat + p) * len + j] = tmp[j] ^ ((u8 *)&sendSet[i][0])[j];
                }
            }
        }

        Hash(rand_S_prefix);

        coproto::sync_wait(socket[1].send(buffer));
        coproto::sync_wait(socket[1].send(rand_S_prefix));
    });

    std::thread recvMaskThr([&]() {
        PRNG prng(sysRandomSeed());
        auto len = d * sizeof(u64);
        auto numBlocks = (len + sizeof(block) - 1) / sizeof(block);
        std::vector<block> rand_S_prefix(n * prefixLenIfmat, ZeroBlock);

        std::vector<u8> buffer(len * prefixLenIfmat * n);
        coproto::sync_wait(socket[0].recv(buffer));
        coproto::sync_wait(socket[0].recv(rand_S_prefix));

        for (int i = 0; i < rand_R.size(); i++) {
            for (int j = 0; j < d; j++) {
                rand_R[i] += low(g[i * d + j]);
            }
        }

        for (auto v : rand_R) {
            auto prefixes = getPrefix(v, prefixLenIfmat);
            for (auto prefix : prefixes) {
                rand_R_prefix.push_back(prefix);
            }
        }

        Hash(rand_R_prefix);

        auto map = google::dense_hash_map<block, u64, NoHash>{};
        map.resize(rand_R_prefix.size());
        map.set_empty_key(oc::ZeroBlock);
        for (auto i = 0; i < rand_R_prefix.size(); i++) {
            map.insert({ rand_R_prefix[i], i });
        }

        for (auto i = 0; i < rand_S_prefix.size(); i++) {
            if (map.find(rand_S_prefix[i]) != map.end()) {
                std::vector<u8> tmp(len);
                prng.SetSeed(rand_S_prefix[i], numBlocks);
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

void bp25High(const oc::CLP &cmd)
{
    u64 n = cmd.getOr("n", 1ull << cmd.getOr("nn", 10));
    size_t d = cmd.getOr("d", 2);
    int delta = cmd.getOr("delta", 2);
    int verbose = cmd.getOr("v", 0);

    int numTry = cmd.getOr("try", 1);

    int interSize = cmd.getOr("nn", 4);

    PRNG prng(sysRandomSeed());
    std::vector<std::vector<u64>> recvSet;
    std::vector<block> recvListKey;
    std::vector<block> recvListVal;
    std::vector<block> rand_R(n, ZeroBlock);
    std::vector<block> g(n * d);

    std::vector<std::vector<u64>> sendSet;
    std::vector<block> sendListKey;
    std::vector<block> sendListVal;
    std::vector<block> rand_S(n, ZeroBlock);
    std::vector<block> r_S(n * d);
    prng.get(r_S.data(), r_S.size());
    for (int i = 0; i < n; i++) {
        block sum = ZeroBlock;
        for (int j = 0; j < d; j++) {
            sum ^= r_S[i * d + j];
        }
        rand_S[i] = sum;
    }

    sampleDataBp25High(sendSet, recvSet, delta, n, d, prng);

    oc::Timer time;
    time.setTimePoint("begin");

    for (size_t i = 0; i < recvSet.size(); i++) {
        for (int j = 0; j < d; j++) {
            recvListKey.push_back(block(j, recvSet[i][j]));
        }
    }

    for (size_t i = 0; i < sendSet.size(); i++) {
        for (int j = 0; j < d; j++) {
            for (int k = -delta; k <= delta; k++) {
                sendListKey.push_back(block(j, sendSet[i][j] + k));
                sendListVal.push_back(r_S[i * d + j]);
            }
        }
    }

    while (sendListKey.size() < n * d * (2 * delta + 1)) {
        sendListKey.push_back(prng.get<block>());
        sendListVal.push_back(prng.get<block>());
    }

    auto s = time.setTimePoint("preprocess done");

    auto socket = coproto::AsioSocket::makePair();

    std::thread recvThr([&]() {
        OpprfRevcer recver(n * d, n * d * (2 * delta + 1));
        recver.setTimer(time);
        recver.recv(recvListKey, g, socket[0]);
    });

    std::thread sendThr([&]() {
        OpprfSender sender(n * d, n * d * (2 * delta + 1));
        sender.setTimer(time);
        sender.send(sendListKey, sendListVal, socket[1]);
    });

    recvThr.join();
    sendThr.join();

    time.setTimePoint("first OPPRF done");

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

        for (int i = 0; i < rand_R.size(); i++) {
            for (int j = 0; j < d; j++) {
                rand_R[i] += g[i * d + j];
            }
        }

        Hash(rand_R);

        auto map = google::dense_hash_map<block, u64, NoHash>{};
        map.resize(rand_R.size());
        map.set_empty_key(oc::ZeroBlock);
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