#include <coproto/Socket/AsioSocket.h>
#include <cryptoTools/Common/Timer.h>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <thread>
#include <vector>
#include "cmp.h"
#include "opprf.h"

void opprf_test()
{
    OpprfSender sender(1024, 1024);
    OpprfRevcer recver(1024, 1024);

    std::vector<oc::block> keys(1024);
    std::vector<oc::block> values(1024);

    PRNG prng(sysRandomSeed());

    prng.get(keys.data(), keys.size());
    prng.get(values.data(), values.size());

    std::vector<oc::block> decodeKeys(keys.begin(), keys.begin() + recver.decodeSize);
    std::vector<oc::block> decodeValues(recver.decodeSize);

    auto socket = coproto::LocalAsyncSocket::makePair();

    auto sendThr = std::thread([&]() { sender.send(keys, values, socket[0]); });

    auto recvThr = std::thread([&]() { recver.recv(decodeKeys, decodeValues, socket[1]); });

    sendThr.join();
    recvThr.join();

    for (auto i = 0; i < recver.decodeSize; ++i) {
        if (decodeValues[i] != values[i]) {
            throw std::runtime_error("opprf test error");
        }
    }
}

void cmp_test()
{
    u64 n = 128 + (1 << 14);

    MillionaireProtocolRecver recver(n, 64, 4);
    MillionaireProtocolSender sender(n, 64, 4);

    auto socket = coproto::AsioSocket::makePair();

    oc::PRNG prng(oc::sysRandomSeed());

    std::vector<u64> data0(n);
    std::vector<u64> data1(n);
    std::vector<u8> outs0(n);
    std::vector<u8> outs1(n);

    for (u64 i = 0; i < n; ++i) {
        data0[i] = prng.get<u64>();
        data1[i] = prng.get<u64>();
    }

    oc::Timer time;
    auto s = time.setTimePoint("begin cmp");

    std::thread recvThr([&]() { recver.compare(outs0.data(), data0.data(), socket[0]); });

    std::thread sendThr([&]() { sender.compare(outs1.data(), data1.data(), socket[1]); });

    recvThr.join();
    sendThr.join();

    auto e = time.setTimePoint("end cmp");

    int correct = 0;
    for (u64 i = 0; i < n; ++i) {
        bool gt = data1[i] > data0[i];
        if (gt == ((outs1[i] ^ outs0[i]) & 1))
            correct++;
    }
    std::cout << "correct: " << correct << " / " << n << std::endl;

    auto comm = socket[0].bytesReceived() + socket[0].bytesSent();
    auto comp = std::chrono::duration_cast<std::chrono::microseconds>(e - s).count();
    std::cout << "average: " << comm / (n * 1.0) << " bytes" << " " << comp / (n * 1.0) << " microseconds" << std::endl;
}