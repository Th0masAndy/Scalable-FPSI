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