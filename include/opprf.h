#pragma once

#include <coproto/Socket/Socket.h>
#include <cryptoTools/Common/block.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cstddef>
#include <vector>
#include "opprf_impl.h"
#include "utils.h"

class LargeOpprfSender : public fpsi::LargeOpprfSender {
public:
    LargeOpprfSender(size_t _decodeSize, size_t _kvSize);

    size_t decodeSize;
    size_t kvSize;
    oc::PRNG prng;

    void send(std::vector<block> &keys, std::vector<block> &values, coproto::Socket &chl);
};

class LargeOpprfRevcer : public fpsi::LargeOpprfReceiver {
public:
    LargeOpprfRevcer(size_t _decodeSize, size_t _kvSize);

    size_t decodeSize;
    size_t kvSize;
    oc::PRNG prng;

    void recv(std::vector<block> &keys, std::vector<block> &outputs, coproto::Socket &chl);
};
