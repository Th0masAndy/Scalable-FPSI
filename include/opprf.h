#pragma once

#include <coproto/Socket/Socket.h>
#include <cryptoTools/Common/block.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cstddef>
#include <vector>
#include "utils.h"
#include "volePSI/RsOpprf.h"

class OpprfSender : public volePSI::RsOpprfSender {
public:
    OpprfSender(size_t _decodeSize, size_t _kvSize);

    size_t decodeSize;
    size_t kvSize;
    oc::PRNG prng;

    void send(std::vector<block> &keys, std::vector<block> &values, coproto::Socket &chl);
};

class OpprfRevcer : public volePSI::RsOpprfReceiver {
public:
    OpprfRevcer(size_t _decodeSize, size_t _kvSize);

    size_t decodeSize;
    size_t kvSize;
    oc::PRNG prng;

    void recv(std::vector<block> &keys, std::vector<block> &outputs, coproto::Socket &chl);
};