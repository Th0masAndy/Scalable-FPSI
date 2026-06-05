#pragma once

#include <coproto/Socket/Socket.h>
#include <cstdint>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtSender.h>
#include <vector>
#include "utils.h"

class MuxSender {
public:
    MuxSender(uint64_t num_, coproto::Socket *socket_);
    ~MuxSender();

    void muxA(std::vector<u8> &b0, std::vector<u64> &v0, std::vector<u64> &res0);

    void muxA(std::vector<u8> &b0, std::vector<u64> &v0, std::vector<u64> &res0, int bitsLen);

    uint64_t num;

private:
    coproto::Socket *socket;
    osuCrypto::SilentOtExtSender *sender;
    osuCrypto::SilentOtExtReceiver *recver;
    osuCrypto::PRNG *prng;
};

class MuxRecver {
public:
    MuxRecver(uint64_t num_, coproto::Socket *socket_);
    ~MuxRecver();

    void muxA(std::vector<u8> &b1, std::vector<u64> &v1, std::vector<u64> &res1);

    void muxA(std::vector<u8> &b1, std::vector<u64> &v1, std::vector<u64> &res1, int bitsLen);

    uint64_t num;

private:
    coproto::Socket *socket;
    osuCrypto::SilentOtExtSender *sender;
    osuCrypto::SilentOtExtReceiver *recver;
    osuCrypto::PRNG *prng;
};