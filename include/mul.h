#pragma once
#include <coproto/Socket/Socket.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cstdint>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtSender.h>
#include <vector>

class MulSender {
public:
    MulSender(uint64_t num_, coproto::Socket *socket_, int bitsLen_ = 64);
    ~MulSender();
    void mul(std::vector<uint64_t> &blk, std::vector<uint64_t> &val);

    uint64_t num;
    int bitsLen;

private:
    coproto::Socket *socket;
    osuCrypto::SilentOtExtSender *sender;
    osuCrypto::PRNG *prng;
};

class MulRecver {
public:
    MulRecver(uint64_t num_, coproto::Socket *socket_, int bitsLen_ = 64);
    ~MulRecver();
    void mul(std::vector<uint64_t> &blk, std::vector<uint64_t> &val);

    uint64_t num;
    int bitsLen;

private:
    coproto::Socket *socket;
    osuCrypto::SilentOtExtReceiver *receiver;
    osuCrypto::PRNG *prng;
};