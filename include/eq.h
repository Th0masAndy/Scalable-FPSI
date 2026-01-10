#pragma once
#include <coproto/Socket/Socket.h>
#include <sys/types.h>
#include <vector>
#include "utils.h"
#include "volePSI/RsPsi.h"

void ssPEQT(u32 idx, std::vector<block> &input, oc::BitVector &out, coproto::Socket &chl, u32 numThreads);

class PEqTSender {
public:
    PEqTSender(uint64_t num_, uint64_t numThreads_, bool noCompress_, coproto::Socket *socket_);
    ~PEqTSender();
    void eq(std::vector<block> &sendSet);

    uint64_t num;
    uint64_t numThreads;
    bool noCompress;

private:
    volePSI::RsPsiSender *sender;
    coproto::Socket *socket;
};

class PEqTRecver {
public:
    PEqTRecver(uint64_t num_, uint64_t numThreads_, bool noCompress_, coproto::Socket *socket_);
    ~PEqTRecver();
    void eq(std::vector<block> &recvSet, std::vector<u64> &intersection);

    uint64_t num;
    uint64_t numThreads;
    bool noCompress;

private:
    volePSI::RsPsiReceiver *recver;
    coproto::Socket *socket;
};