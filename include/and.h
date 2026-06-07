#pragma once

#include <coproto/Socket/Socket.h>
#include <cstdint>
#include <vector>
#include "utils.h"

class AndSender {
public:
    AndSender(uint64_t num_, coproto::Socket *socket_);

    void andBits(std::vector<u8> &lhs, std::vector<u8> &rhs, std::vector<u8> &out);

    uint64_t num;

private:
    coproto::Socket *socket;
};

class AndRecver {
public:
    AndRecver(uint64_t num_, coproto::Socket *socket_);

    void andBits(std::vector<u8> &lhs, std::vector<u8> &rhs, std::vector<u8> &out);

    uint64_t num;

private:
    coproto::Socket *socket;
};
