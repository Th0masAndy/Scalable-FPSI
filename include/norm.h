#pragma once

#include <array>
#include <coproto/Socket/AsioSocket.h>
#include <cryptoTools/Common/Defines.h>
#include <vector>
#include "utils.h"

void normL0(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl);

void normL1(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl);

void normL2(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl);
