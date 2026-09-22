#pragma once
#include <coproto/Socket/Socket.h>
#include <cryptoTools/Common/BitVector.h>
#include <vector>
#include "utils.h"

void ssPEQT(u32 idx, std::vector<block> &input, oc::BitVector &out, coproto::Socket &chl, u32 numThreads);
