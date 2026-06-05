#pragma once
#include <coproto/Socket/Socket.h>
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <random>
#include <vector>
#include "utils.h"

inline std::random_device rd2;
inline std::mt19937 global_built_in_prg2(rd2());

// using namespace oc;
// using namespace std;

using coproto::Socket;

void genPermutation(u32 size, std::vector<u32> &pi);

void permute(std::vector<u32> pi, std::vector<block> &data);

void senderROT(u32 numOTs, oc::PRNG &prng, oc::BitVector &send0, oc::BitVector &send1, coproto::Socket &chl, u32 numThreads);

void receiverROT(u32 numOTs, oc::PRNG &prng, oc::BitVector &bitV, oc::BitVector &recv, coproto::Socket &chl, u32 numThreads);

void perm_cons(u8 nn, std::vector<u32> &in, std::vector<u32> pi, oc::BitVector &switch_network);
// void Opermute(u32 idx, BitVector &network_choice, BitVector &input, BitVector &output, Socket &chl, u32 numThreads);
void Opermute(u32 idx, oc::BitVector &input, oc::BitVector &output, std::vector<u32> &pi, coproto::Socket &chl, u32 numThreads);

void Opermute1(u32 idx, oc::BitVector &input, oc::BitVector &output, std::vector<u32> &pi, coproto::Socket &chl, u32 numThreads);