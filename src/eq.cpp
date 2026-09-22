#include "eq.h"
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Common/block.h>
#include <volePSI/GMW/Circuit.h>
#include <volePSI/GMW/Gmw.h>

using namespace volePSI;
using namespace oc;

void ssPEQT(u32 idx, std::vector<block> &input, oc::BitVector &out, Socket &chl, u32 numThreads)
{
    u32 numBins = input.size();
    u64 keyBitLength = 40 + oc::log2ceil(numBins);
    u64 keyByteLength = oc::divCeil(keyBitLength, 8);
    PRNG prng(sysRandomSeed());

    oc::Matrix<u8> mLabel(numBins, keyByteLength);
    for (u32 i = 0; i < numBins; ++i) {
        memcpy(&mLabel(i, 0), &input[i], keyByteLength);
    }

    // call gmw
    auto cir = volePSI::isZeroCircuit(keyBitLength);

    // volePSI::BetaCircuit cir = volePSI::isZeroCircuit(keyBitLength);
    Timer t;
    volePSI::Gmw cmp;
    cmp.setTimer(t);
    cmp.init(mLabel.rows(), cir, numThreads, idx, prng.get());

    if (idx == 1) {
        cmp.setInput(0, mLabel);
    } else {
        cmp.implSetInput(0, mLabel, mLabel.cols());
    }

    coproto::sync_wait(cmp.run(chl));

    // if (idx == 1) {
    //     std::cout << t << std::endl;
    // }

    oc::Matrix<u8> mOut;
    mOut.resize(numBins, 1);
    cmp.getOutput(0, mOut);

    // get the final output
    out.resize(numBins);
    for (u32 i = 0; i < numBins; ++i) {
        out[i] = mOut(i, 0) & 1;
    }
    return;
}
