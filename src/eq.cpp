#include "eq.h"
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Common/block.h>
#include <vector>
#include <volePSI/GMW/Circuit.h>
#include <volePSI/GMW/Gmw.h>

using namespace volePSI;
using namespace oc;

PEqTSender::PEqTSender(uint64_t num_, uint64_t numThreads_, bool noCompress_, coproto::Socket *socket_)
    : num(num_), numThreads(numThreads_), noCompress(noCompress_), socket(socket_)
{
    sender = new RsPsiSender();
    sender->init(num, num, 40, oc::ZeroBlock, false, numThreads);

    auto type = oc::DefaultMultType;

    sender->setMultType(type);

    if (noCompress) {
        sender->mCompress = false;
        sender->mMaskSize = sizeof(block);
    }
}

void PEqTSender::eq(std::vector<block> &sendSet)
{
    coproto::sync_wait(sender->run(sendSet, *socket));
}

PEqTSender::~PEqTSender()
{
    delete sender;
}

PEqTRecver::PEqTRecver(uint64_t num_, uint64_t numThreads_, bool noCompress_, coproto::Socket *socket_)
    : num(num_), numThreads(numThreads_), noCompress(noCompress_), socket(socket_)
{
    recver = new RsPsiReceiver();
    recver->init(num, num, 40, oc::ZeroBlock, false, numThreads);

    auto type = oc::DefaultMultType;

    recver->setMultType(type);

    if (noCompress) {
        recver->mCompress = false;
        recver->mMaskSize = sizeof(block);
    }
}

void PEqTRecver::eq(std::vector<block> &recvSet, std::vector<u64> &intersection)
{
    coproto::sync_wait(recver->run(recvSet, *socket));
    intersection = recver->mIntersection;
}

PEqTRecver::~PEqTRecver()
{
    delete recver;
}

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
