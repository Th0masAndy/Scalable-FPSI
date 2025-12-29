#include "opprf.h"
#include <coproto/Common/macoro.h>
#include <cryptoTools/Common/block.h>
#include <libOTe/TwoChooseOne/ConfigureCode.h>
#include <vector>
#include <volePSI/Defines.h>
#include <volePSI/Paxos.h>
#include "utils.h"

OpprfSender::OpprfSender(size_t _decodeSize, size_t _kvSize) : decodeSize(_decodeSize), kvSize(_kvSize)
{
    prng = PRNG(oc::sysRandomSeed());
    setMultType(type);
}

void OpprfSender::send(std::vector<block> &keys, std::vector<block> &values, Socket &chl)
{
    coproto::sync_wait(RsOpprfSender::send(decodeSize, keys, values, prng, 1, chl));
}

void OpprfSender::send(std::vector<block> encodings, Socket &chl)
{
    coproto::sync_wait(chl.send(ZeroBlock));

    coproto::sync_wait(mOprfSender.send(decodeSize, prng, chl, 1));

    coproto::sync_wait(chl.send(encodings));
}

void OpprfSender::encode(std::vector<block> &keys, std::vector<block> &values, std::vector<block> &encodings)
{
    mPaxos.init(kvSize, 1 << 14, 3, 40, PaxosParam::GF128, ZeroBlock);

    encodings.resize(mPaxos.size());

    std::vector<block> diff(values.size());

    mOprfSender.mPaxos.init(decodeSize, 1 << 14, 3, 40, PaxosParam::GF128, ZeroBlock);
    mOprfSender.eval(keys, diff, 1);

    mPaxos.solve<block>(keys, diff, encodings, &prng, 1);
}

OpprfRevcer::OpprfRevcer(size_t _decodeSize, size_t _kvSize) : decodeSize(_decodeSize), kvSize(_kvSize)
{
    prng = PRNG(oc::sysRandomSeed());
    setMultType(type);
}

void OpprfRevcer::recv(std::vector<block> &keys, std::vector<block> &outputs, Socket &chl)
{
    coproto::sync_wait(RsOpprfReceiver::receive(kvSize, keys, outputs, prng, 1, chl));
}
