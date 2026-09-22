#include "opprf.h"
#include <coproto/Common/macoro.h>

using namespace oc;
using namespace volePSI;

OpprfSender::OpprfSender(size_t _decodeSize, size_t _kvSize) : decodeSize(_decodeSize), kvSize(_kvSize)
{
    prng = PRNG(oc::sysRandomSeed());
    setMultType(type);
}

void OpprfSender::send(std::vector<block> &keys, std::vector<block> &values, Socket &chl)
{
    coproto::sync_wait(RsOpprfSender::send(decodeSize, keys, values, prng, 1, chl));
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
