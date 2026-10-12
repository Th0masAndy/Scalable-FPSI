#include "opprf.h"
#include <coproto/Common/macoro.h>

using namespace oc;
using namespace volePSI;

LargeOpprfSender::LargeOpprfSender(size_t _decodeSize, size_t _kvSize) : decodeSize(_decodeSize), kvSize(_kvSize)
{
    prng = PRNG(oc::sysRandomSeed());
    setMultType(type);
}

void LargeOpprfSender::send(std::vector<block> &keys, std::vector<block> &values, Socket &chl)
{
    coproto::sync_wait(fpsi::LargeOpprfSender::send(decodeSize, keys, values, prng, 1, chl));
}

LargeOpprfRevcer::LargeOpprfRevcer(size_t _decodeSize, size_t _kvSize) : decodeSize(_decodeSize), kvSize(_kvSize)
{
    prng = PRNG(oc::sysRandomSeed());
    setMultType(type);
}

void LargeOpprfRevcer::recv(std::vector<block> &keys, std::vector<block> &outputs, Socket &chl)
{
    coproto::sync_wait(fpsi::LargeOpprfReceiver::receive(kvSize, keys, outputs, prng, 1, chl));
}
