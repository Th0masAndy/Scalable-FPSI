#include <cryptoTools/Common/block.h>
#include <cstddef>
#include <vector>
#include <volePSI/Defines.h>
#include <volePSI/Paxos.h>
#include "volePSI/RsOpprf.h"

using namespace volePSI;
using namespace oc;

class OpprfSender : public RsOpprfSender {
public:
    OpprfSender(size_t _decodeSize, size_t _kvSize);

    size_t decodeSize;
    size_t kvSize;
    PRNG prng;

    void send(std::vector<block> &keys, std::vector<block> &values, Socket &chl);

    void send(std::vector<block> encodings, Socket &chl); // support offline encodings

    void encode(std::vector<block> &keys, std::vector<block> &values, std::vector<block> &encodings);
};

class OpprfRevcer : public RsOpprfReceiver {
public:
    OpprfRevcer(size_t _decodeSize, size_t _kvSize);

    size_t decodeSize;
    size_t kvSize;
    PRNG prng;

    void recv(std::vector<block> &keys, std::vector<block> &outputs, Socket &chl);
};