#include "mul.h"
#include <vector>
#include "utils.h"

using namespace osuCrypto;

MulSender::MulSender(uint64_t num_, coproto::Socket *socket_, int bitsLen_) : num(num_), bitsLen(bitsLen_), socket(socket_)
{
    sender = new osuCrypto::SilentOtExtSender();
    sender->configure(num * bitsLen);
    sender->mMultType = type;

    prng = new PRNG(ZeroBlock);
}

MulSender::~MulSender()
{
    delete sender;
    delete prng;
}

void MulSender::mul(std::vector<uint64_t> &inputs, std::vector<uint64_t> &val)
{
    coproto::sync_wait(sender->genSilentBaseOts(*prng, *socket));

    u64 numOts = num * 64;
    std::vector<std::array<block, 2>> messages(numOts);

    coproto::sync_wait(sender->send(messages, *prng, *socket));

    std::vector<u64> correctMessages(numOts);

    std::vector<u8> compressedBits;

    for (int i = 0; i < numOts; i++) {
        int shift = i % 64;

        u64 mask = low(messages[i][0]) + inputs[i / 64];
        correctMessages[i] = (low(messages[i][1]) ^ mask) << shift >> shift;
        for (int j = 0; j < 8 - shift / 8; j++) {
            compressedBits.push_back(static_cast<uint8_t>(correctMessages[i] >> (8 * j)) & 0xFF);
        }
    }

    for (int i = 0; i < numOts; i++) {
        int shift = i % 64;
        val[i / 64] += 0 - ((low(messages[i][0]) << shift >> shift) << shift);
    }

    coproto::sync_wait(socket->send(compressedBits));
}

MulRecver::MulRecver(uint64_t num_, coproto::Socket *socket_, int bitsLen_) : num(num_), bitsLen(bitsLen_), socket(socket_)
{
    receiver = new osuCrypto::SilentOtExtReceiver();
    receiver->configure(num * bitsLen);
    receiver->mMultType = type;

    prng = new PRNG(OneBlock);
}

MulRecver::~MulRecver()
{
    delete receiver;
    delete prng;
}

void MulRecver::mul(std::vector<uint64_t> &inputs, std::vector<uint64_t> &val)
{
    coproto::sync_wait(receiver->genSilentBaseOts(*prng, *socket));

    u64 numOts = num * 64;
    std::vector<u8> bytes(numOts / 8);

    for (int i = 0; i < num; i++) {
        u64 lowbits = inputs[i];
        for (int j = 0; j < 8; j++) {
            bytes[i * 8 + j] = (lowbits >> (8 * j)) & 0xFF;
        }
    }

    std::vector<block> messages(numOts);

    BitVector choiceBit(bytes.data(), numOts);

    coproto::sync_wait(receiver->receive(choiceBit, messages, *prng, *socket));

    std::vector<u64> correctMessages(numOts);

    std::vector<u8> compressedBits;

    coproto::sync_wait(socket->recvResize(compressedBits));

    u64 offset = 0;
    for (int i = 0; i < numOts; i++) {
        u64 msg = 0;
        int shift = i % 64;
        for (int j = 0; j < 8 - shift / 8; j++) {
            msg |= u64(compressedBits[offset++]) << (8 * j);
        }
        if (choiceBit[i] & 1) {
            val[i / 64] += ((low(messages[i]) ^ msg) << shift >> shift) << (i % 64);
        } else {
            val[i / 64] += (low(messages[i]) << shift >> shift) << (i % 64);
        }
    }
}