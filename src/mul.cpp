#include "mul.h"
#include <algorithm>
#include <stdexcept>
#include "utils.h"

using namespace osuCrypto;

namespace {
    // Bound both OT working memory and correction messages (at most 72 MiB).
    constexpr u64 maxBatchOts = 1ull << 24;

    void checkBitLength(int bitsLen)
    {
        if (bitsLen < 8 || bitsLen > 64 || bitsLen % 8 != 0) {
            throw std::invalid_argument("multiplication requires a byte-aligned bit length between 8 and 64");
        }
    }
} // namespace

MulSender::MulSender(uint64_t num_, coproto::Socket *socket_, int bitsLen_) : num(num_), bitsLen(bitsLen_), socket(socket_)
{
    checkBitLength(bitsLen);
    sender = new osuCrypto::SilentOtExtSender();
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
    if (inputs.size() < num || val.size() < num) {
        throw std::invalid_argument("multiplication input/output size is smaller than num");
    }

    const u64 bytesLen = bitsLen / 8;
    const u64 batchSize = maxBatchOts / bitsLen;
    const u64 correctionBytes = 4 * bytesLen * (bytesLen + 1);

    for (u64 begin = 0; begin < num;) {
        const u64 count = std::min(batchSize, num - begin);
        const u64 numOts = count * bitsLen;
        sender->configure(numOts);
        coproto::sync_wait(sender->genSilentBaseOts(*prng, *socket));

        std::vector<std::array<block, 2>> messages(numOts);
        coproto::sync_wait(sender->send(messages, *prng, *socket));

        std::vector<u8> compressedBits(count * correctionBytes);
        u64 offset = 0;

        for (u64 i = 0; i < numOts; i++) {
            int shift = i % bitsLen;
            const u64 index = begin + i / bitsLen;

            u64 mask = low(messages[i][0]) + inputs[index];
            u64 correction = (low(messages[i][1]) ^ mask) << shift >> shift;
            for (u64 j = 0; j < bytesLen - shift / 8; j++) {
                compressedBits[offset++] = static_cast<u8>(correction >> (8 * j));
            }
            val[index] += 0 - ((low(messages[i][0]) << shift >> shift) << shift);
        }

        coproto::sync_wait(socket->send(compressedBits));
        // The queued send borrows compressedBits until flush completes.
        coproto::sync_wait(socket->flush());
        begin += count;
    }
}

MulRecver::MulRecver(uint64_t num_, coproto::Socket *socket_, int bitsLen_) : num(num_), bitsLen(bitsLen_), socket(socket_)
{
    checkBitLength(bitsLen);
    receiver = new osuCrypto::SilentOtExtReceiver();
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
    if (inputs.size() < num || val.size() < num) {
        throw std::invalid_argument("multiplication input/output size is smaller than num");
    }

    const u64 bytesLen = bitsLen / 8;
    const u64 batchSize = maxBatchOts / bitsLen;
    const u64 correctionBytes = 4 * bytesLen * (bytesLen + 1);

    for (u64 begin = 0; begin < num;) {
        const u64 count = std::min(batchSize, num - begin);
        const u64 numOts = count * bitsLen;
        receiver->configure(numOts);
        coproto::sync_wait(receiver->genSilentBaseOts(*prng, *socket));

        BitVector choiceBit(numOts);
        for (u64 i = 0; i < count; i++) {
            u64 lowbits = inputs[begin + i];
            for (u64 j = 0; j < bytesLen; j++) {
                choiceBit.data()[i * bytesLen + j] = static_cast<u8>(lowbits >> (8 * j));
            }
        }

        std::vector<block> messages(numOts);
        coproto::sync_wait(receiver->receive(choiceBit, messages, *prng, *socket));

        std::vector<u8> compressedBits(count * correctionBytes);
        coproto::sync_wait(socket->recv(compressedBits));

        u64 offset = 0;
        for (u64 i = 0; i < numOts; i++) {
            u64 msg = 0;
            int shift = i % bitsLen;
            const u64 index = begin + i / bitsLen;
            for (u64 j = 0; j < bytesLen - shift / 8; j++) {
                msg |= u64(compressedBits[offset++]) << (8 * j);
            }
            if (choiceBit[i] & 1) {
                val[index] += ((low(messages[i]) ^ msg) << shift >> shift) << shift;
            } else {
                val[index] += (low(messages[i]) << shift >> shift) << shift;
            }
        }
        begin += count;
    }
}
