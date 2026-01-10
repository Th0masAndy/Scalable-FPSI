#include "mux.h"
#include <coproto/Socket/Socket.h>
#include <volePSI/GMW/Circuit.h>
#include <volePSI/GMW/Gmw.h>
#include <volePSI/Paxos.h>
#include <volePSI/config.h>
#include "utils.h"

using namespace oc;

MuxSender::MuxSender(uint64_t num_, coproto::Socket *socket_) : num(num_), socket(socket_)
{
    sender = new osuCrypto::SilentOtExtSender();
    sender->configure(num);
    sender->mMultType = type;

    recver = new osuCrypto::SilentOtExtReceiver();
    recver->configure(num);
    recver->mMultType = type;

    prng = new PRNG(ZeroBlock);
}

MuxSender::~MuxSender()
{
    delete sender;
    delete recver;
    delete prng;
}

void MuxSender::muxA(std::vector<u8> &choices, std::vector<u64> &v0, std::vector<u64> &res0)
{
    BitVector b0(choices.size());
    for (u64 i = 0; i < choices.size(); ++i) {
        b0[i] = choices[i] & 1;
    }

    coproto::sync_wait(sender->genSilentBaseOts(*prng, *socket));
    coproto::sync_wait(recver->genSilentBaseOts(*prng, *socket));

    std::vector<std::array<block, 2>> messages(num);
    coproto::sync_wait(sender->send(messages, *prng, *socket));

    std::vector<u64> correctMessages(num);

    for (u64 i = 0; i < num; i++) {
        u64 mask = low(messages[i][0]) + (v0[i] - 2 * u64(b0[i]) * v0[i]);
        correctMessages[i] = low(messages[i][1]) ^ mask;
    }

    coproto::sync_wait(socket->send(correctMessages));

    std::vector<block> message(num);
    coproto::sync_wait(recver->receive(b0, message, *prng, *socket));

    std::vector<u64> correctMessages1(num);

    coproto::sync_wait(socket->recv(correctMessages1));

    for (u64 i = 0; i < num; i++) {
        res0[i] = low(message[i]) ^ (b0[i] ? correctMessages1[i] : 0);
        res0[i] = res0[i] + u64(b0[i]) * v0[i];
        res0[i] = res0[i] - low(messages[i][0]);
    }
}

void MuxSender::muxA(std::vector<u8> &choices, std::vector<u64> &v0, std::vector<u64> &res0, int bitsLen)
{
    int bytesLen = divCeil(bitsLen, 8);
    u64 MASK = (u64(1) << bitsLen) - 1;

    BitVector b0(choices.size());
    for (u64 i = 0; i < choices.size(); ++i) {
        b0[i] = choices[i] & 1;
    }

    // coproto::sync_wait(sender->genSilentBaseOts(*prng, *socket));
    // coproto::sync_wait(recver->genSilentBaseOts(*prng, *socket));

    std::vector<std::array<block, 2>> messages(num);
    coproto::sync_wait(sender->send(messages, *prng, *socket));

    std::vector<u8> correctMessages(num * bytesLen);

    for (u64 i = 0; i < num; i++) {
        u64 mask = low(messages[i][0]) + (v0[i] - 2 * u64(b0[i]) * v0[i]);
        u64 val = low(messages[i][1]) ^ mask;
        memcpy(correctMessages.data() + i * bytesLen, &val, bytesLen);
    }

    coproto::sync_wait(socket->send(correctMessages));

    std::vector<block> message(num);
    coproto::sync_wait(recver->receive(b0, message, *prng, *socket));

    std::vector<u8> correctMessages1(num * bytesLen);

    coproto::sync_wait(socket->recv(correctMessages1));

    for (u64 i = 0; i < num; i++) {
        u64 val;
        memcpy(&val, correctMessages1.data() + i * bytesLen, bytesLen);
        res0[i] = low(message[i]) ^ (b0[i] ? val : 0);
        res0[i] = res0[i] + u64(b0[i]) * v0[i];
        res0[i] = res0[i] - low(messages[i][0]);
        res0[i] &= MASK;
    }
}

MuxRecver::MuxRecver(uint64_t num_, coproto::Socket *socket_) : num(num_), socket(socket_)
{
    sender = new osuCrypto::SilentOtExtSender();
    sender->configure(num);
    sender->mMultType = type;

    recver = new osuCrypto::SilentOtExtReceiver();
    recver->configure(num);
    recver->mMultType = type;

    prng = new PRNG(OneBlock);
}

MuxRecver::~MuxRecver()
{
    delete recver;
    delete sender;
    delete prng;
}

void MuxRecver::muxA(std::vector<u8> &choices, std::vector<u64> &v1, std::vector<u64> &res1)
{
    BitVector b1(choices.size());
    for (u64 i = 0; i < choices.size(); ++i) {
        b1[i] = choices[i] & 1;
    }

    coproto::sync_wait(recver->genSilentBaseOts(*prng, *socket));
    coproto::sync_wait(sender->genSilentBaseOts(*prng, *socket));

    std::vector<block> message(num);
    coproto::sync_wait(recver->receive(b1, message, *prng, *socket));

    std::vector<u64> correctMessages1(num);

    coproto::sync_wait(socket->recv(correctMessages1));

    std::vector<std::array<block, 2>> messages(num);
    coproto::sync_wait(sender->send(messages, *prng, *socket));

    std::vector<u64> correctMessages(num);

    for (u64 i = 0; i < num; i++) {
        u64 mask = low(messages[i][0]) + (v1[i] - 2 * u64(b1[i]) * v1[i]);
        correctMessages[i] = low(messages[i][1]) ^ mask;
    }

    coproto::sync_wait(socket->send(correctMessages));

    for (u64 i = 0; i < num; i++) {
        res1[i] = low(message[i]) ^ (b1[i] ? correctMessages1[i] : 0);
        res1[i] = res1[i] + u64(b1[i]) * v1[i];
        res1[i] = res1[i] - low(messages[i][0]);
    }
}

void MuxRecver::muxA(std::vector<u8> &choices, std::vector<u64> &v1, std::vector<u64> &res1, int bitsLen)
{
    int bytesLen = divCeil(bitsLen, 8);
    u64 MASK = (u64(1) << bitsLen) - 1;

    BitVector b1(choices.size());
    for (u64 i = 0; i < choices.size(); ++i) {
        b1[i] = choices[i] & 1;
    }

    // coproto::sync_wait(recver->genSilentBaseOts(*prng, *socket));
    // coproto::sync_wait(sender->genSilentBaseOts(*prng, *socket));

    std::vector<block> message(num);
    coproto::sync_wait(recver->receive(b1, message, *prng, *socket));

    std::vector<u8> correctMessages1(num * bytesLen);

    coproto::sync_wait(socket->recv(correctMessages1));

    std::vector<std::array<block, 2>> messages(num);
    coproto::sync_wait(sender->send(messages, *prng, *socket));

    std::vector<u8> correctMessages(num * bytesLen);

    for (u64 i = 0; i < num; i++) {
        u64 mask = low(messages[i][0]) + (v1[i] - 2 * u64(b1[i]) * v1[i]);
        u64 val = low(messages[i][1]) ^ mask;
        memcpy(correctMessages.data() + i * bytesLen, &val, bytesLen);
    }

    coproto::sync_wait(socket->send(correctMessages));

    for (u64 i = 0; i < num; i++) {
        u64 val;
        memcpy(&val, correctMessages1.data() + i * bytesLen, bytesLen);
        res1[i] = low(message[i]) ^ (b1[i] ? val : 0);
        res1[i] = res1[i] + u64(b1[i]) * v1[i];
        res1[i] = res1[i] - low(messages[i][0]);
        res1[i] &= MASK;
    }
}
