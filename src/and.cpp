#include "and.h"
#include <algorithm>
#include <coproto/Common/macoro.h>
#include <libOTe/Triple/SilentOtTriple/SilentOtTriple.h>
#include <stdexcept>
#include "cmp.h"

using namespace oc;

namespace {
std::vector<u8> paddedBits(const std::vector<u8> &bits, u64 paddedSize)
{
    if (bits.size() > paddedSize) {
        throw std::runtime_error("AND input is larger than padded size");
    }

    std::vector<u8> out(paddedSize, 0);
    std::copy(bits.begin(), bits.end(), out.begin());
    return out;
}
} // namespace

AndSender::AndSender(uint64_t num_, coproto::Socket *socket_)
    : num(num_), socket(socket_)
{
}

void AndSender::andBits(std::vector<u8> &lhs, std::vector<u8> &rhs, std::vector<u8> &out)
{
    if (lhs.size() != num || rhs.size() != num) {
        throw std::runtime_error("AND input size mismatch");
    }

    PRNG prng(sysRandomSeed());
    SilentOtTriple tripleGen;
    u64 numTriples = roundUpTo(num, 128);
    u64 packedLen = numTriples / 8;

    tripleGen.init(1, numTriples);
    coproto::sync_wait(tripleGen.genBaseOts(prng, *socket));

    std::vector<block> A(numTriples / 128);
    std::vector<block> B(numTriples / 128);
    std::vector<block> C(numTriples / 128);
    coproto::sync_wait(tripleGen.expand(A, B, C, prng, *socket));

    auto ai = reinterpret_cast<u8 *>(A.data());
    auto bi = reinterpret_cast<u8 *>(B.data());
    auto ci = reinterpret_cast<u8 *>(C.data());
    auto x = paddedBits(lhs, numTriples);
    auto y = paddedBits(rhs, numTriples);

    std::vector<u8> ei(packedLen);
    std::vector<u8> fi(packedLen);
    std::vector<u8> e(packedLen);
    std::vector<u8> f(packedLen);

    MillionaireProtocolSender::AND_step_1(ei.data(), fi.data(), x.data(), y.data(), ai, bi, static_cast<int>(numTriples));

    coproto::sync_wait(socket->send(oc::span<u8>(ei.data(), packedLen)));
    coproto::sync_wait(socket->send(oc::span<u8>(fi.data(), packedLen)));

    coproto::sync_wait(socket->recv(oc::span<u8>(e.data(), packedLen)));
    coproto::sync_wait(socket->recv(oc::span<u8>(f.data(), packedLen)));

    for (u64 i = 0; i < packedLen; ++i) {
        e[i] ^= ei[i];
        f[i] ^= fi[i];
    }

    std::vector<u8> andRes(numTriples);
    MillionaireProtocolSender::AND_step_2(andRes.data(), e.data(), f.data(), ei.data(), fi.data(), ai, bi, ci, static_cast<int>(numTriples));

    out.assign(andRes.begin(), andRes.begin() + num);
}

AndRecver::AndRecver(uint64_t num_, coproto::Socket *socket_)
    : num(num_), socket(socket_)
{
}

void AndRecver::andBits(std::vector<u8> &lhs, std::vector<u8> &rhs, std::vector<u8> &out)
{
    if (lhs.size() != num || rhs.size() != num) {
        throw std::runtime_error("AND input size mismatch");
    }

    PRNG prng(sysRandomSeed());
    SilentOtTriple tripleGen;
    u64 numTriples = roundUpTo(num, 128);
    u64 packedLen = numTriples / 8;

    tripleGen.init(0, numTriples);
    coproto::sync_wait(tripleGen.genBaseOts(prng, *socket));

    std::vector<block> A(numTriples / 128);
    std::vector<block> B(numTriples / 128);
    std::vector<block> C(numTriples / 128);
    coproto::sync_wait(tripleGen.expand(A, B, C, prng, *socket));

    auto ai = reinterpret_cast<u8 *>(A.data());
    auto bi = reinterpret_cast<u8 *>(B.data());
    auto ci = reinterpret_cast<u8 *>(C.data());
    auto x = paddedBits(lhs, numTriples);
    auto y = paddedBits(rhs, numTriples);

    std::vector<u8> ei(packedLen);
    std::vector<u8> fi(packedLen);
    std::vector<u8> e(packedLen);
    std::vector<u8> f(packedLen);

    MillionaireProtocolRecver::AND_step_1(ei.data(), fi.data(), x.data(), y.data(), ai, bi, static_cast<int>(numTriples));

    coproto::sync_wait(socket->send(oc::span<u8>(ei.data(), packedLen)));
    coproto::sync_wait(socket->send(oc::span<u8>(fi.data(), packedLen)));

    coproto::sync_wait(socket->recv(oc::span<u8>(e.data(), packedLen)));
    coproto::sync_wait(socket->recv(oc::span<u8>(f.data(), packedLen)));

    for (u64 i = 0; i < packedLen; ++i) {
        e[i] ^= ei[i];
        f[i] ^= fi[i];
    }

    std::vector<u8> andRes(numTriples);
    MillionaireProtocolRecver::AND_step_2(andRes.data(), e.data(), f.data(), ei.data(), fi.data(), ai, bi, ci, static_cast<int>(numTriples));

    out.assign(andRes.begin(), andRes.begin() + num);
}
