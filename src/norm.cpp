#include "norm.h"
#include "cmp.h"
#include "mul.h"
#include "mux.h"
#include <stdexcept>
#include <thread>

namespace {
u64 makeMask(int bitsLen)
{
    if (bitsLen <= 0 || bitsLen > 64) {
        throw std::runtime_error("invalid norm bit length");
    }
    return bitsLen == 64 ? ~0ull : ((1ull << bitsLen) - 1);
}
} // namespace

void normL0(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl)
{
    int bitsLen = bytesLen * 8;
    u64 mask = makeMask(bitsLen);

    auto n = x.size() / d;

    std::thread cmpSendThr([&]() {
        MillionaireProtocolSender sender(x.size(), bitsLen);

        std::vector<u8> cmpShare(x.size());
        sender.drelu(cmpShare.data(), x.data(), chl[1]);

        MuxSender mux(x.size(), &chl[1]);

        std::vector<u64> res(x.size());
        std::vector<u64> x_vec(x.begin(), x.end());

        mux.muxA(cmpShare, x_vec, res, bitsLen);

        std::vector<u64> abs(x.size());
        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = 2 * res[i] - x[i];
        }

        std::vector<u64> dis_max(n, 0);

        MillionaireProtocolSender sender2(n, bitsLen);
        MuxSender mux2(n, &chl[1]);

        for (u64 i = 0; i < d; i++) {
            std::vector<u8> compare_res(n);
            std::vector<u64> curr(n, 0);
            std::vector<u64> res(n, 0);
            for (u64 j = 0; j < n; ++j) {
                curr[j] = abs[j * d + i] - dis_max[j];
            }
            sender2.drelu(compare_res.data(), curr.data(), chl[1]);
            mux2.muxA(compare_res, curr, res, bitsLen);
            for (u64 j = 0; j < n; ++j) {
                dis_max[j] += res[j];
            }
        }

        for (u64 i = 0; i < dis_max.size(); ++i) {
            dis_max[i] = delta - dis_max[i];
        }

        sender2.drelu(resBits1.data(), dis_max.data(), chl[1]);
    });

    std::thread cmpRecvThr([&]() {
        MillionaireProtocolRecver recver(y.size(), bitsLen);

        std::vector<u8> cmpShare(y.size());

        recver.drelu(cmpShare.data(), y.data(), chl[0]);

        MuxRecver mux(y.size(), &chl[0]);
        std::vector<u64> res(y.size());
        std::vector<u64> y_vec(y.begin(), y.end());

        mux.muxA(cmpShare, y_vec, res, bitsLen);

        std::vector<u64> abs(y.size());

        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = 2 * res[i] - y[i];
        }

        std::vector<u64> dis_max(n, 0);

        MillionaireProtocolRecver recver2(n, bitsLen);
        MuxRecver mux2(n, &chl[0]);

        for (u64 i = 0; i < d; i++) {
            std::vector<u8> compare_res(n);
            std::vector<u64> curr(n, 0);
            std::vector<u64> res(n, 0);
            for (u64 j = 0; j < n; ++j) {
                curr[j] = abs[j * d + i] - dis_max[j];
            }
            recver2.drelu(compare_res.data(), curr.data(), chl[0]);
            mux2.muxA(compare_res, curr, res, bitsLen);
            for (u64 j = 0; j < n; ++j) {
                dis_max[j] += res[j];
            }
        }

        for (u64 i = 0; i < dis_max.size(); ++i) {
            dis_max[i] = (-dis_max[i]) & mask;
        }

        recver2.drelu(resBits0.data(), dis_max.data(), chl[0]);
    });

    cmpSendThr.join();
    cmpRecvThr.join();
}

void normL1(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl)
{
    int bitsLen = bytesLen * 8;
    u64 mask = makeMask(bitsLen);

    auto n = x.size() / d;

    std::thread cmpSendThr([&]() {
        MillionaireProtocolSender sender(x.size(), bitsLen);

        std::vector<u8> cmpShare(x.size());
        sender.drelu(cmpShare.data(), x.data(), chl[1]);

        MuxSender mux(x.size(), &chl[1]);

        std::vector<u64> res(x.size());
        std::vector<u64> x_vec(x.begin(), x.end());

        mux.muxA(cmpShare, x_vec, res, bitsLen);

        std::vector<u64> abs(x.size());
        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = 2 * res[i] - x[i];
        }

        std::vector<u64> dis(n, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += abs[i * d + j];
            }
            dis[i] = delta - dis[i];
        }

        MillionaireProtocolSender sender2(dis.size(), bitsLen);

        sender2.drelu(resBits1.data(), dis.data(), chl[1]);
    });

    std::thread cmpRecvThr([&]() {
        MillionaireProtocolRecver recver(y.size(), bitsLen);

        std::vector<u8> cmpShare(y.size());
        recver.drelu(cmpShare.data(), y.data(), chl[0]);

        MuxRecver mux(y.size(), &chl[0]);
        std::vector<u64> res(y.size());
        std::vector<u64> y_vec(y.begin(), y.end());

        mux.muxA(cmpShare, y_vec, res, bitsLen);

        std::vector<u64> abs(y.size());

        for (u64 i = 0; i < abs.size(); ++i) {
            abs[i] = 2 * res[i] - y[i];
        }

        std::vector<u64> dis(y.size() / d, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += abs[i * d + j];
            }
            dis[i] = (-dis[i]) & mask;
        }

        MillionaireProtocolRecver recver2(dis.size(), bitsLen);

        recver2.drelu(resBits0.data(), dis.data(), chl[0]);
    });

    cmpSendThr.join();
    cmpRecvThr.join();
}

void normL2(
    oc::span<u64> x,
    oc::span<u64> y,
    std::vector<u8> &resBits0,
    std::vector<u8> &resBits1,
    u64 d,
    int delta,
    int bytesLen,
    std::array<coproto::AsioSocket, 2> &chl)
{
    u64 delta_p = u64(delta) * u64(delta);
    int bitsLen = bytesLen * 8;

    u64 mask = makeMask(bitsLen);

    auto n = x.size() / d;

    std::thread cmpSendThr([&]() {
        MulSender sender(x.size(), &chl[1], bitsLen);

        std::vector<u64> x_vec(x.begin(), x.end());
        std::vector<u64> dots(x.size());
        sender.mul(x_vec, dots);

        std::vector<u64> absSquare(x.size());
        for (u64 i = 0; i < absSquare.size(); ++i) {
            u64 square = (__uint128_t(x[i]) * __uint128_t(x[i])) & mask;
            absSquare[i] = square + 2 * dots[i];
            absSquare[i] = absSquare[i] & mask;
        }

        std::vector<u64> dis(n, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += absSquare[i * d + j];
            }
            dis[i] = (delta_p - dis[i]) & mask;
        }

        MillionaireProtocolSender sender2(dis.size(), bitsLen);

        sender2.drelu(resBits1.data(), dis.data(), chl[1]);
    });

    std::thread cmpRecvThr([&]() {
        MulRecver recver(y.size(), &chl[0], bitsLen);

        std::vector<u64> y_vec(y.begin(), y.end());
        std::vector<u64> dots(y.size());
        recver.mul(y_vec, dots);

        std::vector<u64> absSquare(y.size());

        for (u64 i = 0; i < absSquare.size(); ++i) {
            u64 square = (__uint128_t(y[i]) * __uint128_t(y[i])) & mask;
            absSquare[i] = square + 2 * dots[i];
            absSquare[i] = absSquare[i] & mask;
        }

        std::vector<u64> dis(y.size() / d, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += absSquare[i * d + j];
            }
            dis[i] = (-dis[i]) & mask;
        }

        MillionaireProtocolRecver recver2(dis.size(), bitsLen);

        recver2.drelu(resBits0.data(), dis.data(), chl[0]);
    });

    cmpSendThr.join();
    cmpRecvThr.join();
}
