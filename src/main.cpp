#include <coproto/Common/macoro.h>
#include <coproto/Socket/AsioSocket.h>
#include <coproto/Socket/Socket.h>
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Common/Timer.h>
#include <cstdint>
#include <iostream>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtSender.h>
#include <thread>
#include <vector>
#include <volePSI/RsOprf.h>
#include "cmp.h"
#include "mul.h"
#include "mux.h"
#include "proto.h"

int main(int argc, char **argv)
{
    oc::CLP cmd(argc, argv);

    int lp = cmd.getOr("p", 0);

    const std::pair<const char *, std::function<void()>> handlers[] = {
        { "bp25low",
          [&] {
              (lp ? bp25LowLpPx : bp25LowPx)(cmd);
          } },
        { "low",
          [&] {
              (lp ? fpsiLowLpPx : fpsiLowLpPx)(cmd); // unified framework for lp and l_inf
          } },
        { "high",
          [&] {
              (lp ? fpsiHighLpPx : fpsiHighPx)(cmd);
          } },
        { "bp25high",
          [&] {
              (lp ? bp25HighLp : bp25High)(cmd);
          } },
    };

    for (const auto &[flag, run] : handlers) {
        if (cmd.isSet(flag)) {
            run();
            break;
        }
    }

    return 0;

    auto chl = coproto::AsioSocket::makePair();

    u64 size = 1 << 16;

    int bitLen = 56;
    u64 mask = (1ULL << bitLen) - 1;

    MulSender sender1(size, &chl[1], 56);
    MulRecver recver1(size, &chl[0], 56);

    std::vector<u64> input1(size);
    std::vector<u64> input2(size);
    std::vector<u64> output1(size);
    std::vector<u64> output2(size);
    oc::PRNG prng1(oc::sysRandomSeed());
    prng1.get(input1.data(), size);
    prng1.get(input2.data(), size);

    for (auto i = 0; i < size; ++i) {
        input1[i] = input1[i] & mask;
        input2[i] = input2[i] & mask;
    }

    Timer t;
    t.setTimePoint("begin");

    std::thread sendThr1([&]() { sender1.mul(input1, output1); });

    std::thread recvThr1([&]() { recver1.mul(input2, output2); });

    sendThr1.join();
    recvThr1.join();

    t.setTimePoint("end");

    std::cout << t << std::endl;

    for (auto i = 0; i < size; ++i) {
        __uint128_t res = (__uint128_t(input1[i]) * __uint128_t(input2[i])) & mask;
        __uint128_t val = (__uint128_t(output1[i]) + __uint128_t(output2[i])) & mask;
        if (res != val) {
            std::cout << "error at " << i << ": " << (u64)res << " " << (u64)val << std::endl;
        }
        if (i == 0) {
            std::cout << "check mul: " << input1[i] << " * " << input2[i] << " = " << (u64)res << std::endl;
        }
    }

    auto comm1 = chl[0].bytesReceived() + chl[0].bytesSent();
    std::cout << "comm = " << comm1 * 8 / (double)size << " bits per input" << std::endl;

    return 0;

    PRNG prng(oc::sysRandomSeed());

    u64 r0 = prng.get<u64>();
    u64 r1 = prng.get<u64>();
    u64 x0 = prng.get<u64>();
    u64 x1 = prng.get<u64>();

    // u64 mask = (1ULL << 56) - 1;
    r0 = r0 & mask;
    r1 = r1 & mask;
    x0 = x0 & mask;
    x1 = x1 & mask;
    x1 = (((r0 + r1) & mask) + (1ULL << 56) - x0) & mask;

    std::cout << ((x0 + x1) & mask) << std::endl;
    std::cout << ((r0 + r1) & mask) << std::endl;

    u64 y0 = 2 * r0 - x0;

    u64 y1 = 2 * r1 - x1;

    std::cout << int64_t((y0 + y1) & mask) << std::endl;

    return 0;

    std::vector<u64> x;
    std::vector<u64> y;

    for (u64 i = 0; i < (1 << 7); ++i) {
        auto v = prng.get<u64>();
        x.push_back(v);
        y.push_back(i * 3 - v);
    }

    int d = 1;

    int delta = 128;

    // auto chl = coproto::AsioSocket::makePair();

    int bitsLen = 64;

    auto n = x.size() / d;

    std::vector<u8> choiceBit(n);

    std::vector<u64> abs_s(x.size());
    std::vector<u64> abs_r(x.size());

    std::thread cmpSendThr([&]() {
        MillionaireProtocolSender sender(x.size(), bitsLen);

        std::vector<u8> cmpShare(x.size());
        sender.drelu(cmpShare.data(), x.data(), chl[1]);

        MuxSender mux(x.size(), &chl[1]);

        std::vector<u64> res(x.size());
        std::vector<u64> x_vec(x.begin(), x.end());

        mux.muxA(cmpShare, x_vec, res);

        coproto::sync_wait(chl[1].send(cmpShare));

        for (u64 i = 0; i < abs_s.size(); ++i) {
            abs_s[i] = 2 * res[i] - x[i];
        }

        std::vector<u64> dis(n, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += abs_s[i * d + j];
            }
            dis[i] = delta - dis[i];
        }

        MillionaireProtocolSender sender2(dis.size(), bitsLen);

        std::vector<u8> resBits1(dis.size());

        sender2.drelu(resBits1.data(), dis.data(), chl[1]);

        coproto::sync_wait(chl[1].send(resBits1));
    });

    std::thread cmpRecvThr([&]() {
        MillionaireProtocolRecver recver(y.size(), bitsLen);

        std::vector<u8> cmpShare(y.size());
        recver.drelu(cmpShare.data(), y.data(), chl[0]);

        MuxRecver mux(y.size(), &chl[0]);
        std::vector<u64> res(y.size());
        std::vector<u64> y_vec(y.begin(), y.end());

        mux.muxA(cmpShare, y_vec, res);

        std::vector<u8> cmpShare1(y.size());
        coproto::sync_wait(chl[0].recv(cmpShare1));
        for (u64 i = 0; i < cmpShare1.size(); ++i) {
            cmpShare[i] = cmpShare[i] ^ cmpShare1[i];
            std::cout << (cmpShare[i] & 1) << std::endl;
        }

        for (u64 i = 0; i < abs_r.size(); ++i) {
            abs_r[i] = 2 * res[i] - y[i];
        }

        std::vector<u64> dis(y.size() / d, 0);

        for (u64 i = 0; i < dis.size(); ++i) {
            for (u64 j = 0; j < d; ++j) {
                dis[i] += abs_r[i * d + j];
            }
            dis[i] = -dis[i];
        }

        MillionaireProtocolRecver recver2(dis.size(), bitsLen);

        recver2.drelu(choiceBit.data(), dis.data(), chl[0]);

        std::vector<u8> resBits1(dis.size());
        coproto::sync_wait(chl[0].recv(resBits1));
        for (u64 i = 0; i < dis.size(); ++i) {
            choiceBit[i] = (resBits1[i] ^ choiceBit[i]) & 1;
        }
    });

    cmpSendThr.join();
    cmpRecvThr.join();

    for (u64 i = 0; i < abs_r.size(); ++i) {
        std::cout << (abs_r[i] + abs_s[i]) << std::endl;
    }

    for (u64 i = 0; i < choiceBit.size(); ++i) {
        std::cout << "choiceBit[" << i << "] = " << (u64)choiceBit[i] << std::endl;
    }

    // for (u64 i = 0; i < x.size(); i++) {
    //     std::cout << x[i] + y[i] << std::endl;
    // }

    return 0;

    // u64 n = 1 << 10;

    // int bitsLen = 16;

    // auto chl = coproto::AsioSocket::makePair();

    // MuxSender muxS(n, &chl[1]);
    // MuxRecver muxR(n, &chl[0]);

    // std::vector<u8> choices0(n, 1);
    // std::vector<u8> choices1(n, 1);

    // std::vector<u64> v0(n);
    // std::vector<u64> v1(n);

    // oc::PRNG prng(oc::sysRandomSeed());

    // for (u64 i = 0; i < n; ++i) {
    //     v0[i] = prng.get<u64>() & ((1ull << bitsLen) - 1);
    //     v1[i] = prng.get<u64>() & ((1ull << bitsLen) - 1);
    // }

    // std::vector<u64> res0(n);
    // std::vector<u64> res1(n);

    // std::thread recverThr([&]() { muxR.muxA(choices0, v0, res0, bitsLen); });

    // std::thread senderThr([&]() { muxS.muxA(choices1, v1, res1, bitsLen); });

    // senderThr.join();
    // recverThr.join();

    // for (u64 i = 0; i < n; ++i) {
    //     if (((res0[i] + res1[i]) & ((1ull << bitsLen) - 1)) != ((v0[i] + v1[i]) & ((1ull << bitsLen) - 1))) {
    //         std::cout << "error at " << i << std::endl;
    //     }
    // }

    // return 0;

    // BitVector output0;
    // BitVector output1;

    // std::vector<u8> bit(2);
    // bit[0] = 1;
    // bit[1] = 1;
    // BitVector t(bit.data(), 2 * 8);
    // for (auto v : t) {
    //     std::cout << v << std::endl;
    // }

    // BitVector bits1(8);
    // BitVector bits0(8);

    // for (int i = 0; i < 8; i++) {
    //     bits1[i] = (i == 7);
    //     bits0[i] = 0;
    // }

    // auto sockets = coproto::AsioSocket::makePair();

    // std::thread permuteSendThr([&]() {
    //     std::vector<u32> pi;
    //     Opermute1(0, bits1, output1, pi, sockets[1], 1);
    //     for (auto v : pi) {
    //         std::cout << "pi: " << v << std::endl;
    //     }
    // });

    // std::thread permuteRecvThr([&]() {
    //     std::vector<u32> pi;
    //     Opermute1(1, bits0, output0, pi, sockets[0], 1);
    // });

    // permuteSendThr.join();
    // permuteRecvThr.join();

    // for (u64 i = 0; i < output0.size(); i++) {
    //     std::cout << output0[i] << " " << output1[i] << std::endl;
    // }

    // return 0;

    // volePSI::Baxos paxos;
    // volePSI::u64 n = 1 << 20;
    // u64 c = 128;

    // paxos.init(n, n / 4, 3, 40, volePSI::PaxosParam::Binary, oc::ZeroBlock);

    // std::vector<oc::block> keys(n);
    // // std::vector<oc::block> vals(n);
    // // std::vector<oc::block> encodings(paxos.size());

    // oc::PRNG prng(oc::sysRandomSeed());

    // oc::Matrix<u8> values(n, c), values2(n, c), p(paxos.size(), c);

    // prng.get(keys.data(), keys.size());
    // prng.get(values.data(), values.size());

    // osuCrypto::Timer time;
    // time.setTimePoint("begin encode");

    // // paxos.solve<oc::block>(keys, vals, encodings, &prng, 1);
    // paxos.solve<u8>(keys, values, p);

    // time.setTimePoint("end encode");
    // std::cout << time << std::endl;

    // u64 n = (1 << 20) + 17;

    MillionaireProtocolRecver recver(n, 64, 4);
    MillionaireProtocolSender sender(n, 64, 4);

    auto socket = coproto::AsioSocket::makePair();

    // oc::PRNG prng(oc::sysRandomSeed());

    std::vector<u64> data0(n);
    std::vector<u64> data1(n);
    std::vector<u8> outs0(n);
    std::vector<u8> outs1(n);

    prng.get(data0.data(), n);
    prng.get(data1.data(), n);

    oc::Timer time;
    auto s = time.setTimePoint("begin cmp");

    std::thread recvThr([&]() {
        recver.compare(outs0.data(), data0.data(), socket[0]);
        recver.compare(outs0.data(), data0.data(), socket[0]);
    });

    std::thread sendThr([&]() {
        sender.compare(outs1.data(), data1.data(), socket[1]);
        sender.compare(outs1.data(), data1.data(), socket[1]);
    });

    recvThr.join();
    sendThr.join();

    auto e = time.setTimePoint("end cmp");

    int correct = 0;
    for (u64 i = 0; i < n; ++i) {
        bool gt = data1[i] > data0[i];
        if (gt == ((outs1[i] ^ outs0[i]) & 1))
            correct++;
    }
    std::cout << "correct: " << correct << " / " << n << std::endl;

    auto comm = socket[0].bytesReceived() + socket[0].bytesSent();
    auto comp = std::chrono::duration_cast<std::chrono::microseconds>(e - s).count();
    std::cout << "average: " << comm / (n * 1.0) << " bytes" << " " << comp / (n * 1.0) << " microseconds" << std::endl;

    // u64 n = 1 << 18;

    // SilentOtExtReceiver recv;
    // SilentOtExtSender send;

    // oc::PRNG prng(oc::sysRandomSeed());

    // auto socket = coproto::AsioSocket::makePair();

    // std::thread t0([&]() {
    //     recv.configure(n);
    //     BitVector choices(n);
    //     std::vector<block> messages(n);

    //     coproto::sync_wait(recv.receive(choices, messages, prng, socket[0]));
    // });

    // std::thread t1([&]() {
    //     send.configure(n);
    //     std::vector<std::array<block, 2>> messages(n);
    //     for (u64 i = 0; i < n; ++i) {
    //         messages[i][0] = prng.get<block>();
    //         messages[i][1] = prng.get<block>();
    //     }

    //     coproto::sync_wait(send.send(messages, prng, socket[1]));
    // });
    // t0.join();
    // t1.join();

    // auto comm = socket[0].bytesReceived() + socket[0].bytesSent();
    // // auto comp = std::chrono::duration_cast<std::chrono::microseconds>(e - s).count();
    // std::cout << "average: " << comm << " bytes" << std::endl;
}