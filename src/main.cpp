#include <coproto/Common/macoro.h>
#include <coproto/Socket/AsioSocket.h>
#include <coproto/Socket/Socket.h>
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Common/Timer.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtSender.h>
#include <thread>
#include <vector>
#include "cmp.h"
#include "proto.h"

int main(int argc, char **argv)
{
    oc::CLP cmd(argc, argv);

    int lp = cmd.getOr("p", 0);

    const std::pair<const char *, std::function<void()>> handlers[] = {
        { "bp25",
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
              bp25High(cmd);
          } },
    };

    for (const auto &[flag, run] : handlers) {
        if (cmd.isSet(flag)) {
            run();
            break;
        }
    }

    return 0;

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

    // u64 n = 128 + (1 << 14);

    // MillionaireProtocolRecver recver(n, 64, 4);
    // MillionaireProtocolSender sender(n, 64, 4);

    // auto socket = coproto::AsioSocket::makePair();

    // oc::PRNG prng(oc::sysRandomSeed());

    // u64 data0[n];
    // u64 data1[n];
    // u8 outs0[n];
    // u8 outs1[n];

    // for (u64 i = 0; i < n; ++i) {
    //     data0[i] = prng.get<u64>();
    //     data1[i] = prng.get<u64>();
    // }

    // oc::Timer time;
    // auto s = time.setTimePoint("begin cmp");

    // std::thread t0([&]() { recver.compare(outs0, data0, socket[0]); });

    // std::thread t1([&]() { sender.compare(outs1, data1, socket[1]); });

    // t0.join();
    // t1.join();

    // auto e = time.setTimePoint("end cmp");

    // int correct = 0;
    // for (u64 i = 0; i < n; ++i) {
    //     bool gt = data1[i] > data0[i];
    //     if (gt == ((outs1[i] ^ outs0[i]) & 1))
    //         correct++;
    // }
    // std::cout << "correct: " << correct << " / " << n << std::endl;

    // auto comm = socket[0].bytesReceived() + socket[0].bytesSent();
    // auto comp = std::chrono::duration_cast<std::chrono::microseconds>(e - s).count();
    // std::cout << "average: " << comm / (n * 1.0) << " bytes" << " " << comp / (n * 1.0) << " microseconds" << std::endl;

    u64 n = 1 << 18;

    SilentOtExtReceiver recv;
    SilentOtExtSender send;

    oc::PRNG prng(oc::sysRandomSeed());

    auto socket = coproto::AsioSocket::makePair();

    std::thread t0([&]() {
        recv.configure(n);
        BitVector choices(n);
        std::vector<block> messages(n);

        coproto::sync_wait(recv.receive(choices, messages, prng, socket[0]));
    });

    std::thread t1([&]() {
        send.configure(n);
        std::vector<std::array<block, 2>> messages(n);
        for (u64 i = 0; i < n; ++i) {
            messages[i][0] = prng.get<block>();
            messages[i][1] = prng.get<block>();
        }

        coproto::sync_wait(send.send(messages, prng, socket[1]));
    });
    t0.join();
    t1.join();

    auto comm = socket[0].bytesReceived() + socket[0].bytesSent();
    // auto comp = std::chrono::duration_cast<std::chrono::microseconds>(e - s).count();
    std::cout << "average: " << comm << " bytes" << std::endl;
}