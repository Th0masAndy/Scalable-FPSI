#include <functional>
#include <utility>
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
}