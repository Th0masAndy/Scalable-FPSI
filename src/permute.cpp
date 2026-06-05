#include "permute.h"
#include <cstdint>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtSender.h>
#include <stack>

using namespace oc;

class WaksmanNetwork {
public:
    explicit WaksmanNetwork(std::vector<int32_t> &dest)
    {
        auto n = (int32_t)dest.size();
        auto logN = int32_t(ceil(log2(n)));
        int32_t levels = 2 * logN - 1;
        waksman_perm.resize(n);
        waksman_inv_perm.resize(n);
        waksman_network.resize(levels);
        for (int32_t i = 0; i < levels; ++i) {
            waksman_network[i].resize(n / 2);
            std::fill(waksman_network[i].begin(), waksman_network[i].end(), -1);
        }
        // the input is [0, 1, ..., n)
        std::vector<int32_t> src(n);
        for (int32_t i = 0; i < n; ++i) {
            src[i] = i;
        }
        gen_waksman_route(logN, 0, 0, src, dest);
    }

    std::vector<std::vector<int8_t>> get_waksman_network()
    {
        return waksman_network;
    }

    ~WaksmanNetwork()
    {
        waksman_perm.clear();
        waksman_perm.shrink_to_fit();
        waksman_inv_perm.clear();
        waksman_perm.shrink_to_fit();
        waksman_network.clear();
        waksman_network.shrink_to_fit();
        waksman_path.clear();
        waksman_path.shrink_to_fit();
    }

private:
    /**
     * [N] -> [T]
     */
    std::vector<int32_t> waksman_perm;
    /**
     * [N] <- [T]
     */
    std::vector<int32_t> waksman_inv_perm;
    /**
     * Waksman network
     */
    std::vector<std::vector<int8_t>> waksman_network;
    /**
     * path
     */
    std::vector<int8_t> waksman_path;

    int32_t waksman_right_cycle_shift(int32_t num, int32_t logN)
    {
        return ((num & 1) << (logN - 1)) | (num >> 1);
    }

    /**
     * depth-first search.
     *
     * @param idx switching benes_network index.
     * @param route routh for the given index.
     */
    void waksman_depth_first_search(int32_t idx)
    {
        std::stack<std::pair<int32_t, int8_t>> stack;
        stack.push({ idx, 0 });
        std::pair<int32_t, int8_t> idxRoutePair;
        while (!stack.empty()) {
            idxRoutePair = stack.top();
            stack.pop();
            waksman_path[idxRoutePair.first] = idxRoutePair.second;
            // if the next item in the vertical array is unassigned
            if (waksman_path[idxRoutePair.first ^ 1] < 0) {
                // the next item is always assigned the opposite of this item,
                // unless it was part of path/cycle of previous node
                stack.push({ idxRoutePair.first ^ 1, idxRoutePair.second ^ (int8_t)1 });
            }
            idx = waksman_perm[waksman_inv_perm[idxRoutePair.first] ^ 1];
            if (waksman_path[idx] < 0) {
                stack.push({ idx, idxRoutePair.second ^ (int8_t)1 });
            }
        }
    }

    void waksman_even_depth_first_search()
    {
        assert(waksman_path.size() > 4 && waksman_path.size() % 2 == 0);
        // set the last path to be 0
        int32_t idx = waksman_perm[waksman_path.size() - 1];
        std::stack<std::pair<int32_t, int8_t>> stack;
        stack.push({ idx, 1 });
        std::pair<int32_t, int8_t> idxRoutePair;
        while (!stack.empty()) {
            idxRoutePair = stack.top();
            stack.pop();
            waksman_path[idxRoutePair.first] = idxRoutePair.second;
            // if the next item in the vertical array is unassigned
            if (waksman_path[idxRoutePair.first ^ 1] < 0) {
                // the next item is always assigned the opposite of this item,
                // unless it was part of path/cycle of previous node
                stack.push({ idxRoutePair.first ^ 1, idxRoutePair.second ^ (int8_t)1 });
            }
            idx = waksman_perm[waksman_inv_perm[idxRoutePair.first] ^ 1];
            if (waksman_path[idx] < 0) {
                stack.push({ idx, idxRoutePair.second ^ (int8_t)1 });
            }
        }
    }

    void gen_quadruple_switches(int8_t *switches, const std::vector<int32_t> &src, const std::vector<int32_t> &dest)
    {
        assert(src.size() == 4);
        assert(dest.size() == 4);
        if (dest[0] == src[0]) {
            // [0, 1, 2, 3] -> [0, ?, ?, ?]
            if (dest[1] == src[1]) {
                // [0, 1, 2, 3] -> [0, 1, ?, ?]
                if (dest[2] == src[2]) {
                    /*
                     * [0, 1, 2, 3] -> [0, 1, 2, 3], █ █ █ = 0 0 0
                     *                               █ █ □   0 0
                     */
                    switches[0] = 0, switches[1] = 0, switches[2] = 0, switches[3] = 0, switches[4] = 0;
                } else {
                    assert(dest[2] == src[3]);
                    /*
                     * [0, 1, 2, 3] -> [0, 1, 3, 2], █ █ █ = 0 0 0
                     *                               █ █ □   1 0
                     */
                    switches[0] = 0, switches[1] = 1, switches[2] = 0, switches[3] = 0, switches[4] = 0;
                }
            } else if (dest[1] == src[2]) {
                // [0, 1, 2, 3] -> [0, 2, ?, ?]
                if (dest[2] == src[1]) {
                    /*
                     * [0, 1, 2, 3] -> [0, 2, 1, 3], █ █ █ = 1 1 1
                     *                               █ █ □   0 0
                     */
                    switches[0] = 1, switches[1] = 0, switches[2] = 1, switches[3] = 0, switches[4] = 1;
                } else {
                    assert(dest[2] == src[3]);
                    /*
                     * [0, 1, 2, 3] -> [0, 2, 3, 1], █ █ █ = 0 0 0
                     *                               █ █ □   1 1
                     */
                    switches[0] = 0, switches[1] = 1, switches[2] = 0, switches[3] = 1, switches[4] = 0;
                }
            } else {
                assert(dest[1] == src[3]);
                // [0, 1, 2, 3] -> [0, 3, ?, ?]
                if (dest[2] == src[1]) {
                    /*
                     * [0, 1, 2, 3] -> [0, 3, 1, 2], █ █ █ = 1 1 1
                     *                               █ █ □   1 0
                     */
                    switches[0] = 1, switches[1] = 1, switches[2] = 1, switches[3] = 0, switches[4] = 1;
                } else {
                    assert(dest[2] == src[2]);
                    /*
                     * [0, 1, 2, 3] -> [0, 3, 2, 1], █ █ █ = 0 0 0
                     *                               █ █ □   0 1
                     */
                    switches[0] = 0, switches[1] = 0, switches[2] = 0, switches[3] = 1, switches[4] = 0;
                }
            }
        } else if (dest[0] == src[1]) {
            // [0, 1, 2, 3] -> [1, ?, ?, ?]
            if (dest[1] == src[0]) {
                // [0, 1, 2, 3] -> [1, 0, ?, ?]
                if (dest[2] == src[2]) {
                    /*
                     * [0, 1, 2, 3] -> [1, 0, 2, 3], █ █ █ = 0 0 1
                     *                               █ █ □   0 0
                     */
                    switches[0] = 0, switches[1] = 0, switches[2] = 0, switches[3] = 0, switches[4] = 1;
                } else {
                    assert(dest[2] == src[3]);
                    /*
                     * [0, 1, 2, 3] -> [1, 0, 3, 2], █ █ █ = 0 0 1
                     *                               █ █ □   1 0
                     */
                    switches[0] = 0, switches[1] = 1, switches[2] = 0, switches[3] = 0, switches[4] = 1;
                }
            } else if (dest[1] == src[2]) {
                // [0, 1, 2, 3] -> [1, 2, ?, ?]
                if (dest[2] == src[0]) {
                    /*
                     * [0, 1, 2, 3] -> [1, 2, 0, 3], █ █ █ = 0 1 1
                     *                               █ █ □   0 0
                     */
                    switches[0] = 0, switches[1] = 0, switches[2] = 1, switches[3] = 0, switches[4] = 1;
                } else {
                    assert(dest[2] == src[3]);
                    /*
                     * [0, 1, 2, 3] -> [1, 2, 3, 0], █ █ █ = 1 0 0
                     *                               █ █ □   1 1
                     */
                    switches[0] = 1, switches[1] = 1, switches[2] = 0, switches[3] = 1, switches[4] = 0;
                }
            } else {
                assert(dest[1] == src[3]);
                // [0, 1, 2, 3] -> [1, 3, ?, ?]
                if (dest[2] == src[0]) {
                    /*
                     * [0, 1, 2, 3] -> [1, 3, 0, 2], █ █ █ = 0 1 1
                     *                               █ █ □   1 0
                     */
                    switches[0] = 0, switches[1] = 1, switches[2] = 1, switches[3] = 0, switches[4] = 1;
                } else {
                    assert(dest[2] == src[2]);
                    /*
                     * [0, 1, 2, 3] -> [1, 3, 2, 0], █ █ █ = 1 0 0
                     *                               █ █ □   0 1
                     */
                    switches[0] = 1, switches[1] = 0, switches[2] = 0, switches[3] = 1, switches[4] = 0;
                }
            }
        } else if (dest[0] == src[2]) {
            // [0, 1, 2, 3] -> [2, ?, ?, ?]
            if (dest[1] == src[0]) {
                // [0, 1, 2, 3] -> [2, 0, ?, ?]
                if (dest[2] == src[1]) {
                    /*
                     * [0, 1, 2, 3] -> [2, 0, 1, 3], █ █ █ = 1 1 0
                     *                               █ █ □   0 0
                     */
                    switches[0] = 1, switches[1] = 0, switches[2] = 1, switches[3] = 0, switches[4] = 0;
                } else {
                    assert(dest[2] == src[3]);
                    /*
                     * [0, 1, 2, 3] -> [2, 0, 3, 1], █ █ █ = 0 0 1
                     *                               █ █ □   1 1
                     */
                    switches[0] = 0, switches[1] = 1, switches[2] = 0, switches[3] = 1, switches[4] = 1;
                }
            } else if (dest[1] == src[1]) {
                // [0, 1, 2, 3] -> [2, 1, ?, ?]
                if (dest[2] == src[0]) {
                    /*
                     * [0, 1, 2, 3] -> [2, 1, 0, 3], █ █ █ = 0 1 0
                     *                               █ █ □   0 0
                     */
                    switches[0] = 0, switches[1] = 0, switches[2] = 1, switches[3] = 0, switches[4] = 0;
                } else {
                    assert(dest[2] == src[3]);
                    /*
                     * [0, 1, 2, 3] -> [2, 1, 3, 0], █ █ █ = 1 0 1
                     *                               █ █ □   1 1
                     */
                    switches[0] = 1, switches[1] = 1, switches[2] = 0, switches[3] = 1, switches[4] = 1;
                }
            } else {
                assert(dest[1] == src[3]);
                // [0, 1, 2, 3] -> [2, 3, ?, ?]
                if (dest[2] == src[0]) {
                    /*
                     * [0, 1, 2, 3] -> [2, 3, 0, 1], █ █ █ = 0 1 0
                     *                               █ █ □   0 1
                     */
                    switches[0] = 0, switches[1] = 0, switches[2] = 1, switches[3] = 1, switches[4] = 0;
                } else {
                    assert(dest[2] == src[1]);
                    /*
                     * [0, 1, 2, 3] -> [2, 3, 1, 0], █ █ █ = 1 1 0
                     *                               █ █ □   0 1
                     */
                    switches[0] = 1, switches[1] = 0, switches[2] = 1, switches[3] = 1, switches[4] = 0;
                }
            }
        } else {
            assert(dest[0] == src[3]);
            // [0, 1, 2, 3] -> [3, ?, ?, ?]
            if (dest[1] == src[0]) {
                // [0, 1, 2, 3] -> [3, 0, ?, ?]
                if (dest[2] == src[1]) {
                    /*
                     * [0, 1, 2, 3] -> [3, 0, 1, 2], █ █ █ = 1 1 0
                     *                               █ █ □   1 0
                     */
                    switches[0] = 1, switches[1] = 1, switches[2] = 1, switches[3] = 0, switches[4] = 0;
                } else {
                    assert(dest[2] == src[2]);
                    /*
                     * [0, 1, 2, 3] -> [3, 0, 2, 1], █ █ █ = 0 0 1
                     *                               █ █ □   0 1
                     */
                    switches[0] = 0, switches[1] = 0, switches[2] = 0, switches[3] = 1, switches[4] = 1;
                }
            } else if (dest[1] == src[1]) {
                // [0, 1, 2, 3] -> [3, 1, ?, ?]
                if (dest[2] == src[0]) {
                    /*
                     * [0, 1, 2, 3] -> [3, 1, 0, 2], █ █ █ = 0 1 0
                     *                               █ █ □   1 0
                     */
                    switches[0] = 0, switches[1] = 1, switches[2] = 1, switches[3] = 0, switches[4] = 0;
                } else {
                    assert(dest[2] == src[2]);
                    /*
                     * [0, 1, 2, 3] -> [3, 1, 2, 0], █ █ █ = 1 0 1
                     *                               █ █ □   0 1
                     */
                    switches[0] = 1, switches[1] = 0, switches[2] = 0, switches[3] = 1, switches[4] = 1;
                }
            } else {
                assert(dest[1] == src[2]);
                // [0, 1, 2, 3] -> [3, 2, ?, ?]
                if (dest[2] == src[0]) {
                    /*
                     * [0, 1, 2, 3] -> [3, 2, 0, 1], █ █ █ = 0 1 0
                     *                               █ █ □   1 1
                     */
                    switches[0] = 0, switches[1] = 1, switches[2] = 1, switches[3] = 1, switches[4] = 0;
                } else {
                    assert(dest[2] == src[1]);
                    /*
                     * [0, 1, 2, 3] -> [3, 2, 1, 0], █ █ █ = 1 1 1
                     *                               █ █ □   0 1
                     */
                    switches[0] = 1, switches[1] = 0, switches[2] = 1, switches[3] = 1, switches[4] = 1;
                }
            }
        }
    }

    void gen_waksman_route(int32_t subLogN, int32_t lvl_p, int32_t perm_idx, const std::vector<int32_t> &src, const std::vector<int32_t> &dest)
    {
        auto subN = (int32_t)src.size();
        if (subN == 2) {
            assert(subLogN == 1 || subLogN == 2);
            if (subLogN == 1) {
                // logN == 1, we have 2 * log(N) - 1 = 1 level (█)
                waksman_network[lvl_p][perm_idx] = (int8_t)(src[0] != dest[0]);
            } else {
                // logN == 2，we have 2 * logN - 1 = 3 levels (□ █ □).
                waksman_network[lvl_p][perm_idx] = 2;
                waksman_network[lvl_p + 1][perm_idx] = (int8_t)(src[0] != dest[0]);
                waksman_network[lvl_p + 2][perm_idx] = 2;
            }
        } else if (subN == 3) {
            assert(subLogN == 2);
            if (src[0] == dest[0]) {
                /*
                 * 0 -> 0，1 -> 1，2 -> 2, the benes_network is:
                 * █ □ █ = 0   0
                 * □ █ □     0
                 *
                 * 0 -> 0，1 -> 2，2 -> 1, the benes_network is:
                 * █ □ █ = 0   0
                 * □ █ □     1
                 */
                waksman_network[lvl_p][perm_idx] = (int8_t)0;
                waksman_network[lvl_p + 2][perm_idx] = (int8_t)0;
                if (src[1] == dest[1]) {
                    waksman_network[lvl_p + 1][perm_idx] = (int8_t)0;
                } else {
                    waksman_network[lvl_p + 1][perm_idx] = (int8_t)1;
                }
            } else if (src[0] == dest[1]) {
                /*
                 * 0 -> 1，1 -> 0，2 -> 2, the benes_network is:
                 * █ □ █ = 0   1
                 * □ █ □     0
                 *
                 * 0 -> 1，1 -> 2，2 -> 0, the benes_network is:
                 * █ □ █ = 0   1
                 * □ █ □     1
                 */
                waksman_network[lvl_p][perm_idx] = (int8_t)0;
                waksman_network[lvl_p + 2][perm_idx] = (int8_t)1;
                if (src[1] == dest[0]) {
                    waksman_network[lvl_p + 1][perm_idx] = (int8_t)0;
                } else {
                    waksman_network[lvl_p + 1][perm_idx] = (int8_t)1;
                }
            } else {
                /*
                 * 0 -> 2，1 -> 0，2 -> 1, the benes_network is:
                 * █ □ █ = 1   0
                 * □ █ □     1
                 *
                 * 0 -> 2，1 -> 1，2 -> 0, the benes_network is:
                 * █ □ █ = 1   1
                 * □ █ □     1
                 */
                waksman_network[lvl_p][perm_idx] = (int8_t)1;
                waksman_network[lvl_p + 1][perm_idx] = (int8_t)1;
                if (src[1] == dest[0]) {
                    waksman_network[lvl_p + 2][perm_idx] = (int8_t)0;
                } else {
                    waksman_network[lvl_p + 2][perm_idx] = (int8_t)1;
                }
            }
            return;
        } else if (subN == 4) {
            assert(subLogN == 2 || subLogN == 3);
            auto *switches = new int8_t[5];
            gen_quadruple_switches(switches, src, dest);
            if (subLogN == 2) {
                waksman_network[lvl_p][perm_idx] = switches[0];
                waksman_network[lvl_p][perm_idx + 1] = switches[1];
                waksman_network[lvl_p + 1][perm_idx] = switches[2];
                waksman_network[lvl_p + 1][perm_idx + 1] = switches[3];
                waksman_network[lvl_p + 2][perm_idx] = switches[4];
                waksman_network[lvl_p + 2][perm_idx + 1] = 2;
            } else {
                waksman_network[lvl_p][perm_idx] = switches[0];
                waksman_network[lvl_p][perm_idx + 1] = switches[1];
                waksman_network[lvl_p + 1][perm_idx] = 2;
                waksman_network[lvl_p + 1][perm_idx + 1] = 2;
                waksman_network[lvl_p + 2][perm_idx] = switches[2];
                waksman_network[lvl_p + 2][perm_idx + 1] = switches[3];
                waksman_network[lvl_p + 3][perm_idx] = 2;
                waksman_network[lvl_p + 3][perm_idx + 1] = 2;
                waksman_network[lvl_p + 4][perm_idx] = switches[4];
                waksman_network[lvl_p + 4][perm_idx + 1] = 2;
            }
            delete[] switches;
        } else {
            int32_t i, j, x;
            uint8_t s;
            int32_t subLevel = 2 * subLogN - 1;
            // top subnetwork map, with size Math.floor(n / 2)
            std::vector<int32_t> topSrc(0);
            std::vector<int32_t> topDest(subN / 2);
            // bottom subnetwork map, with size Math.ceil(n / 2)
            std::vector<int32_t> bottomSrc(0);
            std::vector<int32_t> bottomDest(int(ceil(subN * 0.5)));
            // create forward/backward lookup tables
            // subSrcList stores the position map. For example, src = [2, 4, 6], dest = [6, 4, 2].
            // We re-organize the map to the form [0, subN - 1) -> [0, subN - 1)
            for (i = 0; i < subN; ++i) {
                waksman_inv_perm[src[i]] = i;
            }
            for (i = 0; i < subN; ++i) {
                waksman_perm[i] = waksman_inv_perm[dest[i]];
            }
            for (i = 0; i < subN; ++i) {
                waksman_inv_perm[waksman_perm[i]] = i;
            }
            // shorten the array
            waksman_path.resize(subN);
            // path, initialized by -1, we use 2 for empty node
            std::fill(waksman_path.begin(), waksman_path.end(), (int8_t)-1);
            if (subN % 2 == 1) {
                // handling odd n, the last node directly links to the bottom subnetwork.
                waksman_path[subN - 1] = (int8_t)1;
                waksman_path[waksman_perm[subN - 1]] = (int8_t)1;
                // if values - 1 == benes_perm[values - 1], then the last one is also a direct link. Handle other cases.
                if (waksman_perm[subN - 1] != subN - 1) {
                    int32_t idx = waksman_perm[waksman_inv_perm[subN - 1] ^ 1];
                    waksman_depth_first_search(idx);
                }
            } else {
                // handling even n
                waksman_even_depth_first_search();
            }
            // set other switches
            for (i = 0; i < subN; ++i) {
                if (waksman_path[i] < 0) {
                    waksman_depth_first_search(i);
                }
            }
            // create left part of the network.
            for (i = 0; i < subN - 1; i += 2) {
                waksman_network[lvl_p][perm_idx + i / 2] = waksman_path[i];
                for (j = 0; j < 2; ++j) {
                    x = waksman_right_cycle_shift((i | j) ^ waksman_path[i], subLogN);
                    if (x < subN / 2) {
                        topSrc.push_back(src[i | j]);
                    } else {
                        bottomSrc.push_back(src[i | j]);
                    }
                }
            }
            if (subN % 2 == 1) {
                // add one more switch for the odd case.
                bottomSrc.push_back(src[subN - 1]);
            }
            // create right part of the subnetwork.
            for (i = 0; i < subN - 1; i += 2) {
                s = waksman_network[lvl_p + subLevel - 1][perm_idx + i / 2] = waksman_path[waksman_perm[i]];
                for (j = 0; j < 2; ++j) {
                    x = waksman_right_cycle_shift((i | j) ^ s, subLogN);
                    if (x < subN / 2) {
                        topDest[i / 2] = src[waksman_perm[i | j]];
                    } else {
                        bottomDest[i / 2] = src[waksman_perm[i | j]];
                    }
                }
            }
            if (subN % 2 == 1) {
                // add one more switch for the odd case.
                bottomDest[subN / 2] = dest[subN - 1];
            } else {
                // remove one switch for the even case.
                waksman_network[lvl_p + subLevel - 1][perm_idx + subN / 2 - 1] = 2;
            }
            // create top subnetwork, with (log(N) - 1) levels
            gen_waksman_route(subLogN - 1, lvl_p + 1, perm_idx, topSrc, topDest);
            // create bottom subnetwork with (log(N) - 1) levels.
            gen_waksman_route(subLogN - 1, lvl_p + 1, perm_idx + subN / 4, bottomSrc, bottomDest);
        }
    }
};

void genPermutation(u32 size, std::vector<u32> &pi)
{
    pi.resize(size);
    for (size_t i = 0; i < pi.size(); ++i) {
        pi[i] = i;
    }
    std::shuffle(pi.begin(), pi.end(), global_built_in_prg2);
    return;
}

void permute(std::vector<u32> pi, std::vector<block> &data)
{
    std::vector<block> res(data.size());
    for (size_t i = 0; i < pi.size(); ++i) {
        res[i] = data[pi[i]];
    }
    data.assign(res.begin(), res.end());
}

void senderROT(u32 numOTs, PRNG &prng, BitVector &send0, BitVector &send1, Socket &chl, u32 numThreads)
{
    // The protocol receiver acts as the sender in ROT.
    send0.resize(numOTs);
    send1.resize(numOTs);

    SilentOtExtSender sender;
    sender.configure(numOTs, numThreads, 128);

    // sender.configure(numOTs, 40, 1, SilentSecType::SemiHonest);

    std::vector<std::array<block, 2>> sendMsg(numOTs); // Stores the sender output blocks.

    // coproto::sync_wait(sender.genBaseCors(std::nullopt, prng, chl, true));// Run base OT.
    coproto::sync_wait(sender.genBaseOts(prng, chl)); // Run base OT.

    coproto::sync_wait(sender.silentSend(sendMsg, prng, chl)); // Run silent OT.

    for (u32 i = 0; i < numOTs; i++) {
        const unsigned char *bytes = sendMsg[i][0].data();
        bool bit = (bytes[0] >> 0) & 1;
        send0[i] = bit;

        bytes = sendMsg[i][1].data();
        bit = (bytes[0] >> 0) & 1;
        send1[i] = bit;
    }
    BitVector d;
    d.resize(numOTs);

    coproto::sync_wait(chl.recv(d));
    for (u32 i = 0; i < numOTs; i++) {
        if (d[i] == 1) {
            send0[i] ^= send1[i];
            send1[i] ^= send0[i];
            send0[i] ^= send1[i];
        }
    }
}

void receiverROT(u32 numOTs, PRNG &prng, BitVector &bitV, BitVector &recv, Socket &chl, u32 numThreads)
{
    // The protocol sender acts as the receiver in ROT.
    bitV.resize(numOTs);
    recv.resize(numOTs);
    BitVector bitV0 = bitV;
    SilentOtExtReceiver receiver;
    receiver.configure(numOTs, numThreads, 128);
    // receiver.configure(numOTs, 40, 1, SilentSecType::SemiHonest);

    std::vector<block> recvMsg(numOTs); // Stores the receiver output blocks.

    // coproto::sync_wait(receiver.genBaseCors(std::nullopt, prng, chl, true));// Run base OT.
    coproto::sync_wait(receiver.genBaseOts(prng, chl));                                   // Run base OT.
    coproto::sync_wait(receiver.silentReceive(bitV, recvMsg, prng, chl, OTType::Random)); // Run silent OT.

    for (u32 i = 0; i < numOTs; i++) {
        const unsigned char *bytes = recvMsg[i].data();
        bool bit = (bytes[0] >> 0) & 1;
        recv[i] = bit;
    }
    BitVector diff = bitV ^ bitV0;
    coproto::sync_wait(chl.send(diff));
    bitV = bitV0;
}

void perm_cons(u8 nn, std::vector<u32> &in, std::vector<u32> pi, BitVector &switch_bit)
{
    u32 elem_num = 1 << nn;
    u32 switch_num = 1 << (nn - 1);
    std::vector<u32> out;
    out.resize(elem_num);
    for (u32 i = 0; i < elem_num; i++) {
        out[i] = pi[i];
    }
    // std::vector<uint8_t> switch_bit(switch_num);
    switch_bit.resize((2 * nn - 1) * (1 << (nn - 1)));
    u32 switchnum_inner;
    u32 j_num;
    u32 jloop;
    std::vector<u32> all_l(elem_num);
    std::vector<u32> all_r(elem_num);

    std::vector<u32> in1(elem_num);
    std::vector<u32> out1(elem_num);

    u32 midpoint; // midpoint is the position of the current node, with values from 0 to n.
    // nn is the logarithm of the number of elements, pi is the permutation, and switch_bit stores each switch bit.
    for (u8 i = 0; i < nn; i++) {
        switchnum_inner = (1 << (nn - 1 - i));
        j_num = (1 << (nn - i));
        jloop = (1 << i);
        std::vector<bool> allocated(1 << nn);

        // midpoint=in[0];
        // printf("current i=%d, midpoint=%d\n",static_cast<int>(i),midpoint);
        // Compute both all_l and all_r.
        if (i == 0) { // First iteration: initialize all arrays with the default in and out.
            for (int m = 0; m < (elem_num); m++) {
                all_r[out[m]] = m; // all_r gives the position of the right-side external element.
                all_l[in[m]] = m;  // all_l gives the position of the left-side external element.
            }

        } else {
            // for (int q = 0; q < (elem_num); q++)
            // {

            //     in[q]=in1[q];
            //     out[q]=out1[q];
            // }
            memcpy(in.data(), in1.data(), (elem_num) * sizeof(int));
            memcpy(out.data(), out1.data(), (elem_num) * sizeof(int));

            // In later iterations, update all_l and all_r with the in1 and out1 from the previous iteration.
            for (int m = 0; m < (elem_num); m++) {
                all_r[out[m]] = m;
                all_l[in[m]] = m;
            }
        }
        for (int j = 0; j < (jloop); j++) {
            if (i == (nn - 1)) { // Last layer: determine the switch bit by checking equality directly.
                if (in[2 * j] == out[2 * j]) {
                    switch_bit[(nn - 1) * (switch_num) + j] = 0; // This value needs to be changed.
                } else {
                    switch_bit[(nn - 1) * (switch_num) + j] = 1;
                }

            } else {
                bool dump = 0;
                // std::vector<bool> allocated(switchnum_inner);
                // std::vector<bool> allocated(1 << nn);

                for (int k = 0; k < (switchnum_inner); k++) {
                    if (k == 0) { // First iteration: set the first switch to 0 (off) by default.
                        midpoint = in[j * (j_num)];
                        allocated[midpoint] = 1;
                        switch_bit[i * (switch_num) + j * (switchnum_inner)] = 0; // The first switch defaults to off, so set it to 0.
                        in1[all_l[midpoint]] = midpoint;
                        in1[all_l[midpoint] + (switchnum_inner)] = in[all_l[midpoint] + 1]; // Recently changed here.
                        allocated[in[all_l[midpoint] + 1]] = 1;
                        // Here.
                        if (all_r[midpoint] % 2 == 1) {
                            switch_bit[(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] =
                                1; //+j*(1 << (nn-1-i))
                            out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + (all_r[midpoint] % (j_num)) / 2] = midpoint;
                            midpoint = out[all_r[midpoint] - 1];
                            out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;

                        } else {
                            switch_bit[(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] =
                                0; // Here.
                            out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + (all_r[midpoint] % (j_num)) / 2] = midpoint;
                            midpoint = out[all_r[midpoint] + 1];
                            out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;
                        }

                    } else {
                        while (allocated[midpoint] == 1) // This should run after the right side ends, not before the left side starts.
                        {
                            midpoint = in[j * (j_num) + ((all_l[midpoint] - j * (j_num)) + 2) % (j_num)];
                            dump = 1;
                        }

                        if (dump == 1) { // A redundant loop occurred and must be corrected.
                            dump = 0;    // Turn the loop flag off again.
                            allocated[midpoint] = 1;
                            switch_bit[i * (switch_num) + j * (switchnum_inner) + ((all_l[midpoint] % (j_num)) / 2)] = 0;
                            in1[j * (j_num) + (all_l[midpoint] - j * (j_num)) / 2 + (switchnum_inner)] = midpoint;
                            midpoint = in[all_l[midpoint] - 1]; // Logic for this case.
                            in1[j * (j_num) + (all_l[midpoint] - j * (j_num)) / 2] = midpoint;
                            allocated[midpoint] = 1;

                            if (all_r[midpoint] % 2 == 1) {
                                switch_bit
                                    [(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] =
                                        1; //+j*(1 << (nn-1-i))
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + (all_r[midpoint] % (j_num)) / 2] = midpoint;
                                midpoint = out[all_r[midpoint] - 1];
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;

                            } else {
                                switch_bit
                                    [(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] =
                                        0; // Here.
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + (all_r[midpoint] % (j_num)) / 2] = midpoint;
                                midpoint = out[all_r[midpoint] + 1];
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;
                            }

                        } else {
                            if (all_l[midpoint] % 2 == 1) {
                                switch_bit[i * (switch_num) + j * (switchnum_inner) + ((all_l[midpoint] - j * (j_num)) / 2)] = 0;
                                in1[all_l[midpoint] - (all_l[midpoint] % (j_num)) + ((all_l[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint; // zheli
                                allocated[midpoint] = 1;
                                midpoint = in[all_l[midpoint] - 1];
                                allocated[midpoint] = 1;
                                in1[all_l[midpoint] - (all_l[midpoint] % (j_num)) + ((all_l[midpoint] % (j_num)) / 2)] = midpoint;

                            } else {
                                switch_bit[i * (switch_num) + j * (switchnum_inner) + ((all_l[midpoint] - j * (j_num)) / 2)] = 1;
                                in1[all_l[midpoint] - (all_l[midpoint] % (j_num)) + ((all_l[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;
                                allocated[midpoint] = 1;
                                midpoint = in[all_l[midpoint] + 1];
                                allocated[midpoint] = 1;
                                in1[(all_l[midpoint] - (all_l[midpoint] % (j_num)) + ((all_l[midpoint] % (j_num)) / 2))] = midpoint; // There is an r here.
                            }

                            if (all_r[midpoint] % 2 == 0) { // zheli
                                switch_bit
                                    [(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] = 0;
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2)] = midpoint;
                                midpoint = out[all_r[midpoint] + 1];
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;

                            } else {
                                switch_bit
                                    [(2 * nn - 1) * (switch_num) - (i + 1) * (switch_num) + j * (switchnum_inner) + ((all_r[midpoint] - j * (j_num)) / 2)] = 1;
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2)] = midpoint;
                                midpoint = out[all_r[midpoint] - 1];
                                out1[all_r[midpoint] - (all_r[midpoint] % (j_num)) + ((all_r[midpoint] % (j_num)) / 2) + (switchnum_inner)] = midpoint;
                            }
                        }
                    }
                }
            }
        }
        // for (int q = 0; q < (elem_num); q++)
        // {
        //     // in[q]=in1[q];
        //     // out[q]=out1[q];
        // }
    }
}

void Opermute(u32 idx, BitVector &input, BitVector &output, std::vector<u32> &pi, Socket &chl, u32 numThreads)
{
    // Combines the functions above to build the permutation network and run the two-party interaction.
    u32 num = input.size();                   // Number of useful elements in the permutation.
    u8 nn = std::bit_width(input.size() - 1); // nn is the permutation network parameter.
    u32 elem_num = 1 << nn;
    u32 switch_num = (2 * nn - 1) * (1 << (nn - 1)); // Total number of switches.
    input.resize(elem_num);                          // Extend both input and output to 2^nn elements.
    // std::cout<<"forrmal input:"<<idx<<" :"<<input<<std::endl;

    output.resize(elem_num);
    BitVector input0;
    input0.resize(elem_num);
    BitVector d(switch_num);
    if (idx == 0) { // sender
        // Generate the permutation.

        std::vector<u32> norm; // Original order: 0, 1, 2, ...
        BitVector reality;     // Records which elements are useful bits.
        norm.resize(elem_num);
        reality.resize(elem_num);
        pi.resize(elem_num);
        for (u32 i = 0; i < elem_num; i++) { // Permutation input and output.
            norm[i] = i;
            pi[i] = i;
        }

        std::shuffle(pi.begin(), pi.end(), global_built_in_prg2); // Generate a random permutation.
        // std::cout<<"pi0:"<<std::endl;
        // for(u32 i =0;i<elem_num;i++){
        //     std::cout<<pi[i]<<std::endl;

        // }

        for (u32 i = 0; i < elem_num; i++) { // Assign values to reality.
            if (pi[i] < num) {
                reality[i] = 1;
            }
        }
        // std::cout<<"Sreality:"<<reality<<std::endl;

        // Build the network.
        BitVector network_choice;
        Timer timer1;
        timer1.setTimePoint("start");
        perm_cons(nn, norm, pi, network_choice); // Build the network and obtain each switch choice bit.
        timer1.setTimePoint("ours");
        std::cout << timer1 << std::endl;

        // std::cout<<"pi permutation network:"<<network_choice<<std::endl;
        // std::cout<<"pi1:"<<std::endl;
        // for(u32 i =0;i<elem_num;i++){
        //     std::cout<<pi[i]<<std::endl;

        // }

        // Run ROT for each switch.
        bool isSender = (idx == 1); // S is the OT receiver.
        PRNG prng(sysRandomSeed());
        BitVector bitV, recv;
        bitV = network_choice;
        // std::cout<<"sender choice bit:"<<bitV<<std::endl;

        receiverROT(switch_num, prng, bitV, recv, chl, numThreads);
        // std::cout<<"sender rc:"<<recv<<std::endl;

        // Compute the final output locally.
        coproto::sync_wait(chl.recv(d));
        // std::cout<<"sender d:"<<d<<std::endl;

        // std::cout<<"S111"<<std::endl;

        for (u8 i = 0; i < (nn - 1); i++) { // From left to middle.
            u32 group_num = 1 << i;
            for (u32 j = 0; j < group_num; j++) {
                u32 switch_num0 = 1 << (nn - 1 - i);
                u32 elem_num0 = 1 << (nn - i);
                for (u32 k = 0; k < switch_num0; k++) {
                    // std::cout<<"index value"<<(i*(elem_num/2)+j*switch_num0+k)<<std::endl;

                    if (bitV[i * (elem_num / 2) + j * switch_num0 + k] == 0) {
                        input0[j * elem_num0 + k] = input[j * elem_num0 + 2 * k] ^ recv[i * (elem_num / 2) + j * switch_num0 + k];
                        input0[j * elem_num0 + k + switch_num0] =
                            input[j * elem_num0 + 2 * k + 1] ^ recv[i * (elem_num / 2) + j * switch_num0 + k] ^ d[i * (elem_num / 2) + j * switch_num0 + k];
                    } else if (bitV[i * (elem_num / 2) + j * switch_num0 + k] == 1) {
                        input0[j * elem_num0 + k] =
                            input[j * elem_num0 + 2 * k + 1] ^ recv[i * (elem_num / 2) + j * switch_num0 + k] ^ d[i * (elem_num / 2) + j * switch_num0 + k];
                        input0[j * elem_num0 + k + switch_num0] = input[j * elem_num0 + 2 * k] ^ recv[i * (elem_num / 2) + j * switch_num0 + k];
                        // std::cout<<"       1209"<<std::endl;
                    }
                }
            }
            input = input0;
            // std::cout<<"S input state:"<<input<<std::endl;
        }
        // std::cout<<"Sinput0:"<<input<<std::endl;

        for (u32 i = 0; i < (elem_num / 2); i++) { // Middle.
            // std::cout<<"i:"<<i<<std::endl;

            if (bitV[(nn - 1) * (elem_num / 2) + i] == 0) {
                input0[(i << 1)] = input[i << 1] ^ recv[(nn - 1) * (elem_num / 2) + i];
                input0[(i << 1) + 1] = input[(i << 1) + 1] ^ recv[(nn - 1) * (elem_num / 2) + i] ^ d[(nn - 1) * (elem_num / 2) + i];
            } else if (bitV[(nn - 1) * (elem_num / 2) + i] == 1) {
                input0[(i << 1)] = input[(i << 1) + 1] ^ recv[(nn - 1) * (elem_num / 2) + i] ^ d[(nn - 1) * (elem_num / 2) + i];
                input0[(i << 1) + 1] = input[i << 1] ^ recv[(nn - 1) * (elem_num / 2) + i];
            }
        }
        input = input0;
        // std::cout<<"S input state:"<<input<<std::endl;

        // std::cout<<"Sinput1:"<<input<<std::endl;

        for (u8 i = 0; i < (nn - 1); i++) { // From middle to right.
            u32 group_num = 1 << (nn - 2 - i);
            for (u32 j = 0; j < group_num; j++) {
                u32 switch_num0 = 1 << (i + 1);
                u32 elem_num0 = 1 << (i + 2);
                for (u32 k = 0; k < switch_num0; k++) {
                    if (bitV[(nn + i) * (elem_num / 2) + j * switch_num0 + k] == 0) {
                        input0[j * elem_num0 + 2 * k] = input[j * elem_num0 + k] ^ recv[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                        input0[j * elem_num0 + 2 * k + 1] = input[j * elem_num0 + k + switch_num0] ^ recv[(nn + i) * (elem_num / 2) + j * switch_num0 + k] ^
                            d[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                    } else if (bitV[(nn + i) * (elem_num / 2) + j * switch_num0 + k] == 1) {
                        input0[j * elem_num0 + 2 * k] = input[j * elem_num0 + k + switch_num0] ^ recv[(nn + i) * (elem_num / 2) + j * switch_num0 + k] ^
                            d[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                        input0[j * elem_num0 + 2 * k + 1] = input[j * elem_num0 + k] ^ recv[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                    }
                }
            }
            input = input0;
            // std::cout<<"S input state:"<<input<<std::endl;
        }

        // output=input;
        coproto::sync_wait(chl.send(reality));

        // std::cout<<"senderoutput:"<<input<<std::endl;
        for (u32 jelly = 0, i = 0; i < elem_num; i++) {
            if (reality[i] == 1) {
                output[jelly] = input[i];
                jelly++;
            }
        }
        output.resize(num);

        std::vector<u32> pi0(num);
        for (u32 i = 0, j = 0; i < elem_num; i++) { // Resize pi from the original numBin size to the numElements size.
            if (pi[i] < num) {
                pi0[j] = pi[i];
                j++;
            }
        }

        pi.resize(num);

        for (u32 i = 0; i < num; i++) {
            pi[i] = pi0[i];
        }

    } else if (idx == 1) {
        // Run ROT for each switch.
        bool isSender = (idx == 1); // R is the OT sender.
        PRNG prng(sysRandomSeed());
        BitVector send0, send1;
        senderROT(switch_num, prng, send0, send1, chl, numThreads);

        // Compute the final output locally.
        for (u8 i = 0; i < (nn - 1); i++) { // From left to middle.
            u32 group_num = 1 << i;
            for (u32 j = 0; j < group_num; j++) {
                u32 switch_num0 = 1 << (nn - 1 - i);
                u32 elem_num0 = 1 << (nn - i);
                for (u32 k = 0; k < switch_num0; k++) {
                    input0[j * elem_num0 + k] = input[j * elem_num0 + 2 * k] ^ send0[i * (elem_num / 2) + j * switch_num0 + k];
                    input0[j * elem_num0 + k + switch_num0] = input[j * elem_num0 + 2 * k] ^ send1[i * (elem_num / 2) + j * switch_num0 + k];
                    d[i * (elem_num / 2) + j * switch_num0 + k] =
                        input[j * elem_num0 + 2 * k] ^ input[j * elem_num0 + 2 * k + 1] ^ input0[j * elem_num0 + k] ^ input0[j * elem_num0 + k + switch_num0];
                } // Reached this point.
            }
            input = input0;
            // std::cout<<"R input state:"<<input<<std::endl;
        }

        // std::cout<<"Rinput0:"<<input<<std::endl;

        for (u32 i = 0; i < (elem_num / 2); i++) { // Middle.
            input0[i << 1] = input[i << 1] ^ send0[(nn - 1) * (elem_num / 2) + i];
            input0[(i << 1) + 1] = input[i << 1] ^ send1[(nn - 1) * (elem_num / 2) + i];
            d[(nn - 1) * (elem_num / 2) + i] = input0[i << 1] ^ input0[(i << 1) + 1] ^ input[i << 1] ^ input[(i << 1) + 1];
        }
        input = input0;

        for (u8 i = 0; i < (nn - 1); i++) {    // From middle to right.
            u32 group_num = 1 << (nn - 2 - i); //
            for (u32 j = 0; j < group_num; j++) {
                u32 switch_num0 = 1 << (i + 1); //
                u32 elem_num0 = 1 << (i + 2);   //
                for (u32 k = 0; k < switch_num0; k++) {
                    input0[j * elem_num0 + 2 * k] = input[j * elem_num0 + k] ^ send0[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                    input0[j * elem_num0 + 2 * k + 1] = input[j * elem_num0 + k] ^ send1[(nn + i) * (elem_num / 2) + j * switch_num0 + k];
                    d[(nn + i) * (elem_num / 2) + j * switch_num0 + k] =
                        input0[j * elem_num0 + 2 * k] ^ input0[j * elem_num0 + 2 * k + 1] ^ input[j * elem_num0 + k] ^ input[j * elem_num0 + k + switch_num0];
                }
            }
            input = input0;
            // std::cout<<"R input state:"<<input<<std::endl;
        }

        // std::cout<<"receiver d:"<<d<<std::endl;

        coproto::sync_wait(chl.send(d));
        BitVector reality;
        reality.resize(elem_num);
        coproto::sync_wait(chl.recv(reality));
        // std::cout<<"Routput:"<<input<<std::endl;
        // std::cout<<"Rreality:"<<reality<<std::endl;

        for (u32 jelly = 0, i = 0; i < elem_num; i++) {
            if (reality[i] == 1) {
                output[jelly] = input[i];
                jelly++;
            }
        }
        output.resize(num);

        // output=input;
    }
}

void Opermute1(u32 idx, BitVector &input, BitVector &output, std::vector<u32> &pi, Socket &chl, u32 numThreads)
{
    u32 num = input.size();                     // Number of elements in the permutation.
    u8 nn = std::bit_width(input.size() - 1);   // nn is the upper bound of the logarithm of the number of elements.
    u32 layers = 2 * nn - 1;                    // Total number of network layers.
    u32 switch_num_layer = num / 2;             // Number of switches per layer.
    u32 switch_num = layers * switch_num_layer; // Total number of switches.
    BitVector d(switch_num);                    // d values communicated between sender and receiver.

    pi.resize(num);
    output.resize(num);
    for (u32 i = 0; i < num; i++) { // Permutation input and output.
        pi[i] = i;
    }

    std::shuffle(pi.begin(), pi.end(), global_built_in_prg2); // Generate a random permutation.

    std::vector<int32_t> pi1(pi.begin(), pi.end());
    Timer timer2;
    timer2.setTimePoint("start");
    WaksmanNetwork network(pi1); // Generate the network from the random permutation. The sender uses it for routing; the receiver uses it for indicators: -1, 0, 1, 2.
    auto switches = network.get_waksman_network();
    timer2.setTimePoint("alibaba");
    // std::cout << timer2 << std::endl;

    if (idx == 0) {
        // Run the ROT interaction.

        bool isSender = (idx == 1); // S is the OT receiver.
        PRNG prng(sysRandomSeed());
        BitVector recv;
        BitVector bitV(switch_num);
        // bitV.resize(switch_num);
        for (u32 i = 0; i < layers; i++) {
            for (u32 j = 0; j < switch_num_layer; j++) {
                if (switches[i][j] == 1) {
                    bitV[i * switch_num_layer + j] = 1;
                } else {
                    bitV[i * switch_num_layer + j] = 0;
                }
            }
        }

        // bitV=network_choice;
        // std::cout<<"sender choice bit:"<<bitV<<std::endl;
        receiverROT(switch_num, prng, bitV, recv, chl, numThreads);
        // std::cout<<"sender rc:"<<recv<<std::endl;

        coproto::sync_wait(chl.recv(d));
        // Run the local computation.

        u32 loog = 1;
        u32 start = 0;
        std::vector<u32> group_size(switch_num_layer);
        std::vector<u32> switch_start(switch_num_layer);
        BitVector input0(num);

        // group_size.resize(1);
        // switch_start.resize(1);
        group_size[0] = num;
        switch_start[0] = 0;

        // Compute the first 2logN-2 layers.
        for (u8 i = 0; i < nn - 2; i++) { // One layer.
            start = 0;
            u32 iswitchnum = i * switch_num_layer;
            for (u32 j = 0; j < loog; j++) { // One subnetwork in a layer.
                u32 half_groupsize = group_size[j] / 2;
                for (u32 k = 0; k < half_groupsize; k++) { // Element computation in the subnetwork.

                    if (bitV[iswitchnum + switch_start[j] + k] == 1) {
                        input0[start + k] = input[start + k * 2 + 1] ^ recv[iswitchnum + switch_start[j] + k] ^ d[iswitchnum + switch_start[j] + k];
                        input0[start + k + half_groupsize] = input[start + k * 2] ^ recv[iswitchnum + switch_start[j] + k];
                    } else {
                        input0[start + k] = input[start + k * 2] ^ recv[iswitchnum + switch_start[j] + k];
                        input0[start + k + half_groupsize] =
                            input[start + k * 2 + 1] ^ recv[iswitchnum + switch_start[j] + k] ^ d[iswitchnum + switch_start[j] + k];
                    }
                }
                if (group_size[j] % 2 == 1) {
                    input0[start + group_size[j] - 1] = input[start + group_size[j] - 1];
                }
                start += group_size[j];
            }
            // group_size.resize(loog<<1);
            // switch_start.resize(loog<<1);
            for (int ite = loog - 1; ite >= 0; ite--) {
                // Write group_size.
                group_size[(ite << 1) + 1] = (group_size[ite] + 1) / 2;                        // Upper half after division by 2.
                group_size[ite << 1] = group_size[ite] / 2;                                    // Lower half after division by 2.
                switch_start[(ite << 1) + 1] = switch_start[ite] + (group_size[ite << 1] / 2); // Starting switch position of the right subnetwork.
                switch_start[ite << 1] = switch_start[ite];                                    // Starting switch position of the left subnetwork.
            }

            loog = loog << 1;
            input = input0;
            // input   = std::move(input0);
            // std::swap(input0, input);
        }
        start = 0;

        // cout<<"wuwuwu"<<endl;
        // Compute the middle three layers.
        for (u32 i = 0; i < group_size.size(); i++) {
            u32 nn_2layer = switch_num_layer * (nn - 2);
            u32 nn_1layer = switch_num_layer * (nn - 1);
            u32 nn_0layer = switch_num_layer * nn;

            if (group_size[i] == 4) {
                // First layer.
                if (bitV[nn_2layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 1] ^ recv[nn_2layer + switch_start[i]] ^ d[nn_2layer + switch_start[i]];
                    input0[start + 2] = input[start] ^ recv[nn_2layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_2layer + switch_start[i]];
                    input0[start + 2] = input[start + 1] ^ recv[nn_2layer + switch_start[i]] ^ d[nn_2layer + switch_start[i]];
                }

                if (bitV[nn_2layer + switch_start[i] + 1] == 1) {
                    input0[start + 1] = input[start + 3] ^ recv[nn_2layer + switch_start[i] + 1] ^ d[nn_2layer + switch_start[i] + 1];
                    input0[start + 3] = input[start + 2] ^ recv[nn_2layer + switch_start[i] + 1];
                } else {
                    input0[start + 1] = input[start + 2] ^ recv[nn_2layer + switch_start[i] + 1];
                    input0[start + 3] = input[start + 3] ^ recv[nn_2layer + switch_start[i] + 1] ^ d[nn_2layer + switch_start[i] + 1];
                }
                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];

                // input   = std::move(input0);
                // std::swap(input, input0);

                // Second layer.
                if (bitV[nn_1layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 1] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                    input0[start + 1] = input[start] ^ recv[nn_1layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_1layer + switch_start[i]];
                    input0[start + 1] = input[start + 1] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                }

                if (bitV[nn_1layer + switch_start[i] + 1] == 1) {
                    input0[start + 2] = input[start + 3] ^ recv[nn_1layer + switch_start[i] + 1] ^ d[nn_1layer + switch_start[i] + 1];
                    input0[start + 3] = input[start + 2] ^ recv[nn_1layer + switch_start[i] + 1];
                } else {
                    input0[start + 2] = input[start + 2] ^ recv[nn_1layer + switch_start[i] + 1];
                    input0[start + 3] = input[start + 3] ^ recv[nn_1layer + switch_start[i] + 1] ^ d[nn_1layer + switch_start[i] + 1];
                }

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // Third layer.
                if (bitV[nn_0layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 2] ^ recv[nn_0layer + switch_start[i]] ^ d[nn_0layer + switch_start[i]];
                    input0[start + 1] = input[start] ^ recv[nn_0layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_0layer + switch_start[i]];
                    input0[start + 1] = input[start + 2] ^ recv[nn_0layer + switch_start[i]] ^ d[nn_0layer + switch_start[i]];
                }

                input0[start + 2] = input[start + 1];
                input0[start + 3] = input[start + 3]; // In general, this switch can be ignored and handled as a direct connection.

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];
                // input   = std::move(input0);
                // std::swap(input, input0);

            } else if (group_size[i] == 3) {
                // First layer.
                if (bitV[nn_2layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 1] ^ recv[nn_2layer + switch_start[i]] ^ d[nn_2layer + switch_start[i]];
                    input0[start + 1] = input[start] ^ recv[nn_2layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_2layer + switch_start[i]];
                    input0[start + 1] = input[start + 1] ^ recv[nn_2layer + switch_start[i]] ^ d[nn_2layer + switch_start[i]];
                }

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // Second layer.
                if (bitV[nn_1layer + switch_start[i]] == 1) {
                    input0[start + 1] = input[start + 2] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                    input0[start + 2] = input[start + 1] ^ recv[nn_1layer + switch_start[i]];
                } else {
                    input0[start + 1] = input[start + 1] ^ recv[nn_1layer + switch_start[i]];
                    input0[start + 2] = input[start + 2] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                }

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // Third layer.
                if (bitV[nn_0layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 1] ^ recv[nn_0layer + switch_start[i]] ^ d[nn_0layer + switch_start[i]];
                    input0[start + 1] = input[start] ^ recv[nn_0layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_0layer + switch_start[i]];
                    input0[start + 1] = input[start + 1] ^ recv[nn_0layer + switch_start[i]] ^ d[nn_0layer + switch_start[i]];
                }

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

            } else if (group_size[i] == 2) { // Modify the third-layer input directly, while d still needs to be modified in the middle layer.
                if (bitV[nn_1layer + switch_start[i]] == 1) {
                    input0[start] = input[start + 1] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                    input0[start + 1] = input[start] ^ recv[nn_1layer + switch_start[i]];
                } else {
                    input0[start] = input[start] ^ recv[nn_1layer + switch_start[i]];
                    input0[start + 1] = input[start + 1] ^ recv[nn_1layer + switch_start[i]] ^ d[nn_1layer + switch_start[i]];
                }

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                // input   = std::move(input0);
                // std::swap(input, input0);
            }
            start += group_size[i];
        }
        // cout<<"wuwuwu"<<endl;

        // Compute the last 2logN-2 layers.
        loog = loog / 2;
        for (u8 i = 0; i < nn - 2; i++) { // Compute the second half of the network.
            for (int ite = 0; ite < loog; ite++) {
                // Write group_size.
                group_size[ite] = group_size[ite << 1] + group_size[(ite << 1) + 1]; // Upper half after division by 2.
                switch_start[ite] = switch_start[ite << 1];                          // Starting switch position of the left subnetwork.
            }
            // group_size.resize(loog);
            // switch_start.resize(loog);

            start = 0;
            u32 nn_1ilayer = (nn + 1 + i) * switch_num_layer;
            for (u32 j = 0; j < loog; j++) {
                for (u32 k = 0; k < group_size[j] / 2; k++) { // Element computation in the subnetwork.
                    if (bitV[nn_1ilayer + switch_start[j] + k] == 1) {
                        input0[start + 2 * k] =
                            input[start + k + (group_size[j] / 2)] ^ recv[nn_1ilayer + switch_start[j] + k] ^ d[nn_1ilayer + switch_start[j] + k];
                        input0[start + 2 * k + 1] = input[start + k] ^ recv[nn_1ilayer + switch_start[j] + k];
                    } else {
                        input0[start + 2 * k] = input[start + k] ^ recv[nn_1ilayer + switch_start[j] + k];
                        input0[start + 2 * k + 1] =
                            input[start + k + (group_size[j] / 2)] ^ recv[nn_1ilayer + switch_start[j] + k] ^ d[nn_1ilayer + switch_start[j] + k];
                    }
                }
                if (group_size[j] % 2 == 1) {
                    input0[start + group_size[j] - 1] = input[start + group_size[j] - 1];
                }
                start += group_size[j];
            }
            loog = loog / 2;
            input = input0;
            // input   = std::move(input0);
            // std::swap(input, input0);
        }
        output = input;

    } else if (idx == 1) {
        // Run the ROT interaction.
        bool isSender = (idx == 1); // R is the OT sender.
        PRNG prng(sysRandomSeed());
        BitVector send0, send1;
        senderROT(switch_num, prng, send0, send1, chl, numThreads);

        // Run the local computation.

        // Compute the first 2logN-2 layers.
        u32 loog = 1;
        u32 start = 0;
        std::vector<u32> group_size(switch_num_layer);
        std::vector<u32> switch_start(switch_num_layer);
        BitVector input0(num);

        // group_size.resize(1);
        // switch_start.resize(1);
        group_size[0] = num;
        switch_start[0] = 0;
        // cout<<"wuwuwu"<<endl;

        for (u8 i = 0; i < nn - 2; i++) { // One layer.
            start = 0;

            for (u32 j = 0; j < loog; j++) { // One subnetwork in a layer.
                u32 half_group = group_size[j] / 2;
                for (u32 k = 0; k < half_group; k++) {                                                                         // Element computation in the subnetwork.
                    input0[start + k] = input[start + k * 2] ^ send0[switch_num_layer * i + switch_start[j] + k];              // u0 = x0 − r0,
                    input0[start + k + half_group] = input[start + k * 2] ^ send1[switch_num_layer * i + switch_start[j] + k]; // u1 = x0 − r1
                    d[switch_num_layer * i + switch_start[j] + k] = input[start + k * 2] ^ input[start + k * 2 + 1] ^
                        send0[switch_num_layer * i + switch_start[j] + k] ^ send1[switch_num_layer * i + switch_start[j] + k]; // d = (x0 + x1) − (u0 + u1)
                }
                if (group_size[j] % 2 == 1) {
                    input0[start + group_size[j] - 1] = input[start + group_size[j] - 1];
                }
                start += group_size[j];
            }
            // group_size.resize(loog<<1);
            // switch_start.resize(loog<<1);
            for (int ite = loog - 1; ite >= 0; ite--) {
                // Write group_size.
                group_size[(ite << 1) + 1] = (group_size[ite] + 1) / 2;                        // Upper half after division by 2.
                group_size[ite << 1] = group_size[ite] / 2;                                    // Lower half after division by 2.
                switch_start[(ite << 1) + 1] = switch_start[ite] + (group_size[ite << 1] / 2); // Starting switch position of the right subnetwork.
                switch_start[ite << 1] = switch_start[ite];                                    // Starting switch position of the left subnetwork.
            }
            loog = loog << 1;
            input = input0;
            // input   = std::move(input0);
            // std::swap(input, input0);
        }
        start = 0;
        // cout<<"wuwuwu"<<endl;

        // Compute the middle three layers. Remove the outer recursion and write this part in as much detail as possible.
        for (u32 i = 0; i < group_size.size(); i++) {
            u32 nn_2layer = switch_num_layer * (nn - 2);
            u32 nn_1layer = switch_num_layer * (nn - 1);
            u32 nn_0layer = switch_num_layer * nn;
            if (group_size[i] == 4) {
                // First layer.
                input0[start] = input[start] ^ send0[nn_2layer + switch_start[i]];
                input0[start + 1] = input[start + 2] ^ send0[nn_2layer + switch_start[i] + 1];
                input0[start + 2] = input[start] ^ send1[nn_2layer + switch_start[i]];
                input0[start + 3] = input[start + 2] ^ send1[nn_2layer + switch_start[i] + 1];
                d[nn_2layer + switch_start[i]] = input[start] ^ input[start + 1] ^ send0[nn_2layer + switch_start[i]] ^ send1[nn_2layer + switch_start[i]];
                d[nn_2layer + switch_start[i] + 1] =
                    input[start + 2] ^ input[start + 3] ^ send0[nn_2layer + switch_start[i] + 1] ^ send1[nn_2layer + switch_start[i] + 1];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // Second layer.
                input0[start] = input[start] ^ send0[nn_1layer + switch_start[i]];
                input0[start + 1] = input[start] ^ send1[nn_1layer + switch_start[i]];
                d[nn_1layer + switch_start[i]] = input[start] ^ input[start + 1] ^ send0[nn_1layer + switch_start[i]] ^ send1[nn_1layer + switch_start[i]];

                input0[start + 2] = input[start + 2] ^ send0[nn_1layer + switch_start[i] + 1];
                input0[start + 3] = input[start + 2] ^ send1[nn_1layer + switch_start[i] + 1];
                d[nn_1layer + switch_start[i] + 1] =
                    input[start + 2] ^ input[start + 3] ^ send0[nn_1layer + switch_start[i] + 1] ^ send1[nn_1layer + switch_start[i] + 1];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // Third layer.
                input0[start] = input[start] ^ send0[nn_0layer + switch_start[i]];
                input0[start + 1] = input[start] ^ send1[nn_0layer + switch_start[i]];
                d[nn_0layer + switch_start[i]] = input[start] ^ input[start + 2] ^ send0[nn_0layer + switch_start[i]] ^ send1[nn_0layer + switch_start[i]];

                input0[start + 2] = input[start + 1];
                input0[start + 3] = input[start + 3]; // In general, this switch can be ignored and handled as a direct connection.
                // d[switch_num_layer*nn+switch_start[i]+1]=input[start+2]^input[start+3]^send0[switch_num_layer*nn+switch_start[i]+1]^send1[switch_num_layer*nn+switch_start[i]+1];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                input[start + 3] = input0[start + 3];
                // input   = std::move(input0);
                // std::swap(input, input0);

            } else if (group_size[i] == 3) {
                // First layer.
                input0[start] = input[start] ^ send0[nn_2layer + switch_start[i]];
                input0[start + 1] = input[start] ^ send1[nn_2layer + switch_start[i]];
                d[nn_2layer + switch_start[i]] = input[start] ^ input[start + 1] ^ send0[nn_2layer + switch_start[i]] ^ send1[nn_2layer + switch_start[i]];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // Second layer.
                input0[start + 1] = input[start + 1] ^ send0[nn_1layer + switch_start[i]];
                input0[start + 2] = input[start + 1] ^ send1[nn_1layer + switch_start[i]];
                d[nn_1layer + switch_start[i]] = input[start + 1] ^ input[start + 2] ^ send0[nn_1layer + switch_start[i]] ^ send1[nn_1layer + switch_start[i]];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

                // Third layer.
                input0[start] = input[start] ^ send0[nn_0layer + switch_start[i]];
                input0[start + 1] = input[start] ^ send1[nn_0layer + switch_start[i]];
                d[nn_0layer + switch_start[i]] = input[start] ^ input[start + 1] ^ send0[nn_0layer + switch_start[i]] ^ send1[nn_0layer + switch_start[i]];

                // input   = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                input[start + 2] = input0[start + 2];
                // input   = std::move(input0);
                // std::swap(input, input0);

            } else if (group_size[i] == 2) { // Modify the third-layer input directly, while d still needs to be modified in the middle layer.
                input0[start] = input[start] ^ send0[nn_1layer + switch_start[i]];
                input0[start + 1] = input[start] ^ send1[nn_1layer + switch_start[i]];
                d[nn_1layer + switch_start[i]] = input[start] ^ input[start + 1] ^ send0[nn_1layer + switch_start[i]] ^ send1[nn_1layer + switch_start[i]];

                input = input0;
                input[start] = input0[start];
                input[start + 1] = input0[start + 1];
                // input   = std::move(input0);
                // std::swap(input, input0);
            }
            start += group_size[i];
        }
        // cout<<"wuwuwu"<<endl;

        // Compute the last 2logN-2 layers.
        loog = loog / 2;

        for (u8 i = 0; i < nn - 2; i++) {
            u32 nn_1iswitchnum = (nn + 1 + i) * switch_num_layer;

            for (int ite = 0; ite < loog; ite++) {
                // Write group_size.
                group_size[ite] = group_size[ite << 1] + group_size[(ite << 1) + 1];
                switch_start[ite] = switch_start[ite << 1];
            }
            // group_size.resize(loog);
            // switch_start.resize(loog);

            start = 0;
            for (u32 j = 0; j < loog; j++) { // One subnetwork in a layer.

                for (u32 k = 0; k < group_size[j] / 2; k++) {                                               // Element computation in the subnetwork.
                    input0[start + 2 * k] = input[start + k] ^ send0[nn_1iswitchnum + switch_start[j] + k]; // u0 = x0 − r0,
                    input0[start + 2 * k + 1] = input[start + k] ^ send1[nn_1iswitchnum + switch_start[j] + k];
                    d[nn_1iswitchnum + switch_start[j] + k] = input[start + k] ^ input[start + k + (group_size[j] / 2)] ^
                        send0[nn_1iswitchnum + switch_start[j] + k] ^ send1[nn_1iswitchnum + switch_start[j] + k];
                }
                if (group_size[j] % 2 == 1) {
                    input0[start + group_size[j] - 1] = input[start + group_size[j] - 1];
                }

                start += group_size[j];
            }
            loog = loog / 2;
            input = input0;
            // input   = std::move(input0);
            // std::swap(input, input0);
        }
        output = input;
        // Final communication phase.
        coproto::sync_wait(chl.send(d));
    }
}
