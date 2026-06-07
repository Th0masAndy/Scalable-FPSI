#include <coproto/Socket/Socket.h>
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cstdint>
#include <libOTe/Triple/SilentOtTriple/SilentOtTriple.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtReceiver.h>
#include <libOTe/TwoChooseOne/Silent/SilentOtExtSender.h>
#include <utils.h>

using namespace oc;

const u64 M = 4;

#define USE_CHEETAH true

class NcoOTSender {
public:
    NcoOTSender(u64 size)
    {
        prng = new PRNG(oc::sysRandomSeed());
        mNum = size;
        mOt = new SilentOtExtSender();
        mOt->configure(mNum * 4);
    }

    ~NcoOTSender()
    {
        delete prng;
    }

    void send(u8 **messages, Socket &chl);

private:
    PRNG *prng;
    u64 mNum;
    SilentOtExtSender *mOt;
};

class NcoOTRecver {
public:
    NcoOTRecver(u64 size)
    {
        prng = new PRNG(oc::sysRandomSeed());
        mNum = size;
        mOt = new SilentOtExtReceiver();
        mOt->configure(mNum * 4);
    }

    ~NcoOTRecver()
    {
        delete prng;
    }

    void recv(u8 *outs, u8 *choices, Socket &chl);

private:
    PRNG *prng;
    u64 mNum;
    SilentOtExtReceiver *mOt;
};

// this code is from opencheetah (https://github.com/Alibaba-Gemini-Lab/OpenCheetah) and replaced the OT module with libOTe
// Cheetah's variant MillionaireProtocol when USE_CHEETAH=1
class MillionaireProtocolSender {
public:
    NcoOTSender *otpack;
    SilentOtTriple *triple_gen;
    PRNG *prng;
    int party;
    int l, r, log_alpha, beta, beta_pow;
    int num_digits, num_triples_corr, num_triples_std, log_num_digits;
    int num_triples;
    int num_cmps;
    uint8_t mask_beta, mask_r;
    u64 num_triples_round;

    MillionaireProtocolSender(int num_cmps, int bitlength = 32, int radix_base = M);

    void configure(int bitlength, int radix_base = M);

    ~MillionaireProtocolSender();

    void drelu(uint8_t *res, uint64_t *data, osuCrypto::Socket &chl);

    void compare(uint8_t *res, uint64_t *data, osuCrypto::Socket &chl, bool greater_than = true, bool equality = false, int radix_base = M);

    void set_leaf_ot_messages(uint8_t *ot_messages, uint8_t digit, int N, uint8_t mask_cmp, uint8_t mask_eq, bool greater_than, bool eq = true);

    /**************************************************************************************************
     *                         AND computation related functions
     **************************************************************************************************/

    void traverse_and_compute_ANDs(SilentOtTriple &triple_gen, int num_cmps, uint8_t *leaf_res_eq, uint8_t *leaf_res_cmp, osuCrypto::Socket &chl);

    static void AND_step_1(
        uint8_t *ei, // evaluates batch of 8 ANDs
        uint8_t *fi,
        uint8_t *xi,
        uint8_t *yi,
        uint8_t *ai,
        uint8_t *bi,
        int num_ANDs);

    static void AND_step_2(
        uint8_t *zi, // evaluates batch of 8 ANDs
        uint8_t *e,
        uint8_t *f,
        uint8_t *ei,
        uint8_t *fi,
        uint8_t *ai,
        uint8_t *bi,
        uint8_t *ci,
        int num_ANDs);
};

class MillionaireProtocolRecver {
public:
    NcoOTRecver *otpack;
    SilentOtTriple *triple_gen;
    PRNG *prng;
    int party;
    int l, r, log_alpha, beta, beta_pow;
    int num_digits, num_triples_corr, num_triples_std, log_num_digits;
    int num_triples;
    int num_cmps;
    uint8_t mask_beta, mask_r;
    u64 num_triples_round;

    MillionaireProtocolRecver(int num_cmps, int bitlength = 32, int radix_base = M);

    void configure(int bitlength, int radix_base = M);

    ~MillionaireProtocolRecver();

    void drelu(uint8_t *res, uint64_t *data, osuCrypto::Socket &chl);

    void compare(uint8_t *res, uint64_t *data, osuCrypto::Socket &chl, bool greater_than = true, bool equality = false, int radix_base = M);

    void set_leaf_ot_messages(uint8_t *ot_messages, uint8_t digit, int N, uint8_t mask_cmp, uint8_t mask_eq, bool greater_than, bool eq = true);

    /**************************************************************************************************
     *                         AND computation related functions
     **************************************************************************************************/

    void traverse_and_compute_ANDs(SilentOtTriple &triple_gen, int num_cmps, uint8_t *leaf_res_eq, uint8_t *leaf_res_cmp, osuCrypto::Socket &chl);

    static void AND_step_1(
        uint8_t *ei, // evaluates batch of 8 ANDs
        uint8_t *fi,
        uint8_t *xi,
        uint8_t *yi,
        uint8_t *ai,
        uint8_t *bi,
        int num_ANDs);
    static void AND_step_2(
        uint8_t *zi, // evaluates batch of 8 ANDs
        uint8_t *e,
        uint8_t *f,
        uint8_t *ei,
        uint8_t *fi,
        uint8_t *ai,
        uint8_t *bi,
        uint8_t *ci,
        int num_ANDs);
};
