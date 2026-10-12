// © 2022 Visa.
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
// documentation files (the "Software"), to deal in the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all copies or substantial portions of the
// Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
// WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
// COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
// OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// Adapted from upstream volePSI/RsOpprf.cpp (Visa MIT license above).
#include "opprf_impl.h"

namespace fpsi {
    void LargeOpprfSender::oprfEval(span<const block> X, span<u8> out, span<const u8> add, u64 m, u64 numThreads)
    {
        auto n = X.size();
        auto nm = detail::checkedBytes(n, m);

        if (out.size() < nm || add.size() < nm)
            throw RTE_LOC;

        std::unique_ptr<u8[]> tempPtr;
        span<block> oprfOutput;
        if (m == sizeof(block) && !(out.data() == add.data()))
            oprfOutput = span<block>((block *)out.data(), n);
        else {
            tempPtr.reset(new u8[detail::checkedBytes(n, sizeof(block))]);
            oprfOutput = span<block>((block *)tempPtr.get(), n);
        }

        mOprfSender.eval(X, oprfOutput, numThreads);

        if (m <= sizeof(block)) {
            // short string case
            for (u64 i = 0, ij = 0; i < n; ++i) {
                for (u64 j = 0; j < m; ++j, ++ij) {
                    out[ij] = add[ij] ^ oprfOutput[i].get<u8>(j);
                }
            }
        } else if (m > sizeof(block)) {
            // long string case

            std::vector<block> buffer(oc::divCeil(m, sizeof(block)));
            auto buffPtr = (u8 *)buffer.data();
            for (u64 i = 0, ij = 0; i < n; ++i) {
                oc::AES enc(oprfOutput[i]);
                enc.ecbEncCounterMode(0, buffer);
                for (u64 j = 0; j < m; ++j, ++ij)
                    out[ij] = add[ij] ^ buffPtr[j];
            }
        }
    }

    Proto LargeOpprfSender::send(u64 recverSize, span<const block> X, MatrixView<u8> val, PRNG &prng, u64 numThreads, Socket &chl)
    {
        auto n = u64{};
        auto m = u64{};
        auto nm = u64{};
        auto hashingSeed = block{};
        auto type = PaxosParam::Binary;
        auto diffPtr = std::unique_ptr<u8[]>{};
        auto diffU8 = span<u8>{};

        setTimePoint("LargeOpprfSender::send begin");
        n = X.size();
        m = val.cols();
        nm = detail::checkedBytes(n, m);
        detail::checkedBytes(n, sizeof(block));
        detail::checkedBytes(recverSize, sizeof(block));
        detail::chunkElements(mChunkBytes, 1);

        if (val.rows() != n)
            throw std::invalid_argument("OPPRF key/value row counts differ");

        if (!n || !m || !recverSize)
            throw std::invalid_argument("OPPRF requires nonempty inputs and nonzero value width");

        hashingSeed = prng.get();
        co_await (chl.send(std::move(hashingSeed)));

        type = m % sizeof(block) ? PaxosParam::Binary : PaxosParam::GF128;
        mPaxos.init(n, 1 << 14, 3, 40, type, hashingSeed);

        if (mTimer)
            mOprfSender.setTimer(*mTimer);

        co_await (mOprfSender.send(recverSize, prng, chl, numThreads));

        diffPtr.reset(new u8[nm]);
        diffU8 = span<u8>(diffPtr.get(), nm);

        if (diffU8.size() < val.size())
            throw RTE_LOC;
        memset(diffU8.data(), 0, val.size());
        oprfEval(X, diffU8, val, m, numThreads);

        mPaxosByteWidth = m;
        mP.resize(detail::checkedBytes(mPaxos.size(), m));

        setTimePoint("LargeOpprfSender::send eval");

        if (m == 16) {
            // special case for block.
            auto diff = span<block>((block *)diffU8.data(), n);
            auto P = span<block>((block *)mP.data(), mPaxos.size());

            // double check the bounds.
            if ((u8 *)(diff.data() + diff.size()) > diffU8.data() + diffU8.size())
                throw RTE_LOC;
            if ((u8 *)(mP.data() + mP.size()) > mP.data() + mP.size())
                throw RTE_LOC;

            // Vendored Baxos shares this PRNG across encoding workers. Its
            // state is not thread-safe, including direct GF128 backfill draws.
            // Keep randomized encoding serial; other stages retain numThreads.
            mPaxos.solve<block>(X, diff, P, &prng, 1);
        } else {
            // throw RTE_LOC;
            auto diff = MatrixView<u8>(diffPtr.get(), n, m);
            auto P = MatrixView<u8>((u8 *)mP.data(), mPaxos.size(), m);

            // double check the bounds.
            if (diff.data() + diff.size() > diffU8.data() + diffU8.size())
                throw RTE_LOC;
            if (mP.data() + mP.size() > mP.data() + mP.size())
                throw RTE_LOC;

            mPaxos.solve<u8>(X, diff, P, &prng, 1);
        }
        diffPtr.reset();

        setTimePoint("LargeOpprfSender::send paxos solve");

        for (u64 offset = 0; offset < mP.size();) {
            auto count = std::min<u64>(mP.size() - offset, mChunkBytes);
            auto subP = span<u8>(mP.data() + offset, count);
            // A named lvalue send completes before this table can be reused.
            co_await chl.send(subP);
            offset += count;
        }

        setTimePoint("LargeOpprfSender::send end");
    }

    void LargeOpprfSender::eval(span<const block> val, MatrixView<u8> output, u64 numThreads)
    {
        auto m = output.cols();

        if (!m || mPaxosByteWidth != m || val.size() != output.rows())
            throw std::invalid_argument("OPPRF eval output dimensions do not match initialized state");
        detail::checkedBytes(output.rows(), m);
        if (val.empty())
            return;

        MatrixView<u8> P(mP.data(), mPaxos.size(), output.cols());

        setTimePoint("LargeOpprfSender::eval begin");
        if (m == sizeof(block)) {
            // spacial case for block
            auto oo = span<block>((block *)output.data(), output.rows());
            auto pp = span<block>((block *)P.data(), P.rows());
            mPaxos.decode<block>(val, oo, pp);
        } else {
            mPaxos.decode<u8>(val, output, P);
        }

        setTimePoint("LargeOpprfSender::eval paxos decode");

        // output +=  oprf(val)
        auto outU8 = span<u8>(output.data(), output.size());
        oprfEval(val, outU8, outU8, m, numThreads);

        setTimePoint("LargeOpprfSender::eval end");
    }

    Proto LargeOpprfReceiver::receive(u64 senderSize, span<const block> values, MatrixView<u8> outputs, PRNG &prng, u64 numThreads, Socket &chl)
    {
        auto n = u64{ 0 }, m = u64{ 0 };
        auto paxos = Baxos{};
        auto type = PaxosParam::Binary;
        auto temp = oc::AlignedUnVector<block>{};
        auto oprfOutput = span<block>{};
        auto p = Matrix<u8>{};

        setTimePoint("LargeOpprfReceiver::receive begin");

        if (values.empty() || values.size() != outputs.rows() || !outputs.cols() || !senderSize)
            throw std::invalid_argument("OPPRF requires nonempty inputs and matching output dimensions");
        detail::chunkElements(mChunkBytes, 1);
        n = outputs.rows();
        m = outputs.cols();
        detail::checkedBytes(n, m);
        detail::checkedBytes(n, sizeof(block));
        detail::checkedBytes(senderSize, sizeof(block));

        co_await chl.recv(paxos.mSeed);
        type = m % sizeof(block) ? PaxosParam::Binary : PaxosParam::GF128;
        paxos.init(senderSize, 1 << 14, 3, 40, type, paxos.mSeed);

        if (mTimer)
            mOprfReceiver.setTimer(*mTimer);

        auto endOffset = reinterpret_cast<std::uintptr_t>(outputs.data() + outputs.size()) % sizeof(block);
        // Expansion proceeds forwards. An aligned tail scratch area is safe
        // only if writing a row cannot overwrite the next unread OPRF block.
        if (m >= sizeof(block) && m - sizeof(block) >= endOffset) {
            // reuse memory. extra logic to make sure oprfOutput is properly aligned.
            auto ptr = outputs.data() + outputs.size();
            auto offset = (u64)ptr % sizeof(block);
            if (offset)
                ptr -= offset;
            assert((((u64)ptr) % sizeof(block)) == 0);

            oprfOutput = span<block>((block *)(ptr)-n, n);

            assert((void *)(oprfOutput.data() + oprfOutput.size()) <= (void *)(outputs.data() + outputs.size()));
        } else {
            temp.resize(n);
            oprfOutput = { temp.data(), temp.size() };
        }

        co_await mOprfReceiver.receive(values, oprfOutput, prng, chl, numThreads);

        detail::checkedBytes(paxos.size(), m);
        p.resize(paxos.size(), m, oc::AllocType::Uninitialized);

        for (u64 offset = 0; offset < p.size();) {
            auto count = std::min<u64>(p.size() - offset, mChunkBytes);
            auto subP = span<u8>(p.data() + offset, count);
            co_await chl.recv(subP);
            offset += count;
        }

        setTimePoint("LargeOpprfReceiver::receive recv");

        if (m == sizeof(block)) {
            auto v = span<block>((block *)values.data(), n);
            auto o = span<block>((block *)outputs.data(), n);
            auto pp = span<block>((block *)p.data(), paxos.size());

            paxos.mAddToDecode = true;
            paxos.decode<block>(v, o, pp, numThreads);
        } else {
            if (m > sizeof(block)) {
                std::vector<block> buffer(oc::divCeil(m, sizeof(block)));
                auto buffPtr = (u8 *)buffer.data();

                for (u64 i = 0; i < n; ++i) {
                    oc::AES enc(oprfOutput[i]);
                    enc.ecbEncCounterMode(0, buffer);
                    memcpy(&outputs(i * m), buffPtr, m);
                }
            } else {
                for (u64 i = 0; i < n; ++i)
                    memcpy(&outputs(i * m), &oprfOutput[i], m);
            }

            temp = {};
            paxos.mAddToDecode = true;
            paxos.decode<u8>(values, outputs, p, numThreads);
        }
        temp = {};

        setTimePoint("LargeOpprfReceiver::receive end");

        co_return;
    }
} // namespace fpsi
