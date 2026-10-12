#pragma once
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

#include "oprf_impl.h"

namespace fpsi {
    // Same correction-table construction as RsOpprf. Each value has a fixed byte
    // width per invocation. Programmed keys must be distinct. This implementation
    // retains the full OKVS table to support subsequent sender.eval() calls.
    // Randomized OKVS encoding is serial: this Baxos version shares a non-thread-
    // safe PRNG between workers. Other stages honor numThreads as in RsOpprf.
    class LargeOpprfSender : public oc::TimerAdapter {
    public:
        LargeOprfSender mOprfSender;
        oc::AlignedUnVector<u8> mP;
        u64 mPaxosByteWidth = 0;
        Baxos mPaxos;
        u64 mChunkBytes = u64{ 1 } << 30;

        void setMultType(oc::MultType type)
        {
            mOprfSender.setMultType(type);
        }
        Proto send(u64 receiverSize, span<const block> keys, MatrixView<u8> values, PRNG &prng, u64 numThreads, Socket &chl);
        Proto send(u64 receiverSize, span<const block> keys, span<block> values, PRNG &prng, u64 numThreads, Socket &chl)
        {
            return send<block>(receiverSize, keys, values, prng, numThreads, chl);
        }

        template <class T>
        Proto send(u64 receiverSize, span<const block> keys, MatrixView<T> values, PRNG &prng, u64 numThreads, Socket &chl)
        {
            auto width = detail::checkedBytes(values.cols(), sizeof(T));
            detail::checkedBytes(values.rows(), width);
            return send(receiverSize, keys, MatrixView<u8>((u8 *)values.data(), values.rows(), width), prng, numThreads, chl);
        }
        template <class T>
        Proto send(u64 receiverSize, span<const block> keys, span<T> values, PRNG &prng, u64 numThreads, Socket &chl)
        {
            return send(receiverSize, keys, MatrixView<T>(values.data(), values.size(), 1), prng, numThreads, chl);
        }

        void eval(span<const block> keys, MatrixView<u8> outputs, u64 numThreads = 1);
        void eval(span<const block> keys, span<block> outputs, u64 numThreads = 1)
        {
            eval(keys, MatrixView<block>(outputs.data(), outputs.size(), 1), numThreads);
        }
        template <class T>
        void eval(span<const block> keys, MatrixView<T> outputs, u64 numThreads = 1)
        {
            auto width = detail::checkedBytes(outputs.cols(), sizeof(T));
            detail::checkedBytes(outputs.rows(), width);
            eval(keys, MatrixView<u8>((u8 *)outputs.data(), outputs.rows(), width), numThreads);
        }
        template <class T>
        T eval(block key)
        {
            T out;
            eval<T>({ &key, 1 }, MatrixView<T>(&out, 1, 1), 1);
            return out;
        }

    private:
        void oprfEval(span<const block> keys, span<u8> out, span<const u8> add, u64 width, u64 numThreads);
    };

    class LargeOpprfReceiver : public oc::TimerAdapter {
    public:
        LargeOprfReceiver mOprfReceiver;
        u64 mChunkBytes = u64{ 1 } << 30;

        void setMultType(oc::MultType type)
        {
            mOprfReceiver.setMultType(type);
        }
        Proto receive(u64 senderSize, span<const block> keys, MatrixView<u8> outputs, PRNG &prng, u64 numThreads, Socket &chl);
        Proto receive(u64 senderSize, span<const block> keys, span<block> outputs, PRNG &prng, u64 numThreads, Socket &chl)
        {
            return receive<block>(senderSize, keys, outputs, prng, numThreads, chl);
        }
        template <class T>
        Proto receive(u64 senderSize, span<const block> keys, MatrixView<T> outputs, PRNG &prng, u64 numThreads, Socket &chl)
        {
            auto width = detail::checkedBytes(outputs.cols(), sizeof(T));
            detail::checkedBytes(outputs.rows(), width);
            return receive(senderSize, keys, MatrixView<u8>((u8 *)outputs.data(), outputs.rows(), width), prng, numThreads, chl);
        }
        template <class T>
        Proto receive(u64 senderSize, span<const block> keys, span<T> outputs, PRNG &prng, u64 numThreads, Socket &chl)
        {
            return receive(senderSize, keys, MatrixView<T>(outputs.data(), outputs.size(), 1), prng, numThreads, chl);
        }
    };
} // namespace fpsi
